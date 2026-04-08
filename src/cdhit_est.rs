// Native cd-hit-est greedy incremental clustering.
// Based on: Li et al. (2001) Bioinformatics 17(3):282-283
// and the cd-hit-est C++ source (cdhit-est.c++, cdhit-common.c++).
//
// Performance design:
//   • Word table is a Vec<Vec<usize>> of size 4^word_len (direct array index, no HashMap).
//   • k-mer counting and diagonal voting use pre-allocated Vec buffers reset via dirty lists.
//   • All hot-path allocations are eliminated; buffers live for the lifetime of the run.

use anyhow::{Context, Result};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

// Alignment band half-width when scanning offsets (matches cd-hit-est default -b 20).
const BAND_WIDTH: usize = 20;
// Word length for the fast diagonal-vote pass (4^6 = 4096 table entries).
const DIAG_KMER: usize = 6;
const DIAG_TABLE_SIZE: usize = 4096; // 4^DIAG_KMER
// Maximum representative indices stored per k-mer entry (caps ubiquitous k-mers).
const MAX_KMER_REPS: usize = 128;

struct Record {
    id: String,
    seq: Vec<u8>, // A=0, C=1, G=2, T=3, N=4
}

struct ClusterMember {
    orig_idx: usize,
    identity: f64,
    strand: bool, // true = forward, false = reverse complement
}

/// Reusable per-query working buffers — no allocations in the hot path.
struct WorkBufs {
    // Shared k-mer counts: counts[cluster_idx] = number of shared word_len k-mers.
    counts: Vec<usize>,
    counts_dirty: Vec<usize>,

    // Query DIAG_KMER position table: qkpos[kmer_val] = positions in query.
    // Size = 4^DIAG_KMER = 4096.
    qkpos: Vec<Vec<u16>>,
    qkpos_dirty: Vec<usize>,

    // Diagonal vote array: diag[d] = number of DIAG_KMER matches at offset d.
    // Index = rep_pos - query_pos (always >= 0 since rep is always >= query).
    diag: Vec<u32>,
    diag_dirty: Vec<usize>,

    // Scratch for extracted k-mers.
    qkmers_main: Vec<i32>,
    qkmers_diag: Vec<i32>,
}

impl WorkBufs {
    fn new(max_clusters: usize, max_seq_len: usize) -> Self {
        WorkBufs {
            counts: vec![0; max_clusters],
            counts_dirty: Vec::with_capacity(1024),
            qkpos: vec![Vec::new(); DIAG_TABLE_SIZE],
            qkpos_dirty: Vec::with_capacity(512),
            diag: vec![0; max_seq_len + 1],
            diag_dirty: Vec::with_capacity(4096),
            qkmers_main: Vec::with_capacity(max_seq_len),
            qkmers_diag: Vec::with_capacity(max_seq_len),
        }
    }

    fn ensure_counts(&mut self, n_clusters: usize) {
        if self.counts.len() < n_clusters {
            self.counts.resize(n_clusters, 0);
        }
    }

    fn ensure_diag(&mut self, window: usize) {
        if self.diag.len() <= window {
            self.diag.resize(window + 1, 0);
        }
    }
}

/// Cluster nucleotide sequences by identity using the cd-hit-est greedy algorithm.
pub fn cluster(input: &Path, output: &Path, threshold: f64, log_file: &Path) -> Result<()> {
    let log = File::create(log_file)
        .with_context(|| format!("Cannot create log: {}", log_file.display()))?;
    let mut log_w = BufWriter::new(log);

    writeln!(log_w, "BaitBench native cd-hit-est clustering")?;
    writeln!(log_w, "Input:     {}", input.display())?;
    writeln!(log_w, "Threshold: {:.2}", threshold)?;

    let word_len = word_len_for_threshold(threshold);
    let table_size = 4usize.pow(word_len as u32);
    writeln!(log_w, "Word length: {} ({} table entries)", word_len, table_size)?;

    let records = load_fasta(input)?;
    let n = records.len();
    writeln!(log_w, "Sequences: {}", n)?;

    if n == 0 {
        File::create(output)
            .with_context(|| format!("Cannot create: {}", output.display()))?;
        File::create(clstr_path(output))
            .with_context(|| format!("Cannot create .clstr"))?;
        log_w.flush()?;
        return Ok(());
    }

    let max_len = records.iter().map(|r| r.seq.len()).max().unwrap_or(0);

    // Sort indices by decreasing sequence length.
    let mut sorted_idx: Vec<usize> = (0..n).collect();
    sorted_idx.sort_by(|&a, &b| records[b].seq.len().cmp(&records[a].seq.len()));

    let mut clusters: Vec<Vec<ClusterMember>> = Vec::new();
    let mut is_rep: Vec<bool> = vec![false; n];
    let mut rep_orig: Vec<usize> = Vec::new();

    // Word table: direct Vec indexed by k-mer value, each entry is a list of cluster indices.
    let mut word_table: Vec<Vec<usize>> = vec![Vec::new(); table_size];

    let mut buf = WorkBufs::new(n / 10 + 1, max_len * 2 + 1);

    let mut processed = 0usize;
    for &si in &sorted_idx {
        processed += 1;
        if processed % 5000 == 0 {
            log::info!(
                "Clustering: {}/{} sequences, {} clusters",
                processed,
                n,
                clusters.len()
            );
        }

        let qlen = records[si].seq.len();
        buf.ensure_counts(clusters.len());

        let req_base = required_kmers(qlen, threshold, word_len);

        let mut best: Option<(usize, f64, bool)> = None;

        'strand_loop: for is_forward in [true, false] {
            let rc_buf: Vec<u8>;
            let qseq: &[u8] = if is_forward {
                &records[si].seq
            } else {
                rc_buf = reverse_complement(&records[si].seq);
                &rc_buf
            };

            extract_kmers_into(qseq, word_len, &mut buf.qkmers_main);

            // Count invalid k-mers and saturated k-mers to adjust required threshold.
            let mut skip = 0usize;
            let mut saturated = 0usize;
            for &kv in &buf.qkmers_main {
                if kv < 0 {
                    skip += 1;
                } else {
                    let entry = &word_table[kv as usize];
                    if entry.len() >= MAX_KMER_REPS {
                        saturated += 1;
                    }
                }
            }
            let req = req_base.saturating_sub(skip + saturated).max(1);

            // Count shared k-mers per representative.
            for &kv in &buf.qkmers_main {
                if kv < 0 {
                    continue;
                }
                let entry = &word_table[kv as usize];
                if entry.len() >= MAX_KMER_REPS {
                    continue; // saturated; already subtracted from req
                }
                for &cidx in entry {
                    if buf.counts[cidx] == 0 {
                        buf.counts_dirty.push(cidx);
                    }
                    buf.counts[cidx] += 1;
                }
            }

            // Check candidates that passed the k-mer prefilter.
            for &cidx in &buf.counts_dirty {
                let shared = buf.counts[cidx];
                buf.counts[cidx] = 0;

                if shared < req {
                    continue;
                }

                let rep_seq = &records[rep_orig[cidx]].seq;
                let identity = compute_identity(qseq, rep_seq, &mut buf);

                if identity >= threshold {
                    match best {
                        None => best = Some((cidx, identity, is_forward)),
                        Some((bi, _, _)) if cidx < bi => {
                            best = Some((cidx, identity, is_forward))
                        }
                        _ => {}
                    }
                }
            }
            buf.counts_dirty.clear();

            if best.is_some() {
                break 'strand_loop;
            }
        }

        if let Some((cidx, identity, is_fwd)) = best {
            clusters[cidx].push(ClusterMember {
                orig_idx: si,
                identity,
                strand: is_fwd,
            });
        } else {
            // New representative.
            let cidx = clusters.len();
            clusters.push(vec![ClusterMember {
                orig_idx: si,
                identity: 1.0,
                strand: true,
            }]);
            is_rep[si] = true;
            rep_orig.push(si);

            // Index k-mers; cap each entry at MAX_KMER_REPS.
            extract_kmers_into(&records[si].seq, word_len, &mut buf.qkmers_main);
            for kv in &buf.qkmers_main {
                if *kv >= 0 {
                    let entry = &mut word_table[*kv as usize];
                    if entry.len() < MAX_KMER_REPS {
                        entry.push(cidx);
                    }
                }
            }
        }
    }

    writeln!(log_w, "Clusters:  {}", clusters.len())?;
    log_w.flush()?;

    // Write representatives in original input order.
    let out_file = File::create(output)
        .with_context(|| format!("Cannot create: {}", output.display()))?;
    let mut out_w = BufWriter::new(out_file);
    for i in 0..n {
        if is_rep[i] {
            writeln!(out_w, ">{}", records[i].id)?;
            for chunk in records[i].seq.chunks(80) {
                let line: String = chunk.iter().map(|&b| decode_base(b)).collect();
                writeln!(out_w, "{}", line)?;
            }
        }
    }
    out_w.flush()?;

    // Write cluster file.
    let cp = clstr_path(output);
    let clstr_file = File::create(&cp)
        .with_context(|| format!("Cannot create: {}", cp.display()))?;
    let mut clstr_w = BufWriter::new(clstr_file);
    for (cidx, members) in clusters.iter().enumerate() {
        writeln!(clstr_w, ">Cluster {}", cidx)?;
        let rep_short = first_word(&records[members[0].orig_idx].id);
        let rep_len = records[members[0].orig_idx].seq.len();
        writeln!(clstr_w, "0\t{}nt, >{}... *", rep_len, rep_short)?;
        for (mi, member) in members[1..].iter().enumerate() {
            let mem_short = first_word(&records[member.orig_idx].id);
            let mem_len = records[member.orig_idx].seq.len();
            let strand = if member.strand { "+" } else { "-" };
            writeln!(
                clstr_w,
                "{}\t{}nt, >{}... at {}/{:.2}%",
                mi + 1,
                mem_len,
                mem_short,
                strand,
                member.identity * 100.0
            )?;
        }
    }
    clstr_w.flush()?;

    log::info!(
        "Clustering complete: {} sequences → {} clusters",
        n,
        clusters.len()
    );
    Ok(())
}

// ---------------------------------------------------------------------------
// Identity calculation with pre-allocated buffers
// ---------------------------------------------------------------------------

/// Global identity: identical_positions / (min_len − N_count_in_shorter).
/// `query` must not be longer than `rep`.
fn compute_identity(query: &[u8], rep: &[u8], buf: &mut WorkBufs) -> f64 {
    let qlen = query.len();
    let rlen = rep.len();
    debug_assert!(qlen <= rlen);

    let n_count: usize = query.iter().filter(|&&b| b == 4).count();
    let denom = qlen - n_count;
    if denom == 0 {
        return 0.0;
    }

    let offset = if qlen == rlen {
        0
    } else {
        find_best_offset(query, rep, buf)
    };

    let matches = count_matches(query, &rep[offset..offset + qlen]);
    matches as f64 / denom as f64
}

/// Find the starting offset in `rep` that maximises matching positions with `query`.
/// Uses pre-allocated buffers — no per-call allocation.
fn find_best_offset(query: &[u8], rep: &[u8], buf: &mut WorkBufs) -> usize {
    let qlen = query.len();
    let rlen = rep.len();
    let window_range = rlen - qlen; // valid offsets: 0..=window_range

    // Small window: direct scan.
    if window_range <= 2 * BAND_WIDTH {
        return (0..=window_range)
            .max_by_key(|&off| count_matches(query, &rep[off..off + qlen]))
            .unwrap_or(0);
    }

    if qlen < DIAG_KMER {
        let m0 = count_matches(query, &rep[..qlen]);
        let me = count_matches(query, &rep[window_range..]);
        return if m0 >= me { 0 } else { window_range };
    }

    // Build query DIAG_KMER position table using pre-allocated qkpos.
    extract_kmers_into(query, DIAG_KMER, &mut buf.qkmers_diag);
    for (pos, &kv) in buf.qkmers_diag.iter().enumerate() {
        if kv >= 0 {
            let idx = kv as usize;
            if buf.qkpos[idx].is_empty() {
                buf.qkpos_dirty.push(idx);
            }
            if buf.qkpos[idx].len() < 32 {
                buf.qkpos[idx].push(pos as u16);
            }
        }
    }

    // Vote for diagonals using rep DIAG_KMER k-mers.
    buf.ensure_diag(window_range);
    extract_kmers_into(rep, DIAG_KMER, &mut buf.qkmers_main); // reuse scratch
    for (rpos, &rv) in buf.qkmers_main.iter().enumerate() {
        if rv < 0 {
            continue;
        }
        let idx = rv as usize;
        if buf.qkpos[idx].is_empty() {
            continue;
        }
        for &qpos in &buf.qkpos[idx] {
            let qpos = qpos as usize;
            if rpos >= qpos {
                let diag = rpos - qpos;
                if diag <= window_range {
                    if buf.diag[diag] == 0 {
                        buf.diag_dirty.push(diag);
                    }
                    buf.diag[diag] += 1;
                }
            }
        }
    }

    // Find best diagonal.
    let best_diag = buf
        .diag_dirty
        .iter()
        .max_by_key(|&&d| buf.diag[d])
        .copied()
        .unwrap_or(0)
        .min(window_range);

    // Reset diag buffer.
    for &d in &buf.diag_dirty {
        buf.diag[d] = 0;
    }
    buf.diag_dirty.clear();

    // Reset qkpos buffer.
    for &idx in &buf.qkpos_dirty {
        buf.qkpos[idx].clear();
    }
    buf.qkpos_dirty.clear();

    // Scan a band around best_diag.
    let lo = best_diag.saturating_sub(BAND_WIDTH);
    let hi = (best_diag + BAND_WIDTH).min(window_range);
    (lo..=hi)
        .max_by_key(|&off| count_matches(query, &rep[off..off + qlen]))
        .unwrap_or(best_diag)
}

#[inline]
fn count_matches(a: &[u8], b: &[u8]) -> usize {
    a.iter()
        .zip(b.iter())
        .filter(|(&x, &y)| x < 4 && y < 4 && x == y)
        .count()
}

// ---------------------------------------------------------------------------
// K-mer extraction and encoding
// ---------------------------------------------------------------------------

fn word_len_for_threshold(threshold: f64) -> usize {
    if threshold >= 0.90 {
        8
    } else if threshold >= 0.80 {
        6
    } else if threshold >= 0.75 {
        5
    } else {
        4
    }
}

fn required_kmers(len_eff: usize, threshold: f64, word_len: usize) -> usize {
    if len_eff < word_len {
        return 1;
    }
    let total = (len_eff - word_len + 1) as isize;
    let max_mm = ((1.0 - threshold) * len_eff as f64 + 0.5) as isize;
    (total - max_mm * word_len as isize).max(1) as usize
}

/// Extract k-mers into a pre-allocated Vec (cleared first) using a rolling hash.
fn extract_kmers_into(seq: &[u8], k: usize, out: &mut Vec<i32>) {
    out.clear();
    let n = seq.len();
    if n < k {
        return;
    }
    out.reserve(n - k + 1);

    let kpow = 4u32.pow((k - 1) as u32);
    let mut kmer = 0u32;
    let mut n_count = 0usize;

    for &b in &seq[..k] {
        if b == 4 {
            n_count += 1;
        }
        kmer = kmer * 4 + b.min(3) as u32;
    }
    out.push(if n_count == 0 { kmer as i32 } else { -1 });

    for i in k..n {
        let leaving = seq[i - k];
        let entering = seq[i];
        if leaving == 4 {
            n_count -= 1;
        }
        if entering == 4 {
            n_count += 1;
        }
        kmer = (kmer - leaving.min(3) as u32 * kpow) * 4 + entering.min(3) as u32;
        out.push(if n_count == 0 { kmer as i32 } else { -1 });
    }
}

// ---------------------------------------------------------------------------
// Sequence encoding / decoding
// ---------------------------------------------------------------------------

fn encode_base(b: u8) -> u8 {
    match b | 0x20 {
        b'a' => 0,
        b'c' => 1,
        b'g' => 2,
        b't' | b'u' => 3,
        _ => 4,
    }
}

fn decode_base(b: u8) -> char {
    match b {
        0 => 'A',
        1 => 'C',
        2 => 'G',
        3 => 'T',
        _ => 'N',
    }
}

fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            0 => 3,
            1 => 2,
            2 => 1,
            3 => 0,
            _ => 4,
        })
        .collect()
}

// ---------------------------------------------------------------------------
// I/O helpers
// ---------------------------------------------------------------------------

fn load_fasta(input: &Path) -> Result<Vec<Record>> {
    let file = File::open(input)
        .with_context(|| format!("Cannot open FASTA: {}", input.display()))?;
    let reader = BufReader::new(file);

    let mut records: Vec<Record> = Vec::new();
    let mut cur_id: Option<String> = None;
    let mut cur_seq: Vec<u8> = Vec::new();

    for line in reader.lines() {
        let line = line?;
        let t = line.trim_end();
        if let Some(rest) = t.strip_prefix('>') {
            if let Some(id) = cur_id.take() {
                records.push(Record { id, seq: std::mem::take(&mut cur_seq) });
            }
            cur_id = Some(rest.to_string());
        } else if !t.is_empty() {
            cur_seq.extend(t.bytes().map(encode_base));
        }
    }
    if let Some(id) = cur_id {
        records.push(Record { id, seq: cur_seq });
    }
    Ok(records)
}

fn clstr_path(output: &Path) -> std::path::PathBuf {
    let mut s = output.to_string_lossy().into_owned();
    s.push_str(".clstr");
    std::path::PathBuf::from(s)
}

fn first_word(s: &str) -> &str {
    s.split_whitespace().next().unwrap_or(s)
}
