// Syotti-style greedy bait design algorithm.
// Based on: Alanko et al. (2022) "Syotti: scalable bait design for DNA enrichment"
// Bioinformatics 38, i177-i184. https://doi.org/10.1093/bioinformatics/btac226
//
// Original implementation uses an FM-index (SDSL + Divsufsort). Here we replace
// it with a k-mer hash index, which is suitable for the MB-scale inputs typical
// in BaitBench and requires no additional dependencies.

use anyhow::{Context, Result};
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use crate::seq_utils::{hamming_n_mismatch, reverse_complement};

/// Design probes using the Syotti greedy algorithm.
///
/// Scans each input sequence; at every uncovered position, extracts a bait of
/// length `bait_len` and marks all reference windows within `max_mismatches`
/// Hamming distance as covered (searching both forward and reverse complement).
/// N characters never match anything, even other Ns.
///
/// Returns the number of probes written.
pub fn design_probes(
    input: &Path,
    output: &Path,
    bait_len: usize,
    max_mismatches: usize,
    seed_len: usize,
) -> Result<usize> {
    log::info!("Syotti: loading sequences from {}", input.display());
    let (ids, seqs) = load_fasta(input)?;

    if seqs.is_empty() {
        log::warn!("Syotti: no sequences in input, writing empty output");
        File::create(output).with_context(|| format!("Cannot create {}", output.display()))?;
        return Ok(0);
    }

    log::info!(
        "Syotti: {} sequences loaded, total {} bp",
        seqs.len(),
        seqs.iter().map(|s| s.len()).sum::<usize>()
    );

    log::info!("Syotti: building {}-mer index...", seed_len);
    let index = build_kmer_index(&seqs, seed_len);
    log::info!(
        "Syotti: index contains {} distinct {}-mers",
        index.len(),
        seed_len
    );

    let mut covered: Vec<Vec<bool>> = seqs.iter().map(|s| vec![false; s.len()]).collect();

    let out_file =
        File::create(output).with_context(|| format!("Cannot create {}", output.display()))?;
    let mut writer = BufWriter::new(out_file);

    let mut total_probes = 0usize;

    log::info!(
        "Syotti: greedy probe selection (bait_len={}, max_mismatches={}, seed_len={})...",
        bait_len,
        max_mismatches,
        seed_len
    );

    for seq_idx in 0..seqs.len() {
        let seq = &seqs[seq_idx];

        if seq.len() < bait_len {
            log::debug!(
                "Syotti: skipping '{}' (len {} < bait_len {})",
                ids[seq_idx],
                seq.len(),
                bait_len
            );
            continue;
        }

        let mut pos = 0usize;
        let mut probe_num = 0usize;

        while pos < seq.len() {
            if covered[seq_idx][pos] {
                pos += 1;
                continue;
            }

            // Back up if near the end so the bait stays within bounds.
            let start = pos.min(seq.len() - bait_len);
            let bait: Vec<u8> = seq[start..start + bait_len].to_vec();

            // Cover all windows matching the bait (forward) and its reverse complement.
            let fwd_matches =
                find_matches(&index, &seqs, &bait, bait_len, max_mismatches, seed_len);
            let rc_bait = reverse_complement(&bait);
            let rc_matches =
                find_matches(&index, &seqs, &rc_bait, bait_len, max_mismatches, seed_len);

            for &(m_seq, m_pos) in fwd_matches.iter().chain(rc_matches.iter()) {
                let end = (m_pos + bait_len).min(covered[m_seq].len());
                for off in m_pos..end {
                    covered[m_seq][off] = true;
                }
            }

            probe_num += 1;
            total_probes += 1;
            writeln!(writer, ">probe_{}|syotti_{}", ids[seq_idx], probe_num)?;
            writeln!(writer, "{}", std::str::from_utf8(&bait).unwrap_or("INVALID"))?;

            pos = start + bait_len;
        }
    }

    writer.flush()?;
    log::info!("Syotti: designed {} probes total", total_probes);
    Ok(total_probes)
}


/// Build a k-mer hash index: maps each `seed_len`-mer to all (seq_idx, pos)
/// occurrences across the input sequences. Seeds containing N are skipped.
fn build_kmer_index(
    seqs: &[Vec<u8>],
    seed_len: usize,
) -> HashMap<Vec<u8>, Vec<(usize, usize)>> {
    let mut index: HashMap<Vec<u8>, Vec<(usize, usize)>> = HashMap::new();
    for (seq_idx, seq) in seqs.iter().enumerate() {
        if seq.len() < seed_len {
            continue;
        }
        for pos in 0..=(seq.len() - seed_len) {
            let kmer = &seq[pos..pos + seed_len];
            if kmer.contains(&b'N') {
                continue;
            }
            index.entry(kmer.to_vec()).or_default().push((seq_idx, pos));
        }
    }
    index
}

/// Find all reference windows of length `bait_len` that are within `max_mismatches`
/// Hamming distance of `bait`, using seed-and-extend with the k-mer index.
///
/// For each seed position `i` in the bait, the `seed_len`-mer at that offset is
/// looked up in the index; each hit (seq_idx, seq_pos) gives a candidate window
/// starting at `seq_pos - i`. Candidates are then verified with full Hamming distance.
///
/// Using all overlapping seed positions guarantees no false negatives when the
/// number of mismatches is at most `bait_len - seed_len`.
fn find_matches(
    index: &HashMap<Vec<u8>, Vec<(usize, usize)>>,
    seqs: &[Vec<u8>],
    bait: &[u8],
    bait_len: usize,
    max_mismatches: usize,
    seed_len: usize,
) -> Vec<(usize, usize)> {
    let mut candidates: HashSet<(usize, usize)> = HashSet::new();

    let n_seeds = bait_len.saturating_sub(seed_len) + 1;
    for i in 0..n_seeds {
        let seed = &bait[i..i + seed_len];
        if seed.contains(&b'N') {
            continue;
        }
        if let Some(hits) = index.get(seed) {
            for &(seq_idx, seq_pos) in hits {
                // seed is at offset i in the bait and at seq_pos in the sequence,
                // so the candidate window starts at seq_pos - i.
                if seq_pos < i {
                    continue;
                }
                let window_start = seq_pos - i;
                if window_start + bait_len > seqs[seq_idx].len() {
                    continue;
                }
                candidates.insert((seq_idx, window_start));
            }
        }
    }

    candidates
        .into_iter()
        .filter(|&(seq_idx, window_start)| {
            let window = &seqs[seq_idx][window_start..window_start + bait_len];
            hamming_n_mismatch(bait, window) <= max_mismatches
        })
        .collect()
}

/// Load a FASTA file into (ids, sequences). Sequences are stored as uppercase bytes.
fn load_fasta(path: &Path) -> Result<(Vec<String>, Vec<Vec<u8>>)> {
    let file = File::open(path).with_context(|| format!("Cannot open {}", path.display()))?;
    let reader = BufReader::new(file);

    let mut ids: Vec<String> = Vec::new();
    let mut seqs: Vec<Vec<u8>> = Vec::new();
    let mut current_id: Option<String> = None;
    let mut current_seq: Vec<u8> = Vec::new();

    for line in reader.lines() {
        let line = line?;
        let trimmed = line.trim_end();
        if trimmed.starts_with('>') {
            if let Some(id) = current_id.take() {
                ids.push(id);
                seqs.push(current_seq.clone());
                current_seq.clear();
            }
            let id = trimmed[1..]
                .split_whitespace()
                .next()
                .unwrap_or("")
                .to_string();
            current_id = Some(id);
        } else if !trimmed.is_empty() {
            current_seq.extend(trimmed.as_bytes().iter().map(|b| b.to_ascii_uppercase()));
        }
    }
    if let Some(id) = current_id {
        ids.push(id);
        seqs.push(current_seq);
    }

    Ok((ids, seqs))
}
