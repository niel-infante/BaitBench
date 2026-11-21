//! Ungapped local seed-and-extend alignment of probes against a reference.
//!
//! # Why this exists
//!
//! Probe capture is a *local* phenomenon: 60 bp of perfect duplex out of a 120 bp probe
//! hybridizes and captures the fragment perfectly well. Read mappers are built for the
//! opposite problem — placing a read end-to-end — and discard partial matches. On a
//! panel designed against relatives of the sampled organism, measured against
//! *B. burgdorferi* B31, the embedded mapper placed 8 of 228 probes (21 after
//! retuning) where a local aligner placed 74. Probe density drives per-locus capture
//! efficiency, so those missing probes distort the simulation directly.
//!
//! This is **not** a BLAST reimplementation, and deliberately does not claim to be: no
//! gapped extension, no two-hit heuristic, no Karlin-Altschul statistics, no scoring
//! matrices, no masking. It is seed-and-extend with ungapped local extension and
//! identity/coverage thresholds — enough to decide whether a probe forms a duplex long
//! enough to capture, and nothing more.
//!
//! # Direction of indexing
//!
//! The probes are indexed and the reference is streamed past them, not the other way
//! round. Indexing the reference costs one entry per base: ~1.5 GB for a 93 Mb mock
//! community and ~36 GB for a 2.25 GB host genome. Probe panels are small (6,535 x
//! 120 bp is under a megabyte), so inverting the index makes the memory cost a
//! property of the panel rather than of the reference.

use std::collections::HashSet;

use crate::seq_utils::reverse_complement;

/// Tuning for the seed-and-extend search.
#[derive(Debug, Clone, Copy)]
pub struct SeedExtendParams {
    /// Exact-match seed length. Shorter finds more divergent probes and costs more.
    pub seed_len: usize,
    /// Minimum fraction of the probe that must align.
    pub min_coverage: f64,
    /// Minimum identity over the aligned segment.
    pub min_identity: f64,
    /// Mismatch penalty used when choosing the best local segment (match scores +1).
    pub mismatch_penalty: f64,
    /// Cap on reported hits per probe, best-scoring first.
    pub max_hits_per_probe: usize,
    /// Bases of reference scanned per work unit.
    ///
    /// Candidate diagonals are accumulated per chunk before extension, so this bounds
    /// peak memory: roughly `chunk x entries_per_seed x 32 B x threads`. It also sets
    /// the unit of parallelism, which is what lets a single 250 Mbp scaffold use more
    /// than one core. It cannot change the result — the scan is lossless whatever the
    /// chunking — so it is a performance knob, not an experimental parameter, and is
    /// deliberately not exposed on the CLI.
    pub chunk_bases: usize,
    /// Ignore seeds occurring in more than this many probe positions.
    ///
    /// Low-complexity seeds dominate the cost of the reference scan while carrying
    /// almost no information: each one fans out to thousands of candidate diagonals
    /// that are then rejected by extension. Real aligners cap this too (minimap2's
    /// `mid_occ`). Probes that only seed through such k-mers are still recovered via
    /// their other offsets, since every offset is indexed.
    pub max_seed_occ: usize,
}

impl Default for SeedExtendParams {
    fn default() -> Self {
        // Thresholds mirror the BLAST settings used elsewhere for probe placement
        // (-perc_identity 80 -qcov_hsp_perc 50), so the two paths accept the same
        // duplexes. seed_len 10 is short enough to seed a 120-mer at ~80% identity:
        // by pigeonhole 24 mismatches in 120 bp still leave an exact run of ~5 bp on
        // average, and seeds are tried at every offset.
        SeedExtendParams {
            seed_len: 10,
            min_coverage: 0.5,
            min_identity: 0.8,
            mismatch_penalty: 2.0,
            max_hits_per_probe: 1000,
            // 4 Mbp: small enough that eight concurrent chunks stay well inside a
            // laptop's memory on a repetitive genome, large enough that the
            // seed_len-1 overlap between chunks is negligible duplicated work.
            chunk_bases: 4 << 20,
            max_seed_occ: 200,
        }
    }
}

/// One ungapped local alignment of a probe to a reference sequence.
#[derive(Debug, Clone)]
pub struct LocalHit {
    /// Index into the reference sequence list.
    pub seq_idx: usize,
    /// 0-based, inclusive start of the aligned segment in the reference.
    pub ref_start: usize,
    /// 0-based, exclusive end of the aligned segment in the reference.
    pub ref_end: usize,
    /// True when the probe matched as its reverse complement. Not used for scoring —
    /// stacking energy is symmetric under reverse-complementation — but retained
    /// because it is needed to report orientation and to debug placements.
    #[allow(dead_code)]
    pub reverse: bool,
    /// (probe base, reference base) over the aligned segment, in reference orientation.
    pub pairs: Vec<(u8, u8)>,
    /// Identity over the aligned segment.
    pub identity: f64,
}

/// A k-mer index over a set of short queries (probes).
///
/// Seeds are 2-bit encoded and stored CSR-style (bucket offsets plus one flat entry
/// array) rather than in a hash map. At the default seed length there are only
/// 4^10 ≈ 1M possible seeds, so the lookup becomes an array index and the reference
/// scan needs no hashing at all — the scan is the hot loop, running once per
/// reference base.
pub struct ProbeIndex {
    seed_len: usize,
    /// Start offset into `entries` for each encoded seed; length 4^seed_len + 1.
    offsets: Vec<u32>,
    /// (oriented probe index, offset within that probe), grouped by seed.
    entries: Vec<(u32, u32)>,
    /// Oriented probe sequences: forward at `2 * i`, reverse complement at `2 * i + 1`.
    oriented: Vec<Vec<u8>>,
}

/// 2-bit code for a base, or `None` for anything ambiguous.
#[inline]
fn base_code(b: u8) -> Option<u64> {
    match b.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

/// Largest seed length that can be directly indexed without an unreasonable table.
/// 4^12 offsets is 67 MB; beyond that the flat table stops paying for itself.
pub const MAX_SEED_LEN: usize = 12;

impl ProbeIndex {
    /// Index both orientations of every probe.
    pub fn build(probes: &[Vec<u8>], seed_len: usize) -> Self {
        let seed_len = seed_len.clamp(4, MAX_SEED_LEN);
        let n_seeds = 1usize << (2 * seed_len);

        let mut oriented: Vec<Vec<u8>> = Vec::with_capacity(probes.len() * 2);
        for p in probes {
            let up: Vec<u8> = p.iter().map(|b| b.to_ascii_uppercase()).collect();
            let rc = reverse_complement(&up);
            oriented.push(up);
            oriented.push(rc);
        }

        // Two passes: count per seed, prefix-sum to offsets, then fill.
        let mut counts = vec![0u32; n_seeds + 1];
        let visit = |seq: &[u8], mut f: Box<dyn FnMut(usize, usize) + '_>| {
            let mut code: u64 = 0;
            let mut valid = 0usize;
            let mask = (1u64 << (2 * seed_len)) - 1;
            for (i, &b) in seq.iter().enumerate() {
                match base_code(b) {
                    Some(c) => {
                        code = ((code << 2) | c) & mask;
                        valid += 1;
                    }
                    None => {
                        code = 0;
                        valid = 0;
                    }
                }
                if valid >= seed_len {
                    f(code as usize, i + 1 - seed_len);
                }
            }
        };

        for seq in oriented.iter() {
            visit(seq, Box::new(|code, _off| counts[code] += 1));
        }
        let mut offsets = vec![0u32; n_seeds + 1];
        let mut acc = 0u32;
        for i in 0..n_seeds {
            offsets[i] = acc;
            acc = acc.saturating_add(counts[i]);
        }
        offsets[n_seeds] = acc;

        let mut cursor = offsets.clone();
        let mut entries = vec![(0u32, 0u32); acc as usize];
        for (oi, seq) in oriented.iter().enumerate() {
            visit(
                seq,
                Box::new(|code, off| {
                    let slot = cursor[code] as usize;
                    if slot < entries.len() {
                        entries[slot] = (oi as u32, off as u32);
                        cursor[code] += 1;
                    }
                }),
            );
        }

        ProbeIndex {
            seed_len,
            offsets,
            entries,
            oriented,
        }
    }

    #[inline]
    fn lookup(&self, code: usize) -> &[(u32, u32)] {
        let s = self.offsets[code] as usize;
        let e = self.offsets[code + 1] as usize;
        &self.entries[s..e]
    }

    fn probe_of(&self, oriented_idx: usize) -> (usize, bool) {
        (oriented_idx / 2, oriented_idx % 2 == 1)
    }
}

/// Best-scoring contiguous segment under (+1 match, -`penalty` mismatch).
///
/// Kadane's algorithm over the match/mismatch profile: this is the ungapped local
/// extension, and it is what lets a probe be placed on the strength of part of its
/// length rather than all of it.
/// Count positions where two slices hold the same byte, eight at a time.
///
/// Used as a cheap necessary condition before the full extension. It deliberately
/// *over*-counts relative to the scoring rule — which requires the base to be one of
/// ACGT, so `N` against `N` scores as a mismatch but is counted equal here. Over-counting
/// only ever makes the filter more permissive, so it can admit a candidate the full
/// extension then rejects, but never reject one the full extension would have accepted.
fn count_equal_bytes(a: &[u8], b: &[u8]) -> usize {
    let n = a.len().min(b.len());
    let (mut i, mut eq) = (0usize, 0usize);
    while i + 8 <= n {
        let x = u64::from_ne_bytes(a[i..i + 8].try_into().unwrap())
            ^ u64::from_ne_bytes(b[i..i + 8].try_into().unwrap());
        // Classic has-zero-byte test: each equal byte leaves a set high bit.
        let zeros = x.wrapping_sub(0x0101_0101_0101_0101) & !x & 0x8080_8080_8080_8080;
        eq += zeros.count_ones() as usize;
        i += 8;
    }
    while i < n {
        if a[i] == b[i] {
            eq += 1;
        }
        i += 1;
    }
    eq
}

fn best_segment(probe: &[u8], reference: &[u8], penalty: f64) -> Option<(usize, usize, f64)> {
    let n = probe.len().min(reference.len());
    if n == 0 {
        return None;
    }
    let (mut best, mut best_s, mut best_e) = (0.0f64, 0usize, 0usize);
    let (mut cur, mut cur_s) = (0.0f64, 0usize);
    for i in 0..n {
        let a = probe[i].to_ascii_uppercase();
        let b = reference[i].to_ascii_uppercase();
        let score = if a == b && matches!(a, b'A' | b'C' | b'G' | b'T') {
            1.0
        } else {
            -penalty
        };
        if cur <= 0.0 {
            cur = score;
            cur_s = i;
        } else {
            cur += score;
        }
        if cur > best {
            best = cur;
            best_s = cur_s;
            best_e = i + 1;
        }
    }
    if best <= 0.0 {
        None
    } else {
        Some((best_s, best_e, best))
    }
}

/// Find ungapped local alignments of every probe against every reference sequence.
///
/// The reference is streamed; only the probe index is held in memory.
/// `threads` bounds the parallel scan. A *scoped* pool is used rather than rayon's
/// global one: the global pool can only be configured once per process, and silently
/// falling back to it would mean ignoring `--threads` and taking every core on the
/// machine.
/// Record progress through the reference and log at each 10% boundary.
///
/// Sequences finish out of order under `par_iter`, so progress is accumulated as bases
/// consumed and the next threshold is claimed atomically — whichever thread crosses it
/// logs, exactly once.
fn report_progress(
    done: &std::sync::atomic::AtomicU64,
    next: &std::sync::atomic::AtomicU64,
    total: u64,
    just_done: u64,
) {
    use std::sync::atomic::Ordering;
    if total == 0 {
        return;
    }
    let now = done.fetch_add(just_done, Ordering::Relaxed) + just_done;
    let threshold = next.load(Ordering::Relaxed);
    if now >= threshold && threshold > 0 {
        let step = total / 10;
        let bumped = ((now / step.max(1)) + 1) * step.max(1);
        // Only the thread that successfully advances the marker logs this milestone.
        if next
            .compare_exchange(threshold, bumped, Ordering::Relaxed, Ordering::Relaxed)
            .is_ok()
        {
            log::info!(
                "  probe placement: {:.0}% of reference scanned ({:.2} of {:.2} Gbp)",
                100.0 * now as f64 / total as f64,
                now as f64 / 1e9,
                total as f64 / 1e9
            );
        }
    }
}

pub fn find_local_hits(
    index: &ProbeIndex,
    references: &[&[u8]],
    params: &SeedExtendParams,
    threads: usize,
) -> Vec<Vec<LocalHit>> {
    use rayon::prelude::*;

    let n_probes = index.oriented.len() / 2;

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads.max(1))
        .build()
        .unwrap_or_else(|_| rayon::ThreadPoolBuilder::new().build().unwrap());

    // Work is parallelised over **chunks of reference**, not over reference sequences.
    //
    // Per-sequence parallelism leaves a genome like Ixodes — whose largest scaffolds
    // are 253, 218 and 199 Mbp — using a fraction of the machine, because one scaffold
    // is one task: measured ~2.8 of 10 cores. Chunking the scan was already needed to
    // bound memory, and chunks are independent, so they are also the right unit of
    // parallelism. Placements that straddle a boundary are seeded from both sides and
    // deduplicated after the merge.
    //
    // Progress is reported per chunk and counted in bases, since scaffold sizes vary
    // by orders of magnitude; reporting per completed sequence meant a multi-hour scan
    // of one giant contig emitted nothing at all.
    // Funnel counters: how many diagonals the scan proposes versus how many survive
    // extension. The ratio is the headroom for any candidate-reducing optimisation.
    let n_candidates = std::sync::atomic::AtomicU64::new(0);
    let n_extended = std::sync::atomic::AtomicU64::new(0);
    // Lossless lower bound, per oriented probe, on the matches any acceptable window
    // must contain. Computed per probe rather than once: panels are not required to be
    // uniform length, and taking the first probe's length would set too high a bar for
    // every shorter probe and silently drop valid placements.
    let min_window_matches: Vec<usize> = index
        .oriented
        .iter()
        .map(|p| (params.min_identity * params.min_coverage * p.len() as f64).floor() as usize)
        .collect();
    let total_bases: u64 = references.iter().map(|r| r.len() as u64).sum();
    let done_bases = std::sync::atomic::AtomicU64::new(0);
    let next_report = std::sync::atomic::AtomicU64::new(total_bases / 10);

    // (sequence index, first base of the chunk)
    let chunk_bases = params.chunk_bases.max(index.seed_len);
    let work: Vec<(usize, usize)> = references
        .iter()
        .enumerate()
        .flat_map(|(si, r)| {
            let n = r.len();
            (0..n)
                .step_by(chunk_bases)
                .map(move |c| (si, c))
                .collect::<Vec<_>>()
        })
        .collect();

    let per_chunk: Vec<Vec<(usize, LocalHit)>> = pool.install(|| {
        work.par_iter()
            .map(|&(seq_idx, chunk_start)| {
                let reference = references[seq_idx];
                let mut out: Vec<(usize, LocalHit)> = Vec::new();
                if reference.len() < index.seed_len {
                    report_progress(
                        &done_bases,
                        &next_report,
                        total_bases,
                        reference.len() as u64,
                    );
                    return out;
                }

                let seed_hi = (chunk_start + chunk_bases).min(reference.len());
                // Start early enough that a seed beginning just before the chunk
                // boundary still has its full k-mer in view.
                let scan_lo = chunk_start.saturating_sub(index.seed_len - 1);
                let mask = (1u64 << (2 * index.seed_len)) - 1;

                let mut candidates: HashSet<(u32, i64)> = HashSet::new();
                let mut code: u64 = 0;
                let mut valid = 0usize;
                for (i, &b) in reference.iter().enumerate().skip(scan_lo) {
                    match base_code(b) {
                        Some(c) => {
                            code = ((code << 2) | c) & mask;
                            valid += 1;
                        }
                        None => {
                            code = 0;
                            valid = 0;
                            continue;
                        }
                    }
                    if valid < index.seed_len {
                        continue;
                    }
                    let pos = i + 1 - index.seed_len;
                    if pos < chunk_start {
                        continue; // covered by the previous chunk
                    }
                    if pos >= seed_hi {
                        break;
                    }
                    let seed_hits = index.lookup(code as usize);
                    if seed_hits.len() > params.max_seed_occ {
                        continue; // uninformative seed; see max_seed_occ
                    }
                    for &(oi, off) in seed_hits {
                        let start = pos as i64 - off as i64;
                        if start >= 0
                            && (start as usize) + index.oriented[oi as usize].len()
                                <= reference.len()
                        {
                            candidates.insert((oi, start));
                        }
                    }
                }

                n_candidates.fetch_add(
                    candidates.len() as u64,
                    std::sync::atomic::Ordering::Relaxed,
                );
                for (oi, start) in candidates {
                    n_extended.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                    let probe = &index.oriented[oi as usize];
                    let start = start as usize;
                    let window = &reference[start..start + probe.len()];

                    // An accepted placement needs a segment of at least
                    // `min_coverage * probe_len` bases at `min_identity`, so it carries
                    // at least `min_identity * min_coverage * probe_len` matches, and
                    // that segment lies inside this window. A window with fewer matches
                    // cannot produce one, so it is rejected without paying for
                    // extension — this is the whole runtime, since the scan proposes
                    // ~12,000 candidates per surviving placement.
                    if count_equal_bytes(probe, window) < min_window_matches[oi as usize] {
                        continue;
                    }

                    let Some((s, e, _score)) = best_segment(probe, window, params.mismatch_penalty)
                    else {
                        continue;
                    };
                    let seg_len = e - s;
                    if (seg_len as f64) < params.min_coverage * probe.len() as f64 {
                        continue;
                    }
                    let pairs: Vec<(u8, u8)> = probe[s..e]
                        .iter()
                        .zip(window[s..e].iter())
                        .map(|(&p, &r)| (p.to_ascii_uppercase(), r.to_ascii_uppercase()))
                        .collect();
                    let matches = pairs.iter().filter(|(p, r)| p == r).count();
                    let identity = matches as f64 / seg_len as f64;
                    if identity < params.min_identity {
                        continue;
                    }
                    let (probe_idx, reverse) = index.probe_of(oi as usize);
                    out.push((
                        probe_idx,
                        LocalHit {
                            seq_idx,
                            ref_start: start + s,
                            ref_end: start + e,
                            reverse,
                            pairs,
                            identity,
                        },
                    ));
                }

                report_progress(
                    &done_bases,
                    &next_report,
                    total_bases,
                    (seed_hi - chunk_start) as u64,
                );
                out
            })
            .collect()
    });

    // A placement whose seeds straddle a chunk boundary is found by both chunks.
    // Deduplicate once at the merge: the surviving set is small (hundreds to tens of
    // thousands) next to the hundreds of millions of candidates behind it.
    let mut seen: HashSet<(usize, bool, usize, usize, usize)> = HashSet::new();
    let per_seq: Vec<Vec<(usize, LocalHit)>> = vec![per_chunk
        .into_iter()
        .flatten()
        .filter(|(pi, h)| seen.insert((*pi, h.reverse, h.seq_idx, h.ref_start, h.ref_end)))
        .collect()];

    {
        use std::sync::atomic::Ordering;
        let cand = n_candidates.load(Ordering::Relaxed);
        let ext = n_extended.load(Ordering::Relaxed);
        let kept: u64 = per_seq.iter().map(|v| v.len() as u64).sum();
        log::info!(
            "  probe placement funnel: {} candidate diagonals -> {} extended -> {} passed \
             ({:.1}% of extensions discarded)",
            cand,
            ext,
            kept,
            if ext > 0 {
                100.0 * (ext - kept.min(ext)) as f64 / ext as f64
            } else {
                0.0
            }
        );
    }

    let mut per_probe: Vec<Vec<LocalHit>> = vec![Vec::new(); n_probes];
    for chunk in per_seq {
        for (probe_idx, hit) in chunk {
            per_probe[probe_idx].push(hit);
        }
    }

    // Rank by identity then length, and cap, so truncation keeps the best hits rather
    // than an arbitrary subset.
    for hits in per_probe.iter_mut() {
        hits.sort_by(|a, b| {
            b.identity
                .partial_cmp(&a.identity)
                .unwrap_or(std::cmp::Ordering::Equal)
                .then((b.ref_end - b.ref_start).cmp(&(a.ref_end - a.ref_start)))
                .then(a.seq_idx.cmp(&b.seq_idx))
                .then(a.ref_start.cmp(&b.ref_start))
        });
        hits.truncate(params.max_hits_per_probe);
    }
    per_probe
}

#[cfg(test)]
mod tests {
    use super::*;

    fn bytes(s: &str) -> Vec<u8> {
        s.as_bytes().to_vec()
    }

    const REF: &str = "TTTTTTTTTTACGTTGCAAGGCCTTAGCATCGATCGGATCCTAGCATGCAGGGGGGGGGG";

    #[test]
    fn finds_a_perfect_probe() {
        let probe = bytes(&REF[10..50]);
        let idx = ProbeIndex::build(std::slice::from_ref(&probe), 10);
        let hits = find_local_hits(&idx, &[&bytes(REF)], &SeedExtendParams::default(), 1);
        assert_eq!(hits[0].len(), 1);
        let h = &hits[0][0];
        assert_eq!((h.ref_start, h.ref_end), (10, 50));
        assert!(!h.reverse);
        assert!((h.identity - 1.0).abs() < 1e-9);
    }

    #[test]
    fn finds_a_reverse_complement_probe() {
        let probe = reverse_complement(&bytes(&REF[10..50]));
        let idx = ProbeIndex::build(&[probe], 10);
        let hits = find_local_hits(&idx, &[&bytes(REF)], &SeedExtendParams::default(), 1);
        assert_eq!(hits[0].len(), 1);
        assert!(hits[0][0].reverse);
        assert!((hits[0][0].identity - 1.0).abs() < 1e-9);
    }

    /// The point of local extension: a probe matching over only part of its length is
    /// still placed, because that much duplex still captures.
    #[test]
    fn places_a_probe_matching_over_half_its_length() {
        let mut probe = bytes(&REF[10..50]);
        for b in probe.iter_mut().take(14) {
            *b = b'A'; // wreck the first 14 of 40 bases
        }
        let idx = ProbeIndex::build(&[probe], 10);
        let hits = find_local_hits(&idx, &[&bytes(REF)], &SeedExtendParams::default(), 1);
        assert_eq!(hits[0].len(), 1, "partial match should still be placed");
        let h = &hits[0][0];
        assert!(h.ref_end - h.ref_start >= 20, "segment too short: {h:?}");
        assert!(h.ref_start >= 20, "segment should skip the wrecked prefix");
    }

    /// Below the coverage floor a probe is not placed at all.
    #[test]
    fn rejects_a_probe_matching_too_little() {
        let mut probe = bytes(&REF[10..50]);
        for b in probe.iter_mut().take(32) {
            *b = b'A'; // only 8 of 40 bases survive
        }
        let idx = ProbeIndex::build(&[probe], 10);
        let hits = find_local_hits(&idx, &[&bytes(REF)], &SeedExtendParams::default(), 1);
        assert!(hits[0].is_empty());
    }

    #[test]
    fn best_segment_picks_the_strongest_run() {
        let p = bytes("AAAACCCCGGGG");
        let r = bytes("TTTTCCCCGGGG");
        let (s, e, _) = best_segment(&p, &r, 2.0).unwrap();
        assert_eq!((s, e), (4, 12));
    }

    #[test]
    fn seeds_skip_ambiguous_kmers() {
        let probe = bytes("ACGTNACGTNACGTNACGTN");
        let idx = ProbeIndex::build(&[probe], 10);
        assert!(
            idx.entries.is_empty(),
            "N-containing seeds must not be indexed"
        );
    }
}
