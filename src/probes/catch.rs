// Optimization-based probe design following CATCH (Metsky et al. 2019).
//
// Implements single-dataset mode with Hamming distance coverage model.
// CATCH reference: Capturing sequence diversity in metagenomes with comprehensive
// and scalable probe design. Nature Biotechnology 37, 160–168 (2019).
//
// Strand handling (see `Strand`): with `Strand::Forward` the coverage model matches
// CATCH, which scores a probe only against the target as given. `Strand::Both` also
// credits reverse-complement matches, which is the more physical model for capture of
// double-stranded DNA but produces smaller panels that are not comparable to CATCH's.
//
// Architecture: coverage model → candidate generation → MinHash dedup →
//               coverage computation → greedy set cover → output.
// Each stage is independent; see the Planned Extensions section in the plan
// for how LCF-k and blacklisting slot in without changing the other stages.

use anyhow::{Context, Result};
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use crate::seq_utils::{hamming_n_mismatch, reverse_complement};

// ─── Coverage model ───────────────────────────────────────────────────────────

/// How to decide whether a probe hybridizes to a target window.
///
/// Returns `Some((cov_start, cov_end))` (target coordinates) if the probe
/// hybridizes to `target[window_start..window_start+probe_len]`, or `None`
/// otherwise.
///
/// To add LCF-k support: add `LcfK { mismatches, lcf_threshold, extension }`
/// and implement `covers()` for it — the set-cover code is unchanged.
pub enum CoverageModel {
    /// Full probe must align with ≤ `mismatches` Hamming distance.
    /// Equivalent to LCF-k with `lcf_threshold = probe_len` (the CATCH default).
    Hamming { mismatches: usize, extension: usize },
}

/// Which orientations count as coverage when scoring a probe against a target.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Strand {
    /// Score the probe against the target as given, like CATCH. Panel sizes are
    /// directly comparable to `design.py`.
    Forward,
    /// Also credit a probe whose reverse complement matches the target. Physically
    /// appropriate for double-stranded capture, where a bait equal to revcomp(R)
    /// hybridizes to the + strand of a fragment containing R. Not valid for
    /// single-stranded targets (e.g. RNA-virus baits capture one orientation only),
    /// and yields smaller panels that cannot be compared with CATCH's.
    Both,
}

impl CoverageModel {
    fn covers(&self, probe: &[u8], target: &[u8], window_start: usize) -> Option<(usize, usize)> {
        match self {
            CoverageModel::Hamming {
                mismatches,
                extension,
            } => {
                let probe_len = probe.len();
                let window_end = window_start + probe_len;
                if window_end > target.len() {
                    return None;
                }
                let window = &target[window_start..window_end];
                if hamming_n_mismatch(probe, window) <= *mismatches {
                    let cov_start = window_start.saturating_sub(*extension);
                    let cov_end = (window_end + extension).min(target.len());
                    Some((cov_start, cov_end))
                } else {
                    None
                }
            }
        }
    }
}

/// Build a k-mer hash index: maps each `seed_len`-mer to all (seq_idx, pos)
/// occurrences. Seeds containing N are skipped.
fn build_kmer_index(seqs: &[Vec<u8>], seed_len: usize) -> HashMap<Vec<u8>, Vec<(usize, usize)>> {
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

/// Find all (seq_idx, cov_start, cov_end) triples where `probe` hybridizes to a
/// target's forward strand, using seed-and-extend with the coverage model.
///
/// Only the orientation of `probe` as supplied is checked; the caller decides whether
/// to additionally pass the reverse complement (see `Strand`).
fn find_matches_with_model(
    index: &HashMap<Vec<u8>, Vec<(usize, usize)>>,
    seqs: &[Vec<u8>],
    probe: &[u8],
    seed_len: usize,
    model: &CoverageModel,
) -> Vec<(usize, usize, usize)> {
    let probe_len = probe.len();
    let mut candidates: HashSet<(usize, usize)> = HashSet::new();

    let n_seeds = probe_len.saturating_sub(seed_len) + 1;
    for i in 0..n_seeds {
        let seed = &probe[i..i + seed_len];
        if seed.contains(&b'N') {
            continue;
        }
        if let Some(hits) = index.get(seed) {
            for &(seq_idx, seq_pos) in hits {
                if seq_pos < i {
                    continue;
                }
                let window_start = seq_pos - i;
                if window_start + probe_len <= seqs[seq_idx].len() {
                    candidates.insert((seq_idx, window_start));
                }
            }
        }
    }

    candidates
        .into_iter()
        .filter_map(|(seq_idx, window_start)| {
            model
                .covers(probe, &seqs[seq_idx], window_start)
                .map(|(cs, ce)| (seq_idx, cs, ce))
        })
        .collect()
}

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
        if let Some(stripped) = trimmed.strip_prefix('>') {
            if let Some(id) = current_id.take() {
                ids.push(id);
                seqs.push(current_seq.clone());
                current_seq.clear();
            }
            let id = stripped.split_whitespace().next().unwrap_or("").to_string();
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

// ─── Interval utilities ───────────────────────────────────────────────────────

type IntervalList = Vec<(usize, usize)>;

fn interval_total(intervals: &[(usize, usize)]) -> usize {
    intervals.iter().map(|(s, e)| e - s).sum()
}

/// Count bases in `[start, end)` that fall within uncovered intervals.
fn uncovered_gain(intervals: &[(usize, usize)], start: usize, end: usize) -> usize {
    if start >= end {
        return 0;
    }
    let mut gain = 0;
    for &(istart, iend) in intervals {
        if iend <= start {
            continue;
        }
        if istart >= end {
            break;
        }
        let overlap_s = istart.max(start);
        let overlap_e = iend.min(end);
        if overlap_e > overlap_s {
            gain += overlap_e - overlap_s;
        }
    }
    gain
}

/// Remove `[start, end)` from the sorted, non-overlapping interval list.
fn subtract_interval(intervals: &mut IntervalList, start: usize, end: usize) {
    if start >= end {
        return;
    }
    let mut result = Vec::with_capacity(intervals.len() + 1);
    for &(istart, iend) in intervals.iter() {
        if iend <= start || istart >= end {
            result.push((istart, iend));
        } else {
            if istart < start {
                result.push((istart, start));
            }
            if iend > end {
                result.push((end, iend));
            }
        }
    }
    *intervals = result;
}

// ─── MinHash deduplication ────────────────────────────────────────────────────

// 100 hash functions gives ~1% collision probability for typical probe panels.
// k=8 balances specificity vs speed for probe-length (~120 bp) sequences.
// These values are empirical choices; the original CATCH paper (Metsky et al. 2019)
// uses a similar MinHash scheme but does not prescribe specific values.
const MINHASH_NUM_FUNCS: usize = 100;
const MINHASH_KMER_LEN: usize = 8;
// Mersenne prime 2^31-1; fits in u64 arithmetic without overflow during hashing.
const MINHASH_PRIME: u64 = 2_147_483_647;

/// Generate `MINHASH_NUM_FUNCS` pairs (a, b) for polynomial hash functions
/// h(x) = (a*x + b) % MINHASH_PRIME. Uses a fixed seed for determinism.
fn generate_hash_params(seed: u64) -> Vec<(u64, u64)> {
    let mut rng = StdRng::seed_from_u64(seed);
    (0..MINHASH_NUM_FUNCS)
        .map(|_| {
            let a: u64 = rng.gen_range(1..MINHASH_PRIME);
            let b: u64 = rng.gen_range(0..MINHASH_PRIME);
            (a, b)
        })
        .collect()
}

/// Convert a DNA 8-mer to a u64. Returns None if the k-mer contains N or
/// any non-ACGT base.
fn kmer_to_u64(kmer: &[u8]) -> Option<u64> {
    let mut val = 0u64;
    for &b in kmer {
        let base = match b {
            b'A' => 0u64,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => return None,
        };
        val = val.wrapping_mul(4).wrapping_add(base);
    }
    Some(val)
}

fn minhash_signature(seq: &[u8], params: &[(u64, u64)]) -> Vec<u64> {
    let mut sig = vec![u64::MAX; params.len()];
    if seq.len() < MINHASH_KMER_LEN {
        return sig;
    }
    for i in 0..=(seq.len() - MINHASH_KMER_LEN) {
        let kmer = &seq[i..i + MINHASH_KMER_LEN];
        let Some(kval) = kmer_to_u64(kmer) else {
            continue;
        };
        for (j, &(a, b)) in params.iter().enumerate() {
            let h = a.wrapping_mul(kval).wrapping_add(b) % MINHASH_PRIME;
            if h < sig[j] {
                sig[j] = h;
            }
        }
    }
    sig
}

fn jaccard_estimate(sig1: &[u64], sig2: &[u64]) -> f64 {
    let matches = sig1
        .iter()
        .zip(sig2.iter())
        .filter(|(&a, &b)| a == b)
        .count();
    matches as f64 / sig1.len() as f64
}

/// Remove near-duplicate probes using MinHash. Processes candidates in order;
/// the first probe in each near-duplicate group is kept.
fn minhash_dedup(
    candidates: Vec<(Vec<u8>, String)>,
    threshold: f64,
    params: &[(u64, u64)],
) -> Vec<(Vec<u8>, String)> {
    if threshold <= 0.0 {
        return candidates;
    }
    let mut accepted: Vec<(Vec<u8>, String)> = Vec::new();
    let mut accepted_sigs: Vec<Vec<u64>> = Vec::new();

    for (probe, source_id) in candidates {
        let sig = minhash_signature(&probe, params);
        let is_near_dup = accepted_sigs
            .iter()
            .any(|asig| jaccard_estimate(&sig, asig) >= threshold);
        if !is_near_dup {
            accepted_sigs.push(sig);
            accepted.push((probe, source_id));
        }
    }

    accepted
}

// ─── Main entry point ─────────────────────────────────────────────────────────

/// Design probes using the CATCH optimization algorithm.
///
/// Generates candidate probes by stride-based tiling, applies MinHash
/// near-deduplication, then uses a greedy set-cover algorithm to find a
/// minimum-cardinality probe set that covers each target sequence to at
/// least `coverage_frac` of its length.
///
/// Returns the number of probes written.
pub fn design_probes(
    input: &Path,
    output: &Path,
    probe_len: usize,
    stride: usize,
    mismatches: usize,
    extension: usize,
    coverage_frac: f64,
    minhash_threshold: f64,
    strand: Strand,
) -> Result<usize> {
    const SEED_LEN: usize = 10;
    const MINHASH_SEED: u64 = 42;
    // Skip candidates with > 50% N bases
    const MAX_N_FRAC: f64 = 0.5;

    log::info!("CATCH: loading sequences from {}", input.display());
    let (ids, seqs) = load_fasta(input)?;

    if seqs.is_empty() {
        log::warn!("CATCH: no sequences in input, writing empty output");
        File::create(output).with_context(|| format!("Cannot create {}", output.display()))?;
        return Ok(0);
    }

    log::info!(
        "CATCH: {} sequences loaded, {} bp total",
        seqs.len(),
        seqs.iter().map(|s| s.len()).sum::<usize>()
    );

    let model = CoverageModel::Hamming {
        mismatches,
        extension,
    };

    // ── Step 1: Generate candidate probes by tiling ──────────────────────────
    log::info!(
        "CATCH: generating candidates (probe_len={}, stride={})...",
        probe_len,
        stride
    );
    // (probe_sequence, source_seq_id)
    let mut candidates: Vec<(Vec<u8>, String)> = Vec::new();
    let mut seen: HashSet<Vec<u8>> = HashSet::new();
    let max_n = (probe_len as f64 * MAX_N_FRAC) as usize;

    for (seq_idx, seq) in seqs.iter().enumerate() {
        if seq.len() < probe_len {
            continue;
        }
        let source_id = ids[seq_idx].clone();
        let mut pos = 0usize;

        loop {
            let start = if pos + probe_len <= seq.len() {
                pos
            } else {
                seq.len() - probe_len // right-edge probe
            };

            let probe = seq[start..start + probe_len].to_vec();
            let n_count = probe.iter().filter(|&&b| b == b'N').count();
            if n_count <= max_n && seen.insert(probe.clone()) {
                candidates.push((probe, source_id.clone()));
            }

            if pos + probe_len >= seq.len() {
                break;
            }
            pos += stride;
        }
    }
    log::info!(
        "CATCH: {} unique candidate probes generated",
        candidates.len()
    );

    // ── Step 2: [Extension point] Blacklist filter (future) ─────────────────
    // When --catch-blacklist is added: build k-mer index of blacklist sequences,
    // filter candidates here, no other steps change.

    // ── Step 3: MinHash near-deduplication ──────────────────────────────────
    if minhash_threshold > 0.0 {
        log::info!(
            "CATCH: MinHash deduplication (threshold={:.2}, {} candidates)...",
            minhash_threshold,
            candidates.len()
        );
        let params = generate_hash_params(MINHASH_SEED);
        let before = candidates.len();
        candidates = minhash_dedup(candidates, minhash_threshold, &params);
        log::info!(
            "CATCH: {} candidates after dedup (removed {})",
            candidates.len(),
            before - candidates.len()
        );
    }

    if candidates.is_empty() {
        log::warn!("CATCH: no candidates remain after filtering, writing empty output");
        File::create(output).with_context(|| format!("Cannot create {}", output.display()))?;
        return Ok(0);
    }

    // ── Step 4: Build k-mer index of target sequences ────────────────────────
    log::info!(
        "CATCH: building {}-mer index ({} sequences)...",
        SEED_LEN,
        seqs.len()
    );
    let index = build_kmer_index(&seqs, SEED_LEN);

    // ── Step 5: Compute coverage for each candidate ──────────────────────────
    log::info!(
        "CATCH: computing coverage for {} candidates...",
        candidates.len()
    );
    // probe_coverage[i] = Vec<(seq_idx, cov_start, cov_end)>
    let mut probe_coverage: Vec<Vec<(usize, usize, usize)>> = Vec::with_capacity(candidates.len());

    for (probe, _) in &candidates {
        let mut hits = find_matches_with_model(&index, &seqs, probe, SEED_LEN, &model);
        if strand == Strand::Both {
            let rc = reverse_complement(probe);
            hits.extend(find_matches_with_model(
                &index, &seqs, &rc, SEED_LEN, &model,
            ));
        }
        hits.sort_unstable();
        hits.dedup();
        probe_coverage.push(hits);
    }

    // ── Step 6: Greedy set cover ─────────────────────────────────────────────
    log::info!(
        "CATCH: greedy set cover (coverage_frac={:.2}, strand={}, {} candidates, {} targets)...",
        coverage_frac,
        match strand {
            Strand::Forward => "forward",
            Strand::Both => "both",
        },
        candidates.len(),
        seqs.len()
    );

    let seq_lens: Vec<usize> = seqs.iter().map(|s| s.len()).collect();
    let mut uncovered: Vec<IntervalList> = seq_lens.iter().map(|&len| vec![(0, len)]).collect();
    let mut uncovered_bp: Vec<usize> = seq_lens.clone();

    let mut remaining: Vec<usize> = (0..candidates.len()).collect();
    let mut solution: Vec<usize> = Vec::new();
    let mut iter = 0usize;

    loop {
        // Check termination
        let all_covered = (0..seqs.len()).all(|si| {
            if seq_lens[si] == 0 {
                return true;
            }
            let covered = seq_lens[si] - uncovered_bp[si];
            covered as f64 / seq_lens[si] as f64 >= coverage_frac
        });
        if all_covered || remaining.is_empty() {
            break;
        }

        // Find probe with highest gain
        let mut best_pi = usize::MAX;
        let mut best_gain = 0usize;

        for &pi in &remaining {
            let gain: usize = probe_coverage[pi]
                .iter()
                .map(|&(si, cs, ce)| uncovered_gain(&uncovered[si], cs, ce))
                .sum();
            if gain > best_gain {
                best_gain = gain;
                best_pi = pi;
            }
        }

        if best_gain == 0 || best_pi == usize::MAX {
            break; // no more coverage possible
        }

        solution.push(best_pi);
        remaining.retain(|&pi| pi != best_pi);

        for &(si, cs, ce) in &probe_coverage[best_pi] {
            subtract_interval(&mut uncovered[si], cs, ce);
            uncovered_bp[si] = interval_total(&uncovered[si]);
        }

        iter += 1;
        if iter.is_multiple_of(100) {
            let total_covered: usize = seq_lens
                .iter()
                .zip(uncovered_bp.iter())
                .map(|(&total, &unc)| total.saturating_sub(unc))
                .sum();
            let total_bp: usize = seq_lens.iter().sum();
            log::debug!(
                "CATCH: iteration {}, {} probes selected, {}/{} bp covered ({:.1}%)",
                iter,
                solution.len(),
                total_covered,
                total_bp,
                100.0 * total_covered as f64 / total_bp.max(1) as f64
            );
        }
    }

    log::info!("CATCH: set cover selected {} probes", solution.len());

    // Report realized coverage rather than assuming the request was met. The greedy
    // loop also exits when no candidate adds further gain, so `coverage_frac` is a
    // target, not a guarantee.
    let total_bp: usize = seq_lens.iter().sum();
    let covered_bp: usize = seq_lens
        .iter()
        .zip(uncovered_bp.iter())
        .map(|(&total, &unc)| total.saturating_sub(unc))
        .sum();
    let short_targets = (0..seqs.len())
        .filter(|&si| {
            seq_lens[si] > 0
                && ((seq_lens[si] - uncovered_bp[si]) as f64 / seq_lens[si] as f64) < coverage_frac
        })
        .count();
    let realized = if total_bp > 0 {
        covered_bp as f64 / total_bp as f64
    } else {
        1.0
    };
    log::info!(
        "CATCH: realized coverage {:.4} ({}/{} bp); {} of {} targets below the requested {:.2}",
        realized,
        covered_bp,
        total_bp,
        short_targets,
        seqs.len(),
        coverage_frac
    );
    if short_targets > 0 {
        log::warn!(
            "CATCH: {} target(s) did not reach the requested coverage of {:.2} \
             (no remaining candidate probe could cover them). Realized coverage {:.4}.",
            short_targets,
            coverage_frac,
            realized
        );
    }

    // ── Step 7: Write output ──────────────────────────────────────────────────
    let out_file =
        File::create(output).with_context(|| format!("Cannot create {}", output.display()))?;
    let mut writer = BufWriter::new(out_file);
    let mut counter_per_source: HashMap<String, usize> = HashMap::new();

    for &pi in &solution {
        let (probe, source_id) = &candidates[pi];
        let n = counter_per_source.entry(source_id.clone()).or_insert(0);
        *n += 1;
        writeln!(writer, ">probe_{}|catch_{}", source_id, n)?;
        writeln!(
            writer,
            "{}",
            std::str::from_utf8(probe).unwrap_or("INVALID")
        )?;
    }
    writer.flush()?;

    log::info!(
        "CATCH: wrote {} probes to {}",
        solution.len(),
        output.display()
    );
    Ok(solution.len())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    /// Deterministic pseudo-random DNA so tests do not depend on an RNG crate.
    fn make_seq(len: usize, seed: u64) -> String {
        let mut s = String::with_capacity(len);
        let mut x = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
        for _ in 0..len {
            x = x
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            s.push(match (x >> 33) % 4 {
                0 => 'A',
                1 => 'C',
                2 => 'G',
                _ => 'T',
            });
        }
        s
    }

    fn rc_str(s: &str) -> String {
        s.chars()
            .rev()
            .map(|c| match c {
                'A' => 'T',
                'T' => 'A',
                'C' => 'G',
                'G' => 'C',
                other => other,
            })
            .collect()
    }

    fn write_fasta(
        dir: &std::path::Path,
        name: &str,
        recs: &[(String, String)],
    ) -> std::path::PathBuf {
        let p = dir.join(name);
        let mut f = File::create(&p).unwrap();
        for (id, seq) in recs {
            writeln!(f, ">{}", id).unwrap();
            writeln!(f, "{}", seq).unwrap();
        }
        p
    }

    fn count_probes(p: &std::path::Path) -> usize {
        std::fs::read_to_string(p)
            .unwrap()
            .lines()
            .filter(|l| l.starts_with('>'))
            .count()
    }

    /// A target present only as the reverse complement of another target must be left
    /// uncovered under Strand::Forward (CATCH-compatible) and so must receive its own
    /// probes; under Strand::Both it is credited to the sibling and gets none.
    #[test]
    fn strand_forward_does_not_credit_reverse_complement() {
        let dir = tempfile::tempdir().unwrap();
        let a = make_seq(480, 11);
        let b = rc_str(&a); // identical gene, deposited on the opposite strand
        let input = write_fasta(
            dir.path(),
            "in.fa",
            &[("fwd".to_string(), a), ("rev".to_string(), b)],
        );

        let out_fwd = dir.path().join("fwd.fa");
        let n_fwd =
            design_probes(&input, &out_fwd, 120, 120, 5, 0, 1.0, 0.0, Strand::Forward).unwrap();

        let out_both = dir.path().join("both.fa");
        let n_both =
            design_probes(&input, &out_both, 120, 120, 5, 0, 1.0, 0.0, Strand::Both).unwrap();

        assert_eq!(n_fwd, count_probes(&out_fwd));
        assert_eq!(n_both, count_probes(&out_both));
        // Forward must tile both copies; Both can cover one with the other's probes.
        assert!(
            n_fwd > n_both,
            "forward-only should need more probes than strand-agnostic ({} vs {})",
            n_fwd,
            n_both
        );
        assert_eq!(
            n_fwd, 8,
            "two 480 bp targets at stride 120 need 4 probes each"
        );
        assert_eq!(
            n_both, 4,
            "strand-agnostic should cover both copies with 4 probes"
        );
    }

    /// Forward-only design must fully cover every target it is given.
    #[test]
    fn strand_forward_covers_every_target() {
        let dir = tempfile::tempdir().unwrap();
        let recs: Vec<(String, String)> = (0..5)
            .map(|i| (format!("t{}", i), make_seq(360, 100 + i)))
            .collect();
        let input = write_fasta(dir.path(), "in.fa", &recs);
        let out = dir.path().join("p.fa");
        let n = design_probes(&input, &out, 120, 120, 5, 0, 1.0, 0.0, Strand::Forward).unwrap();
        // 5 unrelated targets x 360 bp / 120 bp stride, no extension => 3 probes each.
        assert_eq!(n, 15);
    }

    /// The coverage model must accept an exact match and reject one over the limit.
    #[test]
    fn coverage_model_respects_mismatch_budget() {
        let model = CoverageModel::Hamming {
            mismatches: 5,
            extension: 0,
        };
        let target: Vec<u8> = make_seq(120, 7).into_bytes();
        assert_eq!(model.covers(&target, &target, 0), Some((0, 120)));

        let mut probe = target.clone();
        for p in probe.iter_mut().take(6) {
            *p = if *p == b'A' { b'C' } else { b'A' };
        }
        assert_eq!(
            model.covers(&probe, &target, 0),
            None,
            "6 mismatches must not cover"
        );
    }

    /// Cover extension widens the credited interval and is clamped to the target.
    #[test]
    fn cover_extension_is_clamped_to_target_bounds() {
        let model = CoverageModel::Hamming {
            mismatches: 0,
            extension: 50,
        };
        let target: Vec<u8> = make_seq(200, 3).into_bytes();
        let probe = target[40..160].to_vec();
        assert_eq!(model.covers(&probe, &target, 40), Some((0, 200)));
    }

    /// An N in the target counts as a mismatch, so a window of Ns is uncoverable.
    #[test]
    fn n_bases_count_as_mismatches() {
        let model = CoverageModel::Hamming {
            mismatches: 5,
            extension: 0,
        };
        let mut target: Vec<u8> = make_seq(120, 5).into_bytes();
        let probe = target.clone();
        for b in target.iter_mut().take(6) {
            *b = b'N';
        }
        assert_eq!(model.covers(&probe, &target, 0), None);
    }
}
