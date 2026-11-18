/// Thermodynamic (and simple) probe-biased fragment simulation.
///
/// Implements two-level multinomial sampling:
///   Level 1 — sample a probe uniformly from probes with ≥1 hit
///   Level 2 — sample a hit from that probe, weighted by score
///             (TNN Boltzmann × sequence_weight  or  just sequence_weight in simple mode)
///   Level 3 — fragment center = probe_center ± uniform jitter (±frag_len/4)
///   Level 4 — fragment length from Normal(mean, sd) clamped to [min, max]
///
/// Background fragments are drawn uniformly, weighted by sequence_weight × seq_len.
use anyhow::{bail, Context, Result};
use rand::distributions::WeightedIndex;
use rand::prelude::*;
use rand_distr::Normal;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use crate::alignment::cigar::{cigar_ref_len, expand_cigar, CigarOp};
use crate::sampling::FragmentLengthParams;
use crate::thermodynamics::{boltzmann_score, delta_g, ThermoModel};

/// Simulate mode — controls how probe hit scores are computed.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum SimulateMode {
    /// ΔG Boltzmann × sequence_weight
    Thermodynamic,
    /// Uniform per alignment × sequence_weight (no temperature dependency)
    Simple,
}

/// A single probe alignment hit against the reference.
#[derive(Debug)]
pub struct ProbeHit {
    pub seq_id: String,
    pub start: usize, // 0-based, inclusive
    pub end: usize,   // 0-based, exclusive
    pub score: f64,   // TNN_boltzmann × seq_weight  (or just seq_weight in simple mode)
}

/// Build the (probe base, paired base) list that `delta_g()` consumes, reading the
/// **actual reference sequence** rather than inferring it.
///
/// The second element of each pair is the base on the strand the probe hybridizes to,
/// which is the complement of the reference base: a probe identical to the + strand
/// binds the - strand. For a matching position that is `complement(probe)`, giving a
/// Watson-Crick pair; for a mismatch it is `complement(reference)`, which correctly
/// scores as a mismatch.
///
/// This replaces an earlier MD-tag reconstruction that synthesised `complement(probe)`
/// at every position MD did not flag. That made `is_wc` unconditionally true outside
/// MD mismatches (so a missing MD tag scored every alignment as a perfect duplex), and
/// inserted MD bases without complementing them (so probe A over reference T — a real
/// mismatch — scored as a Watson-Crick pair). Reading the reference removes the
/// dependence on MD entirely, and with it the CIGAR/MD coordinate desynchronisation on
/// insertions and soft clips.
///
/// `query` must be in reference orientation (the SAM writer reverse-complements SEQ for
/// reverse-strand records), and `ref_start` is 0-based.
///
/// Insertions and deletions emit no pair, so stacking runs across a gap as if the
/// flanks were contiguous. That is pre-existing behaviour and a modelling choice rather
/// than a correctness bug; revisit alongside the mismatch-stacking work.
fn aligned_pairs(query: &[u8], reference: &[u8], ref_start: usize, cigar: &str) -> Vec<(u8, u8)> {
    let ops = expand_cigar(cigar);
    let mut pairs: Vec<(u8, u8)> = Vec::with_capacity(query.len());
    let mut qi = 0usize;
    let mut ri = ref_start;

    for op in &ops {
        match op {
            CigarOp::Match | CigarOp::Equal | CigarOp::Mismatch => {
                if qi < query.len() && ri < reference.len() {
                    let probe_b = query[qi].to_ascii_uppercase();
                    let ref_b = reference[ri].to_ascii_uppercase();
                    pairs.push((probe_b, wc_complement_base(ref_b)));
                }
                qi += 1;
                ri += 1;
            }
            CigarOp::Ins | CigarOp::SoftClip => {
                qi += 1;
            }
            CigarOp::Del | CigarOp::Skip => {
                ri += 1;
            }
            CigarOp::HardClip | CigarOp::Pad => {}
        }
    }

    pairs
}

#[inline]
fn wc_complement_base(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        _ => b'N',
    }
}

/// SAM flag bits
const FLAG_UNMAPPED: u16 = 0x4;

/// Place probes with the built-in ungapped local seed-and-extend aligner and score them.
///
/// Runs in-process: no intermediate alignment file, and the probe index (not the
/// reference) is what is held in memory, so cost scales with the panel rather than the
/// genome.
///
/// Unlike the read-mapper path this reports *partial* placements — a probe matching
/// over half its length still forms a duplex long enough to capture — which is the
/// behaviour that matters for panels designed against relatives of the sampled
/// organism.
pub fn probe_hits_seed_extend(
    sequences: &HashMap<String, String>,
    probes_path: &Path,
    weights: &HashMap<String, f64>,
    mode: SimulateMode,
    model: &ThermoModel,
    threads: usize,
) -> Result<HashMap<String, Vec<ProbeHit>>> {
    use crate::alignment::seed_extend::{find_local_hits, ProbeIndex, SeedExtendParams};

    // Only sequences that can actually contribute fragments are worth aligning to.
    let mut ref_ids: Vec<&String> = sequences
        .keys()
        .filter(|id| weights.get(*id).copied().unwrap_or(0.0) > 0.0)
        .collect();
    ref_ids.sort_unstable();
    let ref_seqs: Vec<Vec<u8>> = ref_ids
        .iter()
        .map(|id| sequences[*id].as_bytes().to_vec())
        .collect();

    let probe_map = crate::fasta::parse_fasta(probes_path)?;
    let mut probe_names: Vec<&String> = probe_map.keys().collect();
    probe_names.sort_unstable();
    let probe_seqs: Vec<Vec<u8>> = probe_names
        .iter()
        .map(|n| probe_map[*n].as_bytes().to_vec())
        .collect();

    let params = SeedExtendParams::default();
    let index = ProbeIndex::build(&probe_seqs, params.seed_len);
    let per_probe = find_local_hits(&index, &ref_seqs, &params, threads);

    let mut hits_by_probe: HashMap<String, Vec<ProbeHit>> = HashMap::new();
    let mut n_parsed = 0usize;
    for (pi, hits) in per_probe.iter().enumerate() {
        if hits.is_empty() {
            continue;
        }
        let name = probe_names[pi].to_string();
        let entry = hits_by_probe.entry(name).or_default();
        for h in hits {
            let seq_id = ref_ids[h.seq_idx].to_string();
            let seq_weight = weights.get(&seq_id).copied().unwrap_or(0.0);
            if seq_weight <= 0.0 {
                continue;
            }
            let score = match mode {
                SimulateMode::Thermodynamic => {
                    // Pairs come back as (probe, reference); delta_g wants the base on
                    // the strand the probe binds, i.e. the reference's complement.
                    let pairs: Vec<(u8, u8)> = h
                        .pairs
                        .iter()
                        .map(|&(p, r)| (p, wc_complement_base(r)))
                        .collect();
                    boltzmann_score(delta_g(&pairs, model), model) * seq_weight
                }
                SimulateMode::Simple => seq_weight,
            };
            entry.push(ProbeHit {
                seq_id,
                start: h.ref_start,
                end: h.ref_end,
                score,
            });
            n_parsed += 1;
        }
    }

    for hits in hits_by_probe.values_mut() {
        hits.sort_unstable_by(|a, b| {
            a.seq_id
                .cmp(&b.seq_id)
                .then(a.start.cmp(&b.start))
                .then(a.end.cmp(&b.end))
        });
    }

    log::info!(
        "Placed {} probe hits with seed-and-extend over {} reference sequences",
        n_parsed,
        ref_seqs.len()
    );
    Ok(hits_by_probe)
}

/// Load probe hits from a blastn tabular file produced by `external::blastn::probe_align`.
///
/// Columns: `qseqid sseqid sstart send sstrand qseq sseq`.
///
/// Pairs come straight from the aligned `qseq`/`sseq` strings, so unlike the SAM path
/// there is no CIGAR walk, no MD tag, and no strand handling — BLAST reports both
/// strands already in alignment orientation. Gap columns emit no pair, matching
/// `aligned_pairs`.
pub fn load_probe_hits_blast(
    tsv_path: &Path,
    weights: &HashMap<String, f64>,
    mode: SimulateMode,
    model: &ThermoModel,
) -> Result<HashMap<String, Vec<ProbeHit>>> {
    let file = File::open(tsv_path)
        .with_context(|| format!("Cannot open BLAST output: {}", tsv_path.display()))?;
    let reader = BufReader::new(file);

    let mut hits_by_probe: HashMap<String, Vec<ProbeHit>> = HashMap::new();
    let mut n_parsed = 0usize;
    let mut n_skipped = 0usize;

    for line_result in reader.lines() {
        let line = line_result?;
        if line.trim().is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 7 {
            continue;
        }
        let probe_name = f[0].to_string();
        let rname = f[1];
        let sstart: usize = f[2].parse().unwrap_or(0);
        let send: usize = f[3].parse().unwrap_or(0);
        let qseq = f[5].as_bytes();
        let sseq = f[6].as_bytes();

        let seq_weight = weights.get(rname).copied().unwrap_or(0.0);
        if seq_weight <= 0.0 {
            n_skipped += 1;
            continue;
        }

        // BLAST coordinates are 1-based inclusive and descend on the minus strand.
        let (lo, hi) = if sstart <= send { (sstart, send) } else { (send, sstart) };
        let start = lo.saturating_sub(1);
        let end = hi;

        let score = match mode {
            SimulateMode::Thermodynamic => {
                let mut pairs: Vec<(u8, u8)> = Vec::with_capacity(qseq.len());
                for (&q, &r) in qseq.iter().zip(sseq.iter()) {
                    if q == b'-' || r == b'-' {
                        continue; // gap column: no base pair, as in aligned_pairs
                    }
                    pairs.push((
                        q.to_ascii_uppercase(),
                        wc_complement_base(r.to_ascii_uppercase()),
                    ));
                }
                let dg = delta_g(&pairs, model);
                boltzmann_score(dg, model) * seq_weight
            }
            SimulateMode::Simple => seq_weight,
        };

        hits_by_probe.entry(probe_name).or_default().push(ProbeHit {
            seq_id: rname.to_string(),
            start,
            end,
            score,
        });
        n_parsed += 1;
    }

    for hits in hits_by_probe.values_mut() {
        hits.sort_unstable_by(|a, b| {
            a.seq_id
                .cmp(&b.seq_id)
                .then(a.start.cmp(&b.start))
                .then(a.end.cmp(&b.end))
        });
    }

    log::info!(
        "Loaded {} probe hits ({} skipped) from {}",
        n_parsed,
        n_skipped,
        tsv_path.display()
    );
    Ok(hits_by_probe)
}

/// Parse a probe-vs-reference SAM file and return hits grouped by probe name.
///
/// `sequences` — reference sequences keyed by seq_id (from combined_reference.fa).
/// `weights`   — per-seq_id sampling weights (0.0 = skip).
/// `mode`      — thermodynamic or simple.
/// `model`     — thermodynamic model (temperature + salt); only used in Thermodynamic mode.
pub fn load_probe_hits(
    sam_path: &Path,
    sequences: &HashMap<String, String>,
    weights: &HashMap<String, f64>,
    mode: SimulateMode,
    model: &ThermoModel,
) -> Result<HashMap<String, Vec<ProbeHit>>> {
    let file = File::open(sam_path)
        .with_context(|| format!("Cannot open SAM file: {}", sam_path.display()))?;
    let reader = BufReader::new(file);

    let mut hits_by_probe: HashMap<String, Vec<ProbeHit>> = HashMap::new();
    let mut n_parsed = 0usize;
    let mut n_skipped = 0usize;

    for line_result in reader.lines() {
        let line = line_result?;
        if line.starts_with('@') {
            continue; // SAM header
        }

        let fields: Vec<&str> = line.splitn(12, '\t').collect();
        if fields.len() < 11 {
            continue;
        }

        let probe_name = fields[0].to_string();
        let flag: u16 = fields[1].parse().unwrap_or(0);
        let rname = fields[2];
        let pos: usize = fields[3].parse().unwrap_or(0);
        let cigar = fields[5];
        let seq = fields[9].as_bytes();

        // Skip unmapped reads
        if flag & FLAG_UNMAPPED != 0 || rname == "*" || cigar == "*" {
            n_skipped += 1;
            continue;
        }

        // Skip sequences with zero weight
        let seq_weight = weights.get(rname).copied().unwrap_or(0.0);
        if seq_weight <= 0.0 {
            n_skipped += 1;
            continue;
        }

        // Compute alignment end from CIGAR (reference-consuming operations: M D N = X)
        let ref_len = cigar_ref_len(cigar);
        let start = pos.saturating_sub(1); // SAM is 1-based
        let end = start + ref_len;

        let score = match mode {
            SimulateMode::Thermodynamic => {
                // SEQ is already in reference orientation, so pairs can be read
                // straight off the reference.
                let reference = match sequences.get(rname) {
                    Some(r) => r.as_bytes(),
                    None => {
                        n_skipped += 1;
                        continue;
                    }
                };
                let pairs = aligned_pairs(seq, reference, start, cigar);
                let dg = delta_g(&pairs, model);
                boltzmann_score(dg, model) * seq_weight
            }
            SimulateMode::Simple => seq_weight,
        };

        hits_by_probe
            .entry(probe_name)
            .or_default()
            .push(ProbeHit {
                seq_id: rname.to_string(),
                start,
                end,
                score,
            });

        n_parsed += 1;
    }

    // Probes are read from a HashMap upstream, so SAM record order — and therefore the
    // order of each probe's hit list — varies between runs. Level-2 sampling indexes
    // into these lists, so sort them for reproducibility under a fixed --seed.
    for hits in hits_by_probe.values_mut() {
        hits.sort_unstable_by(|a, b| {
            a.seq_id
                .cmp(&b.seq_id)
                .then(a.start.cmp(&b.start))
                .then(a.end.cmp(&b.end))
        });
    }

    log::info!(
        "Loaded {} probe hits ({} skipped) from {}",
        n_parsed,
        n_skipped,
        sam_path.display()
    );

    Ok(hits_by_probe)
}


/// Sample `n_capture` capture-derived (probe-site-biased) fragments.
///
/// Returns `(header, sequence)` pairs using BaitBench naming:
///   `>{seq_id}_fragment_{counter} start={frag_start} length={frag_len}`
pub fn sample_capture_fragments(
    hits_by_probe: &HashMap<String, Vec<ProbeHit>>,
    sequences: &HashMap<String, String>,
    n_capture: usize,
    len_params: FragmentLengthParams,
    seed: Option<u64>,
    counter_start: usize,
) -> Result<Vec<(String, String)>> {
    let fragment_length_mean = len_params.mean;
    let fragment_length_min = len_params.min;
    let fragment_length_max = len_params.max;
    // Collect probes that have at least one positive-score hit.
    //
    // Sorted: `hits_by_probe` is a HashMap, and `RandomState` reseeds per process, so
    // iteration order varies between runs. Since the RNG indexes into this vector, an
    // unsorted order makes `--seed` select different probes on every invocation.
    let mut probe_names: Vec<&String> = hits_by_probe
        .iter()
        .filter(|(_, hits)| hits.iter().any(|h| h.score > 0.0))
        .map(|(name, _)| name)
        .collect();
    probe_names.sort_unstable();

    if probe_names.is_empty() {
        bail!("No probes with positive-score hits found — cannot sample capture fragments.");
    }

    let mut rng: StdRng = match seed {
        Some(s) => StdRng::seed_from_u64(s.wrapping_add(1)),
        None => StdRng::from_entropy(),
    };

    // Level 1: uniform probe selection
    let probe_dist = WeightedIndex::new(vec![1.0f64; probe_names.len()])?;

    let sd = (fragment_length_max as f64 - fragment_length_min as f64) / 6.0;
    let normal = Normal::new(fragment_length_mean, sd)?;

    let mut results = Vec::with_capacity(n_capture);
    let mut counter = counter_start;
    let mut attempts = 0usize;
    let max_attempts = n_capture * 10;

    while results.len() < n_capture && attempts < max_attempts {
        attempts += 1;

        // Level 1: pick a probe
        let probe_idx = probe_dist.sample(&mut rng);
        let probe_name = probe_names[probe_idx];
        let hits = &hits_by_probe[probe_name];

        // Level 2: pick a hit weighted by score
        let scores: Vec<f64> = hits.iter().map(|h| h.score.max(0.0)).collect();
        if scores.iter().all(|&s| s == 0.0) {
            continue;
        }
        let hit_dist = WeightedIndex::new(&scores)?;
        let hit = &hits[hit_dist.sample(&mut rng)];

        // Get reference sequence
        let seq = match sequences.get(&hit.seq_id) {
            Some(s) => s,
            None => continue,
        };

        // Level 4: fragment length (sample before centering so we can clamp correctly)
        let frag_len_f: f64 = normal.sample(&mut rng);
        let frag_len = (frag_len_f.round() as usize)
            .max(fragment_length_min)
            .min(fragment_length_max)
            .min(seq.len());

        if frag_len == 0 {
            continue;
        }

        // Level 3: fragment center = probe_center ± jitter (±frag_len/4)
        let probe_center = (hit.start + hit.end) / 2;
        let jitter_range = (frag_len / 4).max(1);
        let jitter: i64 = rng.gen_range(-(jitter_range as i64)..=(jitter_range as i64));

        let half = frag_len / 2;
        let center_raw = probe_center as i64 + jitter;
        // Clamp center so fragment fits within sequence
        let center = center_raw
            .max(half as i64)
            .min((seq.len() - frag_len + half) as i64) as usize;

        let frag_start = center.saturating_sub(half);
        let frag_end = (frag_start + frag_len).min(seq.len());
        let actual_len = frag_end - frag_start;

        let fragment = &seq[frag_start..frag_end];

        counter += 1;
        let header = format!(
            ">{}_fragment_{} start={} length={}",
            hit.seq_id, counter, frag_start, actual_len
        );
        results.push((header, fragment.to_string()));
    }

    if results.len() < n_capture {
        log::warn!(
            "Only generated {}/{} capture fragments after {} attempts",
            results.len(),
            n_capture,
            attempts
        );
    }

    Ok(results)
}

/// Sample `n_background` background fragments, weighted by sequence_weight × seq_len.
///
/// Returns (header, sequence) pairs rather than writing to a file directly (so
/// they can be interleaved with capture fragments before writing).
pub fn sample_background_fragments(
    sequences: &HashMap<String, String>,
    weights: &HashMap<String, f64>,
    n_background: usize,
    len_params: FragmentLengthParams,
    seed: Option<u64>,
    counter_start: usize,
) -> Result<Vec<(String, String)>> {
    let fragment_length_mean = len_params.mean;
    let fragment_length_min = len_params.min;
    let fragment_length_max = len_params.max;
    // Filter to sequences with positive weights
    // Sorted for the same reason as the capture sampler: the RNG indexes into this
    // vector, so HashMap iteration order would otherwise leak into a seeded run.
    let mut weighted_seqs: Vec<(&String, &String)> = sequences
        .iter()
        .filter(|(id, _)| weights.get(*id).copied().unwrap_or(0.0) > 0.0)
        .collect();
    weighted_seqs.sort_unstable_by(|a, b| a.0.cmp(b.0));

    if weighted_seqs.is_empty() {
        bail!("No sequences with positive weights found for background sampling!");
    }

    let seq_ids: Vec<&String> = weighted_seqs.iter().map(|(id, _)| *id).collect();
    let seq_weights: Vec<f64> = weighted_seqs
        .iter()
        .map(|(id, seq)| weights[*id] * seq.len() as f64)
        .collect();

    let dist = WeightedIndex::new(&seq_weights)?;
    let sd = (fragment_length_max as f64 - fragment_length_min as f64) / 6.0;
    let normal = Normal::new(fragment_length_mean, sd)?;

    let mut rng: StdRng = match seed {
        Some(s) => StdRng::seed_from_u64(s.wrapping_add(2)),
        None => StdRng::from_entropy(),
    };

    let mut results = Vec::with_capacity(n_background);
    let mut counter = counter_start;

    for _ in 0..n_background {
        let idx = dist.sample(&mut rng);
        let seq_id = seq_ids[idx];
        let seq = weighted_seqs[idx].1;

        let frag_len_f: f64 = normal.sample(&mut rng);
        let frag_len = (frag_len_f.round() as usize)
            .max(fragment_length_min)
            .min(fragment_length_max);

        if seq.len() < frag_len || frag_len == 0 {
            continue;
        }

        let max_start = seq.len() - frag_len;
        let start = rng.gen_range(0..=max_start);
        let fragment = &seq[start..start + frag_len];

        counter += 1;
        let header = format!(
            ">{}_fragment_{} start={} length={}",
            seq_id, counter, start, frag_len
        );
        results.push((header, fragment.to_string()));
    }

    Ok(results)
}

/// Write a list of (header, sequence) fragment pairs to a FASTA file.
pub fn write_fragments(fragments: &[(String, String)], output: &Path) -> Result<()> {
    let file = File::create(output)
        .with_context(|| format!("Cannot create output file: {}", output.display()))?;
    let mut writer = BufWriter::new(file);
    for (header, seq) in fragments {
        writeln!(writer, "{}", header)?;
        writeln!(writer, "{}", seq)?;
    }
    writer.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::seq_utils::reverse_complement;
    use crate::thermodynamics::ThermoModel;

    fn model() -> ThermoModel {
        ThermoModel::new(65.0, 0.05)
    }

    /// 40 bp of mixed-composition reference to align against.
    const REF: &str = "ACGTTGCAAGGCCTTAGCATCGATCGGATCCTAGCATGCA";

    fn pairs_of(query: &str, reference: &str, start: usize, cigar: &str) -> Vec<(u8, u8)> {
        aligned_pairs(query.as_bytes(), reference.as_bytes(), start, cigar)
    }

    /// A probe identical to the reference must pair Watson-Crick at every position:
    /// it binds the opposite strand, so the partner base is complement(reference).
    #[test]
    fn perfect_match_is_all_watson_crick() {
        let q = &REF[0..20];
        let pairs = pairs_of(q, REF, 0, "20M");
        assert_eq!(pairs.len(), 20);
        for (i, (probe_b, partner)) in pairs.iter().enumerate() {
            assert_eq!(*probe_b, q.as_bytes()[i]);
            assert_eq!(*partner, wc_complement_base(*probe_b), "position {i}");
        }
    }

    /// Regression for defect D3. A probe `A` over reference `T` is a real mismatch —
    /// the probe faces the `A` on the opposite strand. The old MD-based builder stored
    /// the uncomplemented reference base, so `is_wc` declared it Watson-Crick and
    /// awarded full stacking energy.
    #[test]
    fn transition_mismatch_is_not_watson_crick() {
        // Reference position 0 is 'A'; make the probe 'T' there.
        let mut q: Vec<u8> = REF.as_bytes()[0..20].to_vec();
        q[0] = b'T';
        let pairs = aligned_pairs(&q, REF.as_bytes(), 0, "20M");

        // Partner is complement(reference 'A') = 'T', same as the probe base: not WC.
        assert_eq!(pairs[0], (b'T', b'T'));

        let perfect = pairs_of(&REF[0..20], REF, 0, "20M");
        let m = model();
        assert!(
            delta_g(&pairs, &m) > delta_g(&perfect, &m),
            "a mismatched probe must be less stable than a perfect one"
        );
    }

    /// The reverse-strand path must agree with the forward one. The SAM writer emits
    /// SEQ in reference orientation, so a probe whose reverse complement matches the
    /// reference arrives here already reverse-complemented — and must score identically
    /// to the equivalent forward probe (defect D1).
    #[test]
    fn reverse_strand_hit_scores_like_its_forward_equivalent() {
        let region = &REF[4..24];
        let probe_rev = reverse_complement(region.as_bytes()); // a minus-strand probe
        let oriented = reverse_complement(&probe_rev); // what the SAM writer emits
        assert_eq!(oriented, region.as_bytes());

        let fwd = pairs_of(region, REF, 4, "20M");
        let rev = aligned_pairs(&oriented, REF.as_bytes(), 4, "20M");
        let m = model();
        assert_eq!(delta_g(&fwd, &m), delta_g(&rev, &m));
    }

    /// Soft clips consume query without pairing, so a clipped alignment must not shift
    /// the remaining bases against the reference (defect D5).
    #[test]
    fn soft_clips_do_not_shift_the_query() {
        // 5 clipped bases, then 15 aligned starting at reference 0.
        let q = format!("GGGGG{}", &REF[0..15]);
        let pairs = pairs_of(&q, REF, 0, "5S15M");
        let ref_bytes = REF.as_bytes();
        assert_eq!(pairs.len(), 15);
        for (i, (probe_b, partner)) in pairs.iter().enumerate() {
            assert_eq!(*probe_b, ref_bytes[i], "clipped prefix shifted the query");
            assert_eq!(*partner, wc_complement_base(ref_bytes[i]));
        }
    }

    /// Insertions consume query only, deletions reference only.
    #[test]
    fn indels_advance_the_correct_strand() {
        // 5M then a 2 bp insertion then 5M: reference advances 10, query 12.
        let q = format!("{}XX{}", &REF[0..5], &REF[5..10]);
        let pairs = pairs_of(&q, REF, 0, "5M2I5M");
        assert_eq!(pairs.len(), 10, "insertion must not pair");
        assert_eq!(pairs[5].0, REF.as_bytes()[5], "reference lost sync after insertion");

        // 5M then a 2 bp deletion then 5M: reference advances 12, query 10.
        let q2 = format!("{}{}", &REF[0..5], &REF[7..12]);
        let pairs2 = pairs_of(&q2, REF, 0, "5M2D5M");
        assert_eq!(pairs2.len(), 10, "deletion must not pair");
        assert_eq!(pairs2[5].0, REF.as_bytes()[7], "query lost sync after deletion");
    }

    /// Lowercase input must score identically to uppercase; ambiguity codes and N
    /// must break the duplex rather than silently pairing.
    #[test]
    fn case_is_normalised_and_ambiguous_bases_break_pairing() {
        let upper = pairs_of(&REF[0..20], REF, 0, "20M");
        let lower = pairs_of(&REF[0..20].to_lowercase(), &REF.to_lowercase(), 0, "20M");
        assert_eq!(upper, lower);

        let mut q: Vec<u8> = REF.as_bytes()[0..20].to_vec();
        q[3] = b'N';
        let with_n = aligned_pairs(&q, REF.as_bytes(), 0, "20M");
        assert_eq!(with_n[3].0, b'N');

        // An IUPAC code in the reference yields an 'N' partner, which is not WC.
        let mut r: Vec<u8> = REF.as_bytes().to_vec();
        r[3] = b'Y';
        let with_iupac = aligned_pairs(&REF.as_bytes()[0..20], &r, 0, "20M");
        assert_eq!(with_iupac[3].1, b'N');
    }

    /// Alignments running past the end of the reference must truncate, not panic.
    #[test]
    fn alignment_past_reference_end_is_truncated() {
        let q = &REF[30..40];
        let pairs = pairs_of(q, REF, 35, "10M");
        assert_eq!(pairs.len(), 5);
    }
}
