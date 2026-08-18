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

/// Parse a CIGAR string into a list of (probe_base, ref_base) pairs suitable
/// for delta_g(). The SEQ field of the SAM record is the probe sequence.
/// We derive the aligned pairs by expanding the CIGAR.
///
/// For matches (M/=), we take the probe base and use 'N' as a placeholder for
/// the reference base — but the reference sequence is not stored in the SAM.
/// Instead, we determine whether a position is a WC match from the MD tag or
/// by knowing the reference from our FASTA. Since we only need the *presence*
/// of a WC pair (to decide whether to accumulate stacking), and the SAM doesn't
/// always carry the reference bases, we use a pragmatic approach:
///
/// - Parse the MD tag (if present) to recover mismatched / deleted reference bases.
/// - Positions covered by M/= with no MD mismatch annotation are treated as matches
///   (top and bottom are WC complements of each other).
/// - Positions with an MD mismatch are given the actual ref base from MD.
///
/// This yields accurate stacking energy for typical probe alignments.
fn cigar_and_md_to_pairs(
    probe_seq: &[u8],
    cigar: &str,
    md: Option<&str>,
) -> Vec<(u8, u8)> {
    // Expand CIGAR into per-position operations
    let ops = expand_cigar(cigar);

    // Build MD-based mismatch map: query_position → ref_base
    let md_mismatches = md
        .map(parse_md_mismatches)
        .unwrap_or_default();

    let mut pairs: Vec<(u8, u8)> = Vec::with_capacity(probe_seq.len());
    let mut qi = 0usize; // index into probe_seq

    for op in &ops {
        match op {
            CigarOp::Match | CigarOp::Equal => {
                if qi < probe_seq.len() {
                    let probe_b = probe_seq[qi].to_ascii_uppercase();
                    let ref_b = if let Some(&rb) = md_mismatches.get(&qi) {
                        rb
                    } else {
                        // Treat as a perfect WC match: give the complement base
                        wc_complement_base(probe_b)
                    };
                    pairs.push((probe_b, ref_b));
                    qi += 1;
                }
            }
            CigarOp::Mismatch => {
                if qi < probe_seq.len() {
                    let probe_b = probe_seq[qi].to_ascii_uppercase();
                    // Try MD for the actual ref base; fall back to 'N'
                    let ref_b = md_mismatches.get(&qi).copied().unwrap_or(b'N');
                    pairs.push((probe_b, ref_b));
                    qi += 1;
                }
            }
            CigarOp::Ins | CigarOp::SoftClip => {
                // Consumes query only; no reference pairing
                qi += 1;
            }
            CigarOp::Del | CigarOp::Skip | CigarOp::HardClip | CigarOp::Pad => {
                // Consumes reference only (or nothing from query)
            }
        }
    }

    pairs
}


/// Parse an MD tag and return a map of query_position → mismatched_ref_base.
/// MD format example: "10A5^TT3" means 10 matches, A mismatch, 5 matches, TT deletion, 3 matches.
/// Deletion bases (after '^') are reference-only and must NOT advance the query position.
fn parse_md_mismatches(md: &str) -> HashMap<usize, u8> {
    let mut map = HashMap::new();
    let mut qi = 0usize;
    let mut num = String::new();
    let mut in_deletion = false;

    for ch in md.chars() {
        if ch.is_ascii_digit() {
            num.push(ch);
            in_deletion = false; // a digit always terminates the deletion run
        } else if ch == '^' {
            // Flush any pending match count before entering the deletion
            if !num.is_empty() {
                qi += num.parse::<usize>().unwrap_or(0);
                num.clear();
            }
            in_deletion = true;
        } else if ch.is_ascii_alphabetic() {
            if in_deletion {
                // Reference-only deletion base — skip without advancing query position
            } else {
                // Flush match count, then record mismatch at current query position
                if !num.is_empty() {
                    qi += num.parse::<usize>().unwrap_or(0);
                    num.clear();
                }
                map.insert(qi, ch.to_ascii_uppercase() as u8);
                qi += 1;
            }
        }
    }

    map
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

/// Parse a probe-vs-reference SAM file and return hits grouped by probe name.
///
/// `sequences` — reference sequences keyed by seq_id (from combined_reference.fa).
/// `weights`   — per-seq_id sampling weights (0.0 = skip).
/// `mode`      — thermodynamic or simple.
/// `model`     — thermodynamic model (temperature + salt); only used in Thermodynamic mode.
pub fn load_probe_hits(
    sam_path: &Path,
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

        // Parse MD tag from optional fields
        let md = fields
            .get(11)
            .and_then(|_| {
                // Re-split the full line to get all optional fields
                line.split('\t')
                    .find(|f| f.starts_with("MD:Z:"))
                    .map(|f| &f[5..])
            });

        let score = match mode {
            SimulateMode::Thermodynamic => {
                let pairs = cigar_and_md_to_pairs(seq, cigar, md);
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
    // Collect probes that have at least one positive-score hit
    let probe_names: Vec<&String> = hits_by_probe
        .iter()
        .filter(|(_, hits)| hits.iter().any(|h| h.score > 0.0))
        .map(|(name, _)| name)
        .collect();

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
    let weighted_seqs: Vec<(&String, &String)> = sequences
        .iter()
        .filter(|(id, _)| weights.get(*id).copied().unwrap_or(0.0) > 0.0)
        .collect();

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
