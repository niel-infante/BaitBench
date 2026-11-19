/// Probe-biased fragment simulation by per-fragment retention.
///
/// Every fragment in the library gets one retention probability,
///
/// ```text
///     p_retain(f) = 1 − (1 − p_ns) · Π_i (1 − θ_i)
/// ```
///
/// over the probes `i` whose alignments overlap `f`, where `θ_i = K/(1+K)` is Langmuir
/// occupancy and `p_ns` is non-specific carryover. A fragment overlapping no probe
/// reduces to `p_ns`; one overlapping probes gets carryover plus each probe's
/// independent chance. Output composition therefore *emerges* from abundance ×
/// efficiency rather than being imposed as a quota, which is what lets a dilution
/// series titrate and what keeps a rare organism from being manufactured into reads it
/// has no template for.
///
/// `p_retain` is piecewise constant in position — it changes only where a probe's
/// overlap window starts or ends — so the per-reference profile is integrated exactly
/// and sampled in two stages (reference ∝ abundance × retention mass, then position ∝
/// `p_retain`). Every draw is a keeper, so no rejection loop is needed even when the
/// library is 99.99% off-target.
///
/// This replaces a fixed-quota two-arm sampler that drew a probe uniformly — giving ΔG
/// no influence on *which* probe captured — and filled a `--capture-fraction` quota,
/// which pinned realized composition and capped the realized distractor fraction at
/// `1 − capture_fraction`.
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
use crate::thermodynamics::{delta_g, langmuir_occupancy, ThermoModel};

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
    /// Duplex free energy in kcal/mol, **without** the sequence weight folded in.
    ///
    /// The capture model turns this into occupancy θ = K/(1+K) at the probe
    /// concentration being tested, while abundance enters separately when drawing which
    /// reference a fragment comes from. Keeping the two apart is what lets the
    /// `--capture-fraction` solver sweep concentration without rescoring any alignment.
    ///
    /// In `Simple` mode this is `SIMPLE_DG` for every hit, so all placed probes share
    /// one affinity and the solver sets its absolute level.
    pub dg: f64,
}

/// ΔG assigned to every placed hit in `SimulateMode::Simple`.
///
/// Zero gives K = [probe] at any temperature, so θ = 0.5 at unit concentration and
/// `Simple` mode stays temperature-independent as documented. The absolute level is
/// arbitrary — the `--capture-fraction` solver rescales concentration to reach the
/// requested on-target fraction — what matters is that every hit shares it.
pub const SIMPLE_DG: f64 = 0.0;

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
            let dg = match mode {
                SimulateMode::Thermodynamic => {
                    // Pairs come back as (probe, reference); delta_g wants the base on
                    // the strand the probe binds, i.e. the reference's complement.
                    let pairs: Vec<(u8, u8)> = h
                        .pairs
                        .iter()
                        .map(|&(p, r)| (p, wc_complement_base(r)))
                        .collect();
                    delta_g(&pairs, model)
                }
                SimulateMode::Simple => SIMPLE_DG,
            };
            entry.push(ProbeHit {
                seq_id,
                start: h.ref_start,
                end: h.ref_end,
                dg,
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

        let dg = match mode {
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
                delta_g(&pairs, model)
            }
            SimulateMode::Simple => SIMPLE_DG,
        };

        hits_by_probe.entry(probe_name).or_default().push(ProbeHit {
            seq_id: rname.to_string(),
            start,
            end,
            dg,
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

        let dg = match mode {
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
                delta_g(&pairs, model)
            }
            SimulateMode::Simple => SIMPLE_DG,
        };

        hits_by_probe
            .entry(probe_name)
            .or_default()
            .push(ProbeHit {
                seq_id: rname.to_string(),
                start,
                end,
                dg,
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

// ─── Retention model ──────────────────────────────────────────────────────────

/// Non-specific carryover: the probability that a fragment overlapping no probe is
/// retained anyway. Washing is imperfect, and in a library that is 99.99% off-target
/// this term, not capture, supplies most of the background reads.
pub const DEFAULT_BACKGROUND_RETENTION: f64 = 1e-3;

/// Fraction of a probe's target site a fragment must carry before that probe can pull
/// it. See `build_retention_profile` for why a threshold is used here.
pub const DEFAULT_MIN_OVERLAP_FRAC: f64 = 0.5;

/// Physical parameters of the capture reaction.
#[derive(Debug, Clone, Copy)]
pub struct RetentionParams {
    /// Probe concentration, in the arbitrary units absorbed into K. Set directly, or
    /// back-solved from a requested on-target fraction by `solve_probe_concentration`.
    pub probe_concentration: f64,
    /// Non-specific carryover probability per fragment.
    pub background_retention: f64,
    /// Minimum fraction of a probe's site a fragment must carry.
    pub min_overlap_frac: f64,
    /// Fragment length the profile is built for.
    pub frag_len: usize,
}

/// A run of fragment start positions over which the retention probability is constant.
#[derive(Debug, Clone)]
pub struct RetentionSegment {
    /// First fragment start position in this run.
    pub start: usize,
    /// One past the last fragment start position in this run.
    pub end: usize,
    /// Probability that a fragment starting anywhere in this run is retained at all.
    pub p_retain: f64,
    /// Probability that it is retained *by a probe* rather than carried over. Always
    /// `<= p_retain`.
    pub p_specific: f64,
}

/// Per-reference retention profile: where capture can happen, and how strongly.
#[derive(Debug, Clone)]
pub struct RetentionProfile {
    pub seq_id: String,
    pub segments: Vec<RetentionSegment>,
    /// Cumulative `len * p_retain` over `segments`, for sampling a start position.
    cum: Vec<f64>,
    /// Total retention mass, `Σ len * p_retain`. The expected number of fragments
    /// retained from this reference per unit of abundance.
    pub mass: f64,
    /// The part of `mass` attributable to probe capture, `Σ len * p_specific`.
    pub mass_specific: f64,
    /// Number of start positions available at this fragment length — the reference's
    /// contribution to the *input* library, against which yield is measured.
    pub n_starts: usize,
    /// Of those, how many a probe can reach at all. The reference's contribution to the
    /// input library's on-target share, which is the denominator of enrichment.
    pub covered_starts: usize,
}

impl RetentionProfile {
    /// Sample a fragment start position proportional to `p_retain`.
    fn sample_start(&self, rng: &mut StdRng) -> Option<usize> {
        if self.mass <= 0.0 || self.segments.is_empty() {
            return None;
        }
        let target = rng.gen_range(0.0..self.mass);
        // Binary search the cumulative mass.
        let idx = match self
            .cum
            .binary_search_by(|v| v.partial_cmp(&target).unwrap_or(std::cmp::Ordering::Less))
        {
            Ok(i) => i,
            Err(i) => i,
        }
        .min(self.segments.len() - 1);
        let seg = &self.segments[idx];
        Some(rng.gen_range(seg.start..seg.end.max(seg.start + 1)))
    }
}

/// Walk one reference's retention profile, calling `emit(from, to, p_retain,
/// p_specific)` once per run of constant probability. Returns the number of available
/// fragment start positions, or 0 if the reference is shorter than a fragment.
///
/// A fragment `[s, s+frag_len)` overlaps a probe hit `[h.start, h.end)` by
/// `min(s+frag_len, h.end) - max(s, h.start)`, so requiring at least `m` bases of
/// overlap puts `s` in `[h.start + m - frag_len, h.end - m]`. Each hit therefore
/// contributes a window of start positions over which it can pull, and
///
/// ```text
///     p_retain(s) = 1 − (1 − p_ns) · Π_{i active at s} (1 − θ_i)
///     p_spec(s)   = 1 − Π_{i active at s} (1 − θ_i)
/// ```
///
/// is constant between window boundaries. The profile is thus piecewise constant with
/// at most `2·n_hits + 1` runs, and its integral — the retention mass — is exact
/// without visiting every base. That is what makes sampling possible without
/// rejection: at the Borrelia LOD's lowest input the library is 99.996% distractor, so
/// accept/reject would discard roughly a thousand candidates per keeper.
///
/// The product is accumulated in log space (`Σ ln(1 − θ)`), since hits enter and leave
/// the active set as the sweep advances and dividing back out would drift. Saturated
/// hits (`θ = 1`) are counted separately because `ln 0` is not finite.
fn sweep_retention<F>(
    seq_len: usize,
    hits: &[(usize, usize, f64)], // (start, end, theta)
    frag_len: usize,
    p_ns: f64,
    min_overlap_frac: f64,
    mut emit: F,
) -> usize
where
    F: FnMut(usize, usize, f64, f64),
{
    if frag_len == 0 || seq_len < frag_len {
        return 0;
    }
    let n_starts = seq_len - frag_len + 1;

    // (position, delta_log, delta_saturated)
    let mut events: Vec<(usize, f64, i32)> = Vec::with_capacity(hits.len() * 2);
    for &(h_start, h_end, theta) in hits {
        let hit_len = h_end.saturating_sub(h_start);
        let min_overlap = ((hit_len as f64) * min_overlap_frac).ceil().max(1.0) as usize;
        if min_overlap > frag_len.min(hit_len) {
            continue; // fragment can never carry enough of this site
        }
        let lo = (h_start + min_overlap).saturating_sub(frag_len);
        let hi = (h_end.saturating_sub(min_overlap) + 1).min(n_starts); // exclusive
        if lo >= hi {
            continue;
        }
        let theta = theta.clamp(0.0, 1.0);
        if theta >= 1.0 {
            events.push((lo, 0.0, 1));
            events.push((hi, 0.0, -1));
        } else if theta > 0.0 {
            let l = (1.0 - theta).ln();
            events.push((lo, l, 0));
            events.push((hi, -l, 0));
        }
    }
    events.sort_unstable_by_key(|e| e.0);

    let mut sum_log = 0.0f64;
    let mut n_sat = 0i32;
    let mut pos = 0usize;
    let mut i = 0usize;

    let run = |from: usize, to: usize, sum_log: f64, n_sat: i32, emit: &mut F| {
        if to <= from {
            return;
        }
        let (p_retain, p_spec) = if n_sat > 0 {
            (1.0, 1.0)
        } else {
            let prod = sum_log.exp();
            (1.0 - (1.0 - p_ns) * prod, 1.0 - prod)
        };
        if p_retain > 0.0 {
            emit(from, to, p_retain, p_spec);
        }
    };

    while i < events.len() {
        let at = events[i].0;
        run(pos, at.min(n_starts), sum_log, n_sat, &mut emit);
        while i < events.len() && events[i].0 == at {
            sum_log += events[i].1;
            n_sat += events[i].2;
            i += 1;
        }
        pos = at;
    }
    run(pos, n_starts, sum_log, n_sat, &mut emit);

    n_starts
}

/// Build the retention profile for one reference.
///
/// `min_overlap_frac` is the fraction of a probe's target site a fragment must carry
/// before that probe can pull it. Without it, a fragment sharing a single base with a
/// hit would be credited the full θ of that hit — but θ was computed from the whole
/// duplex, and a 20 bp sliver of a 120 bp site forms nothing stable at hybridization
/// temperature. This is a threshold, not a model of partial duplexes: the physical
/// treatment is to recompute ΔG over the overlapping portion, which is recorded in
/// `improvements/partial-duplex-capture.md`.
///
/// **Approximations:**
/// * The profile is built for a single `frag_len`. Lengths are drawn per fragment, so a
///   fragment longer than `frag_len` may overlap a hit the profile did not credit. The
///   effect is second-order — hit windows are `frag_len + hit_len` wide, so they are
///   dominated by the fragment length — but making it exact means one profile per
///   length bin.
/// * A probe that clears `min_overlap_frac` is credited its full θ regardless of how
///   much of the site the fragment carries beyond that.
pub fn build_retention_profile(
    seq_id: &str,
    seq_len: usize,
    hits: &[(usize, usize, f64)], // (start, end, theta)
    frag_len: usize,
    p_ns: f64,
    min_overlap_frac: f64,
) -> RetentionProfile {
    let mut segments: Vec<RetentionSegment> = Vec::new();
    let mut cum: Vec<f64> = Vec::new();
    let mut mass = 0.0f64;
    let mut mass_specific = 0.0f64;
    let mut covered_starts = 0usize;

    let n_starts = sweep_retention(
        seq_len,
        hits,
        frag_len,
        p_ns,
        min_overlap_frac,
        |from, to, p_retain, p_spec| {
            let len = (to - from) as f64;
            mass += len * p_retain;
            mass_specific += len * p_spec;
            if p_spec > 0.0 {
                covered_starts += to - from;
            }
            segments.push(RetentionSegment {
                start: from,
                end: to,
                p_retain,
                p_specific: p_spec,
            });
            cum.push(mass);
        },
    );

    RetentionProfile {
        seq_id: seq_id.to_string(),
        segments,
        cum,
        mass,
        mass_specific,
        n_starts,
        covered_starts,
    }
}

/// Probe hits regrouped by reference, carrying ΔG.
///
/// The capture model asks "which probes overlap this fragment", so hits must be indexed
/// by reference rather than by probe. ΔG is carried rather than θ so that the
/// concentration solver can re-derive occupancy without rescoring alignments.
///
/// A `BTreeMap` rather than a `HashMap`: the sampler indexes into a vector built from
/// this, and `RandomState` reseeds per process, so hash order would leak into a seeded
/// run.
pub type HitsBySequence = std::collections::BTreeMap<String, Vec<(usize, usize, f64)>>;

/// Regroup probe hits by reference sequence.
pub fn group_hits_by_sequence(hits_by_probe: &HashMap<String, Vec<ProbeHit>>) -> HitsBySequence {
    let mut by_seq: HitsBySequence = Default::default();
    for hits in hits_by_probe.values() {
        for h in hits {
            by_seq
                .entry(h.seq_id.clone())
                .or_default()
                .push((h.start, h.end, h.dg));
        }
    }
    // Sorted so the event sweep is deterministic for equal positions.
    for v in by_seq.values_mut() {
        v.sort_unstable_by(|a, b| {
            a.0.cmp(&b.0)
                .then(a.1.cmp(&b.1))
                .then(a.2.partial_cmp(&b.2).unwrap_or(std::cmp::Ordering::Equal))
        });
    }
    by_seq
}

/// Occupancies below this contribute nothing measurable but still cost two sweep
/// events, so they are dropped. At 1e-9, a million such hits over one reference shift
/// `p_retain` by under 1e-3 of the carryover floor.
const MIN_THETA: f64 = 1e-9;

/// Convert stored ΔG values to occupancies at a given probe concentration.
fn thetas_at(
    hits: &[(usize, usize, f64)],
    probe_conc: f64,
    model: &ThermoModel,
) -> Vec<(usize, usize, f64)> {
    hits.iter()
        .filter_map(|&(s, e, dg)| {
            let theta = langmuir_occupancy(dg, probe_conc, model);
            if theta > MIN_THETA {
                Some((s, e, theta))
            } else {
                None
            }
        })
        .collect()
}

/// What the capture reaction does to the library as a whole.
#[derive(Debug, Clone, Copy)]
pub struct RetentionStats {
    /// Expected fraction of input fragments retained — the library-wide yield.
    pub yield_fraction: f64,
    /// Expected fraction of *retained* fragments that a probe pulled, rather than
    /// non-specific carryover. This is what `--capture-fraction` requests.
    pub on_target_fraction: f64,
    /// On-target share of the output over the on-target share of the input: the
    /// enrichment factor the capture achieves.
    pub enrichment: f64,
    /// References contributing retention mass.
    pub n_references: usize,
}

/// Aggregate per-reference profiles into library-wide statistics, weighting each
/// reference by its abundance.
pub fn retention_stats(
    profiles: &[RetentionProfile],
    weights: &HashMap<String, f64>,
) -> RetentionStats {
    let mut w_mass = 0.0f64;
    let mut w_mass_spec = 0.0f64;
    let mut w_input = 0.0f64;
    let mut w_covered = 0.0f64;
    for p in profiles {
        let w = weights.get(&p.seq_id).copied().unwrap_or(0.0);
        if w <= 0.0 {
            continue;
        }
        w_mass += w * p.mass;
        w_mass_spec += w * p.mass_specific;
        w_input += w * p.n_starts as f64;
        w_covered += w * p.covered_starts as f64;
    }
    let yield_fraction = if w_input > 0.0 { w_mass / w_input } else { 0.0 };
    let on_target_fraction = if w_mass > 0.0 { w_mass_spec / w_mass } else { 0.0 };
    // Fold enrichment, as a bench protocol reports it: the on-target share of what came
    // out over the on-target share of what went in. "On-target" in the input means
    // positions a probe can reach at all — abundance-weighted, so a probe-covered
    // locus on a rare organism counts for little.
    //
    // Not mass_specific/w_input: that cancels against the numerator and leaves 1/yield,
    // which is not an enrichment factor.
    let input_on_target = if w_input > 0.0 {
        w_covered / w_input
    } else {
        0.0
    };
    let enrichment = if input_on_target > 0.0 {
        on_target_fraction / input_on_target
    } else {
        1.0
    };
    RetentionStats {
        yield_fraction,
        on_target_fraction,
        enrichment,
        n_references: profiles.len(),
    }
}

/// Build retention profiles for every reference with positive abundance.
///
/// References with no probe hits are included: they carry background-only mass, and in
/// a realistic library they are the overwhelming majority of it. Omitting them would
/// make the output pure capture and the on-target fraction meaningless.
pub fn build_profiles(
    sequences: &HashMap<String, String>,
    weights: &HashMap<String, f64>,
    by_seq: &HitsBySequence,
    params: &RetentionParams,
    model: &ThermoModel,
) -> Vec<RetentionProfile> {
    let no_hits: Vec<(usize, usize, f64)> = Vec::new();

    let mut ids: Vec<&String> = sequences
        .keys()
        .filter(|id| weights.get(*id).copied().unwrap_or(0.0) > 0.0)
        .collect();
    ids.sort_unstable();

    ids.into_iter()
        .map(|id| {
            let seq_len = sequences[id].len();
            let hits = by_seq.get(id).unwrap_or(&no_hits);
            let thetas = thetas_at(hits, params.probe_concentration, model);
            // A contig shorter than the nominal fragment length still yields a
            // fragment — a short one. Without this it would get zero mass and drop out
            // of the library entirely, which silently removes short distractor contigs
            // from the background they are supposed to supply.
            let frag_len = params.frag_len.min(seq_len).max(1);
            build_retention_profile(
                id,
                seq_len,
                &thetas,
                frag_len,
                params.background_retention,
                params.min_overlap_frac,
            )
        })
        .filter(|p| p.mass > 0.0)
        .collect()
}

/// The panel's capture efficiency at a trial probe concentration: given that a
/// fragment carries enough of a probe's site, how often is it pulled down?
///
/// Pooled over probe-covered positions only, and deliberately **independent of
/// abundance, of library composition, and of `background_retention`**. That
/// independence is the whole point.
///
/// An earlier version of this solver targeted the realized on-target fraction of the
/// output instead. That silently re-imposed composition: adding distractor made the
/// solver raise the probe concentration to hold on-target at the requested level, so
/// the realized distractor fraction saturated at `1 - capture_fraction` — the same
/// defect as the fixed quota it replaced (reviewer item 6), just arrived at by a
/// different route. Capture efficiency is a property of the panel and the reference,
/// so composition is left free to emerge from abundance x efficiency.
///
/// Monotone increasing in concentration, since every theta is, and it reaches 1 when
/// every probe saturates — so any requested efficiency below 1 is achievable.
fn capture_efficiency_at(
    sequences: &HashMap<String, String>,
    by_seq: &HitsBySequence,
    params: &RetentionParams,
    model: &ThermoModel,
    probe_conc: f64,
) -> f64 {
    let mut covered_capture = 0.0f64;
    let mut covered_len = 0.0f64;
    for (id, hits) in by_seq {
        let seq = match sequences.get(id) {
            Some(s) => s,
            None => continue,
        };
        let thetas = thetas_at(hits, probe_conc, model);
        // p_ns = 0, so p_retain == p_specific and the sweep's own `p_retain > 0` guard
        // drops every position no probe reaches.
        sweep_retention(
            seq.len(),
            &thetas,
            params.frag_len,
            0.0,
            params.min_overlap_frac,
            |from, to, _p_retain, p_spec| {
                let len = (to - from) as f64;
                covered_capture += len * p_spec;
                covered_len += len;
            },
        );
    }
    if covered_len > 0.0 {
        covered_capture / covered_len
    } else {
        0.0
    }
}

/// Back-solve the probe concentration that delivers a requested capture efficiency.
///
/// `--capture-fraction` is a preset over the physical parameter: users know roughly how
/// efficiently their protocol pulls down material it has probes for, but not the probe
/// concentration in the arbitrary units K absorbs. Efficiency is monotone increasing in
/// concentration, so a 1-D bisection on `ln [probe]` inverts it in a fixed number of
/// evaluations, with no iteration in the sampling loop.
///
/// What this does *not* do is target the output's composition. The realized on-target
/// and distractor fractions follow from abundance x efficiency, which is what lets a
/// dilution series titrate instead of sitting on a plateau.
pub fn solve_probe_concentration(
    sequences: &HashMap<String, String>,
    by_seq: &HitsBySequence,
    params: &RetentionParams,
    model: &ThermoModel,
    target_efficiency: f64,
) -> f64 {
    if target_efficiency <= 0.0 {
        return 0.0;
    }
    // Bisect in log space, where occupancy is a logistic and so well conditioned. The
    // bracket is wide enough for any realistic duplex — ln K for a perfect 120-mer at
    // 65 °C is about +135, so ln [probe] = -300 drives every theta to zero and +300
    // saturates every one — and stays inside f64's exponent range.
    let (mut lo, mut hi) = (-300.0f64, 300.0f64);
    let max_reachable = capture_efficiency_at(sequences, by_seq, params, model, hi.exp());
    if max_reachable < target_efficiency {
        log::warn!(
            "--capture-fraction {:.4} is above the {:.4} this panel can reach on this              reference: even with every probe saturated, some probe-covered positions              are reachable only by a fragment that cannot carry enough of the site              (--min-overlap-frac). Using the maximum.",
            target_efficiency,
            max_reachable
        );
        return hi.exp();
    }
    for _ in 0..60 {
        let mid = 0.5 * (lo + hi);
        if capture_efficiency_at(sequences, by_seq, params, model, mid.exp()) < target_efficiency {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    (0.5 * (lo + hi)).exp()
}

/// Write the per-reference retention report.
///
/// This is the model's own account of what it expects to happen, before any sampling
/// noise: for each reference, what fraction of its fragments survive the wash
/// (`efficiency`), how much of that is probe capture rather than carryover, and how
/// strongly its best-covered position is pulled.
///
/// It exists because the aggregate on-target fraction hides exactly what the Borrelia
/// LOD data showed matters: per-locus efficiency varies ~14x across loci in the same
/// panel at the same input, and is separable from abundance. `efficiency` is directly
/// comparable to those measured per-locus ratios, so a regression in relative capture
/// is visible without running the whole pipeline and counting reads.
pub fn write_retention_report(
    profiles: &[RetentionProfile],
    weights: &HashMap<String, f64>,
    output: &Path,
) -> Result<()> {
    let file = File::create(output)
        .with_context(|| format!("Cannot create retention report: {}", output.display()))?;
    let mut w = BufWriter::new(file);
    writeln!(
        w,
        "seq_id\tweight\tn_starts\tcovered_starts\tretention_mass\tspecific_mass\tefficiency\ton_target_fraction\tmax_p_retain\tmax_p_specific"
    )?;
    for p in profiles {
        let weight = weights.get(&p.seq_id).copied().unwrap_or(0.0);
        let efficiency = if p.n_starts > 0 {
            p.mass / p.n_starts as f64
        } else {
            0.0
        };
        let on_target = if p.mass > 0.0 {
            p.mass_specific / p.mass
        } else {
            0.0
        };
        let max_p_retain = p
            .segments
            .iter()
            .map(|s| s.p_retain)
            .fold(0.0f64, f64::max);
        let max_p_specific = p
            .segments
            .iter()
            .map(|s| s.p_specific)
            .fold(0.0f64, f64::max);
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{:.6}\t{:.6}\t{:.6e}\t{:.6}\t{:.6}\t{:.6}",
            p.seq_id,
            weight,
            p.n_starts,
            p.covered_starts,
            p.mass,
            p.mass_specific,
            efficiency,
            on_target,
            max_p_retain,
            max_p_specific
        )?;
    }
    w.flush()?;
    Ok(())
}

/// Draw `n_fragments` fragments from the retention profiles.
///
/// Each reference is drawn with probability proportional to `abundance × retention
/// mass`, then a start position proportional to `p_retain` within it. Every draw is a
/// keeper — there is no accept/reject loop, and no quota split between "capture" and
/// "background" arms. Capture and carryover fragments are therefore interleaved in
/// output order by construction, so a downstream tool reading a prefix of the FASTA
/// sees a representative sample rather than pure capture.
pub fn sample_retention_fragments(
    profiles: &[RetentionProfile],
    sequences: &HashMap<String, String>,
    weights: &HashMap<String, f64>,
    n_fragments: usize,
    len_params: FragmentLengthParams,
    seed: Option<u64>,
    counter_start: usize,
) -> Result<Vec<(String, String)>> {
    if profiles.is_empty() {
        bail!(
            "No reference sequence can yield a fragment: every weighted sequence is \
             either shorter than --fragment-length-min or has zero retention. Check \
             --fragment-length-min against the reference, and --background-retention."
        );
    }

    let seq_weights: Vec<f64> = profiles
        .iter()
        .map(|p| weights.get(&p.seq_id).copied().unwrap_or(0.0) * p.mass)
        .collect();
    if seq_weights.iter().all(|&w| w <= 0.0) {
        bail!("No sequences with positive weight and retention mass — cannot sample fragments.");
    }
    let dist = WeightedIndex::new(&seq_weights)?;

    let sd = (len_params.max as f64 - len_params.min as f64) / 6.0;
    let normal = Normal::new(len_params.mean, sd)?;

    let mut rng: StdRng = match seed {
        Some(s) => StdRng::seed_from_u64(s.wrapping_add(1)),
        None => StdRng::from_entropy(),
    };

    let mut results = Vec::with_capacity(n_fragments);
    let mut counter = counter_start;

    for _ in 0..n_fragments {
        let profile = &profiles[dist.sample(&mut rng)];
        let seq = match sequences.get(&profile.seq_id) {
            Some(s) => s,
            None => continue,
        };
        let start = match profile.sample_start(&mut rng) {
            Some(s) => s,
            None => continue,
        };

        let frag_len_f: f64 = normal.sample(&mut rng);
        let frag_len = (frag_len_f.round() as usize)
            .clamp(len_params.min, len_params.max)
            .min(seq.len());
        if frag_len == 0 {
            continue;
        }

        // The profile was built for one nominal length, so a longer draw can run off
        // the end. Shift it back rather than truncating, to keep the length
        // distribution intact.
        let start = start.min(seq.len() - frag_len);
        let fragment = &seq[start..start + frag_len];

        counter += 1;
        results.push((
            format!(
                ">{}_fragment_{} start={} length={}",
                profile.seq_id, counter, start, frag_len
            ),
            fragment.to_string(),
        ));
    }

    Ok(results)
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

    // ── Capture model: solver and statistics ─────────────────────────────────

    /// One reference, one probe hit, known geometry — small enough to verify the
    /// aggregate statistics against arithmetic done by hand.
    fn one_hit_fixture() -> (HashMap<String, String>, HashMap<String, f64>, HitsBySequence) {
        let mut sequences = HashMap::new();
        sequences.insert("t".to_string(), "A".repeat(1000));
        sequences.insert("d".to_string(), "C".repeat(1000));
        let mut weights = HashMap::new();
        weights.insert("t".to_string(), 1.0);
        weights.insert("d".to_string(), 1.0);
        let mut by_seq: HitsBySequence = Default::default();
        // A strongly favourable duplex (ΔG = −40 kcal/mol) at 400..520.
        by_seq.insert("t".to_string(), vec![(400, 520, -40.0)]);
        (sequences, weights, by_seq)
    }

    fn params(conc: f64) -> RetentionParams {
        RetentionParams {
            probe_concentration: conc,
            background_retention: 1e-3,
            min_overlap_frac: 0.5,
            frag_len: 100,
        }
    }

    /// Capture efficiency must rise with probe concentration and stay a probability.
    /// The solver bisects on this, so monotonicity is what makes it correct.
    #[test]
    fn capture_efficiency_is_monotone_in_concentration() {
        let (sequences, _w, by_seq) = one_hit_fixture();
        let m = model();
        let mut prev = -1.0;
        for exp in [-40.0, -30.0, -20.0, -10.0, 0.0, 10.0] {
            let e = capture_efficiency_at(&sequences, &by_seq, &params(1.0), &m, f64::exp(exp));
            assert!(
                (0.0..=1.0).contains(&e),
                "efficiency must be a probability, got {e} at ln[probe]={exp}"
            );
            assert!(e >= prev, "efficiency fell from {prev} to {e} as concentration rose");
            prev = e;
        }
        assert!(prev > 0.99, "efficiency should saturate at high concentration, got {prev}");
    }

    /// Capture efficiency describes the panel, not the sample: adding distractor must
    /// not change it. This is what stops the solver cancelling out composition, which
    /// is how the quota model pinned the realized distractor fraction.
    #[test]
    fn capture_efficiency_ignores_library_composition() {
        let (mut sequences, _w, by_seq) = one_hit_fixture();
        let m = model();
        let a = capture_efficiency_at(&sequences, &by_seq, &params(1e-6), &m, 1e-6);
        for i in 0..50 {
            sequences.insert(format!("extra_{i}"), "G".repeat(5000));
        }
        let b = capture_efficiency_at(&sequences, &by_seq, &params(1e-6), &m, 1e-6);
        assert_eq!(a, b, "efficiency changed when unprobed sequence was added");
    }

    /// The solver must actually hit the efficiency it was asked for.
    #[test]
    fn solver_reaches_the_requested_efficiency() {
        let (sequences, _w, by_seq) = one_hit_fixture();
        let m = model();
        for target in [0.05, 0.25, 0.5, 0.9] {
            let p = params(1.0);
            let conc = solve_probe_concentration(&sequences, &by_seq, &p, &m, target);
            let got = capture_efficiency_at(&sequences, &by_seq, &p, &m, conc);
            assert!(
                (got - target).abs() < 1e-6,
                "asked for efficiency {target}, solver gave {got} at [probe]={conc:.3e}"
            );
        }
    }

    /// Zero requested efficiency means no probes, not a degenerate concentration.
    #[test]
    fn solver_returns_zero_for_zero_efficiency() {
        let (sequences, _w, by_seq) = one_hit_fixture();
        assert_eq!(
            solve_probe_concentration(&sequences, &by_seq, &params(1.0), &model(), 0.0),
            0.0
        );
    }

    /// Enrichment is the output's on-target share over the input's. With one probe
    /// covering a known number of start positions, both are computable by hand.
    #[test]
    fn enrichment_compares_output_share_to_input_share() {
        let (sequences, weights, by_seq) = one_hit_fixture();
        let m = model();
        let p = params(1.0); // ΔG = −40 at 65 °C saturates, so θ ≈ 1
        let profiles = build_profiles(&sequences, &weights, &by_seq, &p, &m);
        let stats = retention_stats(&profiles, &weights);

        // Hit is 120 bp, min overlap 60, fragment 100 → starts in [400+60-100, 520-60]
        // = [360, 460], so 101 covered positions out of 901 per reference, over two
        // equally weighted references of 901 starts each.
        let covered: usize = profiles.iter().map(|q| q.covered_starts).sum();
        assert_eq!(covered, 101, "expected 101 probe-reachable starts, got {covered}");
        let input_share = 101.0 / (2.0 * 901.0);
        let expected = stats.on_target_fraction / input_share;
        assert!(
            (stats.enrichment - expected).abs() < 1e-9,
            "enrichment {} should equal output share / input share {}",
            stats.enrichment,
            expected
        );
        // At saturation the correct formula and the buggy one coincide, because
        // `mass_specific` approaches `covered_starts`. Below saturation they separate,
        // so that is where the regression is detectable: the earlier version divided by
        // `mass_specific / w_input`, which cancels the numerator and always yields
        // exactly `1 / yield` whatever the occupancy.
        let half = {
            let mut q = params(1.0);
            q.probe_concentration =
                solve_probe_concentration(&sequences, &by_seq, &q, &m, 0.5);
            q
        };
        let profiles = build_profiles(&sequences, &weights, &by_seq, &half, &m);
        let stats = retention_stats(&profiles, &weights);
        let reciprocal_yield = 1.0 / stats.yield_fraction;
        assert!(
            (stats.enrichment - reciprocal_yield).abs() > 1.0,
            "at 50% capture efficiency enrichment ({}) must not equal 1/yield ({})",
            stats.enrichment,
            reciprocal_yield
        );
    }

    /// Unprobed references must still contribute background mass — they are most of a
    /// realistic library, and dropping them would make the on-target fraction
    /// meaningless.
    #[test]
    fn unprobed_references_carry_background_mass() {
        let (sequences, weights, by_seq) = one_hit_fixture();
        let profiles = build_profiles(&sequences, &weights, &by_seq, &params(1.0), &model());
        let d = profiles.iter().find(|p| p.seq_id == "d").expect("distractor profile");
        assert_eq!(d.covered_starts, 0);
        assert!((d.mass - 901.0 * 1e-3).abs() < 1e-9, "mass {} != 901 * p_ns", d.mass);
    }

    /// A contig shorter than the nominal fragment length still belongs in the library.
    /// Dropping it would silently remove short distractor contigs from the background.
    #[test]
    fn short_references_are_not_dropped() {
        let mut sequences = HashMap::new();
        sequences.insert("short".to_string(), "A".repeat(40)); // < frag_len 100
        let mut weights = HashMap::new();
        weights.insert("short".to_string(), 1.0);
        let by_seq: HitsBySequence = Default::default();
        let profiles = build_profiles(&sequences, &weights, &by_seq, &params(1.0), &model());
        let p = profiles
            .iter()
            .find(|p| p.seq_id == "short")
            .expect("short contig dropped from the library");
        assert_eq!(p.n_starts, 1, "a 40 bp contig offers exactly one fragment");
        assert!((p.mass - 1e-3).abs() < 1e-12, "mass {} != p_ns", p.mass);
    }

    /// Zero-weight references are absent from the sample, so they must not appear.
    #[test]
    fn zero_weight_references_are_excluded() {
        let (sequences, mut weights, by_seq) = one_hit_fixture();
        weights.insert("d".to_string(), 0.0);
        let profiles = build_profiles(&sequences, &weights, &by_seq, &params(1.0), &model());
        assert!(profiles.iter().all(|p| p.seq_id != "d"));
    }

    /// Hits must be regrouped by reference, in sorted order, for the sweep to be
    /// deterministic under a fixed seed.
    #[test]
    fn hits_group_by_sequence_in_sorted_order() {
        let mut hits_by_probe: HashMap<String, Vec<ProbeHit>> = HashMap::new();
        hits_by_probe.insert(
            "p1".to_string(),
            vec![
                ProbeHit { seq_id: "b".into(), start: 50, end: 90, dg: -5.0 },
                ProbeHit { seq_id: "a".into(), start: 10, end: 40, dg: -6.0 },
            ],
        );
        hits_by_probe.insert(
            "p2".to_string(),
            vec![ProbeHit { seq_id: "a".into(), start: 5, end: 30, dg: -7.0 }],
        );
        let by_seq = group_hits_by_sequence(&hits_by_probe);
        assert_eq!(by_seq.len(), 2);
        assert_eq!(by_seq["a"], vec![(5, 30, -7.0), (10, 40, -6.0)]);
        assert_eq!(by_seq["b"], vec![(50, 90, -5.0)]);
    }

    // ── Retention profile ────────────────────────────────────────────────────

    /// With no probe hits every start position is retained at the non-specific rate.
    #[test]
    fn profile_without_hits_is_pure_background() {
        let p = build_retention_profile("s", 1000, &[], 100, 0.01, 0.5);
        let n_starts = (1000 - 100 + 1) as f64;
        assert_eq!(p.segments.len(), 1);
        assert!((p.segments[0].p_retain - 0.01).abs() < 1e-12);
        assert!((p.mass - 0.01 * n_starts).abs() < 1e-9, "mass {}", p.mass);
    }

    /// A saturated probe retains every fragment that carries enough of its site. With a
    /// minimum overlap of m the window of start positions is [h0+m-L, h1-m], which
    /// reduces to the any-overlap window [h0-L+1, h1-1] at m = 1.
    #[test]
    fn saturated_hit_gives_certain_capture_over_its_window() {
        let frag = 100usize;
        let (h0, h1) = (500usize, 620usize); // 120 bp hit
        let hit_len = h1 - h0;

        // m = 1: any overlap counts.
        let any = build_retention_profile("s", 2000, &[(h0, h1, 1.0)], frag, 0.0, 1.0 / hit_len as f64);
        let w_any: usize = any.segments.iter().map(|s| s.end - s.start).sum();
        assert_eq!(w_any, frag + hit_len - 1, "any-overlap window");

        // m = 60 (half the site): the window narrows by 2*(m-1).
        let half = build_retention_profile("s", 2000, &[(h0, h1, 1.0)], frag, 0.0, 0.5);
        let w_half: usize = half.segments.iter().map(|s| s.end - s.start).sum();
        assert_eq!(w_half, frag + hit_len - 1 - 2 * (60 - 1), "half-overlap window");
        assert!(w_half < w_any, "a stricter threshold must narrow the window");
        assert!((half.mass - w_half as f64).abs() < 1e-9);
    }

    /// A fragment shorter than the required overlap can never be captured by that probe,
    /// however well the probe matches.
    #[test]
    fn fragment_too_short_for_the_required_overlap_is_never_captured() {
        // 120 bp site, 50% threshold => needs 60 bp, but fragments are 40 bp.
        let p = build_retention_profile("s", 5000, &[(1000, 1120, 1.0)], 40, 0.0, 0.5);
        assert_eq!(p.mass, 0.0, "no capture possible");
    }

    /// Independent probes combine as 1 - prod(1 - theta): two overlapping 0.5 hits give
    /// 0.75 where both are active, not 0.5 or 1.0.
    #[test]
    fn overlapping_hits_combine_independently() {
        let p = build_retention_profile("s", 3000, &[(1000, 1120, 0.5), (1040, 1160, 0.5)], 50, 0.0, 0.25);
        let peak = p.segments.iter().map(|s| s.p_retain).fold(0.0f64, f64::max);
        assert!((peak - 0.75).abs() < 1e-9, "expected 0.75, got {peak}");
    }

    /// Background retention is a floor everywhere; a hit lifts it above that floor.
    #[test]
    fn background_is_a_floor_under_the_profile() {
        let p_ns = 0.02;
        let p = build_retention_profile("s", 3000, &[(1000, 1120, 0.6)], 50, p_ns, 0.25);
        let min = p.segments.iter().map(|s| s.p_retain).fold(f64::MAX, f64::min);
        let max = p.segments.iter().map(|s| s.p_retain).fold(0.0f64, f64::max);
        assert!((min - p_ns).abs() < 1e-12, "floor should be p_ns, got {min}");
        assert!((max - 0.608).abs() < 1e-9, "peak should be 0.608, got {max}");
    }

    /// Probe density raises retention mass - the effect the real LOD data shows.
    #[test]
    fn more_probes_means_more_retention_mass() {
        let sparse = build_retention_profile("s", 10_000, &[(5000, 5120, 0.9)], 200, 0.001, 0.5);
        let dense = build_retention_profile(
            "s", 10_000,
            &[(5000, 5120, 0.9), (5200, 5320, 0.9), (5400, 5520, 0.9)],
            200, 0.001, 0.5,
        );
        assert!(dense.mass > sparse.mass, "{} vs {}", dense.mass, sparse.mass);
    }

    /// Cross-probe affinity now matters, which the previous model could not express.
    #[test]
    fn weaker_probes_carry_less_mass() {
        let strong = build_retention_profile("s", 10_000, &[(5000, 5120, 0.9)], 200, 0.0, 0.5);
        let weak = build_retention_profile("s", 10_000, &[(5000, 5120, 0.1)], 200, 0.0, 0.5);
        assert!(strong.mass > weak.mass, "{} vs {}", strong.mass, weak.mass);
        assert!((strong.mass / weak.mass - 9.0).abs() < 1e-6, "ratio should track theta");
    }

    /// A reference shorter than the fragment cannot contribute.
    #[test]
    fn reference_shorter_than_fragment_has_no_mass() {
        let p = build_retention_profile("s", 50, &[(10, 40, 1.0)], 200, 0.5, 0.5);
        assert_eq!(p.mass, 0.0);
        assert!(p.segments.is_empty());
    }


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
