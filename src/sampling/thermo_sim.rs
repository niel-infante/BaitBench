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
use crate::alignment::probe_cache::Placement;
use crate::sampling::{FragmentLengthParams, FragmentLengthShape};
use crate::thermodynamics::{delta_g, langmuir_occupancy, min_window_delta_g, ThermoModel};

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
/// SAM 0x10: the alignment is to the reverse strand.
const FLAG_REVERSE: u16 = 0x10;

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
    threads: usize,
) -> Result<Vec<Placement>> {
    use crate::alignment::seed_extend::{find_local_hits, ProbeIndex, SeedExtendParams};

    // Only sequences that can actually contribute fragments are worth aligning to.
    let mut ref_ids: Vec<&String> = sequences
        .keys()
        .filter(|id| weights.get(*id).copied().unwrap_or(0.0) > 0.0)
        .collect();
    ref_ids.sort_unstable();
    // Borrowed, not copied: `sequences` already holds every base, and duplicating a
    // multi-gigabyte reference here doubled the peak for no benefit.
    let ref_seqs: Vec<&[u8]> = ref_ids.iter().map(|id| sequences[*id].as_bytes()).collect();

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

    let mut placements: Vec<Placement> = Vec::new();
    // Consumed, not borrowed: each hit owns a pair vector the size of the probe, and
    // with a cap of 1000 hits per probe a large panel holds well over a gigabyte of
    // them. Moving the vectors out instead of collecting new ones avoids holding two
    // copies of every alignment at the peak.
    for (pi, hits) in per_probe.into_iter().enumerate() {
        if hits.is_empty() {
            continue;
        }
        let name = probe_names[pi].to_string();
        for mut h in hits {
            let seq_id = ref_ids[h.seq_idx].to_string();
            if weights.get(&seq_id).copied().unwrap_or(0.0) <= 0.0 {
                continue;
            }
            // seed-extend reports raw (probe, reference) pairs; the SAM and blast paths
            // already complement. Normalise in place so a Placement always holds
            // exactly what delta_g consumes, whatever produced it.
            for pair in h.pairs.iter_mut() {
                pair.1 = wc_complement_base(pair.1);
            }
            placements.push(Placement {
                probe_name: name.clone(),
                seq_id,
                start: h.ref_start,
                end: h.ref_end,
                reverse: h.reverse,
                pairs: std::mem::take(&mut h.pairs),
            });
        }
    }

    log::info!(
        "Placed {} probe hits with seed-and-extend over {} reference sequences",
        placements.len(),
        ref_seqs.len()
    );
    Ok(placements)
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
) -> Result<Vec<Placement>> {
    let file = File::open(tsv_path)
        .with_context(|| format!("Cannot open BLAST output: {}", tsv_path.display()))?;
    let reader = BufReader::new(file);

    let mut placements: Vec<Placement> = Vec::new();
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
        let (lo, hi) = if sstart <= send {
            (sstart, send)
        } else {
            (send, sstart)
        };
        let start = lo.saturating_sub(1);
        let end = hi;

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

        placements.push(Placement {
            probe_name,
            seq_id: rname.to_string(),
            start,
            end,
            reverse: sstart > send,
            pairs,
        });
        n_parsed += 1;
    }

    log::info!(
        "Loaded {} probe hits ({} skipped) from {}",
        n_parsed,
        n_skipped,
        tsv_path.display()
    );
    Ok(placements)
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
) -> Result<Vec<Placement>> {
    let file = File::open(sam_path)
        .with_context(|| format!("Cannot open SAM file: {}", sam_path.display()))?;
    let reader = BufReader::new(file);

    let mut placements: Vec<Placement> = Vec::new();
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

        // SEQ is already in reference orientation, so pairs can be read straight off
        // the reference.
        let reference = match sequences.get(rname) {
            Some(r) => r.as_bytes(),
            None => {
                n_skipped += 1;
                continue;
            }
        };

        placements.push(Placement {
            probe_name,
            seq_id: rname.to_string(),
            start,
            end,
            reverse: flag & FLAG_REVERSE != 0,
            pairs: aligned_pairs(seq, reference, start, cigar),
        });

        n_parsed += 1;
    }

    log::info!(
        "Loaded {} probe hits ({} skipped) from {}",
        n_parsed,
        n_skipped,
        sam_path.display()
    );

    Ok(placements)
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

/// Turn placements into scored hits, grouped by probe.
///
/// Split from the aligners so that the expensive part (finding where probes land) can
/// be cached and reused, while the cheap part (free energy at this temperature and
/// salt) is redone every run. That is what lets one cached alignment serve a whole
/// temperature sweep — and it is why `Placement` stores base pairs rather than a ΔG.
///
/// Hit lists are sorted because the sampler indexes into them, and `HashMap` iteration
/// order varies between processes, which would otherwise leak into a seeded run.
/// `melting_window` selects the free energy that drives occupancy: 0 uses the whole
/// duplex, and any smaller value uses the least stable window of that many pairs. A
/// long duplex frays rather than melting as a unit, so the global sum makes occupancy a
/// step function in identity — at 58 °C a perfect 120-mer and a 9-mismatch one differ
/// by 22.8 kcal/mol, 10^15 in K. See `min_window_delta_g`.
pub fn score_placements(
    placements: &[Placement],
    mode: SimulateMode,
    model: &ThermoModel,
    melting_window: usize,
) -> HashMap<String, Vec<ProbeHit>> {
    let mut hits_by_probe: HashMap<String, Vec<ProbeHit>> = HashMap::new();
    for p in placements {
        let dg = match mode {
            // `pairs` already holds the base the probe binds, so this is the whole of
            // the per-run work.
            SimulateMode::Thermodynamic => {
                if melting_window == 0 {
                    delta_g(&p.pairs, model)
                } else {
                    min_window_delta_g(&p.pairs, melting_window, model)
                }
            }
            SimulateMode::Simple => SIMPLE_DG,
        };
        hits_by_probe
            .entry(p.probe_name.clone())
            .or_default()
            .push(ProbeHit {
                seq_id: p.seq_id.clone(),
                start: p.start,
                end: p.end,
                dg,
            });
    }
    for hits in hits_by_probe.values_mut() {
        hits.sort_unstable_by(|a, b| {
            a.seq_id
                .cmp(&b.seq_id)
                .then(a.start.cmp(&b.start))
                .then(a.end.cmp(&b.end))
        });
    }
    hits_by_probe
}

// ─── Retention model ──────────────────────────────────────────────────────────

/// Non-specific carryover: the probability that a fragment overlapping no probe is
/// retained anyway. Washing is imperfect, and in a library that is 99.99% off-target
/// this term, not capture, supplies most of the background reads.
pub const DEFAULT_BACKGROUND_RETENTION: f64 = 1e-3;

/// Fraction of a probe's target site a fragment must carry before that probe can pull
/// it. See `build_retention_profile` for why a threshold is used here.
pub const DEFAULT_MIN_OVERLAP_FRAC: f64 = 0.5;

/// Probe concentration used when it is not back-solved, in the arbitrary units K
/// absorbs. `e^-2`, calibrated against the Borrelia LOD per-locus efficiencies with a
/// 10 bp melting window and cross-checked on the D6331 mock community.
///
/// Fixed rather than solved per run **on purpose**. Solving it to hit a requested
/// `--capture-fraction` makes the model a rank-threshold device: the scalar re-centres
/// around any change to the affinity basis, so no refinement of ΔG can be observed.
/// Measured — a 10 bp window changed nothing at all while the concentration was
/// solved, and moved per-locus log-RMSE from 2.71 to 1.02 once it was fixed.
///
/// Provisional. It is one constant fitted to eight loci from one organism; it improved
/// the D6331 mock community too (CCC 0.791 -> 0.803), but it has not been fitted
/// jointly across datasets.
pub const DEFAULT_PROBE_CONCENTRATION: f64 = 0.135_34;

/// Pairs in the cooperative melting unit, used when `--melting-window` is not given.
///
/// 10 bp is the canonical DNA nucleation length, and is independently where the fit to
/// measured per-locus capture is best. Using the whole duplex instead makes occupancy a
/// step function in identity: at 58 °C a perfect 120-mer and a 9-mismatch one differ by
/// 22.8 kcal/mol, 10^15 in K, so every probe either saturates or vanishes.
pub const DEFAULT_MELTING_WINDOW: usize = 10;

/// Physical parameters of the capture reaction.
#[derive(Debug, Clone, Copy)]
pub struct RetentionParams {
    /// Probe concentration, in the arbitrary units absorbed into K. Set directly, or
    /// Fixed; see `DEFAULT_PROBE_CONCENTRATION` for why it is not back-solved.
    pub probe_concentration: f64,
    /// Non-specific carryover probability per fragment.
    pub background_retention: f64,
    /// Fraction of probe-bound fragments that survive wash and elution into the
    /// library — `--capture-fraction`. A uniform multiplier on specific capture, so it
    /// sets how good the run is without touching relative capture between loci.
    pub capture_efficiency: f64,
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
    capture_efficiency: f64,
    mut emit: F,
) -> usize
where
    F: FnMut(usize, usize, f64, f64),
{
    if frag_len == 0 || seq_len == 0 {
        return 0;
    }
    // Every base is a potential fragment start, including the last `frag_len - 1`,
    // whose fragments run off the end and are truncated when extracted.
    //
    // This was `seq_len - frag_len + 1`, which counts *start positions* while the input
    // weights encode DNA *mass* (proportional to length). The two agree to 0.1% when
    // `seq_len >> frag_len` and collapse when they are comparable, so a draft assembly
    // was silently down-weighted against a closed genome in the same run. Measured on
    // the D6331 panel: Candida albicans is 13,062 contigs over 12.58 Mb against a
    // 5.8 kb fragment length, and came out at 0.244x its intended input share — the
    // available-starts fraction, to two significant figures.
    //
    // Truncation costs nothing in correctness: a hit is credited to a start only when
    // the start is at or left of `h_end - min_overlap`, and `h_end <= seq_len`, so the
    // part of the fragment that is cut off never contained any of the hit. An end
    // fragment shorter than `min_overlap` is credited nothing, which is right.
    let n_starts = seq_len;

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
        // `p_bound` is the chance at least one probe holds this fragment; multiplying
        // by `capture_efficiency` gives the chance it also survives wash and elution.
        // A uniform multiplier by construction, so it cannot reshape relative capture
        // between loci — which is the whole reason it is applied here rather than by
        // moving the probe concentration.
        let p_bound = if n_sat > 0 { 1.0 } else { 1.0 - sum_log.exp() };
        let p_spec = capture_efficiency * p_bound;
        let p_retain = 1.0 - (1.0 - p_ns) * (1.0 - p_spec);
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
    capture_efficiency: f64,
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
        capture_efficiency,
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
    let on_target_fraction = if w_mass > 0.0 {
        w_mass_spec / w_mass
    } else {
        0.0
    };
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
            // The nominal length is used as-is, even on a contig shorter than it:
            // `sweep_retention` now sweeps every base and truncates at the end, so a
            // short contig yields short fragments with mass proportional to its length.
            // Clamping `frag_len` to `seq_len` here used to give such a contig exactly
            // one start regardless of how long it was.
            build_retention_profile(
                id,
                seq_len,
                &thetas,
                params.frag_len.max(1),
                params.background_retention,
                params.min_overlap_frac,
                params.capture_efficiency,
            )
        })
        .filter(|p| p.mass > 0.0)
        .collect()
}

/// Fraction of weighted input mass sitting on contigs shorter than `frag_len`.
///
/// Those contigs are sampled correctly — mass is proportional to length and the
/// fragments are truncated at the contig end — but every fragment they yield is shorter
/// than requested, which shifts the library's length distribution. A draft assembly can
/// put most of an organism's sequence there without anything in the output saying so:
/// the D6331 panel's *Candida albicans* reference is 13,062 contigs over 12.58 Mb
/// against a 5.8 kb fragment length, which is 74.6% of its bases.
pub fn short_contig_mass_fraction(
    sequences: &HashMap<String, String>,
    weights: &HashMap<String, f64>,
    frag_len: usize,
) -> (f64, usize, usize) {
    let mut short = 0.0f64;
    let mut total = 0.0f64;
    let mut n_short = 0usize;
    let mut n_total = 0usize;
    for (id, seq) in sequences {
        let w = weights.get(id).copied().unwrap_or(0.0);
        if w <= 0.0 {
            continue;
        }
        let m = w * seq.len() as f64;
        total += m;
        n_total += 1;
        if seq.len() < frag_len {
            short += m;
            n_short += 1;
        }
    }
    let frac = if total <= 0.0 { 0.0 } else { short / total };
    (frac, n_short, n_total)
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
        let max_p_retain = p.segments.iter().map(|s| s.p_retain).fold(0.0f64, f64::max);
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

/// Draws fragment lengths from a truncated Normal or log-normal.
///
/// Truncation is by **rejection**, not clamping. Clamping to `[min, max]` puts a point
/// mass on each bound — with the historical `sd = (max - min) / 6` roughly 0.3% of
/// fragments landed exactly on a bound, and at a user-supplied larger `sd` it becomes
/// a substantial spike at both ends of a distribution that should have none.
///
/// The retry budget keeps a bad parameter set from hanging: if `[min, max]` sits far
/// out in the tail, the draw falls back to clamping rather than looping. That is the
/// same wrong answer as before, but only for parameters that are already unusable.
struct FragmentLengthSampler {
    normal: Normal<f64>,
    params: FragmentLengthParams,
}

/// Rejection attempts before falling back to a clamp.
const LENGTH_DRAW_RETRIES: usize = 100;

impl FragmentLengthSampler {
    fn new(params: FragmentLengthParams) -> Result<Self> {
        let sd = params.sd_or_default();
        if sd <= 0.0 {
            bail!(
                "Fragment length standard deviation must be positive, got {}. \
                 Set --fragment-length-sd, or widen --fragment-length-min/--max.",
                sd
            );
        }
        let normal = match params.shape {
            FragmentLengthShape::Normal => Normal::new(params.mean, sd)?,
            // Match the requested mean and sd on the natural scale, rather than
            // treating them as the parameters of the underlying Normal — a user asking
            // for a 5.8 kb mean means 5.8 kb, not exp(5.8 kb).
            FragmentLengthShape::LogNormal => {
                if params.mean <= 0.0 {
                    bail!(
                        "A log-normal fragment length needs a positive mean, got {}.",
                        params.mean
                    );
                }
                let var = sd * sd;
                let mu =
                    (params.mean * params.mean / (var + params.mean * params.mean).sqrt()).ln();
                let sigma = (1.0 + var / (params.mean * params.mean)).ln().sqrt();
                Normal::new(mu, sigma)?
            }
        };
        Ok(Self { normal, params })
    }

    fn draw_raw<R: Rng>(&self, rng: &mut R) -> f64 {
        let x = self.normal.sample(rng);
        match self.params.shape {
            FragmentLengthShape::Normal => x,
            FragmentLengthShape::LogNormal => x.exp(),
        }
    }

    fn sample<R: Rng>(&self, rng: &mut R) -> usize {
        for _ in 0..LENGTH_DRAW_RETRIES {
            let v = self.draw_raw(rng).round();
            if v >= self.params.min as f64 && v <= self.params.max as f64 {
                return v as usize;
            }
        }
        (self.draw_raw(rng).round().max(0.0) as usize).clamp(self.params.min, self.params.max)
    }
}

/// Library amplification: distinct molecules, then reads drawn from them.
///
/// **Duplicate rate is not a parameter — it is an outcome of depth over complexity.**
/// Measured on the Borrelia LOD dilution series (reads mapped to the B. burgdorferi
/// chromosome, duplicates by identical start position and strand):
///
/// | condition | reads | distinct | duplicate | complexity M | N/M |
/// |---|---|---|---|---|---|
/// | s1E2 | 234 | 199 | 15.0% | 702 | 0.33 |
/// | s1E3 | 1,527 | 1,202 | 21.3% | 3,058 | 0.50 |
/// | s1E4 | 26,467 | 9,854 | 62.8% | 10,779 | 2.46 |
/// | s1E5 | 113,239 | 25,304 | 77.7% | 25,612 | 4.42 |
/// | s1E6 | 870,233 | 73,381 | 91.6% | 73,382 | 11.86 |
/// | Pos1 | 739,921 | 79,881 | 89.2% | 79,889 | 9.26 |
///
/// Complexity rises with template input, and the duplicate rate follows `N/M` along the
/// Lander–Waterman curve `1 - (1 - e^(-N/M)) / (N/M)`. So the knob is **how many
/// molecules survive capture**, and duplication falls out. Exposing a duplicate rate
/// directly would invert cause and effect and could not produce jackpotting.
///
/// Note this also corrects the plan's Stage 6 gate, which expected 0.9–7% duplication
/// rising at *low* input. Both halves are wrong for these data: the real rates are
/// 15–92%, and they rise with *high* input, because depth climbs faster than
/// complexity. The 0.9–7% figures came from FastQC on a different experiment entirely.
#[derive(Clone, Copy, Debug, Default)]
pub struct AmplificationParams {
    /// Distinct molecules surviving capture. `None` means unlimited — every read is its
    /// own molecule and there is no duplication, which is the behaviour before Stage 6.
    pub library_molecules: Option<usize>,
    /// Length scale (bp) of an exponential amplification bias toward short molecules:
    /// relative efficiency `exp(-length / scale)`. `None` amplifies every molecule
    /// equally.
    ///
    /// Motivated by a measured effect whose mechanism is **not** established: across all
    /// 18 enriched D6331 runs the library is ~5x shorter than the non-enriched arm
    /// (0.9–1.9 kb against 5.5–6.1 kb). Post-capture PCR favouring short templates is
    /// the usual explanation, but a size-selection step would look identical in this
    /// data. Off by default for that reason.
    pub length_scale: Option<f64>,
}

/// What amplification did, for reporting.
#[derive(Clone, Copy, Debug, Default)]
pub struct AmplificationStats {
    pub reads: usize,
    pub distinct_molecules: usize,
    pub duplicate_fraction: f64,
    /// Reads drawn from the single most-amplified molecule — the jackpotting signal.
    pub max_molecule_reads: usize,
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
    amp: AmplificationParams,
) -> Result<(Vec<(String, String)>, AmplificationStats)> {
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

    let sampler = FragmentLengthSampler::new(len_params)?;

    let mut rng: StdRng = match seed {
        Some(s) => StdRng::seed_from_u64(s.wrapping_add(1)),
        None => StdRng::from_entropy(),
    };

    // How many distinct molecules the library holds. Without a complexity limit this is
    // `n_fragments`, so every read is its own molecule and nothing duplicates.
    //
    // `M > N` is not clamped away: a library can easily hold more molecules than you
    // sequence reads from, and that is the shallow-sequencing regime where duplication
    // is *low*. Clamping to `n_fragments` forced N/M = 1 and a 36.8% floor, against a
    // measured 15.0% at s1E2 (N/M = 0.33) and 21.3% at s1E3 (N/M = 0.50).
    let library_size = amp.library_molecules.map(|m| m.max(1));
    let n_molecules = library_size.unwrap_or(n_fragments);

    // When the library is larger than the read count, materialising every molecule is
    // wasted work — most are never sequenced. Under uniform amplification the reads'
    // molecule choices can be drawn first, and only the distinct ones built. That is
    // O(n_fragments) regardless of how large the library is.
    let uniform_amp = amp.length_scale.is_none_or(|s| s <= 0.0);
    let sparse = library_size.is_some_and(|m| m > n_fragments) && uniform_amp;

    // Phase 1: draw the distinct molecules that survived capture.
    struct Molecule {
        seq_id_idx: usize,
        start: usize,
        len: usize,
    }
    // Under the sparse path, pick which library slots the reads come from first; the
    // number of distinct slots is what has to be built.
    let sparse_counts: Option<Vec<usize>> = if sparse {
        let m = n_molecules;
        let mut hits: HashMap<usize, usize> = HashMap::new();
        for _ in 0..n_fragments {
            *hits.entry(rng.gen_range(0..m)).or_insert(0) += 1;
        }
        let mut counts: Vec<usize> = hits.into_values().collect();
        counts.sort_unstable();
        Some(counts)
    } else {
        None
    };
    let to_build = sparse_counts
        .as_ref()
        .map(|c| c.len())
        .unwrap_or(n_molecules);

    let mut molecules: Vec<Molecule> = Vec::with_capacity(to_build);
    for _ in 0..to_build {
        let idx = dist.sample(&mut rng);
        let profile = &profiles[idx];
        let seq = match sequences.get(&profile.seq_id) {
            Some(s) => s,
            None => continue,
        };
        let start = match profile.sample_start(&mut rng) {
            Some(s) => s,
            None => continue,
        };
        let drawn = sampler.sample(&mut rng);
        // Truncate at the contig end rather than shifting the start back. Shifting
        // moved probability mass off the last `frag_len` bases and piled it onto the
        // final valid start, distorting the positional distribution that the retention
        // profile had just computed and biasing coverage inward at every contig end.
        // Truncation is also what happens physically: the molecule at the end of a
        // contig is simply shorter.
        let end = (start + drawn).min(seq.len());
        let len = end - start;
        if len == 0 {
            continue;
        }
        molecules.push(Molecule {
            seq_id_idx: idx,
            start,
            len,
        });
    }
    if molecules.is_empty() {
        bail!(
            "No molecules could be drawn — every sampled reference yielded a zero-length fragment."
        );
    }

    // Phase 2: sequence the library. With no complexity limit each molecule is read
    // once, which reproduces the pre-Stage-6 output exactly. Otherwise reads are drawn
    // with replacement, so duplicates and jackpotting emerge rather than being imposed.
    let reads_per_molecule: Vec<usize> = if amp.library_molecules.is_none() {
        vec![1; molecules.len()]
    } else if let Some(mut counts) = sparse_counts {
        // One count per built molecule. A molecule that failed to materialise drops its
        // reads, same as the dense path.
        counts.truncate(molecules.len());
        counts.resize(molecules.len(), 0);
        counts
    } else {
        let amp_weights: Vec<f64> = match amp.length_scale {
            Some(scale) if scale > 0.0 => molecules
                .iter()
                .map(|m| (-(m.len as f64) / scale).exp())
                .collect(),
            _ => vec![1.0; molecules.len()],
        };
        let picker = WeightedIndex::new(&amp_weights)?;
        let mut counts = vec![0usize; molecules.len()];
        for _ in 0..n_fragments {
            counts[picker.sample(&mut rng)] += 1;
        }
        counts
    };

    let tag_molecules = amp.library_molecules.is_some();
    let mut results = Vec::with_capacity(n_fragments);
    let mut counter = counter_start;
    let mut distinct = 0usize;
    let mut max_molecule_reads = 0usize;

    for (mol_idx, (m, &n_reads)) in molecules.iter().zip(reads_per_molecule.iter()).enumerate() {
        if n_reads == 0 {
            continue;
        }
        distinct += 1;
        max_molecule_reads = max_molecule_reads.max(n_reads);
        let profile = &profiles[m.seq_id_idx];
        let seq = match sequences.get(&profile.seq_id) {
            Some(s) => s,
            None => continue,
        };
        let fragment = &seq[m.start..m.start + m.len];
        for _ in 0..n_reads {
            counter += 1;
            // `molecule=` identifies PCR duplicates: reads sharing it are copies of one
            // captured molecule. The `{seq_id}_fragment_{n}` id stays unique and keeps
            // its contract, so source extraction downstream is unaffected.
            //
            // Emitted only when amplification is on, so that without
            // `--library-molecules` the header is byte-identical to what every previous
            // version wrote. `fragments.fa` is an output users and gates compare.
            let header = if tag_molecules {
                format!(
                    ">{}_fragment_{} start={} length={} molecule={}",
                    profile.seq_id, counter, m.start, m.len, mol_idx
                )
            } else {
                format!(
                    ">{}_fragment_{} start={} length={}",
                    profile.seq_id, counter, m.start, m.len
                )
            };
            results.push((header, fragment.to_string()));
        }
    }

    let reads = results.len();
    let stats = AmplificationStats {
        reads,
        distinct_molecules: distinct,
        duplicate_fraction: if reads > 0 {
            1.0 - distinct as f64 / reads as f64
        } else {
            0.0
        },
        max_molecule_reads,
    };

    Ok((results, stats))
}

#[cfg(test)]
mod stage6_tests {
    use super::*;
    use crate::sampling::{FragmentLengthParams, FragmentLengthShape};

    fn fixture() -> (
        Vec<RetentionProfile>,
        HashMap<String, String>,
        HashMap<String, f64>,
    ) {
        let mut seqs = HashMap::new();
        let mut w = HashMap::new();
        seqs.insert("a".to_string(), "ACGT".repeat(50_000)); // 200 kb
        w.insert("a".to_string(), 1.0);
        let profiles = vec![build_retention_profile(
            "a",
            200_000,
            &[],
            500,
            1e-3,
            0.5,
            0.7,
        )];
        (profiles, seqs, w)
    }

    fn lp() -> FragmentLengthParams {
        FragmentLengthParams {
            mean: 500.0,
            min: 100,
            max: 1_000,
            sd: Some(150.0),
            shape: FragmentLengthShape::Normal,
        }
    }

    fn draw(n: usize, amp: AmplificationParams) -> (Vec<(String, String)>, AmplificationStats) {
        let (p, s, w) = fixture();
        sample_retention_fragments(&p, &s, &w, n, lp(), Some(42), 0, amp).unwrap()
    }

    /// Without a complexity limit nothing duplicates — this is the pre-Stage-6 path and
    /// it must be untouched.
    #[test]
    fn unlimited_complexity_produces_no_duplicates() {
        let (frags, stats) = draw(2_000, AmplificationParams::default());
        assert_eq!(frags.len(), 2_000);
        assert_eq!(stats.distinct_molecules, 2_000);
        assert!(stats.duplicate_fraction.abs() < 1e-12);
        assert!(
            !frags[0].0.contains("molecule="),
            "no complexity limit should not tag molecules: {}",
            frags[0].0
        );
    }

    /// Duplicate rate must follow Lander-Waterman, `1 - (1 - e^-x) / x` for x = N/M,
    /// rather than being whatever the code happens to produce.
    #[test]
    fn duplicate_rate_follows_lander_waterman() {
        for (n, m) in [(2_000usize, 4_000usize), (2_000, 1_000), (4_000, 500)] {
            let (_, stats) = draw(
                n,
                AmplificationParams {
                    library_molecules: Some(m),
                    length_scale: None,
                },
            );
            // No clamping of M to N: a library larger than the read count is the
            // shallow-sequencing regime, and it is where duplication is lowest.
            let x = n as f64 / m as f64;
            let expected = 1.0 - (1.0 - (-x).exp()) / x;
            assert!(
                (stats.duplicate_fraction - expected).abs() < 0.03,
                "N={n} M={m}: duplicate {:.4} vs Lander-Waterman {:.4}",
                stats.duplicate_fraction,
                expected
            );
        }
    }

    /// Deeper sequencing of the same library duplicates more. This is the direction the
    /// real LOD series shows, and the opposite of the plan's stated expectation.
    #[test]
    fn duplication_rises_with_depth_at_fixed_complexity() {
        let rate = |n: usize| {
            draw(
                n,
                AmplificationParams {
                    library_molecules: Some(500),
                    length_scale: None,
                },
            )
            .1
            .duplicate_fraction
        };
        let (shallow, deep) = (rate(600), rate(5_000));
        assert!(
            deep > shallow + 0.3,
            "duplication should climb with depth: {shallow:.3} -> {deep:.3}"
        );
    }

    /// Reads carry a molecule tag, and the number of distinct tags equals the number of
    /// molecules actually sequenced.
    #[test]
    fn molecule_tags_identify_duplicates() {
        let (frags, stats) = draw(
            1_000,
            AmplificationParams {
                library_molecules: Some(200),
                length_scale: None,
            },
        );
        let tags: std::collections::HashSet<&str> = frags
            .iter()
            .filter_map(|(h, _)| h.split("molecule=").nth(1))
            .collect();
        assert_eq!(tags.len(), stats.distinct_molecules);
        assert!(stats.max_molecule_reads > 1, "expected some jackpotting");
    }

    /// Reads sharing a molecule tag must be byte-identical — they are copies of one
    /// captured molecule, not independent draws from the same region.
    #[test]
    fn duplicate_reads_are_identical_sequence() {
        let (frags, _) = draw(
            800,
            AmplificationParams {
                library_molecules: Some(50),
                length_scale: None,
            },
        );
        let mut by_tag: HashMap<&str, &str> = HashMap::new();
        for (h, seq) in &frags {
            let tag = h.split("molecule=").nth(1).unwrap();
            if let Some(prev) = by_tag.get(tag) {
                assert_eq!(*prev, seq.as_str(), "duplicates of molecule {tag} differ");
            } else {
                by_tag.insert(tag, seq);
            }
        }
        assert!(by_tag.len() > 1);
    }

    /// The length bias shortens the sequenced library, which is what it exists to model.
    #[test]
    fn amplification_length_bias_shortens_the_library() {
        let mean_len = |scale: Option<f64>| {
            let (frags, _) = draw(
                4_000,
                AmplificationParams {
                    library_molecules: Some(2_000),
                    length_scale: scale,
                },
            );
            frags.iter().map(|(_, s)| s.len() as f64).sum::<f64>() / frags.len() as f64
        };
        let unbiased = mean_len(None);
        let biased = mean_len(Some(300.0));
        assert!(
            biased < unbiased * 0.95,
            "length bias should shorten the library: {unbiased:.0} -> {biased:.0}"
        );
    }
}

#[cfg(test)]
mod stage5_tests {
    use super::*;
    use crate::sampling::{FragmentLengthParams, FragmentLengthShape};
    use rand::SeedableRng;

    fn params(
        mean: f64,
        min: usize,
        max: usize,
        sd: Option<f64>,
        shape: FragmentLengthShape,
    ) -> FragmentLengthParams {
        FragmentLengthParams {
            mean,
            min,
            max,
            sd,
            shape,
        }
    }

    /// Retention mass must track DNA mass, not the number of valid start positions.
    ///
    /// One 12 kb contig against twelve 1 kb contigs of the same total length: with no
    /// probes both are pure carryover, so a library drawn from them should split 50/50.
    /// Counting starts (`L - f + 1`) gave the fragmented side 12 x (1000 - 5800 + 1)
    /// -> clamped to 1 start each, and it all but vanished.
    #[test]
    fn mass_is_proportional_to_length_not_to_start_count() {
        let frag_len = 5_800;
        let one_big = build_retention_profile("big", 12_000, &[], frag_len, 1e-3, 0.5, 0.7);
        let many_small: f64 = (0..12)
            .map(|i| {
                build_retention_profile(&format!("s{i}"), 1_000, &[], frag_len, 1e-3, 0.5, 0.7).mass
            })
            .sum();
        let ratio = many_small / one_big.mass;
        assert!(
            (ratio - 1.0).abs() < 1e-9,
            "12 x 1kb should carry the same mass as 1 x 12kb, got ratio {ratio}"
        );
    }

    /// A contig shorter than the fragment length still contributes mass in proportion
    /// to its length, rather than a single start.
    #[test]
    fn short_contig_mass_scales_with_its_length() {
        let frag_len = 5_000;
        let a = build_retention_profile("a", 500, &[], frag_len, 1e-3, 0.5, 0.7).mass;
        let b = build_retention_profile("b", 1_000, &[], frag_len, 1e-3, 0.5, 0.7).mass;
        assert!((b / a - 2.0).abs() < 1e-9, "expected 2x, got {}", b / a);
    }

    /// Lengths are truncated by redrawing, so no fragment lands exactly on a bound more
    /// often than the density there warrants. Clamping produced a visible spike.
    #[test]
    fn lengths_are_rejected_not_clamped() {
        let p = params(175.0, 150, 200, Some(40.0), FragmentLengthShape::Normal);
        let sampler = FragmentLengthSampler::new(p).unwrap();
        let mut rng = StdRng::seed_from_u64(7);
        let draws: Vec<usize> = (0..20_000).map(|_| sampler.sample(&mut rng)).collect();
        assert!(
            draws.iter().all(|&d| (150..=200).contains(&d)),
            "out of bounds draw"
        );
        let at_min = draws.iter().filter(|&&d| d == 150).count();
        let at_max = draws.iter().filter(|&&d| d == 200).count();
        // With sd=40 over a 50 bp window, clamping would pile >25% on each bound.
        assert!(
            at_min < draws.len() / 20,
            "spike at min: {at_min}/{}",
            draws.len()
        );
        assert!(
            at_max < draws.len() / 20,
            "spike at max: {at_max}/{}",
            draws.len()
        );
    }

    /// `--fragment-length-sd` actually changes the spread; the old hardwired
    /// `(max - min) / 6` ignored any such request.
    #[test]
    fn sd_is_honoured() {
        let mut rng = StdRng::seed_from_u64(11);
        let mut spread = |sd: f64| {
            let s = FragmentLengthSampler::new(params(
                5_800.0,
                1_000,
                20_000,
                Some(sd),
                FragmentLengthShape::Normal,
            ))
            .unwrap();
            let d: Vec<f64> = (0..5_000).map(|_| s.sample(&mut rng) as f64).collect();
            let m = d.iter().sum::<f64>() / d.len() as f64;
            (d.iter().map(|x| (x - m) * (x - m)).sum::<f64>() / d.len() as f64).sqrt()
        };
        let narrow = spread(200.0);
        let wide = spread(2_000.0);
        assert!(
            wide > narrow * 4.0,
            "sd ignored: narrow {narrow}, wide {wide}"
        );
    }

    /// The log-normal is parameterised so the requested mean is the mean of the
    /// *lengths*, not of their logarithms, and it is right-skewed.
    #[test]
    fn lognormal_matches_the_requested_mean_and_is_skewed() {
        let s = FragmentLengthSampler::new(params(
            5_800.0,
            1,
            1_000_000,
            Some(2_000.0),
            FragmentLengthShape::LogNormal,
        ))
        .unwrap();
        let mut rng = StdRng::seed_from_u64(3);
        let d: Vec<f64> = (0..50_000).map(|_| s.sample(&mut rng) as f64).collect();
        let mean = d.iter().sum::<f64>() / d.len() as f64;
        assert!(
            (mean / 5_800.0 - 1.0).abs() < 0.05,
            "mean {mean} off target"
        );
        let mut sorted = d.clone();
        sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let median = sorted[sorted.len() / 2];
        assert!(
            median < mean,
            "log-normal should be right-skewed: median {median} mean {mean}"
        );
    }

    /// The short-contig warning fires on a fragmented assembly and stays quiet on a
    /// contiguous one.
    #[test]
    fn short_contig_fraction_is_reported() {
        let mut seqs = HashMap::new();
        let mut w = HashMap::new();
        seqs.insert("closed".to_string(), "A".repeat(10_000));
        w.insert("closed".to_string(), 1.0);
        let (f, n, total) = short_contig_mass_fraction(&seqs, &w, 5_800);
        assert_eq!(n, 0);
        assert_eq!(total, 1);
        assert!(
            f < 1e-12,
            "contiguous assembly should report nothing, got {f}"
        );

        for i in 0..10 {
            seqs.insert(format!("draft{i}"), "A".repeat(1_000));
            w.insert(format!("draft{i}"), 1.0);
        }
        let (f, n, total) = short_contig_mass_fraction(&seqs, &w, 5_800);
        assert_eq!(n, 10);
        assert_eq!(total, 11);
        // 10 x 1000 short against 10000 long = half the mass.
        assert!((f - 0.5).abs() < 1e-9, "expected 0.5, got {f}");
    }
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
    fn one_hit_fixture() -> (
        HashMap<String, String>,
        HashMap<String, f64>,
        HitsBySequence,
    ) {
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
            capture_efficiency: 1.0,
            frag_len: 100,
        }
    }

    /// The whole point of applying `--capture-fraction` as a multiplier: it changes
    /// how much is captured without changing *what* is captured relative to what else.
    /// Solving the probe concentration to hit a target could not promise this, because
    /// occupancy saturates and moving the concentration reshapes the pattern.
    #[test]
    fn capture_efficiency_scales_without_reshaping() {
        let (sequences, weights, by_seq) = one_hit_fixture();
        let m = model();
        let spec = |f: f64| {
            let mut p = params(1.0);
            p.capture_efficiency = f;
            build_profiles(&sequences, &weights, &by_seq, &p, &m)
                .iter()
                .find(|q| q.seq_id == "t")
                .map(|q| q.mass_specific)
                .unwrap_or(0.0)
        };
        let (a, b) = (spec(1.0), spec(0.1));
        assert!(a > 0.0 && b > 0.0);
        assert!(
            ((b / a) - 0.1).abs() < 1e-9,
            "specific mass must scale exactly with capture efficiency: {b}/{a}"
        );
    }

    /// At zero efficiency only non-specific carryover remains, so the library is pure
    /// background rather than empty — a bad run still sequences something.
    #[test]
    fn zero_capture_efficiency_leaves_the_carryover_floor() {
        let (sequences, weights, by_seq) = one_hit_fixture();
        let mut p = params(1.0);
        p.capture_efficiency = 0.0;
        let profiles = build_profiles(&sequences, &weights, &by_seq, &p, &model());
        let t = profiles
            .iter()
            .find(|q| q.seq_id == "t")
            .expect("target profile");
        assert_eq!(t.mass_specific, 0.0, "nothing is specifically captured");
        assert!(
            (t.mass - 1000.0 * 1e-3).abs() < 1e-9,
            "every position should retain at exactly p_ns, got {}",
            t.mass
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
        assert_eq!(
            covered, 101,
            "expected 101 probe-reachable starts, got {covered}"
        );
        let input_share = 101.0 / (2.0 * 1000.0);
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
            q.capture_efficiency = 0.5;
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
        let d = profiles
            .iter()
            .find(|p| p.seq_id == "d")
            .expect("distractor profile");
        assert_eq!(d.covered_starts, 0);
        assert!(
            (d.mass - 1000.0 * 1e-3).abs() < 1e-9,
            "mass {} != 1000 * p_ns",
            d.mass
        );
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
        // Every base is a start, so a short contig carries mass in proportion to its
        // length rather than contributing a single fragment. Counting starts gave it
        // exactly 1, which is what made fragmented assemblies vanish against closed
        // genomes in the same run.
        assert_eq!(p.n_starts, 40, "every base of a 40 bp contig is a start");
        assert!(
            (p.mass - 40.0 * 1e-3).abs() < 1e-12,
            "mass {} != 40 * p_ns",
            p.mass
        );
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
                ProbeHit {
                    seq_id: "b".into(),
                    start: 50,
                    end: 90,
                    dg: -5.0,
                },
                ProbeHit {
                    seq_id: "a".into(),
                    start: 10,
                    end: 40,
                    dg: -6.0,
                },
            ],
        );
        hits_by_probe.insert(
            "p2".to_string(),
            vec![ProbeHit {
                seq_id: "a".into(),
                start: 5,
                end: 30,
                dg: -7.0,
            }],
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
        let p = build_retention_profile("s", 1000, &[], 100, 0.01, 0.5, 1.0);
        let n_starts = 1000.0;
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
        let any = build_retention_profile(
            "s",
            2000,
            &[(h0, h1, 1.0)],
            frag,
            0.0,
            1.0 / hit_len as f64,
            1.0,
        );
        let w_any: usize = any.segments.iter().map(|s| s.end - s.start).sum();
        assert_eq!(w_any, frag + hit_len - 1, "any-overlap window");

        // m = 60 (half the site): the window narrows by 2*(m-1).
        let half = build_retention_profile("s", 2000, &[(h0, h1, 1.0)], frag, 0.0, 0.5, 1.0);
        let w_half: usize = half.segments.iter().map(|s| s.end - s.start).sum();
        assert_eq!(
            w_half,
            frag + hit_len - 1 - 2 * (60 - 1),
            "half-overlap window"
        );
        assert!(
            w_half < w_any,
            "a stricter threshold must narrow the window"
        );
        assert!((half.mass - w_half as f64).abs() < 1e-9);
    }

    /// A fragment shorter than the required overlap can never be captured by that probe,
    /// however well the probe matches.
    #[test]
    fn fragment_too_short_for_the_required_overlap_is_never_captured() {
        // 120 bp site, 50% threshold => needs 60 bp, but fragments are 40 bp.
        let p = build_retention_profile("s", 5000, &[(1000, 1120, 1.0)], 40, 0.0, 0.5, 1.0);
        assert_eq!(p.mass, 0.0, "no capture possible");
    }

    /// Independent probes combine as 1 - prod(1 - theta): two overlapping 0.5 hits give
    /// 0.75 where both are active, not 0.5 or 1.0.
    #[test]
    fn overlapping_hits_combine_independently() {
        let p = build_retention_profile(
            "s",
            3000,
            &[(1000, 1120, 0.5), (1040, 1160, 0.5)],
            50,
            0.0,
            0.25,
            1.0,
        );
        let peak = p.segments.iter().map(|s| s.p_retain).fold(0.0f64, f64::max);
        assert!((peak - 0.75).abs() < 1e-9, "expected 0.75, got {peak}");
    }

    /// Background retention is a floor everywhere; a hit lifts it above that floor.
    #[test]
    fn background_is_a_floor_under_the_profile() {
        let p_ns = 0.02;
        let p = build_retention_profile("s", 3000, &[(1000, 1120, 0.6)], 50, p_ns, 0.25, 1.0);
        let min = p
            .segments
            .iter()
            .map(|s| s.p_retain)
            .fold(f64::MAX, f64::min);
        let max = p.segments.iter().map(|s| s.p_retain).fold(0.0f64, f64::max);
        assert!(
            (min - p_ns).abs() < 1e-12,
            "floor should be p_ns, got {min}"
        );
        assert!(
            (max - 0.608).abs() < 1e-9,
            "peak should be 0.608, got {max}"
        );
    }

    /// Probe density raises retention mass - the effect the real LOD data shows.
    #[test]
    fn more_probes_means_more_retention_mass() {
        let sparse =
            build_retention_profile("s", 10_000, &[(5000, 5120, 0.9)], 200, 0.001, 0.5, 1.0);
        let dense = build_retention_profile(
            "s",
            10_000,
            &[(5000, 5120, 0.9), (5200, 5320, 0.9), (5400, 5520, 0.9)],
            200,
            0.001,
            0.5,
            1.0,
        );
        assert!(
            dense.mass > sparse.mass,
            "{} vs {}",
            dense.mass,
            sparse.mass
        );
    }

    /// Cross-probe affinity now matters, which the previous model could not express.
    #[test]
    fn weaker_probes_carry_less_mass() {
        let strong = build_retention_profile("s", 10_000, &[(5000, 5120, 0.9)], 200, 0.0, 0.5, 1.0);
        let weak = build_retention_profile("s", 10_000, &[(5000, 5120, 0.1)], 200, 0.0, 0.5, 1.0);
        assert!(strong.mass > weak.mass, "{} vs {}", strong.mass, weak.mass);
        assert!(
            (strong.mass / weak.mass - 9.0).abs() < 1e-6,
            "ratio should track theta"
        );
    }

    /// A reference shorter than the fragment still contributes, in proportion to its
    /// length, and its probe window is still scored.
    ///
    /// This replaces `reference_shorter_than_fragment_has_no_mass`, which asserted that
    /// such a reference contributes nothing. That was the behaviour, and it was wrong:
    /// the input weights encode DNA mass, so a short contig dropping to zero mass
    /// silently deleted it from the library. On the D6331 panel it put *Candida
    /// albicans* — 13,062 contigs over 12.58 Mb against a 5.8 kb fragment — at 0.244x
    /// its intended share.
    #[test]
    fn reference_shorter_than_fragment_still_contributes() {
        // 50 bp reference, 30 bp hit at [10, 40), min overlap 15, saturated theta.
        // Starts [0, 26) carry enough of the site; [26, 50) are carryover only.
        let p = build_retention_profile("s", 50, &[(10, 40, 1.0)], 200, 0.5, 0.5, 1.0);
        assert_eq!(p.n_starts, 50, "every base of the contig is a start");
        assert!(!p.segments.is_empty());
        let expected = 26.0 * 1.0 + 24.0 * 0.5;
        assert!(
            (p.mass - expected).abs() < 1e-9,
            "mass {} != {expected}",
            p.mass
        );
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
        assert_eq!(
            pairs[5].0,
            REF.as_bytes()[5],
            "reference lost sync after insertion"
        );

        // 5M then a 2 bp deletion then 5M: reference advances 12, query 10.
        let q2 = format!("{}{}", &REF[0..5], &REF[7..12]);
        let pairs2 = pairs_of(&q2, REF, 0, "5M2D5M");
        assert_eq!(pairs2.len(), 10, "deletion must not pair");
        assert_eq!(
            pairs2[5].0,
            REF.as_bytes()[7],
            "query lost sync after deletion"
        );
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
