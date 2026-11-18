use anyhow::Result;
use std::path::Path;

use crate::external::{blastn, minimap2};
use crate::fasta;
use crate::sampling::{thermo_sim, weights, FragmentLengthParams};
use crate::sampling::thermo_sim::SimulateMode;
use crate::thermodynamics::ThermoModel;


pub struct SimulateArgs<'a> {
    pub reference: &'a Path,
    pub weights: &'a Path,
    pub probes: &'a Path,
    pub num_fragments: usize,
    /// Fraction of fragments from probe binding sites (0.0–1.0)
    pub capture_fraction: f64,
    pub simulate_mode: SimulateMode,
    /// Hybridization temperature in °C (thermodynamic mode only)
    pub hybridization_temperature: f64,
    /// Na+ concentration in molar (e.g. 0.05 for 50 mM); thermodynamic mode only.
    /// SantaLucia (1998) tables assume 1 M; other values apply an entropy correction.
    pub na_conc_m: f64,
    pub seed: Option<u64>,
    pub output: &'a Path,
    pub fragment_length_mean: f64,
    pub fragment_length_min: usize,
    pub fragment_length_max: usize,
    pub threads: usize,
    /// Aligner used to place probes on the reference for scoring.
    pub aligner: crate::cli::ProbeAligner,
}

pub fn execute(args: &SimulateArgs) -> Result<()> {
    if args.fragment_length_mean < args.fragment_length_min as f64
        || args.fragment_length_mean > args.fragment_length_max as f64
    {
        anyhow::bail!(
            "--fragment-length-mean ({}) must be within [--fragment-length-min ({}), --fragment-length-max ({})]",
            args.fragment_length_mean,
            args.fragment_length_min,
            args.fragment_length_max,
        );
    }

    if let Some(seed) = args.seed {
        log::info!("Using random seed: {}", seed);
    }

    log::info!("Parsing FASTA reference...");
    let sequences = fasta::parse_fasta(args.reference)?;

    log::info!("Parsing weights file...");
    let wts = weights::parse_weights(args.weights)?;

    // --- Step 1: Align probes to reference ---
    //
    // The aligner choice matters more here than it looks. Minimizer seeding loses
    // divergent probes, and probe density drives per-locus capture, so a panel built
    // against relatives of the sampled organism can be substantially truncated before
    // the thermodynamic model sees it. `--aligner blast` (blastn-short) recovers those.
    use crate::cli::ProbeAligner;
    let temp_hits = match args.aligner {
        ProbeAligner::Blast => args.output.with_extension("probe_hits.tsv"),
        _ => args.output.with_extension("probe_hits.sam"),
    };
    let temp_log = args.output.with_extension("probe_align.log");

    log::info!(
        "Aligning probes to reference ({})...",
        match args.aligner {
            ProbeAligner::Minimap2 => "rammap, sr preset retuned for divergent probes",
            ProbeAligner::SeedExtend => "built-in ungapped local seed-and-extend",
            ProbeAligner::Blast => "blastn-short",
        }
    );
    match args.aligner {
        ProbeAligner::Blast => {
            blastn::check_available()?;
            blastn::probe_align(
                args.reference,
                args.probes,
                &temp_hits,
                &temp_log,
                args.threads,
            )?;
        }
        ProbeAligner::Minimap2 => {
            minimap2::probe_align(
                "sr",
                args.reference,
                args.probes,
                &temp_hits,
                &temp_log,
                args.threads,
                1000,
            )?;
        }
        ProbeAligner::SeedExtend => {} // done in-process below, no intermediate file
    }

    // --- Step 2: Load and score probe hits ---
    let thermo_model = ThermoModel::new(args.hybridization_temperature, args.na_conc_m);
    log::info!(
        "{}",
        match args.aligner {
            ProbeAligner::Blast => "Loading probe hits from blastn output...",
            ProbeAligner::SeedExtend => "Scoring probe placements...",
            ProbeAligner::Minimap2 => "Loading probe hits from SAM...",
        }
    );
    let hits_by_probe = match args.aligner {
        ProbeAligner::Blast => {
            thermo_sim::load_probe_hits_blast(&temp_hits, &wts, args.simulate_mode, &thermo_model)?
        }
        ProbeAligner::SeedExtend => thermo_sim::probe_hits_seed_extend(
            &sequences,
            args.probes,
            &wts,
            args.simulate_mode,
            &thermo_model,
            args.threads,
        )?,
        ProbeAligner::Minimap2 => thermo_sim::load_probe_hits(
            &temp_hits,
            &sequences,
            &wts,
            args.simulate_mode,
            &thermo_model,
        )?,
    };

    let n_probes_with_hits = hits_by_probe.len();
    let n_total_hits: usize = hits_by_probe.values().map(|v| v.len()).sum();

    // Report the probe alignment RATE, not just the count. A panel designed against
    // relatives of the sampled organism can lose a large share of its probes to the
    // aligner's seeding before the thermodynamic model ever sees them, and probe
    // density is the dominant driver of per-locus capture. Silently simulating with a
    // fraction of the panel produces plausible-looking but wrong coverage.
    let n_probes_total = fasta::parse_fasta(args.probes)
        .map(|p| p.len())
        .unwrap_or(n_probes_with_hits);
    let align_rate = if n_probes_total > 0 {
        n_probes_with_hits as f64 / n_probes_total as f64
    } else {
        0.0
    };
    log::info!(
        "Found {} of {} probes with hits ({:.1}%), {} total hits",
        n_probes_with_hits,
        n_probes_total,
        100.0 * align_rate,
        n_total_hits
    );
    if align_rate < 0.5 && n_probes_total > 0 {
        if matches!(args.aligner, ProbeAligner::Blast | ProbeAligner::SeedExtend) {
            // Already on the sensitive path, so a low rate is informative rather than
            // actionable: a cross-species panel legitimately has probes with no
            // detectable homology to this reference.
            log::info!(
                "{:.1}% of probes aligned. For a panel designed against relatives of \
                 these sequences that is expected; the remainder have no homology at \
                 80% identity over half the probe.",
                100.0 * align_rate
            );
        } else {
            log::warn!(
                "Only {:.1}% of probes aligned to the reference. If the panel targets \
                 relatives of these sequences rather than the sequences themselves, \
                 much of it may be too divergent for minimizer seeding. Try \
                 --probe-aligner blast, which is markedly more sensitive here.",
                100.0 * align_rate
            );
        }
    }

    // --- Step 3: Sample fragments ---
    let capture_fraction = args.capture_fraction.clamp(0.0, 1.0);
    let n_capture = (args.num_fragments as f64 * capture_fraction).round() as usize;
    let n_background = args.num_fragments.saturating_sub(n_capture);

    log::info!(
        "Sampling {} capture fragments (capture_fraction={:.2}) + {} background fragments",
        n_capture,
        capture_fraction,
        n_background
    );

    let mut all_fragments: Vec<(String, String)> = Vec::with_capacity(args.num_fragments);

    let len_params = FragmentLengthParams {
        mean: args.fragment_length_mean,
        min: args.fragment_length_min,
        max: args.fragment_length_max,
    };

    if n_capture > 0 {
        let capture_frags = thermo_sim::sample_capture_fragments(
            &hits_by_probe,
            &sequences,
            n_capture,
            len_params,
            args.seed,
            0,
        )?;
        log::info!("Generated {} capture fragments", capture_frags.len());
        all_fragments.extend(capture_frags);
    }

    if n_background > 0 {
        let background_frags = thermo_sim::sample_background_fragments(
            &sequences,
            &wts,
            n_background,
            len_params,
            args.seed,
            all_fragments.len(),
        )?;
        log::info!("Generated {} background fragments", background_frags.len());
        all_fragments.extend(background_frags);
    }

    // --- Step 4: Write output ---
    log::info!("Writing {} fragments to {}", all_fragments.len(), args.output.display());
    thermo_sim::write_fragments(&all_fragments, args.output)?;

    // --- Cleanup temp files ---
    let _ = std::fs::remove_file(&temp_hits);
    let _ = std::fs::remove_file(&temp_log);

    log::info!(
        "Simulate complete: {} fragments written",
        all_fragments.len()
    );
    Ok(())
}
