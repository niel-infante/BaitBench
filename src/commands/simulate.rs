use anyhow::Result;
use std::path::Path;

use crate::external::minimap2;
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
    let temp_sam = args.output.with_extension("probe_hits.sam");
    let temp_log = args.output.with_extension("probe_align.log");

    log::info!("Aligning probes to reference...");
    minimap2::probe_align(
        "sr",
        args.reference,
        args.probes,
        &temp_sam,
        &temp_log,
        args.threads,
        1000,
    )?;

    // --- Step 2: Load and score probe hits ---
    let thermo_model = ThermoModel::new(args.hybridization_temperature, args.na_conc_m);
    log::info!("Loading probe hits from SAM...");
    let hits_by_probe = thermo_sim::load_probe_hits(
        &temp_sam,
        &wts,
        args.simulate_mode,
        &thermo_model,
    )?;

    let n_probes_with_hits = hits_by_probe.len();
    let n_total_hits: usize = hits_by_probe.values().map(|v| v.len()).sum();
    log::info!(
        "Found {} probes with {} total hits",
        n_probes_with_hits,
        n_total_hits
    );

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
    let _ = std::fs::remove_file(&temp_sam);
    let _ = std::fs::remove_file(&temp_log);

    log::info!(
        "Simulate complete: {} fragments written",
        all_fragments.len()
    );
    Ok(())
}
