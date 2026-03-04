use anyhow::{bail, Context, Result};
use std::collections::{HashMap, HashSet};
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use crate::alignment::coverage;
use crate::commands::{capture, enrich, filter, map_reads, prepare, sequence, simulate};
use crate::external::{minimap2, rscript};
use crate::io_utils;

pub struct CoverageCurveArgs<'a> {
    pub targets: &'a Path,
    pub distractors: &'a [PathBuf],
    pub probes: &'a Path,
    pub sample: &'a HashMap<String, f64>,
    // Sweep dimensions — each has 1+ values
    pub ct_display_values: Vec<f64>,     // original CT values as entered by the user (for display/TSV/dirs)
    pub ct_distractor_fractions: Vec<f64>, // resolved distractor fractions (for pipeline use)
    pub fe_values: Vec<Option<f64>>,     // None = no enrichment
    pub ns_values: Vec<Option<usize>>,   // None = all sequences
    pub swept_params: Vec<String>,       // which params are swept (for dir naming + report)
    // Pipeline params
    pub num_fragments: usize,
    pub read_length: usize,
    pub seed: Option<u64>,
    pub fragment_length_mean: f64,
    pub fragment_length_min: usize,
    pub fragment_length_max: usize,
    pub capture_method: capture::CaptureMethod,
    pub max_mismatches: u32,
    pub min_match_bases: u32,
    pub blast_db: Option<&'a str>,
    pub host_fasta: Option<&'a Path>,
    pub minimap_preset: &'a str,
    pub host_minimap_preset: &'a str,
    pub threads: usize,
    pub outdir: PathBuf,
    pub no_report: bool,
}

struct DepthCurveRow {
    ct: f64,
    fold_enrichment: f64, // 0.0 = no enrichment
    num_sequences: usize, // 0 = all
    reference_id: String,
    depth_threshold: usize,
    pct_covered: f64,
}

/// Build a directory name for a parameter combo, using only swept params.
fn combo_dir_name(
    ct: f64,
    fe: Option<f64>,
    ns: Option<usize>,
    swept: &[String],
) -> String {
    if swept.is_empty() {
        return "run".to_string();
    }
    let mut parts = Vec::new();
    if swept.contains(&"ct".to_string()) {
        parts.push(format!("ct_{}", ct));
    }
    if swept.contains(&"fold_enrichment".to_string()) {
        parts.push(format!(
            "fe_{}",
            fe.map(|v| format!("{}", v)).unwrap_or_else(|| "none".into())
        ));
    }
    if swept.contains(&"num_sequences".to_string()) {
        parts.push(format!(
            "ns_{}",
            ns.map(|v| v.to_string()).unwrap_or_else(|| "all".into())
        ));
    }
    parts.join("_")
}

pub fn execute(args: &CoverageCurveArgs) -> Result<()> {
    // Validate inputs
    if !args.targets.exists() {
        bail!("Targets file not found: {}", args.targets.display());
    }
    if !args.probes.exists() {
        bail!("Probes file not found: {}", args.probes.display());
    }
    if args.sample.is_empty() {
        bail!("Sample is empty — at least one sample target is required");
    }

    minimap2::check_available()?;
    fs::create_dir_all(&args.outdir)?;

    let sample_ids: HashSet<String> = args.sample.keys().cloned().collect();

    let total_combos = args.ct_display_values.len() * args.fe_values.len() * args.ns_values.len();

    log::info!("=============================================");
    log::info!("BaitBench - Coverage Curve Analysis");
    log::info!("=============================================");
    log::info!("Targets        : {}", args.targets.display());
    log::info!("Probes         : {}", args.probes.display());
    log::info!(
        "Sample         : {}",
        io_utils::format_sample_display(args.sample)
    );
    log::info!(
        "Sample targets : {} ({})",
        sample_ids.len(),
        {
            let mut ids: Vec<_> = sample_ids.iter().cloned().collect();
            ids.sort();
            ids.join(", ")
        }
    );
    if args.swept_params.is_empty() {
        log::info!("Mode           : single run (no parameters swept)");
    } else {
        log::info!(
            "Swept params   : {}",
            args.swept_params.join(", ")
        );
    }
    if args.swept_params.contains(&"ct".to_string()) {
        log::info!("CT values      : {:?}", args.ct_display_values);
    }
    if args.swept_params.contains(&"fold_enrichment".to_string()) {
        log::info!("FE values      : {:?}", args.fe_values);
    }
    if args.swept_params.contains(&"num_sequences".to_string()) {
        log::info!("NS values      : {:?}", args.ns_values);
    }
    log::info!("Total combos   : {}", total_combos);
    log::info!("Num fragments  : {}", args.num_fragments);
    log::info!("Read length    : {}", args.read_length);
    log::info!("Output dir     : {}", args.outdir.display());
    log::info!("=============================================");

    let mut all_curves: Vec<DepthCurveRow> = Vec::new();
    let mut combo_idx = 0usize;

    // Outer loop: CT values (affects prepare/simulate/capture)
    for (ct_idx, (&ct_display, &distractor_fraction)) in args
        .ct_display_values
        .iter()
        .zip(args.ct_distractor_fractions.iter())
        .enumerate()
    {
        // Shared prep directory for this CT value
        let prep_dir = args.outdir.join(format!("_prep_ct_{}", ct_idx));
        fs::create_dir_all(&prep_dir)?;

        log::info!(
            "--- Preparing reference (ct={}, distractor_fraction={:.6}) ---",
            ct_display, distractor_fraction
        );
        run_prepare_simulate_capture(&prep_dir, distractor_fraction, args)?;

        // Middle loop: fold-enrichment values (affects enrich step)
        for fe_val in &args.fe_values {
            // Run enrich if needed
            let capture_output = if let Some(fe) = fe_val {
                let enrich_dir = prep_dir.join(format!("_enrich_fe_{}", fe));
                fs::create_dir_all(&enrich_dir)?;
                let enriched_path = enrich_dir.join("enriched.fa");
                if !enriched_path.exists() {
                    log::info!("  Applying {:.1}x fold enrichment...", fe);
                    enrich::execute(&enrich::EnrichArgs {
                        captured: &prep_dir.join("captured.fa"),
                        fragments: &prep_dir.join("fragments.fa"),
                        targets: &prep_dir.join("targets.txt"),
                        distractors: &prep_dir.join("distractors.txt"),
                        fold_enrichment: *fe,
                        seed: args.seed,
                        output: &enriched_path,
                    })?;
                }
                enriched_path
            } else {
                prep_dir.join("captured.fa")
            };

            // Inner loop: num-sequences values (affects sequence step)
            for ns_val in &args.ns_values {
                combo_idx += 1;
                let dir_name =
                    combo_dir_name(ct_display, *fe_val, *ns_val, &args.swept_params);
                let combo_dir = args.outdir.join(&dir_name);
                fs::create_dir_all(&combo_dir)?;

                log::info!(
                    "--- Combo {}/{}: {} ---",
                    combo_idx,
                    total_combos,
                    dir_name
                );

                // Sequence
                log::info!("  Sequencing...");
                sequence::execute(&sequence::SequenceArgs {
                    input: &capture_output,
                    output: &combo_dir.join("reads.fa"),
                    read_length: args.read_length,
                    num_sequences: *ns_val,
                    seed: args.seed,
                })?;

                // Optional host filtering
                let reads_for_mapping = if let Some(host) = args.host_fasta {
                    log::info!("  Filtering host reads...");
                    filter::execute(&filter::FilterArgs {
                        host,
                        reads: &combo_dir.join("reads.fa"),
                        minimap_preset: args.host_minimap_preset,
                        output: &combo_dir.join("filtered.fa"),
                        log_file: &combo_dir.join("host_filter.log"),
                    })?;
                    combo_dir.join("filtered.fa")
                } else {
                    combo_dir.join("reads.fa")
                };

                // Map reads — use reference from prep dir
                log::info!("  Mapping reads...");
                map_reads::execute(&map_reads::MapArgs {
                    reference: &prep_dir.join("combined_reference.fa"),
                    reads: &reads_for_mapping,
                    minimap_preset: args.minimap_preset,
                    output: &combo_dir.join("mapped.sam"),
                    log_file: &combo_dir.join("mapping.log"),
                })?;

                // Compute coverage
                let sam_path = combo_dir.join("mapped.sam");
                let coverage_result = coverage::compute_coverage(&sam_path)?;

                // Extract depth curves for sample targets
                let fe_tsv = fe_val.unwrap_or(0.0);
                let ns_tsv = ns_val.unwrap_or(0);

                for sample_id in &sample_ids {
                    let curves = match coverage_result.coverage.get(sample_id) {
                        Some(depths) => {
                            let ref_len = coverage_result
                                .ref_lengths
                                .get(sample_id)
                                .copied()
                                .unwrap_or(depths.len());
                            compute_depth_curve(depths, ref_len)
                        }
                        None => {
                            log::warn!(
                                "No coverage found for sample target '{}' in combo {}",
                                sample_id,
                                dir_name
                            );
                            vec![(1, 0.0)]
                        }
                    };

                    for (threshold, pct) in curves {
                        all_curves.push(DepthCurveRow {
                            ct: ct_display,
                            fold_enrichment: fe_tsv,
                            num_sequences: ns_tsv,
                            reference_id: sample_id.clone(),
                            depth_threshold: threshold,
                            pct_covered: pct,
                        });
                    }
                }
            }
        }
    }

    // Write aggregated depth curves TSV
    let curves_path = args.outdir.join("coverage_curve_depth_curves.tsv");
    write_depth_curves_tsv(&curves_path, &all_curves)?;
    log::info!("Depth curves written to {}", curves_path.display());

    // Generate report
    if args.no_report {
        log::info!("Skipping report generation (--no-report)");
    } else if rscript::check_available() {
        log::info!("Generating coverage curve report...");
        match generate_report(&curves_path, &sample_ids, &args.swept_params, &args.outdir) {
            Ok(()) => log::info!(
                "Report generated: {}",
                args.outdir.join("coverage_curve_report.html").display()
            ),
            Err(e) => log::warn!("Report generation failed (non-fatal): {}", e),
        }
    } else {
        log::warn!("Rscript not found — skipping HTML report.");
    }

    log::info!("=============================================");
    log::info!("Coverage curve analysis complete!");
    log::info!("Results in {}", args.outdir.display());
    log::info!("=============================================");

    Ok(())
}

/// Run prepare + simulate + capture into a shared prep directory.
fn run_prepare_simulate_capture(
    outdir: &Path,
    distractor_fraction: f64,
    args: &CoverageCurveArgs,
) -> Result<()> {
    // Skip if already done (idempotent for nested loops)
    if outdir.join("captured.fa").exists() {
        log::info!("  Reusing cached prepare/simulate/capture from {}", outdir.display());
        return Ok(());
    }

    // Step 1: Prepare
    log::info!("  Preparing reference...");
    prepare::execute(&prepare::PrepareArgs {
        targets: args.targets,
        distractors: args.distractors,
        sample: Some(args.sample),
        distractor_fraction,
        outdir,
    })?;

    // Step 2: Simulate
    log::info!("  Simulating fragments...");
    simulate::execute(&simulate::SimulateArgs {
        reference: &outdir.join("combined_reference.fa"),
        weights: &outdir.join("weights.txt"),
        num_fragments: args.num_fragments,
        seed: args.seed,
        output: &outdir.join("fragments.fa"),
        fragment_length_mean: args.fragment_length_mean,
        fragment_length_min: args.fragment_length_min,
        fragment_length_max: args.fragment_length_max,
    })?;

    // Step 3: Capture
    log::info!("  Simulating capture...");
    capture::execute(&capture::CaptureArgs {
        method: args.capture_method,
        probes: args.probes,
        fragments: &outdir.join("fragments.fa"),
        max_mismatches: args.max_mismatches,
        min_match_bases: args.min_match_bases,
        blast_db: args.blast_db,
        output: &outdir.join("captured.fa"),
        log_file: &outdir.join("capture.log"),
        threads: args.threads,
    })?;

    Ok(())
}

/// Compute the depth threshold curve for a single reference.
///
/// Returns a vector of (depth_threshold, pct_covered) tuples for thresholds 1..max_depth.
/// Uses a histogram + reverse cumulative sum for O(n + max_depth) efficiency.
fn compute_depth_curve(depths: &[u32], ref_length: usize) -> Vec<(usize, f64)> {
    if ref_length == 0 || depths.is_empty() {
        return vec![(1, 0.0)];
    }

    let max_depth = *depths.iter().max().unwrap_or(&0) as usize;
    if max_depth == 0 {
        return vec![(1, 0.0)];
    }

    // Build histogram: hist[d] = count of positions with depth exactly d
    let mut hist = vec![0usize; max_depth + 1];
    for &d in depths {
        hist[d as usize] += 1;
    }
    // Account for positions beyond the depth vector (0-depth)
    let zero_depth_extra = if ref_length > depths.len() {
        ref_length - depths.len()
    } else {
        0
    };
    hist[0] += zero_depth_extra;

    // Reverse cumulative sum: cumsum[t] = number of positions with depth >= t
    let mut cumsum = vec![0usize; max_depth + 2];
    for t in (0..=max_depth).rev() {
        cumsum[t] = cumsum[t + 1] + hist[t];
    }

    // Build curve for thresholds 1..max_depth
    let mut curve = Vec::with_capacity(max_depth);
    for t in 1..=max_depth {
        let pct = cumsum[t] as f64 / ref_length as f64 * 100.0;
        curve.push((t, pct));
    }

    curve
}

fn write_depth_curves_tsv(path: &Path, rows: &[DepthCurveRow]) -> Result<()> {
    let file = File::create(path)
        .with_context(|| format!("Cannot create depth curves file: {}", path.display()))?;
    let mut w = BufWriter::new(file);

    writeln!(
        w,
        "ct\tfold_enrichment\tnum_sequences\treference_id\tdepth_threshold\tpct_covered"
    )?;
    for row in rows {
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{:.2}",
            row.ct, row.fold_enrichment, row.num_sequences, row.reference_id,
            row.depth_threshold, row.pct_covered
        )?;
    }

    w.flush()?;
    Ok(())
}

fn generate_report(
    sweep_tsv: &Path,
    sample_ids: &HashSet<String>,
    swept_params: &[String],
    outdir: &Path,
) -> Result<()> {
    let r_dir = rscript::find_r_dir()
        .ok_or_else(|| anyhow::anyhow!("Cannot find R scripts directory."))?;

    let script = r_dir.join("coverage_curve.R");
    if !script.exists() {
        bail!(
            "Coverage curve R script not found: {}",
            script.display()
        );
    }

    let sweep_abs = fs::canonicalize(sweep_tsv)?;
    let output_path = if outdir.is_absolute() {
        outdir.join("coverage_curve_report.html")
    } else {
        std::env::current_dir()?
            .join(outdir)
            .join("coverage_curve_report.html")
    };

    let sweep_str = sweep_abs.to_str().unwrap_or("");
    let output_str = output_path.to_str().unwrap_or("");

    let mut ids: Vec<&String> = sample_ids.iter().collect();
    ids.sort();
    let ids_str = ids
        .iter()
        .map(|s| s.as_str())
        .collect::<Vec<_>>()
        .join(",");

    let swept_str = swept_params.join(",");

    rscript::run_rscript(
        &script,
        &[
            "--sweep",
            sweep_str,
            "--sample-ids",
            &ids_str,
            "--swept-params",
            &swept_str,
            "--output",
            output_str,
        ],
    )
}
