use anyhow::{bail, Context, Result};
use std::collections::{HashMap, HashSet};
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use crate::alignment::coverage;
use crate::cleanup;
use crate::cli::ReportMode;
use crate::commands::{filter, map_reads, prepare, sequence, simulate};
use crate::commands::report::{substitute_rmd_params, rmd_output_path};
use crate::external::{minimap2, rscript};
use crate::io_utils;
use crate::io_utils::prefixed_join;
use crate::sampling::thermo_sim::SimulateMode;

pub struct CoverageCurveArgs<'a> {
    pub targets: &'a Path,
    pub genomes: Option<&'a Path>,
    pub distractors: &'a [PathBuf],
    pub probes: &'a Path,
    pub sample: &'a HashMap<String, f64>,
    pub sample_target_map: Option<&'a HashMap<String, Vec<String>>>,
    // Sweep dimensions — each has 1+ values
    pub ct_display_values: Vec<f64>,       // original CT values as entered (for display/dirs)
    pub ct_distractor_fractions: Vec<f64>, // resolved distractor fractions (for pipeline use)
    pub cf_values: Vec<f64>,               // capture fraction values (0.0–1.0)
    pub ns_values: Vec<Option<usize>>,     // None = all sequences
    pub swept_params: Vec<String>,         // which params are swept (for dir naming + report)
    // Simulate params
    pub simulate_mode: SimulateMode,
    pub hybridization_temperature: f64,
    // Pipeline params
    pub num_fragments: usize,
    pub read_length: usize,
    pub seed: Option<u64>,
    pub fragment_length_mean: f64,
    pub fragment_length_min: usize,
    pub fragment_length_max: usize,
    pub host_fasta: Option<&'a Path>,
    pub minimap_preset: &'a str,
    pub host_minimap_preset: &'a str,
    pub threads: usize,
    pub outdir: PathBuf,
    pub output_prefix: String,
    pub report: ReportMode,
    pub cleanup: bool,
}

struct DepthCurveRow {
    ct: f64,
    capture_fraction: f64,
    num_sequences: usize, // 0 = all
    reference_id: String,
    depth_threshold: usize,
    pct_covered: f64,
}

/// Build a directory name for a parameter combo, using only swept params.
fn combo_dir_name(
    ct: f64,
    cf: f64,
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
    if swept.contains(&"capture_fraction".to_string()) {
        parts.push(format!("cf_{:.2}", cf));
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

    let has_genomes = args.genomes.is_some();

    minimap2::check_available()?;
    fs::create_dir_all(&args.outdir)?;

    let sample_ids: HashSet<String> = args.sample.keys().cloned().collect();

    let total_combos = args.ct_display_values.len() * args.cf_values.len() * args.ns_values.len();

    let sim_mode_str = match args.simulate_mode {
        SimulateMode::Thermodynamic => "thermodynamic",
        SimulateMode::Simple => "simple",
    };

    log::info!("=============================================");
    log::info!("BaitBench - Coverage Curve Analysis");
    log::info!("=============================================");
    log::info!("Targets        : {}", args.targets.display());
    if let Some(genomes) = args.genomes {
        log::info!("Genomes        : {}", genomes.display());
    }
    log::info!("Probes         : {}", args.probes.display());
    log::info!(
        "Sample         : {}",
        io_utils::format_sample_display(args.sample)
    );
    log::info!(
        "Sample entries : {} ({})",
        sample_ids.len(),
        {
            let mut ids: Vec<_> = sample_ids.iter().cloned().collect();
            ids.sort();
            ids.join(", ")
        }
    );
    log::info!("Simulate mode  : {}", sim_mode_str);
    if matches!(args.simulate_mode, SimulateMode::Thermodynamic) {
        log::info!("Hybridization temp : {}°C", args.hybridization_temperature);
    }
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
    if args.swept_params.contains(&"capture_fraction".to_string()) {
        log::info!("CF values      : {:?}", args.cf_values);
    }
    if args.swept_params.contains(&"num_sequences".to_string()) {
        log::info!("NS values      : {:?}", args.ns_values);
    }
    log::info!("Total combos   : {}", total_combos);
    log::info!("Num fragments  : {}", args.num_fragments);
    log::info!("Read length    : {}", args.read_length);
    log::info!("Output dir     : {}", args.outdir.display());
    log::info!("=============================================");

    let pfx = &args.output_prefix;
    let mut all_curves: Vec<DepthCurveRow> = Vec::new();
    let mut combo_idx = 0usize;

    // In genome mode, coverage curves track sample *targets* (derived from mapping),
    // not sample genome IDs. Resolved after first prepare step.
    let mut curve_target_ids: Option<HashSet<String>> = None;

    // Outer loop: CT values (affects prepare)
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
        run_prepare(&prep_dir, distractor_fraction, args)?;

        // Resolve curve target IDs on first iteration
        if curve_target_ids.is_none() {
            if has_genomes {
                let map_path = prefixed_join(&prep_dir, pfx, "sample_target_map.txt");
                if map_path.exists() {
                    let stm = io_utils::parse_sample_target_map(&map_path)?;
                    let mut target_ids = HashSet::new();
                    for sample_id in &sample_ids {
                        if let Some(targets) = stm.get(sample_id) {
                            for t in targets {
                                target_ids.insert(t.clone());
                            }
                        }
                    }
                    log::info!(
                        "Coverage curves will track {} sample targets (derived from genome mapping)",
                        target_ids.len()
                    );
                    curve_target_ids = Some(target_ids);
                } else {
                    curve_target_ids = Some(sample_ids.clone());
                }
            } else {
                curve_target_ids = Some(sample_ids.clone());
            }
        }

        let tracking_ids = curve_target_ids.as_ref().unwrap();

        // Determine mapping reference
        let mapping_reference = if has_genomes {
            prefixed_join(&prep_dir, pfx, "mapping_reference.fa")
        } else {
            prefixed_join(&prep_dir, pfx, "combined_reference.fa")
        };

        // Middle loop: capture-fraction values (affects simulate)
        for &cf in &args.cf_values {
            let sim_dir = prep_dir.join(format!("_sim_cf_{:.2}", cf));
            fs::create_dir_all(&sim_dir)?;

            run_simulate(&sim_dir, &prep_dir, cf, args)?;

            let fragments_path = prefixed_join(&sim_dir, pfx, "fragments.fa");

            // Inner loop: num-sequences values (affects sequence step)
            for ns_val in &args.ns_values {
                combo_idx += 1;
                let dir_name =
                    combo_dir_name(ct_display, cf, *ns_val, &args.swept_params);
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
                    input: &fragments_path,
                    output: &prefixed_join(&combo_dir, pfx, "reads.fa"),
                    read_length: args.read_length,
                    num_sequences: *ns_val,
                    seed: args.seed,
                })?;

                // Optional host filtering
                let reads_for_mapping = if let Some(host) = args.host_fasta {
                    log::info!("  Filtering host reads...");
                    filter::execute(&filter::FilterArgs {
                        host,
                        reads: &prefixed_join(&combo_dir, pfx, "reads.fa"),
                        minimap_preset: args.host_minimap_preset,
                        output: &prefixed_join(&combo_dir, pfx, "filtered.fa"),
                        log_file: &prefixed_join(&combo_dir, pfx, "host_filter.log"),
                    })?;
                    prefixed_join(&combo_dir, pfx, "filtered.fa")
                } else {
                    prefixed_join(&combo_dir, pfx, "reads.fa")
                };

                // Map reads
                log::info!("  Mapping reads...");
                map_reads::execute(&map_reads::MapArgs {
                    reference: &mapping_reference,
                    reads: &reads_for_mapping,
                    minimap_preset: args.minimap_preset,
                    output: &prefixed_join(&combo_dir, pfx, "mapped.sam"),
                    log_file: &prefixed_join(&combo_dir, pfx, "mapping.log"),
                })?;

                // Compute coverage
                let sam_path = prefixed_join(&combo_dir, pfx, "mapped.sam");
                let coverage_result = coverage::compute_coverage(&sam_path)?;

                let ns_tsv = ns_val.unwrap_or(0);

                for target_id in tracking_ids {
                    let curves = match coverage_result.coverage.get(target_id) {
                        Some(depths) => {
                            let ref_len = coverage_result
                                .ref_lengths
                                .get(target_id)
                                .copied()
                                .unwrap_or(depths.len());
                            compute_depth_curve(depths, ref_len)
                        }
                        None => {
                            log::warn!(
                                "No coverage found for target '{}' in combo {}",
                                target_id,
                                dir_name
                            );
                            vec![(1, 0.0)]
                        }
                    };

                    for (threshold, pct) in curves {
                        all_curves.push(DepthCurveRow {
                            ct: ct_display,
                            capture_fraction: cf,
                            num_sequences: ns_tsv,
                            reference_id: target_id.clone(),
                            depth_threshold: threshold,
                            pct_covered: pct,
                        });
                    }
                }
            }
        }
    }

    // Write aggregated depth curves TSV
    let curves_path = prefixed_join(&args.outdir, pfx, "coverage_curve_depth_curves.tsv");
    write_depth_curves_tsv(&curves_path, &all_curves)?;
    log::info!("Depth curves written to {}", curves_path.display());

    // Write run_params.tsv
    let params_path = prefixed_join(&args.outdir, pfx, "run_params.tsv");
    write_run_params(&params_path, args)?;

    // Generate report
    let tracking_ids = curve_target_ids.as_ref().unwrap_or(&sample_ids);
    match args.report {
        ReportMode::None => {
            log::info!("Skipping report generation (--report none)");
        }
        ReportMode::Full => {
            if rscript::check_available() {
                let report_path = prefixed_join(&args.outdir, pfx, "coverage_curve_report.html");
                log::info!("Generating coverage curve report...");
                match generate_report(&curves_path, tracking_ids, &args.swept_params, &params_path, &report_path) {
                    Ok(()) => log::info!(
                        "Report generated: {}",
                        report_path.display()
                    ),
                    Err(e) => log::warn!("Report generation failed (non-fatal): {}", e),
                }
            } else {
                log::warn!("Rscript not found — skipping HTML report.");
            }
        }
        ReportMode::Rmd => {
            let report_path = prefixed_join(&args.outdir, pfx, "coverage_curve_report.html");
            log::info!("Generating coverage curve RMarkdown file...");
            match write_coverage_curve_rmd(&curves_path, tracking_ids, &args.swept_params, &params_path, &report_path) {
                Ok(()) => {}
                Err(e) => log::warn!("RMarkdown generation failed (non-fatal): {}", e),
            }
        }
    }

    // Cleanup intermediate directories if requested
    if args.cleanup {
        log::info!("Cleaning up intermediate files...");
        cleanup::cleanup_all_subdirs(&args.outdir);
    }

    log::info!("=============================================");
    log::info!("Coverage curve analysis complete!");
    log::info!("Results in {}", args.outdir.display());
    log::info!("=============================================");

    Ok(())
}

/// Run prepare into a shared prep directory for this CT value.
fn run_prepare(
    outdir: &Path,
    distractor_fraction: f64,
    args: &CoverageCurveArgs,
) -> Result<()> {
    let pfx = &args.output_prefix;

    // Skip if already done (idempotent)
    if prefixed_join(outdir, pfx, "weights.txt").exists() {
        log::info!("  Reusing cached prepare from {}", outdir.display());
        return Ok(());
    }

    log::info!("  Preparing reference...");
    prepare::execute(&prepare::PrepareArgs {
        targets: args.targets,
        genomes: args.genomes,
        distractors: args.distractors,
        sample: Some(args.sample),
        sample_target_map: args.sample_target_map,
        distractor_fraction,
        outdir,
        output_prefix: pfx,
    })?;

    Ok(())
}

/// Run simulate (probe-biased + background) for a specific capture fraction.
/// Reads from the prep dir's reference + weights.
fn run_simulate(
    sim_dir: &Path,
    prep_dir: &Path,
    capture_fraction: f64,
    args: &CoverageCurveArgs,
) -> Result<()> {
    let pfx = &args.output_prefix;

    // Skip if already done (idempotent)
    if prefixed_join(sim_dir, pfx, "fragments.fa").exists() {
        log::info!("  Reusing cached simulate from {}", sim_dir.display());
        return Ok(());
    }

    log::info!("  Simulating fragments (capture_fraction={:.2})...", capture_fraction);
    simulate::execute(&simulate::SimulateArgs {
        reference: &prefixed_join(prep_dir, pfx, "combined_reference.fa"),
        weights: &prefixed_join(prep_dir, pfx, "weights.txt"),
        probes: args.probes,
        num_fragments: args.num_fragments,
        capture_fraction,
        simulate_mode: args.simulate_mode,
        hybridization_temperature: args.hybridization_temperature,
        seed: args.seed,
        output: &prefixed_join(sim_dir, pfx, "fragments.fa"),
        fragment_length_mean: args.fragment_length_mean,
        fragment_length_min: args.fragment_length_min,
        fragment_length_max: args.fragment_length_max,
        threads: args.threads,
    })?;

    Ok(())
}

/// Compute the depth threshold curve for a single reference.
fn compute_depth_curve(depths: &[u32], ref_length: usize) -> Vec<(usize, f64)> {
    if ref_length == 0 || depths.is_empty() {
        return vec![(1, 0.0)];
    }

    let max_depth = *depths.iter().max().unwrap_or(&0) as usize;
    if max_depth == 0 {
        return vec![(1, 0.0)];
    }

    let mut hist = vec![0usize; max_depth + 1];
    for &d in depths {
        hist[d as usize] += 1;
    }
    let zero_depth_extra = if ref_length > depths.len() {
        ref_length - depths.len()
    } else {
        0
    };
    hist[0] += zero_depth_extra;

    let mut cumsum = vec![0usize; max_depth + 2];
    for t in (0..=max_depth).rev() {
        cumsum[t] = cumsum[t + 1] + hist[t];
    }

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
        "ct\tcapture_fraction\tnum_sequences\treference_id\tdepth_threshold\tpct_covered"
    )?;
    for row in rows {
        writeln!(
            w,
            "{}\t{:.4}\t{}\t{}\t{}\t{:.2}",
            row.ct, row.capture_fraction, row.num_sequences, row.reference_id,
            row.depth_threshold, row.pct_covered
        )?;
    }

    w.flush()?;
    Ok(())
}

fn write_run_params(path: &Path, args: &CoverageCurveArgs) -> Result<()> {
    let sim_mode_str = match args.simulate_mode {
        SimulateMode::Thermodynamic => "thermodynamic",
        SimulateMode::Simple => "simple",
    };

    let file = File::create(path)
        .with_context(|| format!("Cannot create params file: {}", path.display()))?;
    let mut w = BufWriter::new(file);

    writeln!(w, "parameter\tflag\tvalue")?;
    writeln!(w, "targets\t--targets\t{}", args.targets.display())?;
    if let Some(genomes) = args.genomes {
        writeln!(w, "genomes\t--genomes\t{}", genomes.display())?;
    }
    for d in args.distractors {
        writeln!(w, "distractors\t--distractors\t{}", d.display())?;
    }
    writeln!(w, "probes\t--probes\t{}", args.probes.display())?;
    writeln!(w, "sample\t--sample\t{}", io_utils::format_sample_display(args.sample))?;
    writeln!(w, "num_fragments\t--num-fragments\t{}", args.num_fragments)?;
    writeln!(w, "read_length\t--read-length\t{}", args.read_length)?;
    writeln!(w, "simulate_mode\t--simulate-mode\t{}", sim_mode_str)?;
    writeln!(w, "hybridization_temperature\t--hybridization-temperature\t{}", args.hybridization_temperature)?;
    writeln!(w, "minimap_preset\t--minimap-preset\t{}", args.minimap_preset)?;
    writeln!(w, "fragment_length_mean\t--fragment-length-mean\t{}", args.fragment_length_mean)?;
    writeln!(w, "fragment_length_min\t--fragment-length-min\t{}", args.fragment_length_min)?;
    writeln!(w, "fragment_length_max\t--fragment-length-max\t{}", args.fragment_length_max)?;
    writeln!(w, "seed\t--seed\t{}", args.seed.map(|s| s.to_string()).unwrap_or_else(|| "none".to_string()))?;
    writeln!(w, "threads\t--threads\t{}", args.threads)?;
    if !args.swept_params.is_empty() {
        writeln!(w, "swept_params\t(info)\t{}", args.swept_params.join(", "))?;
    }
    writeln!(w, "outdir\t--outdir\t{}", args.outdir.display())?;

    w.flush()?;
    Ok(())
}

fn write_coverage_curve_rmd(
    sweep_tsv: &Path,
    sample_ids: &HashSet<String>,
    swept_params: &[String],
    params_path: &Path,
    output_path: &Path,
) -> Result<()> {
    let r_dir = rscript::find_r_dir()
        .ok_or_else(|| anyhow::anyhow!("Cannot find R scripts directory."))?;

    let rmd_template = r_dir.join("coverage_curve.Rmd");
    if !rmd_template.exists() {
        bail!("RMarkdown template not found: {}", rmd_template.display());
    }

    let sweep_abs = fs::canonicalize(sweep_tsv)?;

    let mut ids: Vec<&String> = sample_ids.iter().collect();
    ids.sort();
    let ids_r = format!(
        "!r c({})",
        ids.iter()
            .map(|s| format!("\"{}\"", s))
            .collect::<Vec<_>>()
            .join(", ")
    );

    let swept_r = if swept_params.is_empty() {
        "!r c()".to_string()
    } else {
        format!(
            "!r c({})",
            swept_params
                .iter()
                .map(|s| format!("\"{}\"", s))
                .collect::<Vec<_>>()
                .join(", ")
        )
    };

    let params_abs = fs::canonicalize(params_path)?;

    let params = vec![
        ("sweep_file", sweep_abs.to_str().unwrap_or("")),
        ("sample_ids", &ids_r),
        ("swept_params", &swept_r),
        ("params_file", params_abs.to_str().unwrap_or("")),
    ];

    let template_content = std::fs::read_to_string(&rmd_template)
        .with_context(|| format!("Failed to read template: {}", rmd_template.display()))?;

    let output_content = substitute_rmd_params(&template_content, &params);

    let html_path = if output_path.is_absolute() {
        output_path.to_path_buf()
    } else {
        std::env::current_dir()?.join(output_path)
    };
    let rmd_path = rmd_output_path(&html_path);
    std::fs::write(&rmd_path, &output_content)
        .with_context(|| format!("Failed to write RMarkdown file: {}", rmd_path.display()))?;

    log::info!("RMarkdown file written: {}", rmd_path.display());
    log::info!(
        "Edit and render with: Rscript -e 'rmarkdown::render(\"{}\")'",
        rmd_path.display()
    );
    Ok(())
}

fn generate_report(
    sweep_tsv: &Path,
    sample_ids: &HashSet<String>,
    swept_params: &[String],
    params_path: &Path,
    output_path: &Path,
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
    let params_abs = fs::canonicalize(params_path)?;
    let output_abs = if output_path.is_absolute() {
        output_path.to_path_buf()
    } else {
        std::env::current_dir()?.join(output_path)
    };

    let sweep_str = sweep_abs.to_str().unwrap_or("");
    let params_str = params_abs.to_str().unwrap_or("");
    let output_str = output_abs.to_str().unwrap_or("");

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
            "--params",
            params_str,
            "--output",
            output_str,
        ],
    )
}
