use anyhow::{bail, Context, Result};
use std::collections::HashSet;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

use crate::cleanup;
use crate::cli::ReportMode;
use crate::commands::report::{rmd_output_path, substitute_rmd_params};
use crate::commands::{probe_coverage, xreact};
use crate::external::{minimap2, rscript};
use crate::fasta;
use crate::io_utils::prefixed_join;

pub struct AssessProbesArgs<'a> {
    pub targets: &'a Path,
    pub probes: &'a Path,
    pub genomes: &'a [PathBuf],
    pub threshold: f64,
    pub minimap_preset: &'a str,
    pub proximity: usize,
    pub outdir: &'a Path,
    pub output_prefix: &'a str,
    pub report: ReportMode,
    pub cleanup: bool,
    pub build_stats_file: Option<&'a Path>,
    pub build_params_file: Option<&'a Path>,
    pub refine_threshold: f64,
    pub refine_iterations: Option<usize>,
    pub refine_until_stable: bool,
}

pub fn execute(args: &AssessProbesArgs) -> Result<()> {
    if !args.targets.exists() {
        bail!("Targets file not found: {}", args.targets.display());
    }
    if !args.probes.exists() {
        bail!("Probes file not found: {}", args.probes.display());
    }
    for path in args.genomes {
        if !path.exists() {
            bail!("Genome file not found: {}", path.display());
        }
    }

    fs::create_dir_all(args.outdir)?;
    minimap2::check_available()?;

    let pfx = args.output_prefix;
    let from_build = args.build_stats_file.is_some();

    log::info!("=============================================");
    log::info!("BaitBench - Probe Assessment");
    log::info!("=============================================");
    log::info!("Targets  : {}", args.targets.display());
    log::info!("Probes   : {}", args.probes.display());
    if !args.genomes.is_empty() {
        for g in args.genomes {
            log::info!("Genome   : {}", g.display());
        }
    }
    log::info!("Threshold: {:.1}%", args.threshold);
    log::info!("Preset   : {}", args.minimap_preset);
    log::info!("Proximity: {} bp", args.proximity);
    log::info!("Output   : {}", args.outdir.display());
    if from_build {
        log::info!("Mode     : chained from build-probes");
    }

    // --- Step 1: Run probe coverage analysis ---
    let cov_prefix = format!("{}cov_", pfx);
    log::info!("Running probe coverage analysis...");
    probe_coverage::execute(&probe_coverage::ProbeCoverageArgs {
        targets: args.targets,
        probes: args.probes,
        outdir: args.outdir,
        output_prefix: &cov_prefix,
        minimap_preset: args.minimap_preset,
        proximity: args.proximity,
        report: ReportMode::None,
        cleanup: false,
    })?;

    // --- Step 2: Run cross-reactivity analysis ---
    let xreact_prefix = format!("{}xreact_", pfx);
    log::info!("Running cross-reactivity analysis...");
    xreact::execute(&xreact::XreactArgs {
        probes: args.probes,
        against: args.genomes,
        self_mode: true,
        threshold: args.threshold,
        minimap_preset: args.minimap_preset,
        outdir: args.outdir,
        output_prefix: &xreact_prefix,
        report: ReportMode::None,
        cleanup: false,
    })?;

    // --- Step 3: Write combined run params ---
    let params_path = prefixed_join(args.outdir, pfx, "assess_run_params.tsv");
    write_run_params(&params_path, args)?;

    // --- Step 4: Collect data file paths ---
    let xreact_hits = prefixed_join(args.outdir, &xreact_prefix, "hits.tsv");
    let xreact_summary = prefixed_join(args.outdir, &xreact_prefix, "summary.tsv");
    let cov_summary = prefixed_join(args.outdir, &cov_prefix, "probe_coverage_summary.tsv");
    let cov_depth = prefixed_join(args.outdir, &cov_prefix, "probe_depth.tsv");
    let cov_multi = prefixed_join(args.outdir, &cov_prefix, "multi_mapping_probes.tsv");

    // --- Step 5: Generate report ---
    match args.report {
        ReportMode::None => {
            log::info!("Skipping report generation (--report none)");
        }
        ReportMode::Full => {
            if rscript::check_available() {
                let report_path =
                    prefixed_join(args.outdir, pfx, "assess_probes_report.html");
                log::info!("Generating combined probe assessment report...");
                match generate_assess_report(
                    args.build_stats_file,
                    args.build_params_file,
                    &xreact_hits,
                    &xreact_summary,
                    args.threshold,
                    &cov_summary,
                    &cov_depth,
                    &cov_multi,
                    args.proximity,
                    &params_path,
                    &report_path,
                ) {
                    Ok(()) => log::info!("Report generated: {}", report_path.display()),
                    Err(e) => log::warn!("Report generation failed (non-fatal): {}", e),
                }
            } else {
                log::warn!("Rscript not found -- skipping HTML report.");
            }
        }
        ReportMode::Rmd => {
            let report_path =
                prefixed_join(args.outdir, pfx, "assess_probes_report.html");
            log::info!("Generating probe assessment RMarkdown file...");
            match write_assess_rmd(
                args.build_stats_file,
                args.build_params_file,
                &xreact_hits,
                &xreact_summary,
                args.threshold,
                &cov_summary,
                &cov_depth,
                &cov_multi,
                args.proximity,
                &params_path,
                &report_path,
            ) {
                Ok(()) => {}
                Err(e) => log::warn!("RMarkdown generation failed (non-fatal): {}", e),
            }
        }
    }

    // --- Step 6: Refinement iterations ---
    if args.refine_iterations.is_some() || args.refine_until_stable {
        run_refinement(args, &cov_summary)?;
    }

    // --- Step 7: Cleanup ---
    if args.cleanup {
        log::info!("Cleaning up intermediate files...");
        let cov_intermediates: Vec<String> =
            ["probe_alignment.sam", "probe_alignment.log"]
                .iter()
                .map(|f| format!("{}{}", cov_prefix, f))
                .collect();
        let xreact_intermediates: Vec<String> = ["against.log", "self.log"]
            .iter()
            .map(|f| format!("{}{}", xreact_prefix, f))
            .collect();

        let all_intermediates: Vec<&str> = cov_intermediates
            .iter()
            .chain(xreact_intermediates.iter())
            .map(|s| s.as_str())
            .collect();
        cleanup::cleanup_files(args.outdir, &all_intermediates);
    }

    log::info!("=============================================");
    log::info!("Probe assessment complete!");
    log::info!("Results in {}", args.outdir.display());
    log::info!("=============================================");

    Ok(())
}

fn write_run_params(path: &Path, args: &AssessProbesArgs) -> Result<()> {
    let file = File::create(path)
        .with_context(|| format!("Cannot create params file: {}", path.display()))?;
    let mut w = BufWriter::new(file);

    writeln!(w, "parameter\tflag\tvalue")?;
    writeln!(w, "targets\t--targets\t{}", args.targets.display())?;
    writeln!(w, "probes\t--probes\t{}", args.probes.display())?;
    for g in args.genomes {
        writeln!(w, "genomes\t--genomes\t{}", g.display())?;
    }
    writeln!(w, "threshold\t--threshold\t{:.1}", args.threshold)?;
    writeln!(
        w,
        "minimap_preset\t--minimap-preset\t{}",
        args.minimap_preset
    )?;
    writeln!(w, "proximity\t--proximity\t{}", args.proximity)?;
    writeln!(w, "outdir\t-o\t{}", args.outdir.display())?;

    w.flush()?;
    Ok(())
}

/// Make a path absolute (canonicalize if it exists, otherwise join with cwd).
fn abs_path(p: &Path) -> Result<PathBuf> {
    if p.exists() {
        Ok(std::fs::canonicalize(p)?)
    } else if p.is_absolute() {
        Ok(p.to_path_buf())
    } else {
        Ok(std::env::current_dir()?.join(p))
    }
}

fn abs_path_str(p: &Path) -> Result<String> {
    Ok(abs_path(p)?
        .to_str()
        .unwrap_or("")
        .to_string())
}

fn generate_assess_report(
    build_stats_file: Option<&Path>,
    build_params_file: Option<&Path>,
    xreact_hits: &Path,
    xreact_summary: &Path,
    threshold: f64,
    cov_summary: &Path,
    cov_depth: &Path,
    cov_multi: &Path,
    proximity: usize,
    params_path: &Path,
    output_path: &Path,
) -> Result<()> {
    let r_dir = rscript::find_r_dir()
        .ok_or_else(|| anyhow::anyhow!("Cannot find R scripts directory."))?;

    let script = r_dir.join("assess_probes.R");
    if !script.exists() {
        bail!(
            "Assess probes R script not found: {}",
            script.display()
        );
    }

    let threshold_str = format!("{:.1}", threshold);
    let proximity_str = proximity.to_string();

    let mut r_args: Vec<String> = Vec::new();

    // Required args
    r_args.extend(["--xreact-hits".into(), abs_path_str(xreact_hits)?]);
    r_args.extend(["--xreact-summary".into(), abs_path_str(xreact_summary)?]);
    r_args.extend(["--threshold".into(), threshold_str]);
    r_args.extend(["--cov-summary".into(), abs_path_str(cov_summary)?]);
    r_args.extend(["--cov-depth".into(), abs_path_str(cov_depth)?]);
    r_args.extend(["--proximity".into(), proximity_str]);
    r_args.extend(["--params".into(), abs_path_str(params_path)?]);
    r_args.extend(["--output".into(), abs_path_str(output_path)?]);

    // Optional build-probes data
    if let Some(p) = build_stats_file {
        if p.exists() {
            r_args.extend(["--build-stats".into(), abs_path_str(p)?]);
        }
    }
    if let Some(p) = build_params_file {
        if p.exists() {
            r_args.extend(["--build-params".into(), abs_path_str(p)?]);
        }
    }
    if cov_multi.exists() {
        r_args.extend(["--cov-multi-mapping".into(), abs_path_str(cov_multi)?]);
    }

    let arg_refs: Vec<&str> = r_args.iter().map(|s| s.as_str()).collect();
    rscript::run_rscript(&script, &arg_refs)
}

fn write_assess_rmd(
    build_stats_file: Option<&Path>,
    build_params_file: Option<&Path>,
    xreact_hits: &Path,
    xreact_summary: &Path,
    threshold: f64,
    cov_summary: &Path,
    cov_depth: &Path,
    cov_multi: &Path,
    proximity: usize,
    params_path: &Path,
    output_path: &Path,
) -> Result<()> {
    let r_dir = rscript::find_r_dir()
        .ok_or_else(|| anyhow::anyhow!("Cannot find R scripts directory."))?;

    let rmd_template = r_dir.join("assess_probes.Rmd");
    if !rmd_template.exists() {
        bail!(
            "RMarkdown template not found: {}",
            rmd_template.display()
        );
    }

    let threshold_str = format!("{:.1}", threshold);
    let proximity_str = proximity.to_string();

    let build_stats_str = build_stats_file
        .and_then(|p| if p.exists() { abs_path_str(p).ok() } else { None })
        .unwrap_or_default();
    let build_params_str = build_params_file
        .and_then(|p| if p.exists() { abs_path_str(p).ok() } else { None })
        .unwrap_or_default();
    let cov_multi_str = if cov_multi.exists() {
        abs_path_str(cov_multi)?
    } else {
        String::new()
    };
    let xreact_hits_str = abs_path_str(xreact_hits)?;
    let xreact_summary_str = abs_path_str(xreact_summary)?;
    let cov_summary_str = abs_path_str(cov_summary)?;
    let cov_depth_str = abs_path_str(cov_depth)?;
    let params_path_str = abs_path_str(params_path)?;

    let params = vec![
        ("build_stats_file", build_stats_str.as_str()),
        ("build_params_file", build_params_str.as_str()),
        ("xreact_hits_file", xreact_hits_str.as_str()),
        ("xreact_summary_file", xreact_summary_str.as_str()),
        ("xreact_threshold", &threshold_str),
        ("coverage_summary_file", cov_summary_str.as_str()),
        ("coverage_depth_file", cov_depth_str.as_str()),
        ("coverage_multi_mapping_file", cov_multi_str.as_str()),
        ("coverage_proximity", &proximity_str),
        ("params_file", params_path_str.as_str()),
    ];

    let template_content = std::fs::read_to_string(&rmd_template)
        .with_context(|| format!("Failed to read template: {}", rmd_template.display()))?;

    let output_content = substitute_rmd_params(&template_content, &params);

    let rmd_path = rmd_output_path(output_path);
    std::fs::write(&rmd_path, output_content)
        .with_context(|| format!("Failed to write RMarkdown: {}", rmd_path.display()))?;

    log::info!("RMarkdown file written: {}", rmd_path.display());
    log::info!(
        "Edit and render with: Rscript -e 'rmarkdown::render(\"{}\")'",
        rmd_path.display()
    );
    Ok(())
}

/// Parse a probe_coverage_summary.tsv and return IDs where pct_covered_1x < threshold.
fn read_low_coverage_targets(summary_path: &Path, threshold: f64) -> Result<HashSet<String>> {
    let file = File::open(summary_path)
        .with_context(|| format!("Cannot open coverage summary: {}", summary_path.display()))?;
    let reader = BufReader::new(file);
    let mut low: HashSet<String> = HashSet::new();

    for (i, line) in reader.lines().enumerate() {
        let line = line?;
        if i == 0 {
            continue; // skip header
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 4 {
            continue;
        }
        let pct: f64 = fields[3].parse().unwrap_or(100.0);
        if pct < threshold {
            low.insert(fields[0].to_string());
        }
    }

    Ok(low)
}

/// Run one or more probe-coverage-only iterations on low-coverage targets.
fn run_refinement(args: &AssessProbesArgs, initial_summary: &Path) -> Result<()> {
    let pfx = args.output_prefix;
    let max_iterations = args.refine_iterations.unwrap_or(usize::MAX);

    let mut current_summary = initial_summary.to_path_buf();
    let mut prev_ids: HashSet<String> = HashSet::new();

    for iteration in 1..=max_iterations {
        let low_ids = read_low_coverage_targets(&current_summary, args.refine_threshold)?;

        if low_ids.is_empty() {
            log::info!(
                "Refinement: no targets below {:.1}% 1X coverage — stopping.",
                args.refine_threshold
            );
            break;
        }

        // Stable-check: stop if the low-coverage set hasn't changed since last iteration
        if args.refine_until_stable && low_ids == prev_ids {
            log::info!(
                "Refinement: target set unchanged after iteration {} — stopping.",
                iteration - 1
            );
            break;
        }

        log::info!(
            "Refinement iteration {}: {} targets below {:.1}% 1X coverage",
            iteration,
            low_ids.len(),
            args.refine_threshold
        );

        // Write filtered targets FASTA
        let filtered_targets_path =
            prefixed_join(args.outdir, pfx, &format!("refine_{}_targets.fa", iteration));
        let n_extracted = fasta::extract_by_ids(args.targets, &low_ids, &filtered_targets_path)?;
        log::info!(
            "Refinement iteration {}: wrote {} targets to {}",
            iteration,
            n_extracted,
            filtered_targets_path.display()
        );

        // Run probe coverage on filtered targets
        let cov_prefix = format!("{}refine_{}_cov_", pfx, iteration);
        probe_coverage::execute(&probe_coverage::ProbeCoverageArgs {
            targets: &filtered_targets_path,
            probes: args.probes,
            outdir: args.outdir,
            output_prefix: &cov_prefix,
            minimap_preset: args.minimap_preset,
            proximity: args.proximity,
            report: ReportMode::None,
            cleanup: false,
        })?;

        let iter_summary = prefixed_join(args.outdir, &cov_prefix, "probe_coverage_summary.tsv");
        let iter_depth = prefixed_join(args.outdir, &cov_prefix, "probe_depth.tsv");
        let iter_multi = prefixed_join(args.outdir, &cov_prefix, "multi_mapping_probes.tsv");
        let iter_params = prefixed_join(args.outdir, &cov_prefix, "run_params.tsv");

        // Generate probe-coverage report for this iteration
        let report_path = prefixed_join(
            args.outdir,
            pfx,
            &format!("refine_{}_probe_coverage_report.html", iteration),
        );
        match args.report {
            ReportMode::None => {
                log::info!("Skipping refinement report (--report none)");
            }
            ReportMode::Full => {
                if rscript::check_available() {
                    match probe_coverage::generate_probe_report(
                        &iter_summary,
                        &iter_depth,
                        &iter_multi,
                        &iter_params,
                        &report_path,
                        args.proximity,
                    ) {
                        Ok(()) => log::info!(
                            "Refinement iteration {} report: {}",
                            iteration,
                            report_path.display()
                        ),
                        Err(e) => log::warn!(
                            "Refinement iteration {} report failed (non-fatal): {}",
                            iteration,
                            e
                        ),
                    }
                } else {
                    log::warn!("Rscript not found -- skipping refinement report.");
                }
            }
            ReportMode::Rmd => {
                match probe_coverage::write_probe_coverage_rmd(
                    &iter_summary,
                    &iter_depth,
                    &iter_multi,
                    &iter_params,
                    &report_path,
                    args.proximity,
                ) {
                    Ok(()) => {}
                    Err(e) => log::warn!(
                        "Refinement iteration {} RMarkdown failed (non-fatal): {}",
                        iteration,
                        e
                    ),
                }
            }
        }

        // Cleanup SAM/log intermediates for this iteration if requested
        if args.cleanup {
            let intermediates: Vec<String> = ["probe_alignment.sam", "probe_alignment.log"]
                .iter()
                .map(|f| format!("{}{}", cov_prefix, f))
                .collect();
            let refs: Vec<&str> = intermediates.iter().map(|s| s.as_str()).collect();
            cleanup::cleanup_files(args.outdir, &refs);
            // Also remove the filtered targets FASTA
            if let Err(e) = std::fs::remove_file(&filtered_targets_path) {
                log::warn!("Could not remove {}: {}", filtered_targets_path.display(), e);
            }
        }

        prev_ids = low_ids;
        current_summary = iter_summary;

        // For --refine-until-stable, the loop runs until break conditions above
        // For --refine-iterations, the for-loop bound handles stopping
    }

    Ok(())
}
