use anyhow::{bail, Context, Result};
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::Path;

use crate::cli::ReportMode;
use crate::commands::report::{rmd_output_path, substitute_rmd_params};
use crate::external::rscript;
use crate::io_utils;
use crate::io_utils::prefixed_join;
use crate::target_similarity;

pub struct PanelQcArgs<'a> {
    pub targets: &'a Path,
    pub sample_target_map: &'a Path,
    pub identity_threshold: f64,
    pub minimap_preset: &'a str,
    pub outdir: &'a Path,
    pub output_prefix: &'a str,
    pub report: ReportMode,
    pub cleanup: bool,
}

pub fn execute(args: &PanelQcArgs) -> Result<()> {
    if !args.targets.exists() {
        bail!("Targets file not found: {}", args.targets.display());
    }
    if !args.sample_target_map.exists() {
        bail!(
            "Sample-target-map not found: {}",
            args.sample_target_map.display()
        );
    }

    fs::create_dir_all(args.outdir)?;

    log::info!("=============================================");
    log::info!("BaitBench - Target Panel Discriminability QC");
    log::info!("=============================================");
    log::info!("Targets            : {}", args.targets.display());
    log::info!("Sample-target-map  : {}", args.sample_target_map.display());
    log::info!("Identity threshold : {:.1}%", args.identity_threshold);
    log::info!("Minimap preset     : {}", args.minimap_preset);
    log::info!("Output             : {}", args.outdir.display());

    // Parse genome-to-target mapping
    let genome_to_targets = io_utils::parse_sample_target_map(args.sample_target_map)?;
    let n_species = genome_to_targets.len();
    let n_targets: usize = genome_to_targets.values().map(|v| v.len()).sum();
    log::info!("Species (genomes): {}", n_species);
    log::info!("Total targets: {}", n_targets);

    // Step 1: Compute target-vs-target similarity
    let similarities = target_similarity::compute_target_similarity(
        args.targets,
        args.minimap_preset,
        args.identity_threshold,
        args.outdir,
    )?;

    let pfx = args.output_prefix;

    // Write target similarity TSV
    let sim_path = prefixed_join(args.outdir, pfx, "target_similarity.tsv");
    target_similarity::write_target_similarity(&similarities, &sim_path)?;
    log::info!(
        "Target similarity pairs: {} (written to {})",
        similarities.len(),
        sim_path.display()
    );

    // Step 2: Build context and compute discriminability
    let ctx =
        target_similarity::build_similarity_context(&similarities, &genome_to_targets);
    let discriminability = target_similarity::compute_discriminability(&ctx);

    let disc_path = prefixed_join(args.outdir, pfx, "species_discriminability.tsv");
    target_similarity::write_discriminability(&discriminability, &disc_path)?;

    // Step 3: Build and write confusion matrix
    let (species_ids, matrix) = target_similarity::build_confusion_matrix(&ctx);
    let matrix_path = prefixed_join(args.outdir, pfx, "species_confusion_matrix.tsv");
    target_similarity::write_confusion_matrix(&species_ids, &matrix, &matrix_path)?;

    // Console summary
    let fully_unique = discriminability
        .iter()
        .filter(|d| d.discriminability_score >= 1.0)
        .count();
    let partially_unique = discriminability
        .iter()
        .filter(|d| d.discriminability_score > 0.0 && d.discriminability_score < 1.0)
        .count();
    let no_unique = discriminability
        .iter()
        .filter(|d| d.discriminability_score == 0.0)
        .count();

    log::info!("Species with all unique targets:  {}/{}", fully_unique, n_species);
    log::info!("Species with some unique targets: {}/{}", partially_unique, n_species);
    log::info!(
        "Species with NO unique targets:  {}/{} (cannot discriminate!)",
        no_unique,
        n_species
    );

    if no_unique > 0 {
        log::warn!(
            "{} species have no unique targets and cannot be reliably distinguished. \
             Consider adding more species-specific targets.",
            no_unique
        );
        for d in &discriminability {
            if d.discriminability_score == 0.0 {
                let confusable: Vec<&str> = d
                    .confusable_with
                    .iter()
                    .map(|(s, _)| s.as_str())
                    .collect();
                log::warn!(
                    "  {} ({} targets, all shared) — confusable with: {}",
                    d.species_id,
                    d.total_targets,
                    confusable.join(", ")
                );
            }
        }
    }

    // Write run_params.tsv
    let params_path = prefixed_join(args.outdir, pfx, "run_params.tsv");
    write_run_params(&params_path, args)?;

    // Report generation
    match args.report {
        ReportMode::None => {
            log::info!("Skipping report generation (--report none)");
        }
        ReportMode::Full => {
            if rscript::check_available() {
                let report_path = prefixed_join(args.outdir, pfx, "panel_qc_report.html");
                log::info!("Generating panel QC report...");
                match generate_panel_qc_report(
                    &disc_path,
                    &matrix_path,
                    &sim_path,
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
            let report_path = prefixed_join(args.outdir, pfx, "panel_qc_report.html");
            log::info!("Generating panel QC RMarkdown file...");
            match write_panel_qc_rmd(
                &disc_path,
                &matrix_path,
                &sim_path,
                &params_path,
                &report_path,
            ) {
                Ok(()) => {}
                Err(e) => log::warn!("RMarkdown generation failed (non-fatal): {}", e),
            }
        }
    }

    // Cleanup
    if args.cleanup {
        log::info!("Cleaning up intermediate files...");
        crate::cleanup::cleanup_files(
            args.outdir,
            &["target_vs_target.paf", "target_vs_target.log"],
        );
    }

    log::info!("=============================================");
    log::info!("Panel QC analysis complete!");
    log::info!("Results in {}", args.outdir.display());
    log::info!("=============================================");

    Ok(())
}

fn write_run_params(path: &Path, args: &PanelQcArgs) -> Result<()> {
    let file = File::create(path)
        .with_context(|| format!("Cannot create params file: {}", path.display()))?;
    let mut w = BufWriter::new(file);

    writeln!(w, "parameter\tflag\tvalue")?;
    writeln!(w, "targets\t--targets\t{}", args.targets.display())?;
    writeln!(
        w,
        "sample_target_map\t--sample-target-map\t{}",
        args.sample_target_map.display()
    )?;
    writeln!(
        w,
        "identity_threshold\t--identity-threshold\t{:.1}",
        args.identity_threshold
    )?;
    writeln!(
        w,
        "minimap_preset\t--minimap-preset\t{}",
        args.minimap_preset
    )?;
    writeln!(w, "outdir\t-o\t{}", args.outdir.display())?;

    w.flush()?;
    Ok(())
}

fn generate_panel_qc_report(
    disc_path: &Path,
    matrix_path: &Path,
    sim_path: &Path,
    params_path: &Path,
    output_path: &Path,
) -> Result<()> {
    let r_dir = rscript::find_r_dir()
        .ok_or_else(|| anyhow::anyhow!("Cannot find R scripts directory."))?;

    let script = r_dir.join("panel_qc.R");
    if !script.exists() {
        bail!("Panel QC R script not found: {}", script.display());
    }

    let disc_abs = std::fs::canonicalize(disc_path)?;
    let matrix_abs = std::fs::canonicalize(matrix_path)?;
    let sim_abs = std::fs::canonicalize(sim_path)?;
    let params_abs = std::fs::canonicalize(params_path)?;
    let output_abs = if output_path.is_absolute() {
        output_path.to_path_buf()
    } else {
        std::env::current_dir()?.join(output_path)
    };

    rscript::run_rscript(
        &script,
        &[
            "--discriminability",
            disc_abs.to_str().unwrap_or(""),
            "--matrix",
            matrix_abs.to_str().unwrap_or(""),
            "--similarity",
            sim_abs.to_str().unwrap_or(""),
            "--params",
            params_abs.to_str().unwrap_or(""),
            "--output",
            output_abs.to_str().unwrap_or(""),
        ],
    )
}

fn write_panel_qc_rmd(
    disc_path: &Path,
    matrix_path: &Path,
    sim_path: &Path,
    params_path: &Path,
    output_path: &Path,
) -> Result<()> {
    let r_dir = rscript::find_r_dir()
        .ok_or_else(|| anyhow::anyhow!("Cannot find R scripts directory."))?;

    let rmd_template = r_dir.join("panel_qc.Rmd");
    if !rmd_template.exists() {
        bail!(
            "RMarkdown template not found: {}",
            rmd_template.display()
        );
    }

    let disc_abs = std::fs::canonicalize(disc_path)?;
    let matrix_abs = std::fs::canonicalize(matrix_path)?;
    let sim_abs = std::fs::canonicalize(sim_path)?;
    let params_abs = std::fs::canonicalize(params_path)?;

    let params = vec![
        ("discriminability_file", disc_abs.to_str().unwrap_or("")),
        ("matrix_file", matrix_abs.to_str().unwrap_or("")),
        ("similarity_file", sim_abs.to_str().unwrap_or("")),
        ("params_file", params_abs.to_str().unwrap_or("")),
    ];

    let template_content = std::fs::read_to_string(&rmd_template)
        .with_context(|| format!("Failed to read template: {}", rmd_template.display()))?;

    let output_content = substitute_rmd_params(&template_content, &params);

    let rmd_path = rmd_output_path(output_path);
    std::fs::write(&rmd_path, output_content)
        .with_context(|| format!("Failed to write RMarkdown file: {}", rmd_path.display()))?;

    log::info!("RMarkdown file written: {}", rmd_path.display());
    log::info!(
        "Edit and render with: Rscript -e 'rmarkdown::render(\"{}\")'",
        rmd_path.display()
    );
    Ok(())
}
