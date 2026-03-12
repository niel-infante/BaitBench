use anyhow::{Context, Result, bail};
use std::collections::{HashMap, HashSet};
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use crate::alignment::coverage;
use crate::cli::ReportMode;
use crate::commands::report::{substitute_rmd_params, rmd_output_path};
use crate::external::{minimap2, rscript};
use crate::fasta;

pub struct ProbeCoverageArgs<'a> {
    pub targets: &'a Path,
    pub probes: &'a Path,
    pub outdir: &'a Path,
    pub minimap_preset: &'a str,
    pub proximity: usize,
    pub report: ReportMode,
}

pub fn execute(args: &ProbeCoverageArgs) -> Result<()> {
    if !args.targets.exists() {
        bail!("Targets file not found: {}", args.targets.display());
    }
    if !args.probes.exists() {
        bail!("Probes file not found: {}", args.probes.display());
    }

    fs::create_dir_all(args.outdir)?;

    minimap2::check_available()?;

    log::info!("=============================================");
    log::info!("BaitBench - Probe Coverage Analysis");
    log::info!("=============================================");
    log::info!("Targets  : {}", args.targets.display());
    log::info!("Probes   : {}", args.probes.display());
    log::info!("Preset   : {}", args.minimap_preset);
    log::info!("Proximity: {} bp", args.proximity);
    log::info!("Output   : {}", args.outdir.display());

    let n_targets = fasta::count_sequences(args.targets)?;
    let n_probes = fasta::count_sequences(args.probes)?;
    log::info!("Target sequences: {}", n_targets);
    log::info!("Probe sequences : {}", n_probes);

    // Step 1: Align probes to targets
    let sam_path = args.outdir.join("probe_alignment.sam");
    let log_path = args.outdir.join("probe_alignment.log");
    log::info!("Aligning probes to targets...");
    minimap2::probe_align(
        args.minimap_preset,
        args.targets,
        args.probes,
        &sam_path,
        &log_path,
    )?;

    // Step 2: Compute per-position probe depth (includes secondary alignments)
    log::info!("Computing per-position probe depth...");
    let coverage_result = coverage::compute_probe_coverage(&sam_path)?;

    // Step 3: Calculate per-target probe coverage statistics
    log::info!("Calculating probe coverage statistics...");
    let mut stats: Vec<(String, coverage::ProbeCoverageStats)> = Vec::new();
    for (ref_id, depths) in &coverage_result.coverage {
        let ref_len = coverage_result
            .ref_lengths
            .get(ref_id)
            .copied()
            .unwrap_or(depths.len());
        let s = coverage::calculate_probe_stats(depths, ref_len, args.proximity);
        stats.push((ref_id.clone(), s));
    }
    stats.sort_by(|a, b| a.0.cmp(&b.0));

    // Step 4: Write per-position depth TSV
    let depth_path = args.outdir.join("probe_depth.tsv");
    log::info!("Writing probe depth to {}", depth_path.display());
    coverage::write_coverage_tsv(&depth_path, &coverage_result.coverage)?;

    // Step 5: Write summary statistics TSV
    let summary_path = args.outdir.join("probe_coverage_summary.tsv");
    log::info!(
        "Writing probe coverage summary to {}",
        summary_path.display()
    );
    write_probe_summary(&summary_path, &stats)?;

    // Step 6: Identify probes mapping to multiple targets
    log::info!("Identifying multi-mapping probes...");
    let multi_mapping = find_multi_mapping_probes(&sam_path)?;
    let multi_mapping_path = args.outdir.join("multi_mapping_probes.tsv");
    write_multi_mapping_probes(&multi_mapping_path, &multi_mapping)?;
    log::info!(
        "Multi-mapping probes: {} of {} probes map to multiple targets",
        multi_mapping.len(),
        n_probes
    );

    // Log summary
    let total_targets = stats.len();
    let fully_covered = stats
        .iter()
        .filter(|(_, s)| s.pct_covered_1x >= 100.0)
        .count();
    let well_covered = stats
        .iter()
        .filter(|(_, s)| s.pct_covered_1x >= 90.0)
        .count();
    log::info!("Targets analyzed   : {}", total_targets);
    log::info!("Targets 100% tiled : {}", fully_covered);
    log::info!("Targets >=90% tiled: {}", well_covered);

    // Write run_params.tsv
    let params_path = args.outdir.join("run_params.tsv");
    write_run_params(&params_path, args)?;

    // Step 7: Report
    match args.report {
        ReportMode::None => {
            log::info!("Skipping report generation (--report none)");
        }
        ReportMode::Full => {
            if rscript::check_available() {
                let report_path = args.outdir.join("probe_coverage_report.html");
                log::info!("Generating probe coverage report...");
                match generate_probe_report(&summary_path, &depth_path, &multi_mapping_path, &params_path, &report_path, args.proximity) {
                    Ok(()) => log::info!("Report generated: {}", report_path.display()),
                    Err(e) => log::warn!("Report generation failed (non-fatal): {}", e),
                }
            } else {
                log::warn!("Rscript not found -- skipping HTML report.");
            }
        }
        ReportMode::Rmd => {
            let report_path = args.outdir.join("probe_coverage_report.html");
            log::info!("Generating probe coverage RMarkdown file...");
            match write_probe_coverage_rmd(&summary_path, &depth_path, &multi_mapping_path, &params_path, &report_path, args.proximity) {
                Ok(()) => {}
                Err(e) => log::warn!("RMarkdown generation failed (non-fatal): {}", e),
            }
        }
    }

    log::info!("=============================================");
    log::info!("Probe coverage analysis complete!");
    log::info!("Results in {}", args.outdir.display());
    log::info!("=============================================");

    Ok(())
}

fn write_probe_summary(
    path: &Path,
    stats: &[(String, coverage::ProbeCoverageStats)],
) -> Result<()> {
    let file = File::create(path)
        .with_context(|| format!("Cannot create summary: {}", path.display()))?;
    let mut w = BufWriter::new(file);

    writeln!(
        w,
        "reference_id\tref_length\tcovered_bases\tpct_covered_1x\tmean_depth\tmedian_depth\tpct_covered_2x\tpct_covered_5x\tpct_covered_10x\tmax_gap_length\tnum_gaps\tpct_near_probe"
    )?;

    for (ref_id, s) in stats {
        writeln!(
            w,
            "{}\t{}\t{}\t{:.1}\t{:.2}\t{:.1}\t{:.1}\t{:.1}\t{:.1}\t{}\t{}\t{:.1}",
            ref_id,
            s.ref_length,
            s.covered_bases,
            s.pct_covered_1x,
            s.mean_depth,
            s.median_depth,
            s.pct_covered_2x,
            s.pct_covered_5x,
            s.pct_covered_10x,
            s.max_gap_length,
            s.num_gaps,
            s.pct_near_probe,
        )?;
    }

    w.flush()?;
    Ok(())
}

/// Parse SAM to find probes that map to more than one target.
///
/// Returns a sorted list of (probe_id, sorted Vec<target_ids>).
fn find_multi_mapping_probes(
    sam_path: &Path,
) -> Result<Vec<(String, Vec<String>)>> {
    let file = File::open(sam_path)
        .with_context(|| format!("Cannot open SAM: {}", sam_path.display()))?;
    let reader = BufReader::new(file);

    let mut probe_targets: HashMap<String, HashSet<String>> = HashMap::new();

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('@') {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 3 {
            continue;
        }

        let flag: u16 = fields[1].parse().unwrap_or(0);
        let rname = fields[2];

        // Skip unmapped
        if rname == "*" || flag & 0x4 != 0 {
            continue;
        }

        let probe_id = fields[0];
        probe_targets
            .entry(probe_id.to_string())
            .or_default()
            .insert(rname.to_string());
    }

    let mut multi: Vec<(String, Vec<String>)> = probe_targets
        .into_iter()
        .filter(|(_, targets)| targets.len() > 1)
        .map(|(probe, targets)| {
            let mut t: Vec<String> = targets.into_iter().collect();
            t.sort();
            (probe, t)
        })
        .collect();
    multi.sort_by(|a, b| a.0.cmp(&b.0));

    Ok(multi)
}

fn write_multi_mapping_probes(
    path: &Path,
    probes: &[(String, Vec<String>)],
) -> Result<()> {
    let file = File::create(path)
        .with_context(|| format!("Cannot create multi-mapping file: {}", path.display()))?;
    let mut w = BufWriter::new(file);

    writeln!(w, "probe_id\tnum_targets\ttargets")?;
    for (probe_id, targets) in probes {
        writeln!(w, "{}\t{}\t{}", probe_id, targets.len(), targets.join(","))?;
    }

    w.flush()?;
    Ok(())
}

fn write_run_params(path: &Path, args: &ProbeCoverageArgs) -> Result<()> {
    let file = File::create(path)
        .with_context(|| format!("Cannot create params file: {}", path.display()))?;
    let mut w = BufWriter::new(file);

    writeln!(w, "parameter\tflag\tvalue")?;
    writeln!(w, "targets\t--targets\t{}", args.targets.display())?;
    writeln!(w, "probes\t--probes\t{}", args.probes.display())?;
    writeln!(w, "minimap_preset\t--minimap-preset\t{}", args.minimap_preset)?;
    writeln!(w, "proximity\t--proximity\t{}", args.proximity)?;
    writeln!(w, "outdir\t-o\t{}", args.outdir.display())?;

    w.flush()?;
    Ok(())
}

fn write_probe_coverage_rmd(
    summary_path: &Path,
    depth_path: &Path,
    multi_mapping_path: &Path,
    params_path: &Path,
    output_path: &Path,
    proximity: usize,
) -> Result<()> {
    let r_dir = rscript::find_r_dir()
        .ok_or_else(|| anyhow::anyhow!("Cannot find R scripts directory."))?;

    let rmd_template = r_dir.join("probe_coverage.Rmd");
    if !rmd_template.exists() {
        bail!("RMarkdown template not found: {}", rmd_template.display());
    }

    let summary_abs = std::fs::canonicalize(summary_path)?;
    let depth_abs = std::fs::canonicalize(depth_path)?;
    let multi_abs = std::fs::canonicalize(multi_mapping_path)?;
    let params_abs = std::fs::canonicalize(params_path)?;
    let proximity_str = proximity.to_string();

    let params = vec![
        ("summary_file", summary_abs.to_str().unwrap_or("")),
        ("depth_file", depth_abs.to_str().unwrap_or("")),
        ("multi_mapping_file", multi_abs.to_str().unwrap_or("")),
        ("params_file", params_abs.to_str().unwrap_or("")),
        ("proximity", &proximity_str),
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

fn generate_probe_report(
    summary_path: &Path,
    depth_path: &Path,
    multi_mapping_path: &Path,
    params_path: &Path,
    output_path: &Path,
    proximity: usize,
) -> Result<()> {
    let r_dir = rscript::find_r_dir()
        .ok_or_else(|| anyhow::anyhow!("Cannot find R scripts directory."))?;

    let script = r_dir.join("probe_coverage.R");
    if !script.exists() {
        bail!(
            "Probe coverage R script not found: {}",
            script.display()
        );
    }

    let summary_abs = std::fs::canonicalize(summary_path)?;
    let depth_abs = std::fs::canonicalize(depth_path)?;
    let multi_abs = std::fs::canonicalize(multi_mapping_path)?;
    let params_abs = std::fs::canonicalize(params_path)?;
    let output_abs = if output_path.is_absolute() {
        output_path.to_path_buf()
    } else {
        std::env::current_dir()?.join(output_path)
    };

    let summary_str = summary_abs.to_str().unwrap_or("");
    let depth_str = depth_abs.to_str().unwrap_or("");
    let multi_str = multi_abs.to_str().unwrap_or("");
    let params_str = params_abs.to_str().unwrap_or("");
    let output_str = output_abs.to_str().unwrap_or("");
    let proximity_str = proximity.to_string();

    rscript::run_rscript(
        &script,
        &[
            "--summary",
            summary_str,
            "--depth",
            depth_str,
            "--multi-mapping",
            multi_str,
            "--params",
            params_str,
            "--proximity",
            &proximity_str,
            "--output",
            output_str,
        ],
    )
}
