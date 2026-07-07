use anyhow::{Context, Result, bail};
use std::collections::{HashMap, HashSet};
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

use crate::alignment::coverage;
use crate::cleanup;
use crate::cli::ReportMode;
use crate::commands::report::{substitute_rmd_params, rmd_output_path};
use crate::external::{minimap2, rscript};
use crate::fasta;
use crate::io_utils::prefixed_join;
use crate::sdust;

pub struct ProbeCoverageArgs<'a> {
    pub targets: &'a Path,
    pub probes: &'a Path,
    pub outdir: &'a Path,
    pub output_prefix: &'a str,
    pub minimap_preset: &'a str,
    pub proximity: usize,
    pub report: ReportMode,
    pub cleanup: bool,
}

/// Data returned from a completed probe coverage analysis run.
pub struct CoverageRunData {
    pub coverage: HashMap<String, Vec<u32>>,
    pub ref_lengths: HashMap<String, usize>,
    pub summary_path: PathBuf,
    pub depth_path: PathBuf,
    pub multi_mapping_path: PathBuf,
    pub params_path: PathBuf,
}

/// Result from individual per-target coverage computation.
pub struct IndividualCoverageData {
    pub summary_path: PathBuf,
    pub coverage: HashMap<String, Vec<u32>>,
}

/// Run probe coverage analysis and return all data.
/// Does not generate a report. Use execute() for standalone use with report generation.
pub fn run_probe_coverage(args: &ProbeCoverageArgs) -> Result<CoverageRunData> {
    if !args.targets.exists() {
        bail!("Targets file not found: {}", args.targets.display());
    }
    if !args.probes.exists() {
        bail!("Probes file not found: {}", args.probes.display());
    }

    fs::create_dir_all(args.outdir)?;
    minimap2::check_available()?;

    let n_probes = fasta::count_sequences(args.probes)?;
    let pfx = args.output_prefix;

    let sam_path = prefixed_join(args.outdir, pfx, "probe_alignment.sam");
    let log_path = prefixed_join(args.outdir, pfx, "probe_alignment.log");
    log::info!("Aligning probes to targets...");
    minimap2::probe_align(
        args.minimap_preset,
        args.targets,
        args.probes,
        &sam_path,
        &log_path,
        1,
        1000,
    )?;

    log::info!("Computing per-position probe depth...");
    let coverage_result = coverage::compute_probe_coverage(&sam_path)?;

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

    let depth_path = prefixed_join(args.outdir, pfx, "probe_depth.tsv");
    coverage::write_coverage_intervals(&depth_path, &coverage_result.coverage)?;

    let summary_path = prefixed_join(args.outdir, pfx, "probe_coverage_summary.tsv");
    write_probe_summary(&summary_path, &stats)?;

    log::info!("Identifying multi-mapping probes...");
    let multi_mapping = find_multi_mapping_probes(&sam_path)?;
    let multi_mapping_path = prefixed_join(args.outdir, pfx, "multi_mapping_probes.tsv");
    write_multi_mapping_probes(&multi_mapping_path, &multi_mapping)?;
    log::info!(
        "Multi-mapping probes: {} of {} probes map to multiple targets",
        multi_mapping.len(),
        n_probes
    );

    let total_targets = stats.len();
    let fully_covered = stats.iter().filter(|(_, s)| s.pct_covered_1x >= 100.0).count();
    let well_covered = stats.iter().filter(|(_, s)| s.pct_covered_1x >= 90.0).count();
    log::info!("Targets analyzed   : {}", total_targets);
    log::info!("Targets 100% tiled : {}", fully_covered);
    log::info!("Targets >=90% tiled: {}", well_covered);

    let params_path = prefixed_join(args.outdir, pfx, "run_params.tsv");
    write_run_params(&params_path, args)?;

    if args.cleanup {
        let cleanup_names: Vec<String> = ["probe_alignment.sam", "probe_alignment.log"]
            .iter()
            .map(|f| format!("{}{}", pfx, f))
            .collect();
        let cleanup_refs: Vec<&str> = cleanup_names.iter().map(|s| s.as_str()).collect();
        cleanup::cleanup_files(args.outdir, &cleanup_refs);
    }

    Ok(CoverageRunData {
        coverage: coverage_result.coverage,
        ref_lengths: coverage_result.ref_lengths,
        summary_path,
        depth_path,
        multi_mapping_path,
        params_path,
    })
}

pub fn execute(args: &ProbeCoverageArgs) -> Result<()> {
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

    let data = run_probe_coverage(args)?;
    let pfx = args.output_prefix;

    match args.report {
        ReportMode::None => {
            log::info!("Skipping report generation (--report none)");
        }
        ReportMode::Full => {
            if rscript::check_available() {
                let report_path = prefixed_join(args.outdir, pfx, "probe_coverage_report.html");
                log::info!("Generating probe coverage report...");
                match generate_probe_report(&data.summary_path, &data.depth_path, &data.multi_mapping_path, &data.params_path, &report_path, args.proximity) {
                    Ok(()) => log::info!("Report generated: {}", report_path.display()),
                    Err(e) => log::warn!("Report generation failed (non-fatal): {}", e),
                }
            } else {
                log::warn!("Rscript not found -- skipping HTML report.");
            }
        }
        ReportMode::Rmd => {
            let report_path = prefixed_join(args.outdir, pfx, "probe_coverage_report.html");
            log::info!("Generating probe coverage RMarkdown file...");
            match write_probe_coverage_rmd(&data.summary_path, &data.depth_path, &data.multi_mapping_path, &data.params_path, &report_path, args.proximity) {
                Ok(()) => {}
                Err(e) => log::warn!("RMarkdown generation failed (non-fatal): {}", e),
            }
        }
        ReportMode::BothR => {
            let report_path = prefixed_join(args.outdir, pfx, "probe_coverage_report.html");
            log::info!("Generating probe coverage RMarkdown file...");
            match write_probe_coverage_rmd(&data.summary_path, &data.depth_path, &data.multi_mapping_path, &data.params_path, &report_path, args.proximity) {
                Ok(()) => {}
                Err(e) => log::warn!("RMarkdown generation failed (non-fatal): {}", e),
            }
            if rscript::check_available() {
                log::info!("Generating probe coverage HTML report...");
                match generate_probe_report(&data.summary_path, &data.depth_path, &data.multi_mapping_path, &data.params_path, &report_path, args.proximity) {
                    Ok(()) => log::info!("Report generated: {}", report_path.display()),
                    Err(e) => log::warn!("Report generation failed (non-fatal): {}", e),
                }
            } else {
                log::warn!("Rscript not found — skipping HTML report (Rmd still written).");
            }
        }
    }

    log::info!("=============================================");
    log::info!("Probe coverage analysis complete!");
    log::info!("Results in {}", args.outdir.display());
    log::info!("=============================================");

    Ok(())
}

pub fn write_probe_summary(
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

pub fn write_probe_coverage_rmd(
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

pub fn generate_probe_report(
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

/// Compute per-target probe coverage individually (no probe competition).
///
/// Pre-loads probe sequences once, then for each target writes a temp single-target
/// FASTA and computes coverage in memory (no SAM file per target).
/// Returns the summary TSV path and per-target depth vectors for downstream gap analysis.
pub fn run_individual_coverage(
    targets_path: &Path,
    probes_path: &Path,
    outdir: &Path,
    output_prefix: &str,
    minimap_preset: &str,
    proximity: usize,
    threads: usize,
) -> Result<IndividualCoverageData> {
    let target_ids = fasta::parse_fasta_ids(targets_path)?;
    let n = target_ids.len();
    log::info!("Computing individual coverage for {} targets...", n);

    let probe_seqs = fasta::parse_fasta(probes_path)?;

    let tmp_dir = prefixed_join(outdir, output_prefix, "indiv_tmp");
    fs::create_dir_all(&tmp_dir)?;

    let mut stats: Vec<(String, coverage::ProbeCoverageStats)> = Vec::new();
    let mut all_depths: HashMap<String, Vec<u32>> = HashMap::new();

    let report_interval = (n / 10).max(1);

    for (i, target_id) in target_ids.iter().enumerate() {
        if i > 0 && i % report_interval == 0 {
            log::info!(
                "Individual coverage: {}/{} ({:.0}%)",
                i, n,
                100.0 * i as f64 / n as f64
            );
        }
        log::debug!("Individual coverage [{}/{}]: {}", i + 1, n, target_id);

        let target_fa = tmp_dir.join(format!("{}.fa", i));
        let mut id_set = HashSet::new();
        id_set.insert(target_id.clone());
        let extracted = fasta::extract_by_ids(targets_path, &id_set, &target_fa)?;
        if extracted == 0 {
            log::warn!("Target '{}' not found — skipping", target_id);
            continue;
        }

        let (ref_lengths, cov_map) =
            minimap2::probe_depth_in_memory(minimap_preset, &target_fa, &probe_seqs, 1000, threads)?;

        let _ = fs::remove_file(&target_fa);

        let ref_len = ref_lengths.get(target_id).copied().unwrap_or(0);
        let depths = cov_map.get(target_id).cloned().unwrap_or_else(|| vec![0u32; ref_len]);

        let s = coverage::calculate_probe_stats(&depths, ref_len.max(depths.len()), proximity);
        stats.push((target_id.clone(), s));
        all_depths.insert(target_id.clone(), depths);
    }

    let _ = fs::remove_dir(&tmp_dir);

    stats.sort_by(|a, b| a.0.cmp(&b.0));

    let summary_path = prefixed_join(outdir, output_prefix, "individual_target_coverage_summary.tsv");
    write_probe_summary(&summary_path, &stats)?;
    log::info!("Individual target coverage summary: {}", summary_path.display());

    Ok(IndividualCoverageData {
        summary_path,
        coverage: all_depths,
    })
}

/// Compute the median probe length from the probe FASTA.
/// Returns 120 as a fallback if the file cannot be parsed or is empty.
pub fn compute_median_probe_length(probes_path: &Path) -> usize {
    match fasta::parse_fasta(probes_path) {
        Err(_) => 120,
        Ok(seqs) => {
            if seqs.is_empty() {
                return 120;
            }
            let mut lengths: Vec<usize> = seqs.values().map(|s| s.len()).collect();
            lengths.sort_unstable();
            let mid = lengths.len() / 2;
            if lengths.len() % 2 == 0 {
                (lengths[mid - 1] + lengths[mid]) / 2
            } else {
                lengths[mid]
            }
        }
    }
}

/// Gap detail record for TSV output.
struct GapRecord {
    target_id: String,
    gap_start: usize,  // 1-based
    gap_end: usize,    // 1-based inclusive
    gap_length: usize,
    is_terminal: bool,
    gc_content: f64,
    dust_score: f64,
    gap_type: Option<String>,          // Some if individual_coverage available
    individual_coverage: Option<f64>,  // Some if individual_coverage available
    gap_sequence: String,
}

/// Write gap detail records sorted by gap_length descending.
///
/// If `individual_coverage` is None, the gap_type and individual_coverage columns
/// are written as empty strings (no individual mapping was run).
pub fn compute_gap_details(
    combined_coverage: &HashMap<String, Vec<u32>>,
    combined_ref_lengths: &HashMap<String, usize>,
    individual_coverage: Option<&HashMap<String, Vec<u32>>>,
    targets_path: &Path,
    min_gap_length: usize,
    outdir: &Path,
    prefix: &str,
) -> Result<PathBuf> {
    let target_seqs = fasta::parse_fasta(targets_path)?;

    let mut records: Vec<GapRecord> = Vec::new();

    let mut target_ids: Vec<&String> = combined_coverage.keys().collect();
    target_ids.sort();

    for target_id in target_ids {
        let depths = &combined_coverage[target_id];
        let ref_len = combined_ref_lengths
            .get(target_id)
            .copied()
            .unwrap_or(depths.len());

        let gaps = coverage::collect_gaps(depths, min_gap_length);
        if gaps.is_empty() {
            continue;
        }

        let target_seq = match target_seqs.get(target_id) {
            Some(s) => s.as_str(),
            None => {
                log::warn!("Gap analysis: target '{}' not found in FASTA", target_id);
                continue;
            }
        };

        for gap in gaps {
            let seq_start = gap.start.min(target_seq.len());
            let seq_end = gap.end.min(target_seq.len());
            let gap_sequence = &target_seq[seq_start..seq_end];

            let gc_content = compute_gc_content(gap_sequence);
            let dust_score = sdust::masked_fraction(
                gap_sequence.as_bytes(),
                sdust::DEFAULT_THRESHOLD,
                sdust::DEFAULT_WINDOW,
            );

            let (gap_type, indiv_cov) = if let Some(indiv) = individual_coverage {
                let frac = if let Some(indiv_depths) = indiv.get(target_id) {
                    gap.coverage_fraction(indiv_depths)
                } else {
                    0.0
                };
                let gtype = if frac == 0.0 {
                    "true_gap".to_string()
                } else {
                    "multimapper_gap".to_string()
                };
                (Some(gtype), Some(frac))
            } else {
                (None, None)
            };

            records.push(GapRecord {
                target_id: target_id.clone(),
                gap_start: gap.start + 1,
                gap_end: gap.end,
                gap_length: gap.length(),
                is_terminal: gap.is_terminal(ref_len),
                gc_content,
                dust_score,
                gap_type,
                individual_coverage: indiv_cov,
                gap_sequence: gap_sequence.to_string(),
            });
        }
    }

    // Sort by gap_length descending so the report table shows worst gaps first.
    records.sort_by(|a, b| b.gap_length.cmp(&a.gap_length));

    let gap_path = prefixed_join(outdir, prefix, "gap_details.tsv");
    let file = File::create(&gap_path)
        .with_context(|| format!("Cannot create gap details file: {}", gap_path.display()))?;
    let mut w = BufWriter::new(file);

    writeln!(
        w,
        "target_id\tgap_start\tgap_end\tgap_length\tis_terminal\tgc_content\tdust_score\tgap_type\tindividual_coverage\tgap_sequence"
    )?;

    for r in &records {
        let terminal_str = if r.is_terminal { "TRUE" } else { "FALSE" };
        let gap_type_str = r.gap_type.as_deref().unwrap_or("");
        let indiv_cov_str = r
            .individual_coverage
            .map(|v| format!("{:.4}", v))
            .unwrap_or_default();
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{:.4}\t{:.4}\t{}\t{}\t{}",
            r.target_id,
            r.gap_start,
            r.gap_end,
            r.gap_length,
            terminal_str,
            r.gc_content,
            r.dust_score,
            gap_type_str,
            indiv_cov_str,
            r.gap_sequence,
        )?;
    }
    w.flush()?;

    log::info!(
        "Gap details: {} gaps written to {}",
        records.len(),
        gap_path.display()
    );

    Ok(gap_path)
}

fn compute_gc_content(seq: &str) -> f64 {
    let mut gc = 0usize;
    let mut total = 0usize;
    for b in seq.bytes() {
        match b {
            b'G' | b'g' | b'C' | b'c' => { gc += 1; total += 1; }
            b'A' | b'a' | b'T' | b't' => { total += 1; }
            _ => {}
        }
    }
    if total == 0 { 0.0 } else { gc as f64 / total as f64 }
}
