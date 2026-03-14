use anyhow::{Context, Result};
use serde::Serialize;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use crate::alignment::{coverage, sam};
use crate::fasta;
use crate::io_utils;

pub struct MetricsArgs<'a> {
    pub targets: &'a Path,
    pub distractors: &'a Path,
    pub sample: &'a Path,
    /// Optional sample-target-map file (genome_id → target_id mappings).
    /// When present, enables genome-aware classification where sample IDs are
    /// genome IDs and detection is measured through their mapped target regions.
    pub sample_target_map: Option<&'a Path>,
    pub detected: &'a Path,
    pub fragments: &'a Path,
    pub captured: &'a Path,
    pub sam: &'a Path,
    pub run_name: &'a str,
    pub num_fragments: usize,
    pub seed: &'a str,
    pub output_summary: &'a Path,
    pub output_detail: &'a Path,
    pub output_json: Option<&'a Path>,
    pub output_coverage: Option<&'a Path>,
}

/// Read-level metrics derived from capture and mapping.
struct ReadLevelMetrics {
    /// Captured fragments originating from sample target sequences
    sample_captured: usize,
    /// Captured fragments originating from non-sample target sequences
    nonsample_target_captured: usize,
    /// Captured fragments originating from distractor sequences
    distractor_captured: usize,
    /// Captured fragments from untargeted sample genomes (genomes mode only)
    untargeted_captured: usize,
    /// Mapped reads where source == mapped reference (correct assignment)
    reads_correctly_mapped: usize,
    /// Mapped reads where source != mapped reference (misassignment)
    reads_incorrectly_mapped: usize,
}

/// Classification metrics (3-way or 4-way with untargeted).
struct MetricsResult {
    // Counts
    tp_count: usize,
    fn_count: usize,
    fp_target_count: usize,
    fp_distractor_count: usize,
    fp_total: usize,
    tn_target_count: usize,
    tn_distractor_count: usize,
    tn_total: usize,
    // Rates
    sensitivity: f64,
    specificity: f64,
    precision: f64,
    f1_score: f64,
    // Detail lists
    true_positives: Vec<String>,
    false_negatives: Vec<String>,
    fp_targets: Vec<String>,
    fp_distractors: Vec<String>,
    tn_targets: Vec<String>,
    tn_distractors: Vec<String>,
    unknown_detected: Vec<String>,
    // Untargeted genomes (in sample, no target mapping)
    untargeted_genomes: Vec<String>,
}

/// Resolved context for genome-aware metrics.
/// When a sample-target-map is provided, this holds the derived sets.
struct GenomeContext {
    /// Sample target IDs (derived from sample genomes via mapping)
    sample_targets: HashSet<String>,
    /// All genome IDs (for source classification)
    genome_ids: HashSet<String>,
    /// Sample genome IDs
    sample_genome_ids: HashSet<String>,
    /// Genome-to-target mapping
    genome_to_targets: HashMap<String, Vec<String>>,
    /// Target-to-genome reverse mapping
    target_to_genomes: HashMap<String, Vec<String>>,
    /// Untargeted sample genomes (in sample, no mapping)
    untargeted_genomes: Vec<String>,
}

pub fn execute(args: &MetricsArgs) -> Result<()> {
    // Parse input files
    log::info!("Parsing targets file...");
    let targets = io_utils::parse_id_set(args.targets)?;
    log::info!("  Found {} target references", targets.len());

    log::info!("Parsing distractors file...");
    let distractors = io_utils::parse_id_set(args.distractors)?;
    log::info!("  Found {} distractor references", distractors.len());

    log::info!("Parsing sample file...");
    let sample_ids = io_utils::parse_id_set(args.sample)?;
    log::info!("  Found {} sample references", sample_ids.len());

    // Parse sample-target-map if present (genome-aware mode)
    let genome_ctx = if let Some(map_path) = args.sample_target_map {
        log::info!("Parsing sample-target-map...");
        let genome_to_targets = io_utils::parse_sample_target_map(map_path)?;

        // Build reverse mapping: target → genomes
        let mut target_to_genomes: HashMap<String, Vec<String>> = HashMap::new();
        for (genome_id, target_list) in &genome_to_targets {
            for target_id in target_list {
                target_to_genomes
                    .entry(target_id.clone())
                    .or_default()
                    .push(genome_id.clone());
            }
        }

        // Derive sample targets: targets linked to any sample genome
        let mut sample_targets = HashSet::new();
        for sample_id in &sample_ids {
            if let Some(target_list) = genome_to_targets.get(sample_id) {
                for target_id in target_list {
                    sample_targets.insert(target_id.clone());
                }
            }
        }

        // All genome IDs (from the map keys + sample IDs that may be untargeted)
        let mut genome_ids: HashSet<String> = genome_to_targets.keys().cloned().collect();
        for id in &sample_ids {
            genome_ids.insert(id.clone());
        }

        // Untargeted sample genomes
        let mut untargeted_genomes: Vec<String> = sample_ids
            .iter()
            .filter(|id| !genome_to_targets.contains_key(*id))
            .cloned()
            .collect();
        untargeted_genomes.sort();

        log::info!(
            "  Genome-aware mode: {} sample targets derived, {} untargeted genomes",
            sample_targets.len(),
            untargeted_genomes.len()
        );

        Some(GenomeContext {
            sample_targets,
            genome_ids,
            sample_genome_ids: sample_ids.clone(),
            genome_to_targets,
            target_to_genomes,
            untargeted_genomes,
        })
    } else {
        None
    };

    log::info!("Parsing detection list...");
    let detected = parse_detected(args.detected)?;
    log::info!("  Found {} detected references", detected.len());

    // Count sequences in FASTA files
    log::info!("Counting sequences...");
    let fragments_generated = fasta::count_sequences(args.fragments)?;
    let fragments_captured = fasta::count_sequences(args.captured)?;
    let capture_rate = if fragments_generated > 0 {
        fragments_captured as f64 / fragments_generated as f64
    } else {
        0.0
    };
    log::info!("  Fragments generated: {}", fragments_generated);
    log::info!("  Fragments captured: {}", fragments_captured);
    log::info!("  Capture rate: {:.4}", capture_rate);

    // Per-reference fragment counts (for detail table)
    let generated_per_ref = count_per_source(args.fragments)?;
    let captured_per_ref = count_per_source(args.captured)?;

    // Determine the effective sample set for classification
    // In genome mode: sample targets derived from mapping
    // In standard mode: sample IDs directly
    let effective_sample = if let Some(ctx) = &genome_ctx {
        &ctx.sample_targets
    } else {
        &sample_ids
    };

    // Read-level metrics
    log::info!("Analyzing captured fragments by source...");
    let captured_ids = fasta::parse_fasta_ids(args.captured)?;
    let read_level = compute_read_level_metrics(
        &captured_ids,
        effective_sample,
        &targets,
        &distractors,
        genome_ctx.as_ref(),
        args.sam,
    )?;

    log::info!("  Sample fragments captured: {}", read_level.sample_captured);
    log::info!("  Non-sample target fragments captured: {}", read_level.nonsample_target_captured);
    log::info!("  Distractor fragments captured: {}", read_level.distractor_captured);
    if read_level.untargeted_captured > 0 {
        log::info!("  Untargeted fragments captured: {}", read_level.untargeted_captured);
    }
    log::info!("  Reads correctly mapped: {}", read_level.reads_correctly_mapped);
    log::info!("  Reads incorrectly mapped: {}", read_level.reads_incorrectly_mapped);

    // Compute per-reference coverage from SAM
    log::info!("Computing per-reference coverage...");
    let coverage_result = coverage::compute_coverage(args.sam)?;
    let mut coverage_stats: HashMap<String, coverage::CoverageStats> = HashMap::new();
    for (ref_id, depths) in &coverage_result.coverage {
        let ref_len = coverage_result.ref_lengths.get(ref_id).copied().unwrap_or(depths.len());
        coverage_stats.insert(ref_id.clone(), coverage::calculate_stats(depths, ref_len));
    }

    // Calculate genome-level metrics (classification)
    let metrics = calculate_metrics(effective_sample, &targets, &distractors, &detected, genome_ctx.as_ref());

    log::info!("  True Positives (sample detected): {}", metrics.tp_count);
    log::info!("  False Negatives (sample missed): {}", metrics.fn_count);
    log::info!("  FP targets (non-sample target detected): {}", metrics.fp_target_count);
    log::info!("  FP distractors (distractor detected): {}", metrics.fp_distractor_count);
    log::info!("  TN targets (non-sample target not detected): {}", metrics.tn_target_count);
    log::info!("  TN distractors (distractor not detected): {}", metrics.tn_distractor_count);
    if !metrics.untargeted_genomes.is_empty() {
        log::info!("  Untargeted genomes: {}", metrics.untargeted_genomes.len());
    }
    log::info!("  Sensitivity: {:.4}", metrics.sensitivity);
    log::info!("  Specificity: {:.4}", metrics.specificity);
    log::info!("  Precision: {:.4}", metrics.precision);
    log::info!("  F1 Score: {:.4}", metrics.f1_score);

    let timestamp = chrono::Local::now().format("%Y-%m-%dT%H:%M:%S").to_string();

    // Write summary TSV
    log::info!("Writing summary to {}...", args.output_summary.display());
    write_summary_tsv(
        args.output_summary,
        args.run_name,
        &timestamp,
        args.num_fragments,
        args.seed,
        fragments_generated,
        fragments_captured,
        capture_rate,
        &metrics,
        &read_level,
    )?;

    // Build detail rows (shared between TSV and JSON)
    let detail_rows = build_detail_rows(
        effective_sample, &targets, &distractors, &detected, &metrics,
        &generated_per_ref, &captured_per_ref, &coverage_stats,
        genome_ctx.as_ref(),
    );

    // Write detail TSV
    log::info!("Writing detail to {}...", args.output_detail.display());
    write_detail_tsv(args.output_detail, &detail_rows)?;

    // Write JSON
    if let Some(json_path) = args.output_json {
        log::info!("Writing JSON to {}...", json_path.display());
        write_json(
            json_path,
            args.run_name,
            &timestamp,
            args.num_fragments,
            args.seed,
            fragments_generated,
            fragments_captured,
            capture_rate,
            &metrics,
            &read_level,
            detail_rows,
        )?;
    }

    // Write run-length encoded coverage intervals
    if let Some(cov_path) = args.output_coverage {
        log::info!("Writing coverage profile to {}...", cov_path.display());
        coverage::write_coverage_intervals(cov_path, &coverage_result.coverage)?;
    }

    log::info!("Metrics calculation complete.");
    Ok(())
}

fn compute_read_level_metrics(
    captured_read_names: &[String],
    effective_sample: &HashSet<String>,
    targets: &HashSet<String>,
    distractors: &HashSet<String>,
    genome_ctx: Option<&GenomeContext>,
    sam_path: &Path,
) -> Result<ReadLevelMetrics> {
    let mut sample_captured = 0usize;
    let mut nonsample_target_captured = 0usize;
    let mut distractor_captured = 0usize;
    let mut untargeted_captured = 0usize;

    for name in captured_read_names {
        if let Some(source) = io_utils::extract_source_id(name) {
            if let Some(ctx) = genome_ctx {
                // Genome-aware mode: source is a genome ID
                if ctx.sample_genome_ids.contains(source) {
                    if ctx.genome_to_targets.contains_key(source) {
                        sample_captured += 1;
                    } else {
                        untargeted_captured += 1;
                    }
                } else if ctx.genome_ids.contains(source) {
                    nonsample_target_captured += 1;
                } else if distractors.contains(source) {
                    distractor_captured += 1;
                }
            } else {
                // Standard mode: source is a target ID
                if effective_sample.contains(source) {
                    sample_captured += 1;
                } else if targets.contains(source) {
                    nonsample_target_captured += 1;
                } else if distractors.contains(source) {
                    distractor_captured += 1;
                }
            }
        }
    }

    // Read-level mapping accuracy: compare source to mapped reference
    let mappings = sam::get_read_mappings(sam_path)?;
    let mut reads_correctly_mapped = 0usize;
    let mut reads_incorrectly_mapped = 0usize;

    for (read_name, mapped_ref) in &mappings {
        if let Some(source) = io_utils::extract_source_id(read_name) {
            if let Some(ctx) = genome_ctx {
                // Genome-aware mode: correct if the mapped target is one of the
                // genome's targets
                if let Some(valid_targets) = ctx.genome_to_targets.get(source) {
                    if valid_targets.iter().any(|t| t == mapped_ref) {
                        reads_correctly_mapped += 1;
                    } else {
                        reads_incorrectly_mapped += 1;
                    }
                } else {
                    // Untargeted genome or distractor — any mapping to a target is incorrect
                    reads_incorrectly_mapped += 1;
                }
            } else {
                // Standard mode: source ID should match mapped reference
                if source == mapped_ref {
                    reads_correctly_mapped += 1;
                } else {
                    reads_incorrectly_mapped += 1;
                }
            }
        }
    }

    Ok(ReadLevelMetrics {
        sample_captured,
        nonsample_target_captured,
        distractor_captured,
        untargeted_captured,
        reads_correctly_mapped,
        reads_incorrectly_mapped,
    })
}

/// Count sequences per source genome in a FASTA file.
fn count_per_source(path: &Path) -> Result<HashMap<String, usize>> {
    let ids = fasta::parse_fasta_ids(path)?;
    let mut counts: HashMap<String, usize> = HashMap::new();
    for name in &ids {
        if let Some(source) = io_utils::extract_source_id(name) {
            *counts.entry(source.to_string()).or_insert(0) += 1;
        }
    }
    Ok(counts)
}

fn parse_detected(path: &Path) -> Result<HashMap<String, usize>> {
    let file = File::open(path)
        .with_context(|| format!("Cannot open detection list: {}", path.display()))?;
    let reader = BufReader::new(file);
    let mut detected = HashMap::new();

    for line in reader.lines() {
        let line = line?;
        let line = line.trim().to_string();
        if line.is_empty() {
            continue;
        }
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() >= 2 {
            let ref_id = parts[0].to_string();
            let count: usize = parts[1].parse().unwrap_or(0);
            detected.insert(ref_id, count);
        }
    }

    Ok(detected)
}

/// Classification:
/// - Sample targets: TP if detected, FN if not
/// - Non-sample targets: FP_target if detected, TN_target if not
/// - Distractors: FP_distractor if detected, TN_distractor if not
/// - Untargeted genomes: separate category (not in TP/FP/FN/TN)
///
/// In genome-aware mode, `effective_sample` contains target IDs derived from
/// sample genomes via the mapping. Detection is at the target level.
fn calculate_metrics(
    effective_sample: &HashSet<String>,
    targets: &HashSet<String>,
    distractors: &HashSet<String>,
    detected: &HashMap<String, usize>,
    genome_ctx: Option<&GenomeContext>,
) -> MetricsResult {
    // Sample targets
    let mut true_positives: Vec<String> = effective_sample
        .iter()
        .filter(|id| detected.contains_key(*id))
        .cloned()
        .collect();
    let mut false_negatives: Vec<String> = effective_sample
        .iter()
        .filter(|id| !detected.contains_key(*id))
        .cloned()
        .collect();

    // Non-sample targets (targets NOT in the effective sample)
    let nonsample_targets: HashSet<&String> = targets.iter().filter(|id| !effective_sample.contains(*id)).collect();
    let mut fp_targets: Vec<String> = nonsample_targets
        .iter()
        .filter(|id| detected.contains_key(**id))
        .map(|id| (*id).clone())
        .collect();
    let mut tn_targets: Vec<String> = nonsample_targets
        .iter()
        .filter(|id| !detected.contains_key(**id))
        .map(|id| (*id).clone())
        .collect();

    // Distractors
    let mut fp_distractors: Vec<String> = distractors
        .iter()
        .filter(|id| detected.contains_key(*id))
        .cloned()
        .collect();
    let mut tn_distractors: Vec<String> = distractors
        .iter()
        .filter(|id| !detected.contains_key(*id))
        .cloned()
        .collect();

    // Unknown detected (not in any category)
    let all_known: HashSet<&String> = targets.iter().chain(distractors.iter()).collect();
    let mut unknown_detected: Vec<String> = detected
        .keys()
        .filter(|id| !all_known.contains(id))
        .cloned()
        .collect();

    // Untargeted genomes
    let mut untargeted_genomes: Vec<String> = genome_ctx
        .map(|ctx| ctx.untargeted_genomes.clone())
        .unwrap_or_default();

    let tp = true_positives.len();
    let fn_ = false_negatives.len();
    let fp_target = fp_targets.len();
    let fp_distractor = fp_distractors.len();
    let fp_total = fp_target + fp_distractor;
    let tn_target = tn_targets.len();
    let tn_distractor = tn_distractors.len();
    let tn_total = tn_target + tn_distractor;

    let sensitivity = if tp + fn_ > 0 { tp as f64 / (tp + fn_) as f64 } else { 0.0 };
    let specificity = if tn_total + fp_total > 0 { tn_total as f64 / (tn_total + fp_total) as f64 } else { 0.0 };
    let precision = if tp + fp_total > 0 { tp as f64 / (tp + fp_total) as f64 } else { 0.0 };
    let f1_score = if precision + sensitivity > 0.0 {
        2.0 * (precision * sensitivity) / (precision + sensitivity)
    } else {
        0.0
    };

    true_positives.sort();
    false_negatives.sort();
    fp_targets.sort();
    fp_distractors.sort();
    tn_targets.sort();
    tn_distractors.sort();
    unknown_detected.sort();
    untargeted_genomes.sort();

    MetricsResult {
        tp_count: tp,
        fn_count: fn_,
        fp_target_count: fp_target,
        fp_distractor_count: fp_distractor,
        fp_total,
        tn_target_count: tn_target,
        tn_distractor_count: tn_distractor,
        tn_total,
        sensitivity,
        specificity,
        precision,
        f1_score,
        true_positives,
        false_negatives,
        fp_targets,
        fp_distractors,
        tn_targets,
        tn_distractors,
        unknown_detected,
        untargeted_genomes,
    }
}

fn write_summary_tsv(
    path: &Path,
    run_name: &str,
    timestamp: &str,
    num_fragments: usize,
    seed: &str,
    fragments_generated: usize,
    fragments_captured: usize,
    capture_rate: f64,
    metrics: &MetricsResult,
    read_level: &ReadLevelMetrics,
) -> Result<()> {
    let file = File::create(path)?;
    let mut w = BufWriter::new(file);

    let headers = [
        "run_name", "timestamp", "num_fragments", "seed",
        "fragments_generated", "fragments_captured", "capture_rate",
        "sample_captured", "nonsample_target_captured", "distractor_captured",
        "untargeted_captured",
        "reads_correctly_mapped", "reads_incorrectly_mapped",
        "sample_total", "nonsample_target_total", "distractors_total",
        "tp_count", "fn_count",
        "fp_target_count", "fp_distractor_count", "fp_total",
        "tn_target_count", "tn_distractor_count", "tn_total",
        "sensitivity", "specificity", "precision", "f1_score",
    ];
    writeln!(w, "{}", headers.join("\t"))?;

    let sample_total = metrics.tp_count + metrics.fn_count;
    let nonsample_target_total = metrics.fp_target_count + metrics.tn_target_count;
    let distractors_total = metrics.fp_distractor_count + metrics.tn_distractor_count;

    let values = format!(
        "{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{:.4}\t{:.4}\t{:.4}",
        run_name, timestamp, num_fragments, seed,
        fragments_generated, fragments_captured, capture_rate,
        read_level.sample_captured, read_level.nonsample_target_captured, read_level.distractor_captured,
        read_level.untargeted_captured,
        read_level.reads_correctly_mapped, read_level.reads_incorrectly_mapped,
        sample_total, nonsample_target_total, distractors_total,
        metrics.tp_count, metrics.fn_count,
        metrics.fp_target_count, metrics.fp_distractor_count, metrics.fp_total,
        metrics.tn_target_count, metrics.tn_distractor_count, metrics.tn_total,
        metrics.sensitivity, metrics.specificity, metrics.precision, metrics.f1_score,
    );
    writeln!(w, "{}", values)?;

    w.flush()?;
    Ok(())
}

#[derive(Serialize)]
struct DetailRow {
    reference_id: String,
    category: String,
    expected: String,
    detected: String,
    fragments_generated: usize,
    fragments_captured: usize,
    reads_assigned: usize,
    classification: String,
    ref_length: usize,
    avg_coverage: f64,
    pct_covered_5x: f64,
    pct_covered_20x: f64,
}

fn write_detail_tsv(path: &Path, rows: &[DetailRow]) -> Result<()> {
    let file = File::create(path)?;
    let mut w = BufWriter::new(file);

    writeln!(w, "reference_id\tcategory\texpected\tdetected\tfragments_generated\tfragments_captured\treads_assigned\tclassification\tref_length\tavg_coverage\tpct_covered_5x\tpct_covered_20x")?;

    for row in rows {
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.2}\t{:.1}\t{:.1}",
            row.reference_id, row.category, row.expected, row.detected,
            row.fragments_generated, row.fragments_captured, row.reads_assigned,
            row.classification, row.ref_length, row.avg_coverage,
            row.pct_covered_5x, row.pct_covered_20x
        )?;
    }

    w.flush()?;
    Ok(())
}

#[derive(Serialize)]
struct JsonOutput {
    run_info: RunInfo,
    capture_stats: CaptureStats,
    read_level: ReadLevelStats,
    metrics: JsonMetrics,
    details: JsonDetails,
}

#[derive(Serialize)]
struct RunInfo {
    run_name: String,
    timestamp: String,
    num_fragments: usize,
    seed: String,
}

#[derive(Serialize)]
struct CaptureStats {
    fragments_generated: usize,
    fragments_captured: usize,
    capture_rate: f64,
    sample_captured: usize,
    nonsample_target_captured: usize,
    distractor_captured: usize,
    untargeted_captured: usize,
}

#[derive(Serialize)]
struct ReadLevelStats {
    reads_correctly_mapped: usize,
    reads_incorrectly_mapped: usize,
}

#[derive(Serialize)]
struct JsonMetrics {
    tp_count: usize,
    fn_count: usize,
    fp_target_count: usize,
    fp_distractor_count: usize,
    fp_total: usize,
    tn_target_count: usize,
    tn_distractor_count: usize,
    tn_total: usize,
    sensitivity: f64,
    specificity: f64,
    precision: f64,
    f1_score: f64,
}

#[derive(Serialize)]
struct JsonDetails {
    true_positives: Vec<String>,
    false_negatives: Vec<String>,
    fp_targets: Vec<String>,
    fp_distractors: Vec<String>,
    tn_targets: Vec<String>,
    tn_distractors: Vec<String>,
    unknown_detected: Vec<String>,
    untargeted_genomes: Vec<String>,
    detail_rows: Vec<DetailRow>,
}

fn build_detail_rows(
    effective_sample: &HashSet<String>,
    targets: &HashSet<String>,
    distractors: &HashSet<String>,
    detected: &HashMap<String, usize>,
    metrics: &MetricsResult,
    generated_per_ref: &HashMap<String, usize>,
    captured_per_ref: &HashMap<String, usize>,
    coverage_stats: &HashMap<String, coverage::CoverageStats>,
    genome_ctx: Option<&GenomeContext>,
) -> Vec<DetailRow> {
    let mut rows: Vec<DetailRow> = Vec::new();

    for (ref_id, &count) in detected {
        let (category, expected, classification) = if effective_sample.contains(ref_id) {
            ("sample", "true", "TP")
        } else if targets.contains(ref_id) {
            ("target", "false", "FP_target")
        } else if distractors.contains(ref_id) {
            ("distractor", "false", "FP_distractor")
        } else {
            ("unknown", "false", "UNKNOWN")
        };

        // In genome mode, fragment counts are per genome source, but detected
        // is per target. Use target-to-genome mapping to aggregate fragment counts.
        let (frag_gen, frag_cap) = if let Some(ctx) = genome_ctx {
            if let Some(genome_list) = ctx.target_to_genomes.get(ref_id) {
                let gen: usize = genome_list
                    .iter()
                    .map(|g| generated_per_ref.get(g).copied().unwrap_or(0))
                    .sum();
                let cap: usize = genome_list
                    .iter()
                    .map(|g| captured_per_ref.get(g).copied().unwrap_or(0))
                    .sum();
                (gen, cap)
            } else {
                (
                    generated_per_ref.get(ref_id).copied().unwrap_or(0),
                    captured_per_ref.get(ref_id).copied().unwrap_or(0),
                )
            }
        } else {
            (
                generated_per_ref.get(ref_id).copied().unwrap_or(0),
                captured_per_ref.get(ref_id).copied().unwrap_or(0),
            )
        };

        let cov = coverage_stats.get(ref_id);
        rows.push(DetailRow {
            reference_id: ref_id.clone(),
            category: category.to_string(),
            expected: expected.to_string(),
            detected: "true".to_string(),
            fragments_generated: frag_gen,
            fragments_captured: frag_cap,
            reads_assigned: count,
            classification: classification.to_string(),
            ref_length: cov.map(|c| c.ref_length).unwrap_or(0),
            avg_coverage: cov.map(|c| c.avg_coverage).unwrap_or(0.0),
            pct_covered_5x: cov.map(|c| c.pct_covered_5x).unwrap_or(0.0),
            pct_covered_20x: cov.map(|c| c.pct_covered_20x).unwrap_or(0.0),
        });
    }

    for ref_id in &metrics.false_negatives {
        let cov = coverage_stats.get(ref_id);
        rows.push(DetailRow {
            reference_id: ref_id.clone(),
            category: "sample".to_string(),
            expected: "true".to_string(),
            detected: "false".to_string(),
            fragments_generated: generated_per_ref.get(ref_id).copied().unwrap_or(0),
            fragments_captured: captured_per_ref.get(ref_id).copied().unwrap_or(0),
            reads_assigned: 0,
            classification: "FN".to_string(),
            ref_length: cov.map(|c| c.ref_length).unwrap_or(0),
            avg_coverage: cov.map(|c| c.avg_coverage).unwrap_or(0.0),
            pct_covered_5x: cov.map(|c| c.pct_covered_5x).unwrap_or(0.0),
            pct_covered_20x: cov.map(|c| c.pct_covered_20x).unwrap_or(0.0),
        });
    }

    // Add untargeted genomes to detail rows
    if let Some(ctx) = genome_ctx {
        for genome_id in &ctx.untargeted_genomes {
            let frag_gen = generated_per_ref.get(genome_id).copied().unwrap_or(0);
            let frag_cap = captured_per_ref.get(genome_id).copied().unwrap_or(0);
            rows.push(DetailRow {
                reference_id: genome_id.clone(),
                category: "untargeted".to_string(),
                expected: "false".to_string(),
                detected: if frag_cap > 0 { "true" } else { "false" }.to_string(),
                fragments_generated: frag_gen,
                fragments_captured: frag_cap,
                reads_assigned: 0, // untargeted genomes have no target to detect reads against
                classification: "untargeted".to_string(),
                ref_length: 0,
                avg_coverage: 0.0,
                pct_covered_5x: 0.0,
                pct_covered_20x: 0.0,
            });
        }
    }

    rows.sort_by(|a, b| {
        let order = |c: &str| match c {
            "TP" => 0,
            "FP_target" => 1,
            "FP_distractor" => 2,
            "UNKNOWN" => 3,
            "FN" => 4,
            "untargeted" => 5,
            _ => 99,
        };
        order(&a.classification)
            .cmp(&order(&b.classification))
            .then(b.reads_assigned.cmp(&a.reads_assigned))
    });

    rows
}

fn write_json(
    path: &Path,
    run_name: &str,
    timestamp: &str,
    num_fragments: usize,
    seed: &str,
    fragments_generated: usize,
    fragments_captured: usize,
    capture_rate: f64,
    metrics: &MetricsResult,
    read_level: &ReadLevelMetrics,
    detail_rows: Vec<DetailRow>,
) -> Result<()> {
    let output = JsonOutput {
        run_info: RunInfo {
            run_name: run_name.to_string(),
            timestamp: timestamp.to_string(),
            num_fragments,
            seed: seed.to_string(),
        },
        capture_stats: CaptureStats {
            fragments_generated,
            fragments_captured,
            capture_rate,
            sample_captured: read_level.sample_captured,
            nonsample_target_captured: read_level.nonsample_target_captured,
            distractor_captured: read_level.distractor_captured,
            untargeted_captured: read_level.untargeted_captured,
        },
        read_level: ReadLevelStats {
            reads_correctly_mapped: read_level.reads_correctly_mapped,
            reads_incorrectly_mapped: read_level.reads_incorrectly_mapped,
        },
        metrics: JsonMetrics {
            tp_count: metrics.tp_count,
            fn_count: metrics.fn_count,
            fp_target_count: metrics.fp_target_count,
            fp_distractor_count: metrics.fp_distractor_count,
            fp_total: metrics.fp_total,
            tn_target_count: metrics.tn_target_count,
            tn_distractor_count: metrics.tn_distractor_count,
            tn_total: metrics.tn_total,
            sensitivity: metrics.sensitivity,
            specificity: metrics.specificity,
            precision: metrics.precision,
            f1_score: metrics.f1_score,
        },
        details: JsonDetails {
            true_positives: metrics.true_positives.clone(),
            false_negatives: metrics.false_negatives.clone(),
            fp_targets: metrics.fp_targets.clone(),
            fp_distractors: metrics.fp_distractors.clone(),
            tn_targets: metrics.tn_targets.clone(),
            tn_distractors: metrics.tn_distractors.clone(),
            unknown_detected: metrics.unknown_detected.clone(),
            untargeted_genomes: metrics.untargeted_genomes.clone(),
            detail_rows,
        },
    };

    let file = File::create(path)?;
    serde_json::to_writer_pretty(file, &output)?;
    Ok(())
}
