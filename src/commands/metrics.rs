use anyhow::{Context, Result};
use serde::Serialize;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use crate::alignment::sam;
use crate::fasta;
use crate::io_utils;

pub struct MetricsArgs<'a> {
    pub targets: &'a Path,
    pub distractors: &'a Path,
    pub sample: &'a Path,
    pub detected: &'a Path,
    pub reads: &'a Path,
    pub captured: &'a Path,
    pub sam: &'a Path,
    pub run_name: &'a str,
    pub num_reads: usize,
    pub seed: &'a str,
    pub output_summary: &'a Path,
    pub output_detail: &'a Path,
    pub output_json: Option<&'a Path>,
}

/// Read-level metrics derived from capture and mapping.
struct ReadLevelMetrics {
    /// Captured reads originating from sample target sequences
    sample_captured: usize,
    /// Captured reads originating from non-sample target sequences
    nonsample_target_captured: usize,
    /// Captured reads originating from distractor sequences
    distractor_captured: usize,
    /// Mapped reads where source == mapped reference (correct assignment)
    reads_correctly_mapped: usize,
    /// Mapped reads where source != mapped reference (misassignment)
    reads_incorrectly_mapped: usize,
}

/// 3-way classification metrics.
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
    let sample = io_utils::parse_id_set(args.sample)?;
    log::info!("  Found {} sample references", sample.len());

    log::info!("Parsing detection list...");
    let detected = parse_detected(args.detected)?;
    log::info!("  Found {} detected references", detected.len());

    // Count sequences in FASTA files
    log::info!("Counting sequences...");
    let reads_generated = fasta::count_sequences(args.reads)?;
    let reads_captured = fasta::count_sequences(args.captured)?;
    let capture_rate = if reads_generated > 0 {
        reads_captured as f64 / reads_generated as f64
    } else {
        0.0
    };
    log::info!("  Reads generated: {}", reads_generated);
    log::info!("  Reads captured: {}", reads_captured);
    log::info!("  Capture rate: {:.4}", capture_rate);

    // Read-level metrics: count captured reads by source type
    log::info!("Analyzing captured reads by source...");
    let captured_ids = fasta::parse_fasta_ids(args.captured)?;
    let read_level = compute_read_level_metrics(
        &captured_ids,
        &sample,
        &targets,
        &distractors,
        args.sam,
    )?;

    log::info!("  Sample reads captured: {}", read_level.sample_captured);
    log::info!("  Non-sample target reads captured: {}", read_level.nonsample_target_captured);
    log::info!("  Distractor reads captured: {}", read_level.distractor_captured);
    log::info!("  Reads correctly mapped: {}", read_level.reads_correctly_mapped);
    log::info!("  Reads incorrectly mapped: {}", read_level.reads_incorrectly_mapped);

    // Calculate genome-level metrics (3-way classification)
    let metrics = calculate_metrics(&sample, &targets, &distractors, &detected);

    log::info!("  True Positives (sample detected): {}", metrics.tp_count);
    log::info!("  False Negatives (sample missed): {}", metrics.fn_count);
    log::info!("  FP targets (non-sample target detected): {}", metrics.fp_target_count);
    log::info!("  FP distractors (distractor detected): {}", metrics.fp_distractor_count);
    log::info!("  TN targets (non-sample target not detected): {}", metrics.tn_target_count);
    log::info!("  TN distractors (distractor not detected): {}", metrics.tn_distractor_count);
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
        args.num_reads,
        args.seed,
        reads_generated,
        reads_captured,
        capture_rate,
        &metrics,
        &read_level,
    )?;

    // Write detail TSV
    log::info!("Writing detail to {}...", args.output_detail.display());
    write_detail_tsv(args.output_detail, &sample, &targets, &distractors, &detected, &metrics)?;

    // Write JSON
    if let Some(json_path) = args.output_json {
        log::info!("Writing JSON to {}...", json_path.display());
        write_json(
            json_path,
            args.run_name,
            &timestamp,
            args.num_reads,
            args.seed,
            reads_generated,
            reads_captured,
            capture_rate,
            &metrics,
            &read_level,
        )?;
    }

    log::info!("Metrics calculation complete.");
    Ok(())
}

fn compute_read_level_metrics(
    captured_read_names: &[String],
    sample: &HashSet<String>,
    targets: &HashSet<String>,
    distractors: &HashSet<String>,
    sam_path: &Path,
) -> Result<ReadLevelMetrics> {
    let mut sample_captured = 0usize;
    let mut nonsample_target_captured = 0usize;
    let mut distractor_captured = 0usize;

    for name in captured_read_names {
        if let Some(source) = io_utils::extract_source_id(name) {
            if sample.contains(source) {
                sample_captured += 1;
            } else if targets.contains(source) {
                nonsample_target_captured += 1;
            } else if distractors.contains(source) {
                distractor_captured += 1;
            }
        }
    }

    // Read-level mapping accuracy: compare source to mapped reference
    let mappings = sam::get_read_mappings(sam_path)?;
    let mut reads_correctly_mapped = 0usize;
    let mut reads_incorrectly_mapped = 0usize;

    for (read_name, mapped_ref) in &mappings {
        if let Some(source) = io_utils::extract_source_id(read_name) {
            if source == mapped_ref {
                reads_correctly_mapped += 1;
            } else {
                reads_incorrectly_mapped += 1;
            }
        }
    }

    Ok(ReadLevelMetrics {
        sample_captured,
        nonsample_target_captured,
        distractor_captured,
        reads_correctly_mapped,
        reads_incorrectly_mapped,
    })
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

/// 3-way genome-level classification:
/// - Sample targets: TP if detected, FN if not
/// - Non-sample targets: FP_target if detected, TN_target if not
/// - Distractors: FP_distractor if detected, TN_distractor if not
fn calculate_metrics(
    sample: &HashSet<String>,
    targets: &HashSet<String>,
    distractors: &HashSet<String>,
    detected: &HashMap<String, usize>,
) -> MetricsResult {
    // Sample targets
    let mut true_positives: Vec<String> = sample
        .iter()
        .filter(|id| detected.contains_key(*id))
        .cloned()
        .collect();
    let mut false_negatives: Vec<String> = sample
        .iter()
        .filter(|id| !detected.contains_key(*id))
        .cloned()
        .collect();

    // Non-sample targets (targets that are NOT in the sample)
    let nonsample_targets: HashSet<&String> = targets.iter().filter(|id| !sample.contains(*id)).collect();
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
    }
}

fn write_summary_tsv(
    path: &Path,
    run_name: &str,
    timestamp: &str,
    num_reads: usize,
    seed: &str,
    reads_generated: usize,
    reads_captured: usize,
    capture_rate: f64,
    metrics: &MetricsResult,
    read_level: &ReadLevelMetrics,
) -> Result<()> {
    let file = File::create(path)?;
    let mut w = BufWriter::new(file);

    let headers = [
        "run_name", "timestamp", "num_reads", "seed",
        "reads_generated", "reads_captured", "capture_rate",
        "sample_captured", "nonsample_target_captured", "distractor_captured",
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
        "{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{:.4}\t{:.4}\t{:.4}",
        run_name, timestamp, num_reads, seed,
        reads_generated, reads_captured, capture_rate,
        read_level.sample_captured, read_level.nonsample_target_captured, read_level.distractor_captured,
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
    read_count: usize,
    classification: String,
}

fn write_detail_tsv(
    path: &Path,
    sample: &HashSet<String>,
    targets: &HashSet<String>,
    distractors: &HashSet<String>,
    detected: &HashMap<String, usize>,
    metrics: &MetricsResult,
) -> Result<()> {
    let file = File::create(path)?;
    let mut w = BufWriter::new(file);

    writeln!(w, "reference_id\tcategory\texpected\tdetected\tread_count\tclassification")?;

    let mut rows: Vec<DetailRow> = Vec::new();

    // Detected references
    for (ref_id, &count) in detected {
        let (category, expected, classification) = if sample.contains(ref_id) {
            ("sample", "true", "TP")
        } else if targets.contains(ref_id) {
            ("target", "false", "FP_target")
        } else if distractors.contains(ref_id) {
            ("distractor", "false", "FP_distractor")
        } else {
            ("unknown", "false", "UNKNOWN")
        };

        rows.push(DetailRow {
            reference_id: ref_id.clone(),
            category: category.to_string(),
            expected: expected.to_string(),
            detected: "true".to_string(),
            read_count: count,
            classification: classification.to_string(),
        });
    }

    // False negatives (sample targets not detected)
    for ref_id in &metrics.false_negatives {
        rows.push(DetailRow {
            reference_id: ref_id.clone(),
            category: "sample".to_string(),
            expected: "true".to_string(),
            detected: "false".to_string(),
            read_count: 0,
            classification: "FN".to_string(),
        });
    }

    // Sort: TP=0, FP_target=1, FP_distractor=2, UNKNOWN=3, FN=4, then by read_count descending
    rows.sort_by(|a, b| {
        let order = |c: &str| match c {
            "TP" => 0,
            "FP_target" => 1,
            "FP_distractor" => 2,
            "UNKNOWN" => 3,
            "FN" => 4,
            _ => 99,
        };
        order(&a.classification)
            .cmp(&order(&b.classification))
            .then(b.read_count.cmp(&a.read_count))
    });

    for row in &rows {
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}",
            row.reference_id, row.category, row.expected, row.detected, row.read_count, row.classification
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
    num_reads: usize,
    seed: String,
}

#[derive(Serialize)]
struct CaptureStats {
    reads_generated: usize,
    reads_captured: usize,
    capture_rate: f64,
    sample_captured: usize,
    nonsample_target_captured: usize,
    distractor_captured: usize,
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
}

fn write_json(
    path: &Path,
    run_name: &str,
    timestamp: &str,
    num_reads: usize,
    seed: &str,
    reads_generated: usize,
    reads_captured: usize,
    capture_rate: f64,
    metrics: &MetricsResult,
    read_level: &ReadLevelMetrics,
) -> Result<()> {
    let output = JsonOutput {
        run_info: RunInfo {
            run_name: run_name.to_string(),
            timestamp: timestamp.to_string(),
            num_reads,
            seed: seed.to_string(),
        },
        capture_stats: CaptureStats {
            reads_generated,
            reads_captured,
            capture_rate,
            sample_captured: read_level.sample_captured,
            nonsample_target_captured: read_level.nonsample_target_captured,
            distractor_captured: read_level.distractor_captured,
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
        },
    };

    let file = File::create(path)?;
    serde_json::to_writer_pretty(file, &output)?;
    Ok(())
}
