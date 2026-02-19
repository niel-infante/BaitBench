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
#[allow(dead_code)]
struct ReadLevelMetrics {
    /// Captured reads originating from target sequences
    target_captured: usize,
    /// Captured reads originating from distractor sequences
    distractor_captured: usize,
    /// Captured reads with unknown source
    unknown_captured: usize,
    /// Mapped reads where source == mapped reference (correct assignment)
    reads_correctly_mapped: usize,
    /// Mapped reads where source != mapped reference (misassignment)
    reads_incorrectly_mapped: usize,
    /// Mapped reads with unknown source
    reads_unknown_source: usize,
}

#[derive(Serialize)]
struct MetricsResult {
    tp_count: usize,
    fp_count: usize,
    fn_count: usize,
    tn_count: usize,
    sensitivity: f64,
    specificity: f64,
    precision: f64,
    f1_score: f64,
    fnr: f64,
    fpr: f64,
    true_positives: Vec<String>,
    false_positives: Vec<String>,
    false_negatives: Vec<String>,
    true_negatives: Vec<String>,
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
        &targets,
        &distractors,
        args.sam,
    )?;

    log::info!("  Target reads captured: {}", read_level.target_captured);
    log::info!("  Distractor reads captured: {}", read_level.distractor_captured);
    log::info!("  Reads correctly mapped: {}", read_level.reads_correctly_mapped);
    log::info!("  Reads incorrectly mapped: {}", read_level.reads_incorrectly_mapped);

    // Calculate genome-level metrics
    let metrics = calculate_metrics(&targets, &distractors, &detected);

    log::info!("  True Positives: {}", metrics.tp_count);
    log::info!("  False Positives: {}", metrics.fp_count);
    log::info!("  False Negatives: {}", metrics.fn_count);
    log::info!("  True Negatives: {}", metrics.tn_count);
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
    write_detail_tsv(args.output_detail, &targets, &distractors, &detected, &metrics)?;

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
    targets: &HashSet<String>,
    distractors: &HashSet<String>,
    sam_path: &Path,
) -> Result<ReadLevelMetrics> {
    // Count captured reads by source type
    let mut target_captured = 0usize;
    let mut distractor_captured = 0usize;
    let mut unknown_captured = 0usize;

    for name in captured_read_names {
        if let Some(source) = io_utils::extract_source_id(name) {
            if targets.contains(source) {
                target_captured += 1;
            } else if distractors.contains(source) {
                distractor_captured += 1;
            } else {
                unknown_captured += 1;
            }
        } else {
            unknown_captured += 1;
        }
    }

    // Read-level mapping accuracy: compare source to mapped reference
    let mappings = sam::get_read_mappings(sam_path)?;
    let mut reads_correctly_mapped = 0usize;
    let mut reads_incorrectly_mapped = 0usize;
    let mut reads_unknown_source = 0usize;

    for (read_name, mapped_ref) in &mappings {
        if let Some(source) = io_utils::extract_source_id(read_name) {
            if source == mapped_ref {
                reads_correctly_mapped += 1;
            } else {
                reads_incorrectly_mapped += 1;
            }
        } else {
            reads_unknown_source += 1;
        }
    }

    Ok(ReadLevelMetrics {
        target_captured,
        distractor_captured,
        unknown_captured,
        reads_correctly_mapped,
        reads_incorrectly_mapped,
        reads_unknown_source,
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

fn calculate_metrics(
    targets: &HashSet<String>,
    distractors: &HashSet<String>,
    detected: &HashMap<String, usize>,
) -> MetricsResult {
    let detected_set: HashSet<&String> = detected.keys().collect();
    let targets_ref: HashSet<&String> = targets.iter().collect();
    let distractors_ref: HashSet<&String> = distractors.iter().collect();

    let true_positives: Vec<String> = targets
        .iter()
        .filter(|id| detected.contains_key(*id))
        .cloned()
        .collect();
    let false_negatives: Vec<String> = targets
        .iter()
        .filter(|id| !detected.contains_key(*id))
        .cloned()
        .collect();
    let false_positives: Vec<String> = distractors
        .iter()
        .filter(|id| detected.contains_key(*id))
        .cloned()
        .collect();
    let true_negatives: Vec<String> = distractors
        .iter()
        .filter(|id| !detected.contains_key(*id))
        .cloned()
        .collect();

    let known: HashSet<&String> = targets_ref.union(&distractors_ref).copied().collect();
    let unknown_detected: Vec<String> = detected_set
        .iter()
        .filter(|id| !known.contains(**id))
        .map(|id| (*id).clone())
        .collect();

    let tp = true_positives.len();
    let fp = false_positives.len();
    let fn_ = false_negatives.len();
    let tn = true_negatives.len();

    let sensitivity = if tp + fn_ > 0 { tp as f64 / (tp + fn_) as f64 } else { 0.0 };
    let specificity = if tn + fp > 0 { tn as f64 / (tn + fp) as f64 } else { 0.0 };
    let precision = if tp + fp > 0 { tp as f64 / (tp + fp) as f64 } else { 0.0 };
    let fnr = if tp + fn_ > 0 { fn_ as f64 / (tp + fn_) as f64 } else { 0.0 };
    let fpr = if fp + tn > 0 { fp as f64 / (fp + tn) as f64 } else { 0.0 };
    let f1_score = if precision + sensitivity > 0.0 {
        2.0 * (precision * sensitivity) / (precision + sensitivity)
    } else {
        0.0
    };

    let mut tp_sorted = true_positives;
    let mut fp_sorted = false_positives;
    let mut fn_sorted = false_negatives;
    let mut tn_sorted = true_negatives;
    let mut unk_sorted = unknown_detected;
    tp_sorted.sort();
    fp_sorted.sort();
    fn_sorted.sort();
    tn_sorted.sort();
    unk_sorted.sort();

    MetricsResult {
        tp_count: tp,
        fp_count: fp,
        fn_count: fn_,
        tn_count: tn,
        sensitivity,
        specificity,
        precision,
        f1_score,
        fnr,
        fpr,
        true_positives: tp_sorted,
        false_positives: fp_sorted,
        false_negatives: fn_sorted,
        true_negatives: tn_sorted,
        unknown_detected: unk_sorted,
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
        "target_captured", "distractor_captured",
        "reads_correctly_mapped", "reads_incorrectly_mapped",
        "targets_total", "distractors_total",
        "tp_count", "fp_count", "fn_count", "tn_count",
        "sensitivity", "specificity", "precision", "f1_score",
    ];
    writeln!(w, "{}", headers.join("\t"))?;

    let targets_total = metrics.tp_count + metrics.fn_count;
    let distractors_total = metrics.tn_count + metrics.fp_count;

    let values = format!(
        "{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{:.4}\t{:.4}\t{:.4}",
        run_name, timestamp, num_reads, seed,
        reads_generated, reads_captured, capture_rate,
        read_level.target_captured, read_level.distractor_captured,
        read_level.reads_correctly_mapped, read_level.reads_incorrectly_mapped,
        targets_total, distractors_total,
        metrics.tp_count, metrics.fp_count, metrics.fn_count, metrics.tn_count,
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
        let (category, classification) = if targets.contains(ref_id) {
            ("target", "TP")
        } else if distractors.contains(ref_id) {
            ("distractor", "FP")
        } else {
            ("unknown", "UNKNOWN")
        };

        rows.push(DetailRow {
            reference_id: ref_id.clone(),
            category: category.to_string(),
            expected: if targets.contains(ref_id) { "true" } else { "false" }.to_string(),
            detected: "true".to_string(),
            read_count: count,
            classification: classification.to_string(),
        });
    }

    // False negatives (targets not detected)
    for ref_id in &metrics.false_negatives {
        rows.push(DetailRow {
            reference_id: ref_id.clone(),
            category: "target".to_string(),
            expected: "true".to_string(),
            detected: "false".to_string(),
            read_count: 0,
            classification: "FN".to_string(),
        });
    }

    // Sort: TP=0, FP=1, UNKNOWN=2, FN=3, then by read_count descending
    rows.sort_by(|a, b| {
        let order = |c: &str| match c {
            "TP" => 0,
            "FP" => 1,
            "UNKNOWN" => 2,
            "FN" => 3,
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
    target_captured: usize,
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
    fp_count: usize,
    fn_count: usize,
    tn_count: usize,
    sensitivity: f64,
    specificity: f64,
    precision: f64,
    f1_score: f64,
    fnr: f64,
    fpr: f64,
}

#[derive(Serialize)]
struct JsonDetails {
    true_positives: Vec<String>,
    false_positives: Vec<String>,
    false_negatives: Vec<String>,
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
            target_captured: read_level.target_captured,
            distractor_captured: read_level.distractor_captured,
        },
        read_level: ReadLevelStats {
            reads_correctly_mapped: read_level.reads_correctly_mapped,
            reads_incorrectly_mapped: read_level.reads_incorrectly_mapped,
        },
        metrics: JsonMetrics {
            tp_count: metrics.tp_count,
            fp_count: metrics.fp_count,
            fn_count: metrics.fn_count,
            tn_count: metrics.tn_count,
            sensitivity: metrics.sensitivity,
            specificity: metrics.specificity,
            precision: metrics.precision,
            f1_score: metrics.f1_score,
            fnr: metrics.fnr,
            fpr: metrics.fpr,
        },
        details: JsonDetails {
            true_positives: metrics.true_positives.clone(),
            false_positives: metrics.false_positives.clone(),
            false_negatives: metrics.false_negatives.clone(),
            unknown_detected: metrics.unknown_detected.clone(),
        },
    };

    let file = File::create(path)?;
    serde_json::to_writer_pretty(file, &output)?;
    Ok(())
}
