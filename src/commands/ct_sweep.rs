use anyhow::{bail, Context, Result};
use std::collections::HashSet;
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use crate::alignment::coverage;
use crate::commands::{capture, enrich, filter, map_reads, prepare, sequence, simulate};
use crate::external::{minimap2, rscript};
use crate::io_utils;

pub struct CtSweepArgs<'a> {
    pub targets: &'a Path,
    pub distractors: &'a [PathBuf],
    pub probes: &'a Path,
    pub sample: &'a Path,
    pub ct_values: &'a [f64],
    pub ct_baseline: f64,
    pub ct_baseline_fraction: f64,
    pub num_fragments: usize,
    pub read_length: usize,
    pub num_sequences: Option<usize>,
    pub seed: Option<u64>,
    pub fragment_length_mean: f64,
    pub fragment_length_min: usize,
    pub fragment_length_max: usize,
    pub capture_method: capture::CaptureMethod,
    pub max_mismatches: u32,
    pub min_match_bases: u32,
    pub blast_db: Option<&'a str>,
    pub fold_enrichment: Option<f64>,
    pub host_fasta: Option<&'a Path>,
    pub minimap_preset: &'a str,
    pub host_minimap_preset: &'a str,
    pub threads: usize,
    pub outdir: PathBuf,
    pub no_report: bool,
}

struct DepthCurveRow {
    ct: f64,
    reference_id: String,
    depth_threshold: usize,
    pct_covered: f64,
}

pub fn execute(args: &CtSweepArgs) -> Result<()> {
    // Validate inputs
    if !args.targets.exists() {
        bail!("Targets file not found: {}", args.targets.display());
    }
    if !args.probes.exists() {
        bail!("Probes file not found: {}", args.probes.display());
    }
    if !args.sample.exists() {
        bail!("Sample manifest not found: {}", args.sample.display());
    }
    if args.ct_values.is_empty() {
        bail!("At least one CT value is required");
    }

    minimap2::check_available()?;
    fs::create_dir_all(&args.outdir)?;

    // Parse sample manifest to get sample target IDs
    let sample_manifest = io_utils::parse_sample_manifest(args.sample)?;
    let sample_ids: HashSet<String> = sample_manifest.keys().cloned().collect();
    if sample_ids.is_empty() {
        bail!("Sample manifest is empty");
    }

    let mut sorted_ct_values = args.ct_values.to_vec();
    sorted_ct_values.sort_by(|a, b| a.partial_cmp(b).unwrap());

    log::info!("=============================================");
    log::info!("BaitBench - CT Sweep Analysis");
    log::info!("=============================================");
    log::info!("Targets        : {}", args.targets.display());
    log::info!("Probes         : {}", args.probes.display());
    log::info!("Sample         : {}", args.sample.display());
    log::info!(
        "Sample targets : {} ({})",
        sample_ids.len(),
        sample_ids.iter().cloned().collect::<Vec<_>>().join(", ")
    );
    log::info!(
        "CT values      : {:?}",
        sorted_ct_values
    );
    log::info!("Num fragments  : {}", args.num_fragments);
    log::info!("Read length    : {}", args.read_length);
    log::info!(
        "Num sequences  : {}",
        args.num_sequences
            .map(|s| s.to_string())
            .unwrap_or_else(|| "all".to_string())
    );
    log::info!("Output dir     : {}", args.outdir.display());
    log::info!("=============================================");

    // Run pipeline for each CT value and collect depth curves
    let mut all_curves: Vec<DepthCurveRow> = Vec::new();

    for (i, &ct_val) in sorted_ct_values.iter().enumerate() {
        log::info!(
            "--- CT {:.1} ({}/{}) ---",
            ct_val,
            i + 1,
            sorted_ct_values.len()
        );

        let distractor_fraction =
            crate::ct_to_distractor_fraction(ct_val, args.ct_baseline, args.ct_baseline_fraction)?;

        let ct_dir = args.outdir.join(format!("ct_{}", ct_val));
        fs::create_dir_all(&ct_dir)?;

        // Run pipeline: prepare → simulate → capture → [enrich] → sequence → [filter] → map
        run_pipeline(&ct_dir, distractor_fraction, args)?;

        // Compute coverage from SAM
        let sam_path = ct_dir.join("mapped.sam");
        let coverage_result = coverage::compute_coverage(&sam_path)?;

        // Compute depth curves for sample targets only
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
                    // No coverage for this target at this CT — check if we know the ref length
                    log::warn!(
                        "No coverage found for sample target '{}' at CT {:.1}",
                        sample_id,
                        ct_val
                    );
                    vec![(1, 0.0)]
                }
            };

            for (threshold, pct) in curves {
                all_curves.push(DepthCurveRow {
                    ct: ct_val,
                    reference_id: sample_id.clone(),
                    depth_threshold: threshold,
                    pct_covered: pct,
                });
            }
        }
    }

    // Write aggregated depth curves TSV
    let curves_path = args.outdir.join("ct_sweep_depth_curves.tsv");
    write_depth_curves_tsv(&curves_path, &all_curves)?;
    log::info!("Depth curves written to {}", curves_path.display());

    // Generate report
    if args.no_report {
        log::info!("Skipping report generation (--no-report)");
    } else if rscript::check_available() {
        log::info!("Generating CT sweep report...");
        match generate_ct_sweep_report(&curves_path, &sample_ids, &args.outdir) {
            Ok(()) => log::info!(
                "Report generated: {}",
                args.outdir.join("ct_sweep_report.html").display()
            ),
            Err(e) => log::warn!("Report generation failed (non-fatal): {}", e),
        }
    } else {
        log::warn!("Rscript not found — skipping HTML report.");
    }

    log::info!("=============================================");
    log::info!("CT sweep analysis complete!");
    log::info!("Results in {}", args.outdir.display());
    log::info!("=============================================");

    Ok(())
}

fn run_pipeline(outdir: &Path, distractor_fraction: f64, args: &CtSweepArgs) -> Result<()> {
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

    // Step 3b: Optional enrichment
    let capture_output = if let Some(fe) = args.fold_enrichment {
        log::info!("  Applying {:.1}x fold enrichment...", fe);
        enrich::execute(&enrich::EnrichArgs {
            captured: &outdir.join("captured.fa"),
            fragments: &outdir.join("fragments.fa"),
            targets: &outdir.join("targets.txt"),
            distractors: &outdir.join("distractors.txt"),
            fold_enrichment: fe,
            seed: args.seed,
            output: &outdir.join("enriched.fa"),
        })?;
        outdir.join("enriched.fa")
    } else {
        outdir.join("captured.fa")
    };

    // Step 4: Sequence
    log::info!("  Sequencing...");
    sequence::execute(&sequence::SequenceArgs {
        input: &capture_output,
        output: &outdir.join("reads.fa"),
        read_length: args.read_length,
        num_sequences: args.num_sequences,
        seed: args.seed,
    })?;

    // Step 5: Optional host filtering
    let reads_for_mapping = if let Some(host) = args.host_fasta {
        log::info!("  Filtering host reads...");
        filter::execute(&filter::FilterArgs {
            host,
            reads: &outdir.join("reads.fa"),
            minimap_preset: args.host_minimap_preset,
            output: &outdir.join("filtered.fa"),
            log_file: &outdir.join("host_filter.log"),
        })?;
        outdir.join("filtered.fa")
    } else {
        outdir.join("reads.fa")
    };

    // Step 6: Map reads
    log::info!("  Mapping reads...");
    map_reads::execute(&map_reads::MapArgs {
        reference: &outdir.join("combined_reference.fa"),
        reads: &reads_for_mapping,
        minimap_preset: args.minimap_preset,
        output: &outdir.join("mapped.sam"),
        log_file: &outdir.join("mapping.log"),
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

    writeln!(w, "ct\treference_id\tdepth_threshold\tpct_covered")?;
    for row in rows {
        writeln!(
            w,
            "{}\t{}\t{}\t{:.2}",
            row.ct, row.reference_id, row.depth_threshold, row.pct_covered
        )?;
    }

    w.flush()?;
    Ok(())
}

fn generate_ct_sweep_report(
    sweep_tsv: &Path,
    sample_ids: &HashSet<String>,
    outdir: &Path,
) -> Result<()> {
    let r_dir = rscript::find_r_dir()
        .ok_or_else(|| anyhow::anyhow!("Cannot find R scripts directory."))?;

    let script = r_dir.join("ct_sweep.R");
    if !script.exists() {
        bail!(
            "CT sweep R script not found: {}",
            script.display()
        );
    }

    let sweep_abs = fs::canonicalize(sweep_tsv)?;
    let output_path = if outdir.is_absolute() {
        outdir.join("ct_sweep_report.html")
    } else {
        std::env::current_dir()?.join(outdir).join("ct_sweep_report.html")
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

    rscript::run_rscript(
        &script,
        &["--sweep", sweep_str, "--sample-ids", &ids_str, "--output", output_str],
    )
}
