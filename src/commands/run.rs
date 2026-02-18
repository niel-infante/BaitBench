use anyhow::Result;
use std::fs;
use std::path::{Path, PathBuf};

use crate::commands::{capture, filter, generate_list, map_reads, metrics, prepare, report, simulate};
use crate::external::rscript;

pub struct RunArgs<'a> {
    pub targets: &'a Path,
    pub distractors: &'a Path,
    pub probes: &'a Path,
    pub host_fasta: Option<&'a Path>,
    pub run_name: String,
    pub num_reads: usize,
    pub distractor_fraction: f64,
    pub seed: Option<u64>,
    pub capture_method: capture::CaptureMethod,
    pub max_mismatches: u32,
    pub min_match_bases: u32,
    pub blast_db: Option<String>,
    pub minimap_preset: String,
    pub host_minimap_preset: String,
    pub read_length_mean: f64,
    pub read_length_min: usize,
    pub read_length_max: usize,
    pub outdir: PathBuf,
    pub threads: usize,
}

pub fn execute(args: &RunArgs) -> Result<()> {
    let outdir = &args.outdir;
    fs::create_dir_all(outdir)?;

    log::info!("=============================================");
    log::info!("BaitBench - Probe Capture Testing");
    log::info!("=============================================");
    log::info!("Run name            : {}", args.run_name);
    log::info!("Targets FASTA       : {}", args.targets.display());
    log::info!("Distractors FASTA   : {}", args.distractors.display());
    log::info!("Probes FASTA        : {}", args.probes.display());
    log::info!(
        "Capture method      : {}",
        if args.capture_method == capture::CaptureMethod::Blast { "blast" } else { "minimap2" }
    );
    log::info!(
        "Host FASTA          : {}",
        args.host_fasta
            .map(|p| p.display().to_string())
            .unwrap_or_else(|| "none (skip host filtering)".to_string())
    );
    log::info!("Num reads           : {}", args.num_reads);
    log::info!("Distractor fraction : {}", args.distractor_fraction);
    log::info!("Max mismatches      : {}", args.max_mismatches);
    log::info!("Min match bases     : {}", args.min_match_bases);
    log::info!(
        "Seed                : {}",
        args.seed.map(|s| s.to_string()).unwrap_or_else(|| "random".to_string())
    );
    log::info!("Output dir          : {}", outdir.display());
    log::info!("=============================================");

    // Step 1: Prepare reference
    log::info!("Step 1/7: Preparing reference...");
    prepare::execute(&prepare::PrepareArgs {
        targets: args.targets,
        distractors: args.distractors,
        distractor_fraction: args.distractor_fraction,
        outdir,
    })?;

    // Step 2: Simulate reads
    log::info!("Step 2/7: Simulating reads...");
    simulate::execute(&simulate::SimulateArgs {
        reference: &outdir.join("combined_reference.fa"),
        weights: &outdir.join("weights.txt"),
        num_reads: args.num_reads,
        seed: args.seed,
        output: &outdir.join("reads.fa"),
        read_length_mean: args.read_length_mean,
        read_length_min: args.read_length_min,
        read_length_max: args.read_length_max,
    })?;

    // Step 3: Capture
    log::info!("Step 3/7: Simulating capture...");
    capture::execute(&capture::CaptureArgs {
        method: args.capture_method,
        probes: args.probes,
        reads: &outdir.join("reads.fa"),
        max_mismatches: args.max_mismatches,
        min_match_bases: args.min_match_bases,
        blast_db: args.blast_db.as_deref(),
        output: &outdir.join("captured.fa"),
        log_file: &outdir.join("capture.log"),
        threads: args.threads,
    })?;

    // Step 4: Optional host filtering
    let reads_for_mapping = if let Some(host) = args.host_fasta {
        log::info!("Step 4/7: Filtering host reads...");
        filter::execute(&filter::FilterArgs {
            host,
            reads: &outdir.join("captured.fa"),
            minimap_preset: &args.host_minimap_preset,
            output: &outdir.join("filtered.fa"),
            log_file: &outdir.join("host_filter.log"),
        })?;
        outdir.join("filtered.fa")
    } else {
        log::info!("Step 4/7: Skipping host filtering (no host genome provided)");
        outdir.join("captured.fa")
    };

    // Step 5: Map reads
    log::info!("Step 5/7: Mapping reads to reference...");
    map_reads::execute(&map_reads::MapArgs {
        reference: &outdir.join("combined_reference.fa"),
        reads: &reads_for_mapping,
        minimap_preset: &args.minimap_preset,
        output: &outdir.join("mapped.sam"),
        log_file: &outdir.join("mapping.log"),
    })?;

    // Step 6: Generate detection list
    log::info!("Step 6/7: Generating detection list...");
    generate_list::execute(&generate_list::ListArgs {
        sam: &outdir.join("mapped.sam"),
        output: &outdir.join("detected.list"),
    })?;

    // Step 7: Calculate metrics
    log::info!("Step 7/7: Calculating metrics...");
    let seed_str = args
        .seed
        .map(|s| s.to_string())
        .unwrap_or_else(|| "NA".to_string());
    metrics::execute(&metrics::MetricsArgs {
        targets: &outdir.join("targets.txt"),
        distractors: &outdir.join("distractors.txt"),
        detected: &outdir.join("detected.list"),
        reads: &outdir.join("reads.fa"),
        captured: &outdir.join("captured.fa"),
        run_name: &args.run_name,
        num_reads: args.num_reads,
        seed: &seed_str,
        output_summary: &outdir.join("results.tsv"),
        output_detail: &outdir.join("detected_detail.tsv"),
        output_json: Some(&outdir.join("results.json")),
    })?;

    // Generate report (optional — skip if R not available)
    if rscript::check_available() {
        log::info!("Generating report...");
        match report::execute(&report::ReportArgs {
            summary: &outdir.join("results.tsv"),
            detail: &outdir.join("detected_detail.tsv"),
            run_name: &args.run_name,
            output: &outdir.join("report.html"),
        }) {
            Ok(()) => {}
            Err(e) => log::warn!("Report generation failed (non-fatal): {}", e),
        }
    } else {
        log::warn!(
            "Rscript not found — skipping HTML report. Install R to enable report generation."
        );
    }

    log::info!("=============================================");
    log::info!("Pipeline complete!");
    log::info!("Results in {}", outdir.display());
    log::info!("=============================================");

    Ok(())
}
