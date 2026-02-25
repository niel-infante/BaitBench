use anyhow::Result;
use std::fs;
use std::path::{Path, PathBuf};

use crate::commands::{capture, filter, generate_list, map_reads, metrics, prepare, report, sequence, simulate};
use crate::external::rscript;

pub struct RunArgs<'a> {
    pub targets: &'a Path,
    pub distractors: &'a [PathBuf],
    pub probes: &'a Path,
    pub sample: Option<&'a Path>,
    pub host_fasta: Option<&'a Path>,
    pub run_name: String,
    pub num_fragments: usize,
    pub distractor_fraction: f64,
    pub seed: Option<u64>,
    pub capture_method: capture::CaptureMethod,
    pub max_mismatches: u32,
    pub min_match_bases: u32,
    pub blast_db: Option<String>,
    pub minimap_preset: String,
    pub host_minimap_preset: String,
    pub fragment_length_mean: f64,
    pub fragment_length_min: usize,
    pub fragment_length_max: usize,
    pub read_length: usize,
    pub outdir: PathBuf,
    pub threads: usize,
    pub no_report: bool,
}

pub fn execute(args: &RunArgs) -> Result<()> {
    let outdir = &args.outdir;
    fs::create_dir_all(outdir)?;

    log::info!("=============================================");
    log::info!("BaitBench - Probe Capture Testing");
    log::info!("=============================================");
    log::info!("Run name            : {}", args.run_name);
    log::info!("Targets FASTA       : {}", args.targets.display());
    for (i, d) in args.distractors.iter().enumerate() {
        log::info!("Distractors FASTA {} : {}", i + 1, d.display());
    }
    log::info!("Probes FASTA        : {}", args.probes.display());
    log::info!(
        "Sample manifest     : {}",
        args.sample
            .map(|p| p.display().to_string())
            .unwrap_or_else(|| "none (all targets)".to_string())
    );
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
    log::info!("Num fragments       : {}", args.num_fragments);
    log::info!("Fragment length     : mean={}, min={}, max={}", args.fragment_length_mean, args.fragment_length_min, args.fragment_length_max);
    log::info!("Read length         : {}", args.read_length);
    log::info!("Distractor fraction : {}", args.distractor_fraction);
    log::info!("Max mismatches      : {}", args.max_mismatches);
    log::info!("Min match bases     : {}", args.min_match_bases);
    log::info!(
        "Seed                : {}",
        args.seed.map(|s| s.to_string()).unwrap_or_else(|| "none".to_string())
    );
    log::info!("Output dir          : {}", outdir.display());
    log::info!("=============================================");

    // Write run parameters file (parameter, CLI flag, value)
    // The flag column enables the report to reconstruct the command generically.
    // When new parameters are added, add them here and the report picks them up
    // automatically — no need to update the R template.
    {
        use std::io::Write;
        let params_path = outdir.join("run_params.tsv");
        let mut f = std::io::BufWriter::new(fs::File::create(&params_path)?);
        writeln!(f, "parameter\tflag\tvalue")?;
        writeln!(f, "targets\t--targets\t{}", args.targets.display())?;
        for d in args.distractors {
            writeln!(f, "distractors\t--distractors\t{}", d.display())?;
        }
        writeln!(f, "probes\t--probes\t{}", args.probes.display())?;
        writeln!(f, "sample\t--sample\t{}", args.sample.map(|p| p.display().to_string()).unwrap_or_else(|| "none".to_string()))?;
        writeln!(f, "host_fasta\t--host-fasta\t{}", args.host_fasta.map(|p| p.display().to_string()).unwrap_or_else(|| "none".to_string()))?;
        writeln!(f, "num_fragments\t--num-fragments\t{}", args.num_fragments)?;
        writeln!(f, "distractor_fraction\t--distractor-fraction\t{}", args.distractor_fraction)?;
        writeln!(f, "capture_method\t--capture-method\t{}", if args.capture_method == capture::CaptureMethod::Blast { "blast" } else { "minimap2" })?;
        writeln!(f, "max_mismatches\t--max-mismatches\t{}", args.max_mismatches)?;
        writeln!(f, "min_match_bases\t--min-match-bases\t{}", args.min_match_bases)?;
        writeln!(f, "blast_db\t--blast-db\t{}", args.blast_db.as_deref().unwrap_or("none"))?;
        writeln!(f, "minimap_preset\t--minimap-preset\t{}", args.minimap_preset)?;
        writeln!(f, "host_minimap_preset\t--host-minimap-preset\t{}", args.host_minimap_preset)?;
        writeln!(f, "fragment_length_mean\t--fragment-length-mean\t{}", args.fragment_length_mean)?;
        writeln!(f, "fragment_length_min\t--fragment-length-min\t{}", args.fragment_length_min)?;
        writeln!(f, "fragment_length_max\t--fragment-length-max\t{}", args.fragment_length_max)?;
        writeln!(f, "read_length\t--read-length\t{}", args.read_length)?;
        writeln!(f, "threads\t--threads\t{}", args.threads)?;
        writeln!(f, "seed\t--seed\t{}", args.seed.map(|s| s.to_string()).unwrap_or_else(|| "none".to_string()))?;
        writeln!(f, "outdir\t--outdir\t{}", outdir.display())?;
        f.flush()?;
    }

    // Step 1: Prepare reference
    log::info!("Step 1/8: Preparing reference...");
    prepare::execute(&prepare::PrepareArgs {
        targets: args.targets,
        distractors: args.distractors,
        sample: args.sample,
        distractor_fraction: args.distractor_fraction,
        outdir,
    })?;

    // Step 2: Simulate fragments
    log::info!("Step 2/8: Simulating fragments...");
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
    log::info!("Step 3/8: Simulating capture...");
    capture::execute(&capture::CaptureArgs {
        method: args.capture_method,
        probes: args.probes,
        fragments: &outdir.join("fragments.fa"),
        max_mismatches: args.max_mismatches,
        min_match_bases: args.min_match_bases,
        blast_db: args.blast_db.as_deref(),
        output: &outdir.join("captured.fa"),
        log_file: &outdir.join("capture.log"),
        threads: args.threads,
    })?;

    // Step 4: Sequence captured fragments into reads
    log::info!("Step 4/8: Sequencing captured fragments...");
    sequence::execute(&sequence::SequenceArgs {
        input: &outdir.join("captured.fa"),
        output: &outdir.join("reads.fa"),
        read_length: args.read_length,
    })?;

    // Step 5: Optional host filtering
    let reads_for_mapping = if let Some(host) = args.host_fasta {
        log::info!("Step 5/8: Filtering host reads...");
        filter::execute(&filter::FilterArgs {
            host,
            reads: &outdir.join("reads.fa"),
            minimap_preset: &args.host_minimap_preset,
            output: &outdir.join("filtered.fa"),
            log_file: &outdir.join("host_filter.log"),
        })?;
        outdir.join("filtered.fa")
    } else {
        log::info!("Step 5/8: Skipping host filtering (no host genome provided)");
        outdir.join("reads.fa")
    };

    // Step 6: Map reads
    log::info!("Step 6/8: Mapping reads to reference...");
    map_reads::execute(&map_reads::MapArgs {
        reference: &outdir.join("combined_reference.fa"),
        reads: &reads_for_mapping,
        minimap_preset: &args.minimap_preset,
        output: &outdir.join("mapped.sam"),
        log_file: &outdir.join("mapping.log"),
    })?;

    // Step 7: Generate detection list
    log::info!("Step 7/8: Generating detection list...");
    generate_list::execute(&generate_list::ListArgs {
        sam: &outdir.join("mapped.sam"),
        output: &outdir.join("detected.list"),
    })?;

    // Step 8: Calculate metrics
    log::info!("Step 8/8: Calculating metrics...");
    let seed_str = args
        .seed
        .map(|s| s.to_string())
        .unwrap_or_else(|| "NA".to_string());
    metrics::execute(&metrics::MetricsArgs {
        targets: &outdir.join("targets.txt"),
        distractors: &outdir.join("distractors.txt"),
        sample: &outdir.join("sample.txt"),
        detected: &outdir.join("detected.list"),
        fragments: &outdir.join("fragments.fa"),
        captured: &outdir.join("captured.fa"),
        sam: &outdir.join("mapped.sam"),
        run_name: &args.run_name,
        num_fragments: args.num_fragments,
        seed: &seed_str,
        output_summary: &outdir.join("results.tsv"),
        output_detail: &outdir.join("detected_detail.tsv"),
        output_json: Some(&outdir.join("results.json")),
    })?;

    // Generate report (optional — skip if R not available or --no-report)
    if args.no_report {
        log::info!("Skipping report generation (--no-report)");
    } else if rscript::check_available() {
        log::info!("Generating report...");
        match report::execute(&report::ReportArgs {
            summary: &outdir.join("results.tsv"),
            detail: &outdir.join("detected_detail.tsv"),
            params: &outdir.join("run_params.tsv"),
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
