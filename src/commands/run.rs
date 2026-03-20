use anyhow::Result;
use std::collections::HashMap;
use std::fs;
use std::path::{Path, PathBuf};

use crate::cleanup;
use crate::cli::ReportMode;
use crate::commands::{capture, enrich, filter, generate_list, identify, map_reads, metrics, prepare, report, sequence, simulate};
use crate::external::rscript;
use crate::fasta;
use crate::io_utils;
use crate::io_utils::prefixed_join;

pub struct RunArgs<'a> {
    pub targets: &'a Path,
    pub genomes: Option<&'a Path>,
    pub distractors: &'a [PathBuf],
    pub probes: &'a Path,
    pub sample: Option<&'a HashMap<String, f64>>,
    pub sample_target_map: Option<&'a HashMap<String, Vec<String>>>,
    pub host_fasta: Option<&'a Path>,
    pub run_name: String,
    pub num_fragments: usize,
    pub distractor_fraction: f64,
    pub ct: Option<f64>,
    pub ct_baseline: f64,
    pub ct_baseline_fraction: f64,
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
    pub num_sequences: Option<usize>,
    pub outdir: PathBuf,
    pub threads: usize,
    pub fold_enrichment: Option<f64>,
    pub identify: bool,
    pub identity_threshold: f64,
    pub min_unique_targets: usize,
    pub report: ReportMode,
    pub cleanup: bool,
    pub output_prefix: String,
}

pub fn execute(args: &RunArgs) -> Result<()> {
    let outdir = &args.outdir;
    fs::create_dir_all(outdir)?;

    let has_genomes = args.genomes.is_some();

    log::info!("=============================================");
    log::info!("BaitBench - Probe Capture Testing");
    log::info!("=============================================");
    log::info!("Run name            : {}", args.run_name);
    log::info!("Targets FASTA       : {}", args.targets.display());
    if let Some(genomes) = args.genomes {
        log::info!("Genomes FASTA       : {}", genomes.display());
    }
    for (i, d) in args.distractors.iter().enumerate() {
        log::info!("Distractors FASTA {} : {}", i + 1, d.display());
    }
    log::info!("Probes FASTA        : {}", args.probes.display());
    log::info!(
        "Sample              : {}",
        args.sample
            .map(|s| io_utils::format_sample_display(s))
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
    log::info!(
        "Num sequences       : {}",
        args.num_sequences.map(|s| s.to_string()).unwrap_or_else(|| "all".to_string())
    );
    log::info!("Distractor fraction : {}", args.distractor_fraction);
    if let Some(ct) = args.ct {
        log::info!("CT score            : {}", ct);
        log::info!("CT baseline         : {} (fraction {})", args.ct_baseline, args.ct_baseline_fraction);
    }
    log::info!("Max mismatches      : {}", args.max_mismatches);
    log::info!("Min match bases     : {}", args.min_match_bases);
    log::info!(
        "Fold enrichment     : {}",
        args.fold_enrichment.map(|f| format!("{:.1}x", f)).unwrap_or_else(|| "none (binary capture)".to_string())
    );
    log::info!(
        "Seed                : {}",
        args.seed.map(|s| s.to_string()).unwrap_or_else(|| "none".to_string())
    );
    log::info!("Output dir          : {}", outdir.display());
    log::info!("=============================================");

    let pfx = &args.output_prefix;

    // Write run parameters file (parameter, CLI flag, value)
    {
        use std::io::Write;
        let params_path = prefixed_join(outdir, pfx, "run_params.tsv");
        let mut f = std::io::BufWriter::new(fs::File::create(&params_path)?);
        writeln!(f, "parameter\tflag\tvalue")?;
        writeln!(f, "targets\t--targets\t{}", args.targets.display())?;
        if let Some(genomes) = args.genomes {
            writeln!(f, "genomes\t--genomes\t{}", genomes.display())?;
        }
        for d in args.distractors {
            writeln!(f, "distractors\t--distractors\t{}", d.display())?;
        }
        writeln!(f, "probes\t--probes\t{}", args.probes.display())?;
        writeln!(f, "sample\t--sample\t{}", args.sample.map(|s| io_utils::format_sample_display(s)).unwrap_or_else(|| "none".to_string()))?;
        writeln!(f, "host_fasta\t--host-fasta\t{}", args.host_fasta.map(|p| p.display().to_string()).unwrap_or_else(|| "none".to_string()))?;
        writeln!(f, "num_fragments\t--num-fragments\t{}", args.num_fragments)?;
        writeln!(f, "distractor_fraction\t--distractor-fraction\t{}", args.distractor_fraction)?;
        writeln!(f, "ct\t--ct\t{}", args.ct.map(|v| v.to_string()).unwrap_or_else(|| "none".to_string()))?;
        writeln!(f, "ct_baseline\t--ct-baseline\t{}", args.ct_baseline)?;
        writeln!(f, "ct_baseline_fraction\t--ct-baseline-fraction\t{}", args.ct_baseline_fraction)?;
        writeln!(f, "capture_method\t--capture-method\t{}", if args.capture_method == capture::CaptureMethod::Blast { "blast" } else { "minimap2" })?;
        writeln!(f, "max_mismatches\t--max-mismatches\t{}", args.max_mismatches)?;
        writeln!(f, "min_match_bases\t--min-match-bases\t{}", args.min_match_bases)?;
        writeln!(f, "fold_enrichment\t--fold-enrichment\t{}", args.fold_enrichment.map(|f| f.to_string()).unwrap_or_else(|| "none".to_string()))?;
        writeln!(f, "blast_db\t--blast-db\t{}", args.blast_db.as_deref().unwrap_or("none"))?;
        writeln!(f, "minimap_preset\t--minimap-preset\t{}", args.minimap_preset)?;
        writeln!(f, "host_minimap_preset\t--host-minimap-preset\t{}", args.host_minimap_preset)?;
        writeln!(f, "fragment_length_mean\t--fragment-length-mean\t{}", args.fragment_length_mean)?;
        writeln!(f, "fragment_length_min\t--fragment-length-min\t{}", args.fragment_length_min)?;
        writeln!(f, "fragment_length_max\t--fragment-length-max\t{}", args.fragment_length_max)?;
        writeln!(f, "read_length\t--read-length\t{}", args.read_length)?;
        writeln!(f, "num_sequences\t--num-sequences\t{}", args.num_sequences.map(|s| s.to_string()).unwrap_or_else(|| "none".to_string()))?;
        writeln!(f, "threads\t--threads\t{}", args.threads)?;
        writeln!(f, "seed\t--seed\t{}", args.seed.map(|s| s.to_string()).unwrap_or_else(|| "none".to_string()))?;
        writeln!(f, "outdir\t--outdir\t{}", outdir.display())?;
        f.flush()?;
    }

    // Step 1: Prepare reference
    log::info!("Step 1/8: Preparing reference...");
    prepare::execute(&prepare::PrepareArgs {
        targets: args.targets,
        genomes: args.genomes,
        distractors: args.distractors,
        sample: args.sample,
        sample_target_map: args.sample_target_map,
        distractor_fraction: args.distractor_fraction,
        outdir,
        output_prefix: pfx,
    })?;

    // Step 2: Simulate fragments
    log::info!("Step 2/8: Simulating fragments...");
    simulate::execute(&simulate::SimulateArgs {
        reference: &prefixed_join(outdir, pfx, "combined_reference.fa"),
        weights: &prefixed_join(outdir, pfx, "weights.txt"),
        num_fragments: args.num_fragments,
        seed: args.seed,
        output: &prefixed_join(outdir, pfx, "fragments.fa"),
        fragment_length_mean: args.fragment_length_mean,
        fragment_length_min: args.fragment_length_min,
        fragment_length_max: args.fragment_length_max,
    })?;

    // Step 3: Capture
    log::info!("Step 3/8: Simulating capture...");
    capture::execute(&capture::CaptureArgs {
        method: args.capture_method,
        probes: args.probes,
        fragments: &prefixed_join(outdir, pfx, "fragments.fa"),
        max_mismatches: args.max_mismatches,
        min_match_bases: args.min_match_bases,
        blast_db: args.blast_db.as_deref(),
        output: &prefixed_join(outdir, pfx, "captured.fa"),
        log_file: &prefixed_join(outdir, pfx, "capture.log"),
        threads: args.threads,
    })?;

    // Step 3b: Optional fold enrichment adjustment
    // In genomes mode, enrich classifies by genome IDs (genomes.txt), not target IDs
    let enrich_targets_file = if has_genomes {
        prefixed_join(outdir, pfx, "genomes.txt")
    } else {
        prefixed_join(outdir, pfx, "targets.txt")
    };

    let capture_output = if let Some(fe) = args.fold_enrichment {
        log::info!("Step 3b: Applying {:.1}x fold enrichment...", fe);
        enrich::execute(&enrich::EnrichArgs {
            captured: &prefixed_join(outdir, pfx, "captured.fa"),
            fragments: &prefixed_join(outdir, pfx, "fragments.fa"),
            targets: &enrich_targets_file,
            distractors: &prefixed_join(outdir, pfx, "distractors.txt"),
            fold_enrichment: fe,
            seed: args.seed,
            output: &prefixed_join(outdir, pfx, "enriched.fa"),
        })?;
        prefixed_join(outdir, pfx, "enriched.fa")
    } else {
        prefixed_join(outdir, pfx, "captured.fa")
    };

    // Step 4: Sequence captured fragments into reads
    log::info!("Step 4/8: Sequencing captured fragments...");
    sequence::execute(&sequence::SequenceArgs {
        input: &capture_output,
        output: &prefixed_join(outdir, pfx, "reads.fa"),
        read_length: args.read_length,
        num_sequences: args.num_sequences,
        seed: args.seed,
    })?;
    let reads_sequenced = fasta::count_sequences(&prefixed_join(outdir, pfx, "reads.fa"))?;
    log::info!("  Reads after sequencing: {}", reads_sequenced);

    // Step 5: Optional host filtering
    let (reads_for_mapping, reads_after_filter) = if let Some(host) = args.host_fasta {
        log::info!("Step 5/8: Filtering host reads...");
        filter::execute(&filter::FilterArgs {
            host,
            reads: &prefixed_join(outdir, pfx, "reads.fa"),
            minimap_preset: &args.host_minimap_preset,
            output: &prefixed_join(outdir, pfx, "filtered.fa"),
            log_file: &prefixed_join(outdir, pfx, "host_filter.log"),
        })?;
        let count = fasta::count_sequences(&prefixed_join(outdir, pfx, "filtered.fa"))?;
        log::info!("  Reads after host filtering: {}", count);
        (prefixed_join(outdir, pfx, "filtered.fa"), Some(count))
    } else {
        log::info!("Step 5/8: Skipping host filtering (no host genome provided)");
        (prefixed_join(outdir, pfx, "reads.fa"), None)
    };

    // Step 6: Map reads
    // In genomes mode, map to mapping_reference.fa (targets + distractors)
    // In standard mode, map to combined_reference.fa (targets + distractors, same thing)
    let mapping_reference = if has_genomes {
        prefixed_join(outdir, pfx, "mapping_reference.fa")
    } else {
        prefixed_join(outdir, pfx, "combined_reference.fa")
    };

    log::info!("Step 6/8: Mapping reads to reference...");
    map_reads::execute(&map_reads::MapArgs {
        reference: &mapping_reference,
        reads: &reads_for_mapping,
        minimap_preset: &args.minimap_preset,
        output: &prefixed_join(outdir, pfx, "mapped.sam"),
        log_file: &prefixed_join(outdir, pfx, "mapping.log"),
    })?;

    // Step 7: Generate detection list
    log::info!("Step 7/8: Generating detection list...");
    generate_list::execute(&generate_list::ListArgs {
        sam: &prefixed_join(outdir, pfx, "mapped.sam"),
        output: &prefixed_join(outdir, pfx, "detected.list"),
    })?;

    // Step 8: Calculate metrics and coverage
    log::info!("Step 8/8: Calculating metrics and coverage...");
    let seed_str = args
        .seed
        .map(|s| s.to_string())
        .unwrap_or_else(|| "NA".to_string());

    let sample_target_map_path = if has_genomes {
        Some(prefixed_join(outdir, pfx, "sample_target_map.txt"))
    } else {
        None
    };

    metrics::execute(&metrics::MetricsArgs {
        targets: &prefixed_join(outdir, pfx, "targets.txt"),
        distractors: &prefixed_join(outdir, pfx, "distractors.txt"),
        sample: &prefixed_join(outdir, pfx, "sample.txt"),
        sample_target_map: sample_target_map_path.as_deref(),
        detected: &prefixed_join(outdir, pfx, "detected.list"),
        fragments: &prefixed_join(outdir, pfx, "fragments.fa"),
        captured: &capture_output,
        sam: &prefixed_join(outdir, pfx, "mapped.sam"),
        run_name: &args.run_name,
        num_fragments: args.num_fragments,
        seed: &seed_str,
        output_summary: &prefixed_join(outdir, pfx, "results.tsv"),
        output_detail: &prefixed_join(outdir, pfx, "detected_detail.tsv"),
        output_json: Some(&prefixed_join(outdir, pfx, "results.json")),
        output_coverage: Some(&prefixed_join(outdir, pfx, "coverage.tsv")),
        reads_sequenced: Some(reads_sequenced),
        reads_after_filter,
    })?;

    // Optional: Species identification
    if args.identify {
        if !has_genomes || args.sample_target_map.is_none() {
            log::warn!(
                "--identify requires genome mode (--genomes) and --sample-target-map. Skipping species identification."
            );
        } else {
            log::info!("Step 9: Species identification...");
            match identify::execute(&identify::IdentifyArgs {
                detected_detail: &prefixed_join(outdir, pfx, "detected_detail.tsv"),
                sample_target_map: &prefixed_join(outdir, pfx, "sample_target_map.txt"),
                target_similarity_file: None,
                targets_fasta: Some(args.targets),
                identity_threshold: args.identity_threshold,
                minimap_preset: &args.minimap_preset,
                min_unique_targets: args.min_unique_targets,
                outdir,
                output_prefix: pfx,
            }) {
                Ok(()) => {}
                Err(e) => log::warn!("Species identification failed (non-fatal): {}", e),
            }
        }
    }

    // Determine species calls file path (if identify step ran)
    let species_calls_path = prefixed_join(outdir, pfx, "species_calls.tsv");
    let species_calls_opt = if species_calls_path.exists() {
        Some(species_calls_path.as_path())
    } else {
        None
    };

    // Generate report
    match args.report {
        ReportMode::None => {
            log::info!("Skipping report generation (--report none)");
        }
        ReportMode::Full => {
            if rscript::check_available() {
                log::info!("Generating report...");
                match report::execute(&report::ReportArgs {
                    summary: &prefixed_join(outdir, pfx, "results.tsv"),
                    detail: &prefixed_join(outdir, pfx, "detected_detail.tsv"),
                    params: &prefixed_join(outdir, pfx, "run_params.tsv"),
                    coverage: Some(&prefixed_join(outdir, pfx, "coverage.tsv")),
                    species_calls: species_calls_opt,
                    run_name: &args.run_name,
                    output: &prefixed_join(outdir, pfx, "report.html"),
                    report: ReportMode::Full,
                }) {
                    Ok(()) => {}
                    Err(e) => log::warn!("Report generation failed (non-fatal): {}", e),
                }
            } else {
                log::warn!(
                    "Rscript not found — skipping HTML report. Install R to enable report generation."
                );
            }
        }
        ReportMode::Rmd => {
            log::info!("Generating RMarkdown file...");
            match report::execute(&report::ReportArgs {
                summary: &prefixed_join(outdir, pfx, "results.tsv"),
                detail: &prefixed_join(outdir, pfx, "detected_detail.tsv"),
                params: &prefixed_join(outdir, pfx, "run_params.tsv"),
                coverage: Some(&prefixed_join(outdir, pfx, "coverage.tsv")),
                species_calls: species_calls_opt,
                run_name: &args.run_name,
                output: &prefixed_join(outdir, pfx, "report.html"),
                report: ReportMode::Rmd,
            }) {
                Ok(()) => {}
                Err(e) => log::warn!("RMarkdown generation failed (non-fatal): {}", e),
            }
        }
    }

    // Cleanup intermediate files if requested
    if args.cleanup {
        log::info!("Cleaning up intermediate files...");
        let cleanup_names: Vec<String> = [
            "combined_reference.fa",
            "mapping_reference.fa",
            "weights.txt",
            "targets.txt",
            "distractors.txt",
            "sample.txt",
            "genomes.txt",
            "sample_target_map.txt",
            "fragments.fa",
            "captured.fa",
            "enriched.fa",
            "reads.fa",
            "filtered.fa",
            "mapped.sam",
            "detected.list",
            "capture.log",
            "mapping.log",
            "host_filter.log",
        ]
        .iter()
        .map(|f| format!("{}{}", pfx, f))
        .collect();
        let cleanup_refs: Vec<&str> = cleanup_names.iter().map(|s| s.as_str()).collect();
        cleanup::cleanup_files(outdir, &cleanup_refs);
    }

    log::info!("=============================================");
    log::info!("Pipeline complete!");
    log::info!("Results in {}", outdir.display());
    log::info!("=============================================");

    Ok(())
}
