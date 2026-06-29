use anyhow::Result;
use std::collections::HashMap;
use std::fs;
use std::path::{Path, PathBuf};

use crate::cleanup;
use crate::cli::ReportMode;
use crate::commands::{filter, generate_list, identify, map_reads, metrics, prepare, report, sequence, simulate};
use crate::commands::sequence::ReadSimulator;
use crate::external::rscript;
use crate::fasta;
use crate::io_utils;
use crate::io_utils::prefixed_join;
use crate::sampling::thermo_sim::SimulateMode;

pub struct RunArgs<'a> {
    pub targets: &'a Path,
    pub genomes: Option<&'a Path>,
    pub distractors: &'a [PathBuf],
    pub probes: &'a Path,
    pub sample: Option<&'a HashMap<String, f64>>,
    pub sample_target_map: Option<&'a HashMap<String, Vec<String>>>,
    /// Optional target groups file (seq_id → group_name TSV).
    pub groups: Option<&'a Path>,
    /// Optional distractor groups file (contig_id → group_name TSV).
    pub distractor_groups: Option<&'a Path>,
    pub host_fasta: Option<&'a Path>,
    pub run_name: String,
    pub num_fragments: usize,
    pub distractor_fraction: f64,
    pub ct: Option<f64>,
    pub ct_baseline: f64,
    pub ct_baseline_fraction: f64,
    pub ct_efficiency: f64,
    pub ct_calibration: Option<Vec<String>>,
    pub seed: Option<u64>,
    pub simulate_mode: SimulateMode,
    pub hybridization_temperature: f64,
    /// Na+ concentration in molar (converted from CLI mM)
    pub na_conc_m: f64,
    pub capture_fraction: f64,
    pub minimap_preset: String,
    pub host_minimap_preset: String,
    pub fragment_length_mean: f64,
    pub fragment_length_min: usize,
    pub fragment_length_max: usize,
    pub read_length: usize,
    pub num_sequences: Option<usize>,
    pub simulator: ReadSimulator,
    pub sequencer_profile: String,
    pub coverage_depth: f64,
    pub paired_end: bool,
    pub pe_frag_len_mean: usize,
    pub pe_frag_len_sd: usize,
    pub outdir: PathBuf,
    pub threads: usize,
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

    let sim_mode_str = match args.simulate_mode {
        SimulateMode::Thermodynamic => "thermodynamic",
        SimulateMode::Simple => "simple",
    };

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
    log::info!("Simulate mode       : {}", sim_mode_str);
    if matches!(args.simulate_mode, SimulateMode::Thermodynamic) {
        log::info!("Hybridization temp  : {}°C", args.hybridization_temperature);
        log::info!("Salt concentration  : {} mM Na+", args.na_conc_m * 1000.0);
    }
    log::info!("Capture fraction    : {}", args.capture_fraction);
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
        if let Some(cal) = &args.ct_calibration {
            log::info!("CT calibration      : {} / {}", cal[0], cal[1]);
        } else {
            log::info!("CT baseline         : {} (fraction {})", args.ct_baseline, args.ct_baseline_fraction);
            log::info!("CT efficiency       : {}", args.ct_efficiency);
        }
    }
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
        if let Some(cal) = &args.ct_calibration {
            writeln!(f, "ct_calibration_1\t--ct-calibration\t{}", cal[0])?;
            writeln!(f, "ct_calibration_2\t--ct-calibration\t{}", cal[1])?;
        } else {
            writeln!(f, "ct_baseline\t--ct-baseline\t{}", args.ct_baseline)?;
            writeln!(f, "ct_baseline_fraction\t--ct-baseline-fraction\t{}", args.ct_baseline_fraction)?;
            writeln!(f, "ct_efficiency\t--ct-efficiency\t{}", args.ct_efficiency)?;
        }
        writeln!(f, "simulate_mode\t--simulate-mode\t{}", sim_mode_str)?;
        writeln!(f, "hybridization_temperature\t--hybridization-temperature\t{}", args.hybridization_temperature)?;
        writeln!(f, "salt_concentration\t--salt-concentration\t{}", args.na_conc_m * 1000.0)?;
        writeln!(f, "capture_fraction\t--capture-fraction\t{}", args.capture_fraction)?;
        writeln!(f, "minimap_preset\t--minimap-preset\t{}", args.minimap_preset)?;
        writeln!(f, "host_minimap_preset\t--host-minimap-preset\t{}", args.host_minimap_preset)?;
        writeln!(f, "fragment_length_mean\t--fragment-length-mean\t{}", args.fragment_length_mean)?;
        writeln!(f, "fragment_length_min\t--fragment-length-min\t{}", args.fragment_length_min)?;
        writeln!(f, "fragment_length_max\t--fragment-length-max\t{}", args.fragment_length_max)?;
        writeln!(f, "read_length\t--read-length\t{}", args.read_length)?;
        writeln!(f, "num_sequences\t--num-sequences\t{}", args.num_sequences.map(|s| s.to_string()).unwrap_or_else(|| "none".to_string()))?;
        let sim_name = match args.simulator {
            ReadSimulator::Perfect => "perfect",
            ReadSimulator::Art     => "art",
            ReadSimulator::Badread => "badread",
        };
        writeln!(f, "read_simulator\t--read-simulator\t{}", sim_name)?;
        writeln!(f, "sequencer_profile\t--sequencer-profile\t{}", args.sequencer_profile)?;
        writeln!(f, "coverage_depth\t--coverage-depth\t{}", args.coverage_depth)?;
        writeln!(f, "paired_end\t--paired-end\t{}", args.paired_end)?;
        writeln!(f, "pe_frag_len_mean\t--pe-frag-len-mean\t{}", args.pe_frag_len_mean)?;
        writeln!(f, "pe_frag_len_sd\t--pe-frag-len-sd\t{}", args.pe_frag_len_sd)?;
        writeln!(f, "threads\t--threads\t{}", args.threads)?;
        writeln!(f, "seed\t--seed\t{}", args.seed.map(|s| s.to_string()).unwrap_or_else(|| "none".to_string()))?;
        writeln!(f, "outdir\t--outdir\t{}", outdir.display())?;
        f.flush()?;
    }

    // Step 1: Prepare reference
    log::info!("Step 1/7: Preparing reference...");
    prepare::execute(&prepare::PrepareArgs {
        targets: args.targets,
        genomes: args.genomes,
        distractors: args.distractors,
        sample: args.sample,
        sample_target_map: args.sample_target_map,
        groups: args.groups,
        distractor_groups: args.distractor_groups,
        distractor_fraction: args.distractor_fraction,
        outdir,
        output_prefix: pfx,
    })?;

    // Step 2: Simulate fragments (probe-biased + background)
    log::info!("Step 2/7: Simulating fragments (mode={}, capture_fraction={:.2})...",
        sim_mode_str, args.capture_fraction);
    simulate::execute(&simulate::SimulateArgs {
        reference: &prefixed_join(outdir, pfx, "combined_reference.fa"),
        weights: &prefixed_join(outdir, pfx, "weights.txt"),
        probes: args.probes,
        num_fragments: args.num_fragments,
        capture_fraction: args.capture_fraction,
        simulate_mode: args.simulate_mode,
        hybridization_temperature: args.hybridization_temperature,
        na_conc_m: args.na_conc_m,
        seed: args.seed,
        output: &prefixed_join(outdir, pfx, "fragments.fa"),
        fragment_length_mean: args.fragment_length_mean,
        fragment_length_min: args.fragment_length_min,
        fragment_length_max: args.fragment_length_max,
        threads: args.threads,
    })?;

    // Step 3: Sequence fragments into reads
    log::info!("Step 3/7: Sequencing fragments (simulator={:?})...", args.simulator);
    let reads_r1_path = prefixed_join(outdir, pfx, "reads.fa");
    let reads_r2_path = if args.paired_end {
        Some(prefixed_join(outdir, pfx, "reads_R2.fa"))
    } else {
        None
    };
    sequence::execute(&sequence::SequenceArgs {
        input: &prefixed_join(outdir, pfx, "fragments.fa"),
        output: &reads_r1_path,
        output_r2: reads_r2_path.as_deref(),
        read_length: args.read_length,
        num_sequences: args.num_sequences,
        seed: args.seed,
        simulator: args.simulator,
        sequencer_profile: args.sequencer_profile.clone(),
        coverage_depth: args.coverage_depth,
        paired_end: args.paired_end,
        pe_frag_len_mean: args.pe_frag_len_mean,
        pe_frag_len_sd: args.pe_frag_len_sd,
    })?;
    let reads_sequenced = fasta::count_sequences(&reads_r1_path)?;
    log::info!("  Reads after sequencing: {}", reads_sequenced);

    // Step 4: Optional host filtering
    let filtered_r2_path = reads_r2_path.as_ref().map(|_| prefixed_join(outdir, pfx, "filtered_R2.fa"));
    let (reads_for_mapping, reads_for_mapping_r2, reads_after_filter) = if let Some(host) = args.host_fasta {
        log::info!("Step 4/7: Filtering host reads...");
        filter::execute(&filter::FilterArgs {
            host,
            reads: &reads_r1_path,
            reads_r2: reads_r2_path.as_deref(),
            minimap_preset: &args.host_minimap_preset,
            output: &prefixed_join(outdir, pfx, "filtered.fa"),
            output_r2: filtered_r2_path.as_deref(),
            log_file: &prefixed_join(outdir, pfx, "host_filter.log"),
        })?;
        let count = fasta::count_sequences(&prefixed_join(outdir, pfx, "filtered.fa"))?;
        log::info!("  Reads after host filtering: {}", count);
        (prefixed_join(outdir, pfx, "filtered.fa"), filtered_r2_path, Some(count))
    } else {
        log::info!("Step 4/7: Skipping host filtering (no host genome provided)");
        (reads_r1_path.clone(), reads_r2_path.clone(), None)
    };

    // Step 5: Map reads
    // In genomes mode, map to mapping_reference.fa (targets + distractors)
    // In standard mode, map to combined_reference.fa (targets + distractors, same thing)
    let mapping_reference = if has_genomes {
        prefixed_join(outdir, pfx, "mapping_reference.fa")
    } else {
        prefixed_join(outdir, pfx, "combined_reference.fa")
    };

    log::info!("Step 5/7: Mapping reads to reference...");
    map_reads::execute(&map_reads::MapArgs {
        reference: &mapping_reference,
        reads: &reads_for_mapping,
        reads_r2: reads_for_mapping_r2.as_deref(),
        minimap_preset: &args.minimap_preset,
        output: &prefixed_join(outdir, pfx, "mapped.sam"),
        log_file: &prefixed_join(outdir, pfx, "mapping.log"),
    })?;

    // Step 6: Generate detection list
    log::info!("Step 6/7: Generating detection list...");
    generate_list::execute(&generate_list::ListArgs {
        sam: &prefixed_join(outdir, pfx, "mapped.sam"),
        output: &prefixed_join(outdir, pfx, "detected.list"),
    })?;

    // Step 7: Calculate metrics and coverage
    log::info!("Step 7/7: Calculating metrics and coverage...");
    let seed_str = args
        .seed
        .map(|s| s.to_string())
        .unwrap_or_else(|| "NA".to_string());

    let sample_target_map_path = if has_genomes {
        Some(prefixed_join(outdir, pfx, "sample_target_map.txt"))
    } else {
        None
    };

    // Target groups file written by prepare (only if --groups was provided)
    let target_groups_path = prefixed_join(outdir, pfx, "target_groups.tsv");
    let target_groups_opt = if target_groups_path.exists() {
        Some(target_groups_path)
    } else {
        None
    };

    // Distractor groups file is always written by prepare
    let distractor_groups_path = prefixed_join(outdir, pfx, "distractor_groups.tsv");
    let distractor_groups_opt = if distractor_groups_path.exists() {
        Some(distractor_groups_path)
    } else {
        None
    };

    // Write group detail only when any grouping is active
    let group_detail_path = prefixed_join(outdir, pfx, "group_detail.tsv");
    let group_detail_opt = if target_groups_opt.is_some() || distractor_groups_opt.is_some() {
        Some(group_detail_path)
    } else {
        None
    };

    let fragments_path = prefixed_join(outdir, pfx, "fragments.fa");

    metrics::execute(&metrics::MetricsArgs {
        targets: &prefixed_join(outdir, pfx, "targets.txt"),
        distractors: &prefixed_join(outdir, pfx, "distractors.txt"),
        sample: &prefixed_join(outdir, pfx, "sample.txt"),
        sample_target_map: sample_target_map_path.as_deref(),
        target_groups: target_groups_opt.as_deref(),
        distractor_groups: distractor_groups_opt.as_deref(),
        detected: &prefixed_join(outdir, pfx, "detected.list"),
        fragments: &fragments_path,
        captured: &fragments_path, // fragments IS the post-capture pool in the new pipeline
        sam: &prefixed_join(outdir, pfx, "mapped.sam"),
        run_name: &args.run_name,
        num_fragments: args.num_fragments,
        seed: &seed_str,
        output_summary: &prefixed_join(outdir, pfx, "results.tsv"),
        output_detail: &prefixed_join(outdir, pfx, "detected_detail.tsv"),
        output_group_detail: group_detail_opt.as_deref(),
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
            log::info!("Step 8: Species identification...");
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
                    group_detail: group_detail_opt.as_deref(),
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
                group_detail: group_detail_opt.as_deref(),
                run_name: &args.run_name,
                output: &prefixed_join(outdir, pfx, "report.html"),
                report: ReportMode::Rmd,
            }) {
                Ok(()) => {}
                Err(e) => log::warn!("RMarkdown generation failed (non-fatal): {}", e),
            }
        }
        ReportMode::BothR => {
            log::info!("Generating RMarkdown file and HTML report...");
            match report::execute(&report::ReportArgs {
                summary: &prefixed_join(outdir, pfx, "results.tsv"),
                detail: &prefixed_join(outdir, pfx, "detected_detail.tsv"),
                params: &prefixed_join(outdir, pfx, "run_params.tsv"),
                coverage: Some(&prefixed_join(outdir, pfx, "coverage.tsv")),
                species_calls: species_calls_opt,
                group_detail: group_detail_opt.as_deref(),
                run_name: &args.run_name,
                output: &prefixed_join(outdir, pfx, "report.html"),
                report: ReportMode::BothR,
            }) {
                Ok(()) => {}
                Err(e) => log::warn!("Report generation failed (non-fatal): {}", e),
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
            "reads.fa",
            "reads_R2.fa",
            "filtered.fa",
            "filtered_R2.fa",
            "mapped.sam",
            "detected.list",
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
