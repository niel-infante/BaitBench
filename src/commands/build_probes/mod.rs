mod filters;
mod n_fix;
mod report;
mod tiling;

use anyhow::Result;
use std::fs;
use std::path::{Path, PathBuf};

use crate::cli::{Aligner, ProbeMethod, ReportMode};
use crate::commands::assess_probes;
use crate::external::{cdhit, rscript};
use crate::external::catch as external_catch;
use crate::io_utils::prefixed_join;
use crate::probes::{syotti, catch, probetools};

use filters::{count_fasta_stats, filter_complexity, filter_gc, filter_n_content, filter_short_sequences};
use n_fix::fix_n_bases;
use tiling::tile_probes;
use report::{generate_report, write_rmd, write_run_params, write_stats_tsv};

pub struct BuildProbesArgs<'a> {
    pub targets: &'a Path,
    pub method: ProbeMethod,
    pub probe_length: usize,
    pub step: i64,
    pub catch_probe_stride: usize,
    pub catch_mismatches: usize,
    pub catch_extension: usize,
    pub catch_coverage: f64,
    pub catch_minhash_threshold: f64,
    pub min_gc: f64,
    pub max_gc: f64,
    pub dust_threshold: f64,
    pub dust_window: usize,
    pub max_masked_frac: f64,
    pub max_n_frac: f64,
    pub collapse_threshold: f64,
    pub dedup_threshold: f64,
    pub threads: usize,
    pub genomes: &'a [PathBuf],
    pub threshold: f64,
    pub aligner: Aligner,
    pub minimap_preset: &'a str,
    pub proximity: usize,
    pub skip_assess: bool,
    pub outdir: &'a Path,
    pub output_prefix: &'a str,
    pub report: ReportMode,
    pub cleanup: bool,
    pub refine_threshold: f64,
    pub refine_iterations: Option<usize>,
    pub refine_until_stable: bool,
    pub syotti_mismatches: usize,
    pub syotti_seed_len: usize,
    pub pt_step: usize,
    pub pt_identity: f64,
    pub pt_coverage: f64,
    pub pt_batch_size: usize,
    pub pt_max_panel_size: Option<usize>,
    pub pt_min_depth: u32,
    pub pt_max_iterations: usize,
    pub pt_min_coverage_gain: f64,
    pub no_n_in_probes: bool,
}

pub(super) struct StepStats {
    pub step: String,
    pub sequences: usize,
    pub bases: usize,
}

pub fn execute(args: &BuildProbesArgs) -> Result<()> {
    log::info!("=============================================");
    log::info!("BaitBench - Build Probes");
    log::info!("=============================================");

    if !args.targets.exists() {
        anyhow::bail!("Targets file not found: {}", args.targets.display());
    }

    let cdhit_available = cdhit::is_available();
    if !cdhit_available {
        log::warn!("cd-hit-est not found — target collapse (step 2) and probe deduplication (step 7) will be skipped. The probeset may contain redundant sequences.");
    }

    fs::create_dir_all(args.outdir)?;

    let mut stats: Vec<StepStats> = Vec::new();

    // --- Step 0: Count input ---
    log::info!("Counting input sequences...");
    let (input_seqs, input_bases) = count_fasta_stats(args.targets)?;
    log::info!("  Input: {} sequences, {} bases", input_seqs, input_bases);
    stats.push(StepStats { step: "input".to_string(), sequences: input_seqs, bases: input_bases });

    // --- Step 1: Filter by N content ---
    log::info!(
        "Step 1: Filtering targets by N content (max {:.0}%)...",
        args.max_n_frac * 100.0
    );
    let targets_clean_path = prefixed_join(args.outdir, args.output_prefix, "targets_clean.fa");
    filter_n_content(args.targets, &targets_clean_path, args.max_n_frac)?;
    let (clean_seqs, clean_bases) = count_fasta_stats(&targets_clean_path)?;
    log::info!(
        "  After N filter: {} sequences, {} bases (removed {})",
        clean_seqs, clean_bases, input_seqs - clean_seqs
    );
    stats.push(StepStats { step: "n_filtered".to_string(), sequences: clean_seqs, bases: clean_bases });

    // --- Step 2: Collapse with cd-hit-est ---
    let collapsed_path = prefixed_join(args.outdir, args.output_prefix, "collapsed.fa");
    let (collapsed_seqs, collapsed_bases) = if cdhit_available {
        log::info!(
            "Step 2: Collapsing targets (cd-hit-est, threshold={:.2})...",
            args.collapse_threshold
        );
        let collapse_log = prefixed_join(args.outdir, args.output_prefix, "cdhit_collapse.log");
        cdhit::cluster(&targets_clean_path, &collapsed_path, args.collapse_threshold, args.threads, &collapse_log)?;
        let s = count_fasta_stats(&collapsed_path)?;
        log::info!("  Collapsed: {} sequences, {} bases (removed {})", s.0, s.1, clean_seqs - s.0);
        s
    } else {
        log::warn!("Step 2: Skipping target collapse (cd-hit-est not available).");
        fs::copy(&targets_clean_path, &collapsed_path)?;
        (clean_seqs, clean_bases)
    };
    stats.push(StepStats { step: "collapsed".to_string(), sequences: collapsed_seqs, bases: collapsed_bases });

    // --- Step 3: Filter short sequences ---
    log::info!(
        "Step 3: Filtering sequences shorter than probe length ({} bp)...",
        args.probe_length
    );
    let length_filtered_path = prefixed_join(args.outdir, args.output_prefix, "length_filtered.fa");
    filter_short_sequences(&collapsed_path, &length_filtered_path, args.probe_length)?;
    let (length_filtered_seqs, length_filtered_bases) = count_fasta_stats(&length_filtered_path)?;
    log::info!(
        "  After length filter: {} sequences, {} bases (removed {})",
        length_filtered_seqs, length_filtered_bases, collapsed_seqs - length_filtered_seqs
    );
    stats.push(StepStats {
        step: "length_filtered".to_string(),
        sequences: length_filtered_seqs,
        bases: length_filtered_bases,
    });

    // --- Step 4: Build probes ---
    let probes_raw_path = prefixed_join(args.outdir, args.output_prefix, "probes_raw.fa");
    match args.method {
        ProbeMethod::Tile => {
            let stride = args.probe_length as i64 + args.step;
            if stride <= 0 {
                anyhow::bail!(
                    "Invalid tiling parameters: probe_length ({}) + step ({}) = {} (must be > 0)",
                    args.probe_length, args.step, stride
                );
            }
            log::info!(
                "Step 4: Building probes (tile, length={}, step={}, stride={})...",
                args.probe_length, args.step, stride
            );
            tile_probes(&length_filtered_path, &probes_raw_path, args.probe_length, args.step)?;
        }
        ProbeMethod::CatchLite => {
            log::info!(
                "Step 4: Building probes (catch-lite, length={}, stride={}, mismatches={}, \
                 extension={}, coverage={:.2}, minhash={:.2})...",
                args.probe_length, args.catch_probe_stride, args.catch_mismatches,
                args.catch_extension, args.catch_coverage, args.catch_minhash_threshold,
            );
            catch::design_probes(
                &length_filtered_path, &probes_raw_path, args.probe_length,
                args.catch_probe_stride, args.catch_mismatches, args.catch_extension,
                args.catch_coverage, args.catch_minhash_threshold,
            )?;
        }
        ProbeMethod::SyottiLite => {
            log::info!(
                "Step 4: Building probes (syotti-lite, length={}, mismatches={}, seed_len={})...",
                args.probe_length, args.syotti_mismatches, args.syotti_seed_len
            );
            syotti::design_probes(
                &length_filtered_path, &probes_raw_path, args.probe_length,
                args.syotti_mismatches, args.syotti_seed_len,
            )?;
        }
        ProbeMethod::Catch => {
            log::info!(
                "Step 4: Building probes (external CATCH, length={}, stride={}, mismatches={}, \
                 extension={}, coverage={:.2}, minhash={:.2})...",
                args.probe_length, args.catch_probe_stride, args.catch_mismatches,
                args.catch_extension, args.catch_coverage, args.catch_minhash_threshold,
            );
            external_catch::check_available()?;
            external_catch::design_probes(
                &length_filtered_path, &probes_raw_path, args.probe_length,
                args.catch_probe_stride, args.catch_mismatches, args.catch_extension,
                args.catch_coverage, args.catch_minhash_threshold,
            )?;
        }
        ProbeMethod::ProbeToolsLite => {
            log::info!(
                "Step 4: Building probes (probetools-lite, length={}, step={}, identity={:.2}, \
                 coverage_goal={:.2}, batch={}, max_iterations={})...",
                args.probe_length, args.pt_step, args.pt_identity,
                args.pt_coverage, args.pt_batch_size, args.pt_max_iterations,
            );
            cdhit::check_available()?;
            probetools::design_probes(
                &length_filtered_path, &probes_raw_path, args.probe_length,
                args.pt_step, args.pt_identity, args.pt_coverage, args.pt_batch_size,
                args.pt_max_panel_size, args.pt_min_depth, args.pt_max_iterations,
                args.pt_min_coverage_gain, args.minimap_preset, args.threads, args.outdir,
            )?;
        }
    }
    let (tiled_seqs, tiled_bases) = count_fasta_stats(&probes_raw_path)?;
    log::info!("  Built: {} probes, {} bases", tiled_seqs, tiled_bases);
    stats.push(StepStats { step: "built".to_string(), sequences: tiled_seqs, bases: tiled_bases });

    // --- Step 4b: Fix N bases ---
    let gc_filter_input = if args.no_n_in_probes {
        let probes_n_fixed_path = prefixed_join(args.outdir, args.output_prefix, "probes_n_fixed.fa");
        log::info!("Step 4b: Fixing N bases in probes (--no-n-in-probes)...");
        let (probes_with_n, n_replaced) = fix_n_bases(&probes_raw_path, &probes_n_fixed_path)?;
        log::info!("  Replaced {} N bases across {} probes", n_replaced, probes_with_n);
        probes_n_fixed_path
    } else {
        probes_raw_path.clone()
    };

    // --- Step 5: Filter by GC content ---
    log::info!(
        "Step 5: Filtering probes by GC content ({:.0}%-{:.0}%)...",
        args.min_gc * 100.0, args.max_gc * 100.0
    );
    let probes_gc_path = prefixed_join(args.outdir, args.output_prefix, "probes_gc.fa");
    filter_gc(&gc_filter_input, &probes_gc_path, args.min_gc, args.max_gc)?;
    let (gc_seqs, gc_bases) = count_fasta_stats(&probes_gc_path)?;
    log::info!(
        "  After GC filter: {} probes, {} bases (removed {})",
        gc_seqs, gc_bases, tiled_seqs - gc_seqs
    );
    stats.push(StepStats { step: "gc_filtered".to_string(), sequences: gc_seqs, bases: gc_bases });

    // --- Step 6: Filter by complexity (sDUST) ---
    let (complexity_input, complexity_seqs, complexity_bases);
    if args.max_masked_frac >= 1.0 {
        log::info!("Step 6: Skipping complexity filter (--max-masked-frac >= 1.0)");
        complexity_input = probes_gc_path.clone();
        complexity_seqs = gc_seqs;
        complexity_bases = gc_bases;
    } else {
        log::info!(
            "Step 6: Filtering probes by complexity (sDUST, T={}, W={}, max masked={:.0}%)...",
            args.dust_threshold, args.dust_window, args.max_masked_frac * 100.0
        );
        let probes_complexity_path =
            prefixed_join(args.outdir, args.output_prefix, "probes_complexity.fa");
        filter_complexity(
            &probes_gc_path, &probes_complexity_path,
            args.dust_threshold, args.dust_window, args.max_masked_frac,
        )?;
        let (cs, cb) = count_fasta_stats(&probes_complexity_path)?;
        log::info!(
            "  After complexity filter: {} probes, {} bases (removed {})",
            cs, cb, gc_seqs - cs
        );
        complexity_input = probes_complexity_path;
        complexity_seqs = cs;
        complexity_bases = cb;
    }
    stats.push(StepStats {
        step: "complexity_filtered".to_string(),
        sequences: complexity_seqs,
        bases: complexity_bases,
    });

    // --- Step 7: Deduplicate with cd-hit-est ---
    let probes_final_path = prefixed_join(args.outdir, args.output_prefix, "probes_final.fa");
    let (final_seqs, final_bases) = if cdhit_available {
        log::info!(
            "Step 7: Deduplicating probes (cd-hit-est, threshold={:.2})...",
            args.dedup_threshold
        );
        let dedup_log = prefixed_join(args.outdir, args.output_prefix, "cdhit_dedup.log");
        cdhit::cluster(&complexity_input, &probes_final_path, args.dedup_threshold, args.threads, &dedup_log)?;
        let s = count_fasta_stats(&probes_final_path)?;
        log::info!("  Final: {} probes, {} bases (removed {})", s.0, s.1, complexity_seqs - s.0);
        s
    } else {
        log::warn!("Step 7: Skipping probe deduplication (cd-hit-est not available).");
        fs::copy(&complexity_input, &probes_final_path)?;
        (complexity_seqs, complexity_bases)
    };
    stats.push(StepStats { step: "deduplicated".to_string(), sequences: final_seqs, bases: final_bases });

    // --- Write stats and params ---
    let stats_path = prefixed_join(args.outdir, args.output_prefix, "build_probes_stats.tsv");
    write_stats_tsv(&stats_path, &stats)?;
    log::info!("Pipeline stats written to {}", stats_path.display());

    let params_path = prefixed_join(args.outdir, args.output_prefix, "run_params.tsv");
    write_run_params(&params_path, args, cdhit_available)?;

    // --- Assessment and/or report ---
    if !args.skip_assess {
        log::info!("Chaining into probe assessment...");
        assess_probes::execute(&assess_probes::AssessProbesArgs {
            targets: args.targets,
            probes: &probes_final_path,
            genomes: args.genomes,
            threshold: args.threshold,
            aligner: args.aligner,
            minimap_preset: args.minimap_preset,
            proximity: args.proximity,
            outdir: args.outdir,
            output_prefix: args.output_prefix,
            report: args.report,
            cleanup: args.cleanup,
            build_stats_file: Some(&stats_path),
            build_params_file: Some(&params_path),
            refine_threshold: args.refine_threshold,
            refine_iterations: args.refine_iterations,
            refine_until_stable: args.refine_until_stable,
            no_individual_targets: false,
            gap_min_length: None,
            threads: args.threads,
        })?;
    } else {
        match args.report {
            ReportMode::None => {
                log::info!("Skipping report generation (--report none)");
            }
            ReportMode::Full => {
                if rscript::check_available() {
                    let report_path =
                        prefixed_join(args.outdir, args.output_prefix, "build_probes_report.html");
                    match generate_report(&stats_path, &params_path, &report_path) {
                        Ok(()) => log::info!("Report generated: {}", report_path.display()),
                        Err(e) => log::warn!("Report generation failed (non-fatal): {}", e),
                    }
                } else {
                    log::warn!("Rscript not found -- skipping HTML report.");
                }
            }
            ReportMode::Rmd => {
                let report_path =
                    prefixed_join(args.outdir, args.output_prefix, "build_probes_report.html");
                let rmd_path = crate::commands::report::rmd_output_path(&report_path);
                match write_rmd(&stats_path, &params_path, &rmd_path) {
                    Ok(()) => log::info!("RMarkdown written: {}", rmd_path.display()),
                    Err(e) => log::warn!("RMarkdown generation failed (non-fatal): {}", e),
                }
            }
            ReportMode::BothR => {
                let report_path =
                    prefixed_join(args.outdir, args.output_prefix, "build_probes_report.html");
                let rmd_path = crate::commands::report::rmd_output_path(&report_path);
                match write_rmd(&stats_path, &params_path, &rmd_path) {
                    Ok(()) => log::info!("RMarkdown written: {}", rmd_path.display()),
                    Err(e) => log::warn!("RMarkdown generation failed (non-fatal): {}", e),
                }
                if rscript::check_available() {
                    match generate_report(&stats_path, &params_path, &report_path) {
                        Ok(()) => log::info!("Report generated: {}", report_path.display()),
                        Err(e) => log::warn!("Report generation failed (non-fatal): {}", e),
                    }
                } else {
                    log::warn!("Rscript not found — skipping HTML report (Rmd still written).");
                }
            }
        }
    }

    // --- Cleanup ---
    if args.cleanup {
        log::info!("Cleaning up build intermediate files...");
        let intermediates = [
            "targets_clean.fa",
            "collapsed.fa",
            "collapsed.fa.clstr",
            "length_filtered.fa",
            "probes_raw.fa",
            "probes_n_fixed.fa",
            "probes_gc.fa",
            "probes_complexity.fa",
            "probes_final.fa.clstr",
            "cdhit_collapse.log",
            "cdhit_dedup.log",
        ];
        for name in &intermediates {
            let path = prefixed_join(args.outdir, args.output_prefix, name);
            if path.exists() {
                let _ = fs::remove_file(&path);
            }
        }
        let pt_work = args.outdir.join("probetools_work");
        if pt_work.exists() {
            let _ = fs::remove_dir_all(&pt_work);
        }
    }

    log::info!("=============================================");
    log::info!("Build probes complete!");
    log::info!("  Final probeset: {}", probes_final_path.display());
    log::info!("  {} probes, {} bases", final_seqs, final_bases);
    log::info!("  Results in {}", args.outdir.display());
    log::info!("=============================================");

    Ok(())
}
