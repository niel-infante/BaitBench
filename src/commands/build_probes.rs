use anyhow::{Context, Result};
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

use crate::cli::{ProbeMethod, ReportMode};
use crate::commands::assess_probes;
use crate::external::{cdhit, rscript};
use crate::external::catch as external_catch;
use crate::io_utils::prefixed_join;
use crate::sdust;
use crate::syotti;
use crate::catch;

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
}

/// Statistics for a single pipeline step.
struct StepStats {
    step: String,
    sequences: usize,
    bases: usize,
}

pub fn execute(args: &BuildProbesArgs) -> Result<()> {
    log::info!("=============================================");
    log::info!("BaitBench - Build Probes");
    log::info!("=============================================");

    // Validate inputs
    if !args.targets.exists() {
        anyhow::bail!("Targets file not found: {}", args.targets.display());
    }

    // Check external tools
    let cdhit_available = cdhit::is_available();
    if !cdhit_available {
        log::warn!("cd-hit-est not found — target collapse (step 2) and probe deduplication (step 7) will be skipped. The probeset may contain redundant sequences.");
    }

    // Create output directory
    fs::create_dir_all(args.outdir)?;

    let mut stats: Vec<StepStats> = Vec::new();

    // --- Step 0: Count input ---
    log::info!("Counting input sequences...");
    let (input_seqs, input_bases) = count_fasta_stats(args.targets)?;
    log::info!("  Input: {} sequences, {} bases", input_seqs, input_bases);
    stats.push(StepStats {
        step: "input".to_string(),
        sequences: input_seqs,
        bases: input_bases,
    });

    // --- Step 1: Filter targets by N content ---
    log::info!(
        "Step 1: Filtering targets by N content (max {:.0}%)...",
        args.max_n_frac * 100.0
    );
    let targets_clean_path =
        prefixed_join(args.outdir, args.output_prefix, "targets_clean.fa");
    filter_n_content(args.targets, &targets_clean_path, args.max_n_frac)?;
    let (clean_seqs, clean_bases) = count_fasta_stats(&targets_clean_path)?;
    log::info!(
        "  After N filter: {} sequences, {} bases (removed {})",
        clean_seqs,
        clean_bases,
        input_seqs - clean_seqs
    );
    stats.push(StepStats {
        step: "n_filtered".to_string(),
        sequences: clean_seqs,
        bases: clean_bases,
    });

    // --- Step 2: Collapse targets with cd-hit-est ---
    let collapsed_path = prefixed_join(args.outdir, args.output_prefix, "collapsed.fa");
    let (collapsed_seqs, collapsed_bases) = if cdhit_available {
        log::info!(
            "Step 2: Collapsing targets (cd-hit-est, threshold={:.2})...",
            args.collapse_threshold
        );
        let collapse_log = prefixed_join(args.outdir, args.output_prefix, "cdhit_collapse.log");
        cdhit::cluster(
            &targets_clean_path,
            &collapsed_path,
            args.collapse_threshold,
            args.threads,
            &collapse_log,
        )?;
        let stats = count_fasta_stats(&collapsed_path)?;
        log::info!(
            "  Collapsed: {} sequences, {} bases (removed {})",
            stats.0, stats.1, clean_seqs - stats.0
        );
        stats
    } else {
        log::warn!("Step 2: Skipping target collapse (cd-hit-est not available).");
        fs::copy(&targets_clean_path, &collapsed_path)?;
        (clean_seqs, clean_bases)
    };
    stats.push(StepStats {
        step: "collapsed".to_string(),
        sequences: collapsed_seqs,
        bases: collapsed_bases,
    });

    // --- Step 3: Filter short sequences (< probe length) ---
    log::info!(
        "Step 3: Filtering sequences shorter than probe length ({} bp)...",
        args.probe_length
    );
    let length_filtered_path =
        prefixed_join(args.outdir, args.output_prefix, "length_filtered.fa");
    filter_short_sequences(&collapsed_path, &length_filtered_path, args.probe_length)?;
    let (length_filtered_seqs, length_filtered_bases) = count_fasta_stats(&length_filtered_path)?;
    log::info!(
        "  After length filter: {} sequences, {} bases (removed {})",
        length_filtered_seqs,
        length_filtered_bases,
        collapsed_seqs - length_filtered_seqs
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
                args.probe_length,
                args.catch_probe_stride,
                args.catch_mismatches,
                args.catch_extension,
                args.catch_coverage,
                args.catch_minhash_threshold,
            );
            catch::design_probes(
                &length_filtered_path,
                &probes_raw_path,
                args.probe_length,
                args.catch_probe_stride,
                args.catch_mismatches,
                args.catch_extension,
                args.catch_coverage,
                args.catch_minhash_threshold,
            )?;
        }
        ProbeMethod::SyottiLite => {
            log::info!(
                "Step 4: Building probes (syotti-lite, length={}, mismatches={}, seed_len={})...",
                args.probe_length, args.syotti_mismatches, args.syotti_seed_len
            );
            syotti::design_probes(
                &length_filtered_path,
                &probes_raw_path,
                args.probe_length,
                args.syotti_mismatches,
                args.syotti_seed_len,
            )?;
        }
        ProbeMethod::Catch => {
            log::info!(
                "Step 4: Building probes (external CATCH, length={}, stride={}, mismatches={}, \
                 extension={}, coverage={:.2}, minhash={:.2})...",
                args.probe_length,
                args.catch_probe_stride,
                args.catch_mismatches,
                args.catch_extension,
                args.catch_coverage,
                args.catch_minhash_threshold,
            );
            external_catch::check_available()?;
            external_catch::design_probes(
                &length_filtered_path,
                &probes_raw_path,
                args.probe_length,
                args.catch_probe_stride,
                args.catch_mismatches,
                args.catch_extension,
                args.catch_coverage,
                args.catch_minhash_threshold,
            )?;
        }
    }
    let (tiled_seqs, tiled_bases) = count_fasta_stats(&probes_raw_path)?;
    log::info!("  Built: {} probes, {} bases", tiled_seqs, tiled_bases);
    stats.push(StepStats {
        step: "built".to_string(),
        sequences: tiled_seqs,
        bases: tiled_bases,
    });

    // --- Step 5: Filter probes by GC content ---
    log::info!(
        "Step 5: Filtering probes by GC content ({:.0}%-{:.0}%)...",
        args.min_gc * 100.0,
        args.max_gc * 100.0
    );
    let probes_gc_path =
        prefixed_join(args.outdir, args.output_prefix, "probes_gc.fa");
    filter_gc(
        &probes_raw_path,
        &probes_gc_path,
        args.min_gc,
        args.max_gc,
    )?;
    let (gc_seqs, gc_bases) = count_fasta_stats(&probes_gc_path)?;
    log::info!(
        "  After GC filter: {} probes, {} bases (removed {})",
        gc_seqs,
        gc_bases,
        tiled_seqs - gc_seqs
    );
    stats.push(StepStats {
        step: "gc_filtered".to_string(),
        sequences: gc_seqs,
        bases: gc_bases,
    });

    // --- Step 6: Filter probes by complexity (sDUST) ---
    let (complexity_input, complexity_seqs, complexity_bases);
    if args.max_masked_frac >= 1.0 {
        log::info!("Step 6: Skipping complexity filter (--max-masked-frac >= 1.0)");
        complexity_input = probes_gc_path.clone();
        complexity_seqs = gc_seqs;
        complexity_bases = gc_bases;
    } else {
        log::info!(
            "Step 6: Filtering probes by complexity (sDUST, T={}, W={}, max masked={:.0}%)...",
            args.dust_threshold,
            args.dust_window,
            args.max_masked_frac * 100.0
        );
        let probes_complexity_path =
            prefixed_join(args.outdir, args.output_prefix, "probes_complexity.fa");
        filter_complexity(
            &probes_gc_path,
            &probes_complexity_path,
            args.dust_threshold,
            args.dust_window,
            args.max_masked_frac,
        )?;
        let (cs, cb) = count_fasta_stats(&probes_complexity_path)?;
        log::info!(
            "  After complexity filter: {} probes, {} bases (removed {})",
            cs,
            cb,
            gc_seqs - cs
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

    // --- Step 7: Deduplicate probes with cd-hit-est ---
    let probes_final_path = prefixed_join(args.outdir, args.output_prefix, "probes_final.fa");
    let (final_seqs, final_bases) = if cdhit_available {
        log::info!(
            "Step 7: Deduplicating probes (cd-hit-est, threshold={:.2})...",
            args.dedup_threshold
        );
        let dedup_log = prefixed_join(args.outdir, args.output_prefix, "cdhit_dedup.log");
        cdhit::cluster(
            &complexity_input,
            &probes_final_path,
            args.dedup_threshold,
            args.threads,
            &dedup_log,
        )?;
        let s = count_fasta_stats(&probes_final_path)?;
        log::info!(
            "  Final: {} probes, {} bases (removed {})",
            s.0, s.1, complexity_seqs - s.0
        );
        s
    } else {
        log::warn!("Step 7: Skipping probe deduplication (cd-hit-est not available).");
        fs::copy(&complexity_input, &probes_final_path)?;
        (complexity_seqs, complexity_bases)
    };
    stats.push(StepStats {
        step: "deduplicated".to_string(),
        sequences: final_seqs,
        bases: final_bases,
    });

    // --- Write stats TSV ---
    let stats_path = prefixed_join(args.outdir, args.output_prefix, "build_probes_stats.tsv");
    write_stats_tsv(&stats_path, &stats)?;
    log::info!("Pipeline stats written to {}", stats_path.display());

    // --- Write run params ---
    let params_path = prefixed_join(args.outdir, args.output_prefix, "run_params.tsv");
    write_run_params(&params_path, args, cdhit_available)?;

    // --- Assessment and/or Report ---
    if !args.skip_assess {
        // Chain into assess_probes: runs coverage + xreact + combined report
        log::info!("Chaining into probe assessment...");
        assess_probes::execute(&assess_probes::AssessProbesArgs {
            targets: args.targets,
            probes: &probes_final_path,
            genomes: args.genomes,
            threshold: args.threshold,
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
            all_individual_targets: false,
        })?;
    } else {
        // Original build-probes-only report when assessment is skipped
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
    }

    log::info!("=============================================");
    log::info!("Build probes complete!");
    log::info!("  Final probeset: {}", probes_final_path.display());
    log::info!("  {} probes, {} bases", final_seqs, final_bases);
    log::info!("  Results in {}", args.outdir.display());
    log::info!("=============================================");

    Ok(())
}

/// Count sequences and total bases in a FASTA file (streaming).
fn count_fasta_stats(path: &Path) -> Result<(usize, usize)> {
    let file =
        File::open(path).with_context(|| format!("Cannot open FASTA: {}", path.display()))?;
    let reader = BufReader::new(file);

    let mut num_seqs = 0;
    let mut total_bases = 0;

    for line in reader.lines() {
        let line = line?;
        let trimmed = line.trim_end();
        if trimmed.starts_with('>') {
            num_seqs += 1;
        } else if !trimmed.is_empty() {
            total_bases += trimmed.len();
        }
    }

    Ok((num_seqs, total_bases))
}

/// Tile probes across each sequence using a sliding window.
///
/// For each sequence, probes of `probe_length` bp are generated starting from
/// position 0 and advancing by `stride = probe_length + step` each time.
/// - `step < 0`: probes overlap (e.g., -60 with length 120 → stride 60, 50% overlap)
/// - `step = 0`: probes are perfectly tiled (no overlap, no gap)
/// - `step > 0`: gaps between probes
///
/// A final probe is always placed at the end of the sequence (last `probe_length` bp),
/// regardless of overlap with the previous probe. Sequences shorter than `probe_length`
/// produce a single probe of whatever length is available.
fn tile_probes(input: &Path, output: &Path, probe_length: usize, step: i64) -> Result<()> {
    let stride = (probe_length as i64 + step) as usize;

    let file =
        File::open(input).with_context(|| format!("Cannot open FASTA: {}", input.display()))?;
    let reader = BufReader::new(file);

    let out_file =
        File::create(output).with_context(|| format!("Cannot create: {}", output.display()))?;
    let mut writer = BufWriter::new(out_file);

    let mut current_id: Option<String> = None;
    let mut current_seq = String::new();

    let emit_tiles =
        |id: &str, seq: &str, writer: &mut BufWriter<File>| -> Result<()> {
            let seq_len = seq.len();
            if seq_len == 0 {
                return Ok(());
            }

            if seq_len <= probe_length {
                // Sequence shorter than or equal to probe length: emit as single probe
                writeln!(writer, ">probe_{}|tile_1", id)?;
                writeln!(writer, "{}", seq)?;
                return Ok(());
            }

            let mut tile_num = 0usize;
            let mut start = 0usize;

            // Generate probes at regular stride intervals
            while start + probe_length <= seq_len {
                tile_num += 1;
                writeln!(writer, ">probe_{}|tile_{}", id, tile_num)?;
                writeln!(writer, "{}", &seq[start..start + probe_length])?;
                start += stride;
            }

            // Final probe: anchored to the end of the sequence.
            // Only emit if it doesn't duplicate the last probe exactly.
            let final_start = seq_len - probe_length;
            let last_emitted_start = start - stride; // start of the last probe we wrote
            if final_start != last_emitted_start {
                tile_num += 1;
                writeln!(writer, ">probe_{}|tile_{}", id, tile_num)?;
                writeln!(writer, "{}", &seq[final_start..seq_len])?;
            }

            Ok(())
        };

    for line in reader.lines() {
        let line = line?;
        let trimmed = line.trim_end();
        if trimmed.starts_with('>') {
            if let Some(ref id) = current_id {
                emit_tiles(id, &current_seq, &mut writer)?;
            }
            current_id = Some(trimmed[1..].split_whitespace().next().unwrap_or("").to_string());
            current_seq.clear();
        } else if !trimmed.is_empty() {
            current_seq.push_str(&trimmed.to_uppercase());
        }
    }

    if let Some(ref id) = current_id {
        emit_tiles(id, &current_seq, &mut writer)?;
    }

    writer.flush()?;
    Ok(())
}

/// Compute fraction of N (ambiguous) bases in a sequence.
fn compute_n_frac(seq: &str) -> f64 {
    if seq.is_empty() {
        return 0.0;
    }
    let n_count = seq
        .chars()
        .filter(|c| !matches!(c, 'A' | 'C' | 'G' | 'T' | 'a' | 'c' | 'g' | 't'))
        .count();
    n_count as f64 / seq.len() as f64
}

/// Filter FASTA sequences by N (ambiguous base) content (streaming).
///
/// Keeps sequences whose N fraction is at most `max_n_frac`.
fn filter_n_content(input: &Path, output: &Path, max_n_frac: f64) -> Result<()> {
    let file =
        File::open(input).with_context(|| format!("Cannot open FASTA: {}", input.display()))?;
    let reader = BufReader::new(file);

    let out_file =
        File::create(output).with_context(|| format!("Cannot create: {}", output.display()))?;
    let mut writer = BufWriter::new(out_file);

    let mut current_id: Option<String> = None;
    let mut current_seq = String::new();

    for line in reader.lines() {
        let line = line?;
        let trimmed = line.trim_end();
        if trimmed.starts_with('>') {
            if let Some(ref id) = current_id {
                let frac = compute_n_frac(&current_seq);
                if frac <= max_n_frac {
                    writeln!(writer, ">{}", id)?;
                    writeln!(writer, "{}", current_seq)?;
                }
            }
            current_id = Some(trimmed[1..].to_string());
            current_seq.clear();
        } else if !trimmed.is_empty() {
            current_seq.push_str(&trimmed.to_uppercase());
        }
    }

    if let Some(ref id) = current_id {
        let frac = compute_n_frac(&current_seq);
        if frac <= max_n_frac {
            writeln!(writer, ">{}", id)?;
            writeln!(writer, "{}", current_seq)?;
        }
    }

    writer.flush()?;
    Ok(())
}

/// Filter FASTA sequences shorter than `min_length` (streaming).
fn filter_short_sequences(input: &Path, output: &Path, min_length: usize) -> Result<()> {
    let file =
        File::open(input).with_context(|| format!("Cannot open FASTA: {}", input.display()))?;
    let reader = BufReader::new(file);

    let out_file =
        File::create(output).with_context(|| format!("Cannot create: {}", output.display()))?;
    let mut writer = BufWriter::new(out_file);

    let mut current_id: Option<String> = None;
    let mut current_seq = String::new();

    for line in reader.lines() {
        let line = line?;
        let trimmed = line.trim_end();
        if trimmed.starts_with('>') {
            if let Some(ref id) = current_id {
                if current_seq.len() >= min_length {
                    writeln!(writer, ">{}", id)?;
                    writeln!(writer, "{}", current_seq)?;
                }
            }
            current_id = Some(trimmed[1..].to_string());
            current_seq.clear();
        } else if !trimmed.is_empty() {
            current_seq.push_str(&trimmed.to_uppercase());
        }
    }

    if let Some(ref id) = current_id {
        if current_seq.len() >= min_length {
            writeln!(writer, ">{}", id)?;
            writeln!(writer, "{}", current_seq)?;
        }
    }

    writer.flush()?;
    Ok(())
}

/// Compute GC fraction of a sequence.
fn compute_gc(seq: &str) -> f64 {
    if seq.is_empty() {
        return 0.0;
    }
    let gc_count = seq
        .chars()
        .filter(|c| matches!(c, 'G' | 'C' | 'g' | 'c'))
        .count();
    gc_count as f64 / seq.len() as f64
}

/// Filter FASTA sequences by GC content (streaming).
fn filter_gc(input: &Path, output: &Path, min_gc: f64, max_gc: f64) -> Result<()> {
    let file =
        File::open(input).with_context(|| format!("Cannot open FASTA: {}", input.display()))?;
    let reader = BufReader::new(file);

    let out_file =
        File::create(output).with_context(|| format!("Cannot create: {}", output.display()))?;
    let mut writer = BufWriter::new(out_file);

    let mut current_id: Option<String> = None;
    let mut current_seq = String::new();

    for line in reader.lines() {
        let line = line?;
        let trimmed = line.trim_end();
        if trimmed.starts_with('>') {
            // Check and write previous sequence
            if let Some(ref id) = current_id {
                let gc = compute_gc(&current_seq);
                if gc >= min_gc && gc <= max_gc {
                    writeln!(writer, ">{}", id)?;
                    writeln!(writer, "{}", current_seq)?;
                }
            }
            current_id = Some(trimmed[1..].to_string());
            current_seq.clear();
        } else if !trimmed.is_empty() {
            current_seq.push_str(&trimmed.to_uppercase());
        }
    }

    // Check and write last sequence
    if let Some(ref id) = current_id {
        let gc = compute_gc(&current_seq);
        if gc >= min_gc && gc <= max_gc {
            writeln!(writer, ">{}", id)?;
            writeln!(writer, "{}", current_seq)?;
        }
    }

    writer.flush()?;
    Ok(())
}

/// Filter FASTA sequences by low-complexity content using sDUST (streaming).
///
/// Keeps sequences whose masked fraction is at most `max_masked_frac`.
fn filter_complexity(
    input: &Path,
    output: &Path,
    dust_threshold: f64,
    dust_window: usize,
    max_masked_frac: f64,
) -> Result<()> {
    let file =
        File::open(input).with_context(|| format!("Cannot open FASTA: {}", input.display()))?;
    let reader = BufReader::new(file);

    let out_file =
        File::create(output).with_context(|| format!("Cannot create: {}", output.display()))?;
    let mut writer = BufWriter::new(out_file);

    let mut current_id: Option<String> = None;
    let mut current_seq = String::new();

    for line in reader.lines() {
        let line = line?;
        let trimmed = line.trim_end();
        if trimmed.starts_with('>') {
            if let Some(ref id) = current_id {
                let frac = sdust::masked_fraction(current_seq.as_bytes(), dust_threshold, dust_window);
                if frac <= max_masked_frac {
                    writeln!(writer, ">{}", id)?;
                    writeln!(writer, "{}", current_seq)?;
                }
            }
            current_id = Some(trimmed[1..].to_string());
            current_seq.clear();
        } else if !trimmed.is_empty() {
            current_seq.push_str(&trimmed.to_uppercase());
        }
    }

    if let Some(ref id) = current_id {
        let frac = sdust::masked_fraction(current_seq.as_bytes(), dust_threshold, dust_window);
        if frac <= max_masked_frac {
            writeln!(writer, ">{}", id)?;
            writeln!(writer, "{}", current_seq)?;
        }
    }

    writer.flush()?;
    Ok(())
}

/// Write pipeline statistics TSV.
fn write_stats_tsv(path: &Path, stats: &[StepStats]) -> Result<()> {
    let file = File::create(path).with_context(|| format!("Cannot create: {}", path.display()))?;
    let mut w = BufWriter::new(file);

    writeln!(w, "step\tsequences\tbases")?;
    for s in stats {
        writeln!(w, "{}\t{}\t{}", s.step, s.sequences, s.bases)?;
    }

    w.flush()?;
    Ok(())
}

/// Write run parameters TSV for report reconstruction.
fn write_run_params(path: &Path, args: &BuildProbesArgs, cdhit_available: bool) -> Result<()> {
    let file = File::create(path).with_context(|| format!("Cannot create: {}", path.display()))?;
    let mut w = BufWriter::new(file);

    writeln!(w, "parameter\tflag\tvalue")?;
    writeln!(w, "targets\t--targets\t{}", args.targets.display())?;
    writeln!(w, "method\t--method\t{:?}", args.method)?;
    writeln!(w, "probe_length\t--probe-length\t{}", args.probe_length)?;
    writeln!(w, "step\t--step\t{}", args.step)?;
    writeln!(w, "catch_probe_stride\t--catch-stride\t{}", args.catch_probe_stride)?;
    writeln!(w, "catch_mismatches\t--catch-mismatches\t{}", args.catch_mismatches)?;
    writeln!(w, "catch_extension\t--catch-extension\t{}", args.catch_extension)?;
    writeln!(w, "catch_coverage\t--catch-coverage\t{:.2}", args.catch_coverage)?;
    writeln!(w, "catch_minhash_threshold\t--catch-minhash-threshold\t{:.2}", args.catch_minhash_threshold)?;
    writeln!(w, "min_gc\t--min-gc\t{:.2}", args.min_gc)?;
    writeln!(w, "max_gc\t--max-gc\t{:.2}", args.max_gc)?;
    writeln!(w, "max_n_frac\t--max-n-frac\t{:.2}", args.max_n_frac)?;
    writeln!(w, "dust_threshold\t--dust-threshold\t{}", args.dust_threshold)?;
    writeln!(w, "dust_window\t--dust-window\t{}", args.dust_window)?;
    writeln!(
        w,
        "max_masked_frac\t--max-masked-frac\t{:.2}",
        args.max_masked_frac
    )?;
    writeln!(
        w,
        "collapse_threshold\t--collapse-threshold\t{:.2}",
        args.collapse_threshold
    )?;
    writeln!(
        w,
        "dedup_threshold\t--dedup-threshold\t{:.2}",
        args.dedup_threshold
    )?;
    writeln!(w, "threads\t--threads\t{}", args.threads)?;
    writeln!(w, "outdir\t--outdir\t{}", args.outdir.display())?;
    writeln!(w, "cdhit_available\t--n/a\t{}", cdhit_available)?;

    w.flush()?;
    Ok(())
}

/// Generate the HTML report via Rscript.
fn generate_report(stats_path: &Path, params_path: &Path, output_path: &Path) -> Result<()> {
    let r_dir =
        rscript::find_r_dir().ok_or_else(|| anyhow::anyhow!("Cannot find R scripts directory."))?;

    let script = r_dir.join("build_probes.R");
    if !script.exists() {
        anyhow::bail!("Build probes R script not found: {}", script.display());
    }

    let stats_abs = std::fs::canonicalize(stats_path)?;
    let params_abs = std::fs::canonicalize(params_path)?;
    let output_abs = if output_path.is_absolute() {
        output_path.to_path_buf()
    } else {
        std::env::current_dir()?.join(output_path)
    };

    rscript::run_rscript(
        &script,
        &[
            "--stats",
            stats_abs.to_str().unwrap_or(""),
            "--params",
            params_abs.to_str().unwrap_or(""),
            "--output",
            output_abs.to_str().unwrap_or(""),
        ],
    )
}

/// Write a parameterized RMarkdown file for manual rendering.
fn write_rmd(stats_path: &Path, params_path: &Path, rmd_path: &Path) -> Result<()> {
    let r_dir =
        rscript::find_r_dir().ok_or_else(|| anyhow::anyhow!("Cannot find R scripts directory."))?;

    let template_path = r_dir.join("build_probes.Rmd");
    if !template_path.exists() {
        anyhow::bail!(
            "Build probes Rmd template not found: {}",
            template_path.display()
        );
    }

    let template = std::fs::read_to_string(&template_path)?;

    let stats_abs = std::fs::canonicalize(stats_path)?;
    let params_abs = std::fs::canonicalize(params_path)?;

    let substituted = crate::commands::report::substitute_rmd_params(
        &template,
        &[
            ("stats_file", stats_abs.to_str().unwrap_or("")),
            ("params_file", params_abs.to_str().unwrap_or("")),
        ],
    );

    std::fs::write(rmd_path, substituted)?;
    Ok(())
}
