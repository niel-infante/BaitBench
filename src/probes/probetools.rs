use anyhow::{Context, Result};
use std::collections::HashMap;
use std::fs::{self, File, OpenOptions};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

use crate::alignment::coverage;
use crate::external::{cdhit, minimap2};

/// Design probes using the ProbeTools-Lite algorithm:
/// iterative k-mer clustering + coverage-driven greedy selection
/// (reimplementation of Kuchinski et al. 2022, BMC Genomics).
///
/// Substitutions vs. original: cd-hit-est replaces VSEARCH for clustering;
/// minimap2 replaces BLAST for coverage assessment.
pub fn design_probes(
    input: &Path,
    output: &Path,
    probe_len: usize,
    step: usize,
    identity: f64,
    coverage_goal: f64,
    batch_size: usize,
    max_panel_size: Option<usize>,
    min_depth: u32,
    max_iterations: usize,
    min_coverage_gain: f64,
    minimap_preset: &str,
    threads: usize,
    workdir: &Path,
) -> Result<usize> {
    let targets = load_fasta(input)?;
    if targets.is_empty() {
        log::warn!("probetools-lite: no input sequences — writing empty output.");
        File::create(output)
            .with_context(|| format!("Cannot create output: {}", output.display()))?;
        return Ok(0);
    }

    let pt_work = workdir.join("probetools_work");
    fs::create_dir_all(&pt_work)?;

    // Create (truncate) output panel file
    File::create(output)
        .with_context(|| format!("Cannot create panel file: {}", output.display()))?;

    let mut panel_size = 0usize;
    let mut probe_counter = 0usize;
    let mut iteration = 0usize;
    let mut prev_coverage = -1.0f64;
    let mut current_seqs: Vec<(String, String)> = targets.clone();

    loop {
        // Condition 5: no sequences left to enumerate k-mers from
        if current_seqs.is_empty() {
            log::info!("  [probetools-lite] No under-covered regions remain. Done.");
            break;
        }

        log::info!(
            "  [probetools-lite] Iteration {}: {} input region(s)",
            iteration,
            current_seqs.len()
        );

        // Step A: enumerate k-mers from current sequences
        let kmers_path = pt_work.join(format!("iter_{}_kmers.fa", iteration));
        let n_kmers = enumerate_kmers(&current_seqs, probe_len, step, &kmers_path)?;

        // Condition 6: no k-mers (regions too short or all-N)
        if n_kmers == 0 {
            log::info!("  [probetools-lite] No k-mers generated from under-covered regions. Done.");
            break;
        }

        // Step B: cluster k-mers with cd-hit-est
        let centroids_path = pt_work.join(format!("iter_{}_centroids.fa", iteration));
        let cdhit_log = pt_work.join(format!("iter_{}_cdhit.log", iteration));
        cdhit::cluster(&kmers_path, &centroids_path, identity, threads, &cdhit_log)?;

        let clstr_path = clstr_path_for(&centroids_path);

        // Step C: parse cluster sizes and rank centroids
        let size_map = parse_clstr(&clstr_path)?;
        let ranked = rank_centroids(&centroids_path, &size_map)?;

        if ranked.is_empty() {
            log::info!("  [probetools-lite] No centroids after clustering. Done.");
            break;
        }

        // Step D: select probes for this batch
        let n_select = {
            let capacity = if let Some(max) = max_panel_size {
                let remaining = max.saturating_sub(panel_size);
                if remaining == 0 {
                    // Condition 2: panel already full
                    log::info!(
                        "  [probetools-lite] Panel capacity reached ({} probes). Done.",
                        max
                    );
                    break;
                }
                remaining
            } else {
                usize::MAX
            };
            batch_size.min(ranked.len()).min(capacity)
        };

        // Step E: append selected probes to output panel
        {
            let out_file = OpenOptions::new()
                .create(true)
                .append(true)
                .open(output)
                .with_context(|| format!("Cannot open panel for append: {}", output.display()))?;
            let mut writer = BufWriter::new(out_file);
            for (_name, seq, _size) in ranked.iter().take(n_select) {
                probe_counter += 1;
                writeln!(writer, ">probe_probetools_{}|probetools_{}", probe_counter, probe_counter)?;
                writeln!(writer, "{}", seq)?;
            }
            writer.flush()?;
        }
        panel_size += n_select;
        log::info!(
            "  [probetools-lite] Added {} probes (panel total: {})",
            n_select,
            panel_size
        );

        // Condition 2 (post-add): panel now at max
        if let Some(max) = max_panel_size {
            if panel_size >= max {
                log::info!(
                    "  [probetools-lite] Panel size limit reached ({} probes). Done.",
                    max
                );
                break;
            }
        }

        // Step F: align full panel to original targets
        let sam_path = pt_work.join("panel_align.sam");
        let sam_log = pt_work.join("panel_align.log");
        minimap2::probe_align(
            minimap_preset,
            input,
            output,
            &sam_path,
            &sam_log,
            threads,
            1000,
        )?;

        // Step G: compute per-position coverage
        let cov_result = coverage::compute_probe_coverage(&sam_path)?;

        // Step H: compute 10th-percentile coverage fraction across all targets
        let current_coverage =
            percentile_coverage_for_targets(&cov_result.coverage, &targets, min_depth, 10.0);
        log::info!(
            "  [probetools-lite] 10th-pct coverage: {:.2}% (goal: {:.2}%)",
            current_coverage * 100.0,
            coverage_goal * 100.0,
        );

        // Condition 1: coverage goal reached
        if current_coverage >= coverage_goal {
            log::info!("  [probetools-lite] Coverage goal reached. Done.");
            break;
        }

        // Condition 3: max iterations
        if iteration + 1 >= max_iterations {
            log::warn!(
                "  [probetools-lite] Max iterations ({}) reached. Coverage: {:.2}% / goal {:.2}%.",
                max_iterations,
                current_coverage * 100.0,
                coverage_goal * 100.0,
            );
            break;
        }

        // Condition 4: stagnation guard
        if prev_coverage >= 0.0 {
            let gain = current_coverage - prev_coverage;
            if gain < min_coverage_gain {
                log::warn!(
                    "  [probetools-lite] Coverage gain {:.4}% < threshold {:.4}% — stagnating. Stopping.",
                    gain * 100.0,
                    min_coverage_gain * 100.0,
                );
                break;
            }
        }
        prev_coverage = current_coverage;

        // Step I: extract under-covered regions for next iteration
        let low_cov_path = pt_work.join(format!("iter_{}_lowcov.fa", iteration));
        let n_low_cov = extract_low_coverage_regions(
            &cov_result.coverage,
            &targets,
            min_depth,
            probe_len,
            &low_cov_path,
        )?;

        // Condition 5 (next iteration)
        if n_low_cov == 0 {
            log::info!("  [probetools-lite] No under-covered regions remain. Done.");
            break;
        }

        current_seqs = load_fasta(&low_cov_path)?;
        iteration += 1;
    }

    log::info!(
        "  [probetools-lite] Complete: {} probes, {} iteration(s).",
        panel_size,
        iteration + 1,
    );
    Ok(panel_size)
}

// ---------------------------------------------------------------------------
// K-mer enumeration
// ---------------------------------------------------------------------------

/// Enumerate overlapping k-mers of length `probe_len` at intervals of `step`.
/// Filters out k-mers with >50% N bases. Writes to `output` FASTA.
/// Returns the number of k-mers written.
fn enumerate_kmers(
    seqs: &[(String, String)],
    probe_len: usize,
    step: usize,
    output: &Path,
) -> Result<usize> {
    let file = File::create(output)
        .with_context(|| format!("Cannot create k-mer file: {}", output.display()))?;
    let mut writer = BufWriter::new(file);

    let step = step.max(1);
    let mut global_idx = 0usize;
    let max_n = probe_len / 2; // >50% N threshold

    for (_seq_id, seq) in seqs {
        let seq_bytes = seq.as_bytes();
        let seq_len = seq_bytes.len();
        if seq_len < probe_len {
            continue;
        }
        let mut pos = 0usize;
        while pos + probe_len <= seq_len {
            let kmer = &seq_bytes[pos..pos + probe_len];
            let n_count = kmer.iter().filter(|&&b| b == b'N').count();
            if n_count <= max_n {
                writeln!(writer, ">km{}", global_idx)?;
                writer.write_all(kmer)?;
                writeln!(writer)?;
                global_idx += 1;
            }
            pos += step;
        }
        // Emit k-mer anchored to end of sequence if not already covered
        let last_start = seq_len - probe_len;
        let last_step_start = if seq_len >= probe_len {
            let steps = (last_start) / step;
            steps * step
        } else {
            0
        };
        if last_start != last_step_start {
            let kmer = &seq_bytes[last_start..seq_len];
            let n_count = kmer.iter().filter(|&&b| b == b'N').count();
            if n_count <= max_n {
                writeln!(writer, ">km{}", global_idx)?;
                writer.write_all(kmer)?;
                writeln!(writer)?;
                global_idx += 1;
            }
        }
    }

    writer.flush()?;
    Ok(global_idx)
}

// ---------------------------------------------------------------------------
// cd-hit-est .clstr parsing
// ---------------------------------------------------------------------------

/// Returns path to the .clstr file produced by cd-hit-est for the given output FASTA.
fn clstr_path_for(centroids_path: &Path) -> PathBuf {
    let name = centroids_path
        .file_name()
        .unwrap_or_default()
        .to_string_lossy()
        .to_string();
    centroids_path
        .parent()
        .unwrap_or(Path::new("."))
        .join(format!("{}.clstr", name))
}

/// Parse a cd-hit-est `.clstr` file.
/// Returns a map from centroid sequence name → cluster member count.
fn parse_clstr(clstr_path: &Path) -> Result<HashMap<String, usize>> {
    let file = File::open(clstr_path)
        .with_context(|| format!("Cannot open .clstr: {}", clstr_path.display()))?;
    let reader = BufReader::new(file);

    let mut size_map: HashMap<String, usize> = HashMap::new();
    let mut current_size = 0usize;
    let mut current_rep: Option<String> = None;

    for line in reader.lines() {
        let line = line?;
        let trimmed = line.trim();

        if trimmed.starts_with(">Cluster") {
            if let Some(rep) = current_rep.take() {
                size_map.insert(rep, current_size);
            }
            current_size = 0;
        } else if !trimmed.is_empty() {
            current_size += 1;
            if trimmed.ends_with(" *") {
                if let Some(name) = extract_clstr_name(trimmed) {
                    current_rep = Some(name);
                }
            }
        }
    }
    // Flush final cluster
    if let Some(rep) = current_rep {
        size_map.insert(rep, current_size);
    }

    Ok(size_map)
}

/// Extract the sequence name from a .clstr member line.
///
/// Format: `0\t120nt, >name... *`  (representative)
///      or `1\t120nt, >name... at +/95.00%`
/// cd-hit-est always appends `...` after the (possibly truncated) name.
fn extract_clstr_name(line: &str) -> Option<String> {
    // Split on ", >" to find the name part
    let after_gt = line.split(", >").nth(1)?;
    // Name ends at first "..."
    let name = if let Some(idx) = after_gt.find("...") {
        &after_gt[..idx]
    } else {
        // Fallback: take until first whitespace
        after_gt.split_whitespace().next()?
    };
    Some(name.to_string())
}

// ---------------------------------------------------------------------------
// Centroid ranking
// ---------------------------------------------------------------------------

/// Load centroids FASTA and sort by cluster size (descending).
/// Returns Vec<(name, sequence, cluster_size)>.
fn rank_centroids(
    centroids_path: &Path,
    size_map: &HashMap<String, usize>,
) -> Result<Vec<(String, String, usize)>> {
    let seqs = load_fasta(centroids_path)?;
    let mut ranked: Vec<(String, String, usize)> = seqs
        .into_iter()
        .map(|(name, seq)| {
            let size = size_map.get(&name).copied().unwrap_or(1);
            (name, seq, size)
        })
        .collect();
    ranked.sort_by(|a, b| b.2.cmp(&a.2)); // descending by cluster size
    Ok(ranked)
}

// ---------------------------------------------------------------------------
// Coverage computation
// ---------------------------------------------------------------------------

/// Compute the Nth-percentile of per-target coverage fractions across all targets.
///
/// For each target in `targets`, the coverage fraction is the proportion of
/// positions with depth >= `min_depth`. Targets absent from `depth_map` are
/// treated as 0% covered. Returns a value in [0.0, 1.0].
fn percentile_coverage_for_targets(
    depth_map: &HashMap<String, Vec<u32>>,
    targets: &[(String, String)],
    min_depth: u32,
    percentile: f64,
) -> f64 {
    if targets.is_empty() {
        return 0.0;
    }

    let mut fractions: Vec<f64> = targets
        .iter()
        .map(|(id, seq)| {
            match depth_map.get(id) {
                Some(depths) if !depths.is_empty() => {
                    let covered = depths.iter().filter(|&&d| d >= min_depth).count();
                    covered as f64 / depths.len() as f64
                }
                _ => {
                    // Absent from SAM @SQ headers or empty depths: 0% covered
                    if seq.is_empty() { 1.0 } else { 0.0 }
                }
            }
        })
        .collect();

    fractions.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

    let n = fractions.len();
    let idx = ((n as f64 * percentile / 100.0).floor() as usize).min(n.saturating_sub(1));
    fractions[idx]
}

// ---------------------------------------------------------------------------
// Low-coverage region extraction
// ---------------------------------------------------------------------------

/// Write under-covered sub-sequences of each target to `output` FASTA.
///
/// For each target, finds consecutive runs of positions where depth < `min_depth`,
/// expands runs shorter than `probe_len` bidirectionally, merges overlapping
/// expanded regions, then extracts and writes the corresponding sub-sequences.
///
/// Returns the number of regions written.
fn extract_low_coverage_regions(
    depth_map: &HashMap<String, Vec<u32>>,
    targets: &[(String, String)],
    min_depth: u32,
    probe_len: usize,
    output: &Path,
) -> Result<usize> {
    let file = File::create(output)
        .with_context(|| format!("Cannot create low-cov file: {}", output.display()))?;
    let mut writer = BufWriter::new(file);
    let mut count = 0usize;

    for (seq_id, seq) in targets {
        let seq_len = seq.len();
        if seq_len < probe_len {
            continue;
        }

        let regions = match depth_map.get(seq_id) {
            None => {
                // Entirely uncovered: extract the whole sequence
                vec![(0usize, seq_len)]
            }
            Some(depths) => {
                find_low_coverage_regions(depths, min_depth, probe_len, seq_len)
            }
        };

        for (lcr_idx, (start, end)) in regions.iter().enumerate() {
            let region = &seq[*start..*end];
            if region.len() >= probe_len {
                writeln!(writer, ">{}_lcr{}", seq_id, lcr_idx)?;
                writeln!(writer, "{}", region)?;
                count += 1;
            }
        }
    }

    writer.flush()?;
    Ok(count)
}

/// Identify low-coverage runs in a depth vector, expand short runs to `probe_len`,
/// and merge overlapping regions. Returns (start, end) pairs (0-based, exclusive end).
fn find_low_coverage_regions(
    depths: &[u32],
    min_depth: u32,
    probe_len: usize,
    seq_len: usize,
) -> Vec<(usize, usize)> {
    if seq_len == 0 {
        return Vec::new();
    }

    // Find raw low-coverage runs
    let depth_len = depths.len().min(seq_len);
    let mut raw_runs: Vec<(usize, usize)> = Vec::new();
    let mut run_start: Option<usize> = None;

    for (i, &depth) in depths.iter().enumerate().take(depth_len) {
        if depth < min_depth {
            if run_start.is_none() {
                run_start = Some(i);
            }
        } else if let Some(s) = run_start.take() {
            raw_runs.push((s, i));
        }
    }
    if let Some(s) = run_start {
        raw_runs.push((s, depth_len));
    }
    // Positions beyond depth vector are also uncovered
    if depth_len < seq_len {
        raw_runs.push((depth_len, seq_len));
    }

    if raw_runs.is_empty() {
        return Vec::new();
    }

    // Expand runs shorter than probe_len and merge overlapping regions
    let mut expanded: Vec<(usize, usize)> = raw_runs
        .iter()
        .map(|&(s, e)| expand_run(s, e, probe_len, seq_len))
        .collect();

    expanded.sort_by_key(|&(s, _)| s);

    let mut merged: Vec<(usize, usize)> = Vec::new();
    for (s, e) in expanded {
        match merged.last_mut() {
            Some(last) if s < last.1 => last.1 = last.1.max(e),
            _ => merged.push((s, e)),
        }
    }

    merged
}

/// Expand a [start, end) run to at least `probe_len` bases, clipped to [0, seq_len).
fn expand_run(start: usize, end: usize, probe_len: usize, seq_len: usize) -> (usize, usize) {
    let len = end - start;
    if len >= probe_len {
        return (start, end);
    }
    // Centre the window on the middle of the run
    let center = (start + end) / 2;
    let half = probe_len / 2;
    let new_start = center.saturating_sub(half);
    let new_end = (new_start + probe_len).min(seq_len);
    // Shift left if we hit the right edge
    let new_start = if new_end == seq_len && seq_len >= probe_len {
        seq_len - probe_len
    } else {
        new_start
    };
    (new_start, new_end)
}

// ---------------------------------------------------------------------------
// FASTA utilities
// ---------------------------------------------------------------------------

/// Load all sequences from a FASTA file (preserves insertion order).
fn load_fasta(path: &Path) -> Result<Vec<(String, String)>> {
    let file = File::open(path)
        .with_context(|| format!("Cannot open FASTA: {}", path.display()))?;
    let reader = BufReader::new(file);

    let mut seqs: Vec<(String, String)> = Vec::new();
    let mut current_id: Option<String> = None;
    let mut current_seq = String::new();

    for line in reader.lines() {
        let line = line?;
        let trimmed = line.trim_end();
        if let Some(stripped) = trimmed.strip_prefix('>') {
            if let Some(id) = current_id.take() {
                if !current_seq.is_empty() {
                    seqs.push((id, std::mem::take(&mut current_seq)));
                }
            }
            current_id = Some(
                stripped
                    .split_whitespace()
                    .next()
                    .unwrap_or("")
                    .to_string(),
            );
        } else if !trimmed.is_empty() {
            current_seq.push_str(&trimmed.to_uppercase());
        }
    }
    if let Some(id) = current_id {
        if !current_seq.is_empty() {
            seqs.push((id, current_seq));
        }
    }

    Ok(seqs)
}
