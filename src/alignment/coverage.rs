use anyhow::{Context, Result};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

/// Coverage statistics for a single reference sequence.
pub struct CoverageStats {
    pub ref_length: usize,
    pub avg_coverage: f64,
    pub pct_covered_5x: f64,
    pub pct_covered_20x: f64,
}

/// Result of computing coverage across all references.
pub struct CoverageResult {
    /// Reference name → sequence length (from @SQ headers)
    pub ref_lengths: HashMap<String, usize>,
    /// Reference name → per-position depth vector (0-indexed, length = ref_length)
    pub coverage: HashMap<String, Vec<u32>>,
    /// Reference name → count of primary alignments mapping to it
    pub reads_per_ref: HashMap<String, usize>,
}

/// Parse CIGAR string into (length, operation) pairs.
///
/// E.g. "50M2I3D10M" -> [(50, 'M'), (2, 'I'), (3, 'D'), (10, 'M')]
fn parse_cigar(cigar: &str) -> Vec<(u32, char)> {
    let mut ops = Vec::new();
    let mut num_start = 0;
    for (i, c) in cigar.char_indices() {
        if c.is_ascii_alphabetic() || c == '=' {
            if let Ok(len) = cigar[num_start..i].parse::<u32>() {
                ops.push((len, c));
            }
            num_start = i + c.len_utf8();
        }
    }
    ops
}

/// Compute per-position coverage for all references from a SAM file.
///
/// Single-pass algorithm:
/// 1. Parse @SQ headers to get reference lengths and pre-allocate coverage vectors
/// 2. For each primary alignment, parse POS and CIGAR
/// 3. Walk through CIGAR ops, incrementing depth for M/=/X positions
///
/// Skips: unmapped (RNAME=*), secondary (FLAG & 0x100), supplementary (FLAG & 0x800).
pub fn compute_coverage(sam_path: &Path) -> Result<CoverageResult> {
    let file = File::open(sam_path)
        .with_context(|| format!("Cannot open SAM: {}", sam_path.display()))?;
    let reader = BufReader::new(file);

    let mut ref_lengths: HashMap<String, usize> = HashMap::new();
    let mut coverage: HashMap<String, Vec<u32>> = HashMap::new();
    let mut reads_per_ref: HashMap<String, usize> = HashMap::new();

    for line in reader.lines() {
        let line = line?;

        // Parse @SQ headers for reference lengths
        if line.starts_with('@') {
            if line.starts_with("@SQ") {
                let mut name = None;
                let mut length = None;
                for field in line.split('\t') {
                    if let Some(sn) = field.strip_prefix("SN:") {
                        name = Some(sn.to_string());
                    } else if let Some(ln) = field.strip_prefix("LN:") {
                        length = ln.parse::<usize>().ok();
                    }
                }
                if let (Some(n), Some(l)) = (name, length) {
                    ref_lengths.insert(n.clone(), l);
                    coverage.insert(n, vec![0u32; l]);
                }
            }
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 6 {
            continue;
        }

        let flag: u16 = fields[1].parse().unwrap_or(0);
        let rname = fields[2];
        let pos: usize = fields[3].parse().unwrap_or(0); // 1-based
        let cigar = fields[5];

        // Skip unmapped, secondary, supplementary
        if rname == "*" || cigar == "*" {
            continue;
        }
        if flag & 0x100 != 0 || flag & 0x800 != 0 {
            continue;
        }

        *reads_per_ref.entry(rname.to_string()).or_insert(0) += 1;

        // Get coverage vector for this reference
        let cov = match coverage.get_mut(rname) {
            Some(c) => c,
            None => continue, // reference not in @SQ headers
        };

        // Walk CIGAR and increment coverage for M/=/X ops
        let ops = parse_cigar(cigar);
        let mut ref_pos = pos.saturating_sub(1); // convert to 0-based

        for (len, op) in ops {
            match op {
                'M' | '=' | 'X' => {
                    // Alignment match/mismatch: read covers these reference positions
                    for _ in 0..len {
                        if ref_pos < cov.len() {
                            cov[ref_pos] += 1;
                        }
                        ref_pos += 1;
                    }
                }
                'D' | 'N' => {
                    // Deletion/skip: consume reference but no coverage
                    ref_pos += len as usize;
                }
                'I' | 'S' | 'H' | 'P' => {
                    // Insertion/clip/pad: do not consume reference
                }
                _ => {}
            }
        }
    }

    Ok(CoverageResult {
        ref_lengths,
        coverage,
        reads_per_ref,
    })
}

/// Calculate coverage statistics from a per-position depth vector.
pub fn calculate_stats(coverage: &[u32], ref_length: usize) -> CoverageStats {
    if ref_length == 0 {
        return CoverageStats {
            ref_length: 0,
            avg_coverage: 0.0,
            pct_covered_5x: 0.0,
            pct_covered_20x: 0.0,
        };
    }

    let total_depth: u64 = coverage.iter().map(|&d| d as u64).sum();
    let avg_coverage = total_depth as f64 / ref_length as f64;

    let above_5x = coverage.iter().filter(|&&d| d >= 5).count();
    let above_20x = coverage.iter().filter(|&&d| d >= 20).count();

    CoverageStats {
        ref_length,
        avg_coverage,
        pct_covered_5x: above_5x as f64 / ref_length as f64 * 100.0,
        pct_covered_20x: above_20x as f64 / ref_length as f64 * 100.0,
    }
}

// ---------------------------------------------------------------------------
// Probe tiling coverage (separate from read/fragment coverage above)
// ---------------------------------------------------------------------------

/// Extended coverage statistics for probe tiling QC.
///
/// These metrics evaluate how well probes physically tile target sequences,
/// which is distinct from read/fragment coverage in the simulation pipeline.
pub struct ProbeCoverageStats {
    pub ref_length: usize,
    pub covered_bases: usize,
    pub pct_covered_1x: f64,
    pub mean_depth: f64,
    pub median_depth: f64,
    pub pct_covered_2x: f64,
    pub pct_covered_5x: f64,
    pub pct_covered_10x: f64,
    pub max_gap_length: usize,
    pub num_gaps: usize,
    pub pct_near_probe: f64,
}

/// Compute per-position probe depth from a SAM file.
///
/// Same algorithm as `compute_coverage` but **includes secondary alignments**
/// (does not skip FLAG & 0x100). This is necessary because a single probe can
/// legitimately map to conserved regions across multiple target sequences.
///
/// Still skips: unmapped (RNAME=*), supplementary (FLAG & 0x800).
pub fn compute_probe_coverage(sam_path: &Path) -> Result<CoverageResult> {
    let file = File::open(sam_path)
        .with_context(|| format!("Cannot open SAM: {}", sam_path.display()))?;
    let reader = BufReader::new(file);

    let mut ref_lengths: HashMap<String, usize> = HashMap::new();
    let mut coverage: HashMap<String, Vec<u32>> = HashMap::new();

    for line in reader.lines() {
        let line = line?;

        // Parse @SQ headers for reference lengths
        if line.starts_with('@') {
            if line.starts_with("@SQ") {
                let mut name = None;
                let mut length = None;
                for field in line.split('\t') {
                    if let Some(sn) = field.strip_prefix("SN:") {
                        name = Some(sn.to_string());
                    } else if let Some(ln) = field.strip_prefix("LN:") {
                        length = ln.parse::<usize>().ok();
                    }
                }
                if let (Some(n), Some(l)) = (name, length) {
                    ref_lengths.insert(n.clone(), l);
                    coverage.insert(n, vec![0u32; l]);
                }
            }
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 6 {
            continue;
        }

        let flag: u16 = fields[1].parse().unwrap_or(0);
        let rname = fields[2];
        let pos: usize = fields[3].parse().unwrap_or(0); // 1-based
        let cigar = fields[5];

        // Skip unmapped and supplementary (but NOT secondary)
        if rname == "*" || cigar == "*" {
            continue;
        }
        if flag & 0x800 != 0 {
            continue;
        }

        // Get coverage vector for this reference
        let cov = match coverage.get_mut(rname) {
            Some(c) => c,
            None => continue,
        };

        // Walk CIGAR and increment depth for M/=/X ops
        let ops = parse_cigar(cigar);
        let mut ref_pos = pos.saturating_sub(1); // convert to 0-based

        for (len, op) in ops {
            match op {
                'M' | '=' | 'X' => {
                    for _ in 0..len {
                        if ref_pos < cov.len() {
                            cov[ref_pos] += 1;
                        }
                        ref_pos += 1;
                    }
                }
                'D' | 'N' => {
                    ref_pos += len as usize;
                }
                'I' | 'S' | 'H' | 'P' => {}
                _ => {}
            }
        }
    }

    Ok(CoverageResult {
        ref_lengths,
        coverage,
        reads_per_ref: HashMap::new(),
    })
}

/// Calculate probe tiling quality statistics from a per-position depth vector.
///
/// `proximity_bp`: bases on each side of a covered position that count as
/// "near a probe" (the capture pull-down zone).
pub fn calculate_probe_stats(
    coverage: &[u32],
    ref_length: usize,
    proximity_bp: usize,
) -> ProbeCoverageStats {
    if ref_length == 0 {
        return ProbeCoverageStats {
            ref_length: 0,
            covered_bases: 0,
            pct_covered_1x: 0.0,
            mean_depth: 0.0,
            median_depth: 0.0,
            pct_covered_2x: 0.0,
            pct_covered_5x: 0.0,
            pct_covered_10x: 0.0,
            max_gap_length: 0,
            num_gaps: 0,
            pct_near_probe: 0.0,
        };
    }

    let total_depth: u64 = coverage.iter().map(|&d| d as u64).sum();
    let mean_depth = total_depth as f64 / ref_length as f64;

    // Median: sort depths, take middle value
    let mut sorted = coverage.to_vec();
    sorted.sort_unstable();
    let median_depth = if sorted.len() % 2 == 0 {
        (sorted[sorted.len() / 2 - 1] as f64 + sorted[sorted.len() / 2] as f64) / 2.0
    } else {
        sorted[sorted.len() / 2] as f64
    };

    // Tiered coverage
    let covered_bases = coverage.iter().filter(|&&d| d >= 1).count();
    let above_2x = coverage.iter().filter(|&&d| d >= 2).count();
    let above_5x = coverage.iter().filter(|&&d| d >= 5).count();
    let above_10x = coverage.iter().filter(|&&d| d >= 10).count();

    // Gap analysis: consecutive stretches of depth == 0
    let mut max_gap_length = 0usize;
    let mut num_gaps = 0usize;
    let mut current_gap = 0usize;
    for &d in coverage {
        if d == 0 {
            current_gap += 1;
        } else {
            if current_gap > 0 {
                num_gaps += 1;
                if current_gap > max_gap_length {
                    max_gap_length = current_gap;
                }
                current_gap = 0;
            }
        }
    }
    if current_gap > 0 {
        num_gaps += 1;
        if current_gap > max_gap_length {
            max_gap_length = current_gap;
        }
    }

    // Proximity metric: expand each covered position by ±proximity_bp
    let mut near_probe = vec![false; ref_length];
    for (i, &d) in coverage.iter().enumerate() {
        if d > 0 {
            let start = i.saturating_sub(proximity_bp);
            let end = std::cmp::min(i + proximity_bp + 1, ref_length);
            for pos in start..end {
                near_probe[pos] = true;
            }
        }
    }
    let near_count = near_probe.iter().filter(|&&v| v).count();

    ProbeCoverageStats {
        ref_length,
        covered_bases,
        pct_covered_1x: covered_bases as f64 / ref_length as f64 * 100.0,
        mean_depth,
        median_depth,
        pct_covered_2x: above_2x as f64 / ref_length as f64 * 100.0,
        pct_covered_5x: above_5x as f64 / ref_length as f64 * 100.0,
        pct_covered_10x: above_10x as f64 / ref_length as f64 * 100.0,
        max_gap_length,
        num_gaps,
        pct_near_probe: near_count as f64 / ref_length as f64 * 100.0,
    }
}

/// Write run-length encoded coverage intervals to a TSV file.
///
/// Format: reference_id\tstart\tend\tdepth (1-based inclusive coordinates)
/// Consecutive positions with identical depth are collapsed into a single interval.
/// Writes all references, including those with zero coverage (as a single 0-depth interval).
pub fn write_coverage_intervals(
    path: &Path,
    coverage: &HashMap<String, Vec<u32>>,
) -> Result<()> {
    let file = File::create(path)
        .with_context(|| format!("Cannot create coverage TSV: {}", path.display()))?;
    let mut w = BufWriter::new(file);

    writeln!(w, "reference_id\tstart\tend\tdepth")?;

    // Sort reference names for deterministic output
    let mut ref_names: Vec<&String> = coverage.keys().collect();
    ref_names.sort();

    for ref_name in ref_names {
        let depths = &coverage[ref_name];
        if depths.is_empty() {
            continue;
        }

        let mut run_start: usize = 0; // 0-based index
        let mut run_depth: u32 = depths[0];

        for i in 1..depths.len() {
            if depths[i] != run_depth {
                // Emit the previous run (1-based inclusive)
                writeln!(w, "{}\t{}\t{}\t{}", ref_name, run_start + 1, i, run_depth)?;
                run_start = i;
                run_depth = depths[i];
            }
        }
        // Emit the final run
        writeln!(
            w,
            "{}\t{}\t{}\t{}",
            ref_name,
            run_start + 1,
            depths.len(),
            run_depth
        )?;
    }

    w.flush()?;
    Ok(())
}
