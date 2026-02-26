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

/// Write per-position coverage to a TSV file.
///
/// Format: reference_id\tposition\tdepth (1-based positions)
/// Only writes entries for references that have at least one position with depth > 0.
pub fn write_coverage_tsv(
    path: &Path,
    coverage: &HashMap<String, Vec<u32>>,
) -> Result<()> {
    let file = File::create(path)
        .with_context(|| format!("Cannot create coverage TSV: {}", path.display()))?;
    let mut w = BufWriter::new(file);

    writeln!(w, "reference_id\tposition\tdepth")?;

    // Sort reference names for deterministic output
    let mut ref_names: Vec<&String> = coverage.keys().collect();
    ref_names.sort();

    for ref_name in ref_names {
        let depths = &coverage[ref_name];
        let has_coverage = depths.iter().any(|&d| d > 0);
        if !has_coverage {
            continue;
        }
        for (pos, &depth) in depths.iter().enumerate() {
            writeln!(w, "{}\t{}\t{}", ref_name, pos + 1, depth)?; // 1-based position
        }
    }

    w.flush()?;
    Ok(())
}
