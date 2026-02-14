use anyhow::{Context, Result};
use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

/// Filter PAF records and return the set of passing read IDs.
///
/// Filters based on:
/// - No indels (I or D in cg:Z: CIGAR tag)
/// - Mismatches (NM:i: tag) <= max_mismatches
/// - Matching bases (column 10) >= min_match_bases
pub fn filter_paf(
    paf_path: &Path,
    max_mismatches: u32,
    min_match_bases: u32,
) -> Result<HashSet<String>> {
    let file = File::open(paf_path)
        .with_context(|| format!("Cannot open PAF: {}", paf_path.display()))?;
    let reader = BufReader::new(file);

    let mut passing = HashSet::new();
    let mut total_mappings = 0u64;

    for line in reader.lines() {
        let line = line?;
        let line = line.trim();
        if line.is_empty() {
            continue;
        }

        total_mappings += 1;
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 12 {
            continue;
        }

        let query_name = fields[0];
        let matching_bases: u32 = fields[9].parse().unwrap_or(0);

        if matching_bases < min_match_bases {
            continue;
        }

        // Scan optional tags (fields 12+)
        let mut nm: u32 = 0;
        let mut has_indel = false;

        for field in &fields[12..] {
            if let Some(val) = field.strip_prefix("NM:i:") {
                nm = val.parse().unwrap_or(0);
            } else if let Some(cigar) = field.strip_prefix("cg:Z:") {
                if cigar.contains('I') || cigar.contains('D') {
                    has_indel = true;
                }
            }
        }

        if !has_indel && nm <= max_mismatches {
            passing.insert(query_name.to_string());
        }
    }

    log::info!(
        "PAF filtering: {} total mappings, {} reads passing",
        total_mappings,
        passing.len()
    );

    Ok(passing)
}
