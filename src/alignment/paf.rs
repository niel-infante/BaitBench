use anyhow::{Context, Result};
use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

/// A parsed PAF alignment record with fields needed for cross-reactivity analysis.
#[allow(dead_code)]
pub struct PafRecord {
    pub query_name: String,
    pub query_length: u32,
    pub query_start: u32,
    pub query_end: u32,
    pub target_name: String,
    pub target_length: u32,
    pub target_start: u32,
    pub target_end: u32,
    pub matching_bases: u32,
    pub block_length: u32,
    pub mapq: u32,
}

/// Parse a PAF file into structured records.
///
/// Returns all records without filtering. Skips malformed lines (< 12 columns).
pub fn parse_paf_records(paf_path: &Path) -> Result<Vec<PafRecord>> {
    let file = File::open(paf_path)
        .with_context(|| format!("Cannot open PAF: {}", paf_path.display()))?;
    let reader = BufReader::new(file);

    let mut records = Vec::new();

    for line in reader.lines() {
        let line = line?;
        let line = line.trim();
        if line.is_empty() {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 12 {
            log::debug!("Skipping malformed PAF line ({} columns)", fields.len());
            continue;
        }

        records.push(PafRecord {
            query_name: fields[0].to_string(),
            query_length: fields[1].parse().unwrap_or(0),
            query_start: fields[2].parse().unwrap_or(0),
            query_end: fields[3].parse().unwrap_or(0),
            // fields[4] is strand, skipped
            target_name: fields[5].to_string(),
            target_length: fields[6].parse().unwrap_or(0),
            target_start: fields[7].parse().unwrap_or(0),
            target_end: fields[8].parse().unwrap_or(0),
            matching_bases: fields[9].parse().unwrap_or(0),
            block_length: fields[10].parse().unwrap_or(0),
            mapq: fields[11].parse().unwrap_or(0),
        });
    }

    Ok(records)
}

/// Filter PAF records and return the set of passing read IDs.
///
/// Filters based on:
/// - No indels (I or D in cg:Z: CIGAR tag)
/// - Mismatches (NM:i: tag) <= max_mismatches
/// - Matching bases (column 10) >= min_match_bases
#[allow(dead_code)]
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
        let mut nm: Option<u32> = None;
        let mut has_indel = false;

        for field in &fields[12..] {
            if let Some(val) = field.strip_prefix("NM:i:") {
                nm = val.parse().ok();
            } else if let Some(cigar) = field.strip_prefix("cg:Z:") {
                if cigar.contains('I') || cigar.contains('D') {
                    has_indel = true;
                }
            }
        }

        // NM absent: skip mismatch filter (alignment passes on match-bases alone)
        let nm_passes = nm.map_or(true, |n| n <= max_mismatches);
        if !has_indel && nm_passes {
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
