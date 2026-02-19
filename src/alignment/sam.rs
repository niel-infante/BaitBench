use anyhow::{Context, Result};
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

/// Parse a SAM file and return a map of reference_id -> read_count.
///
/// Skips:
/// - Header lines (starting with @)
/// - Unmapped reads (RNAME == *)
/// - Secondary alignments (FLAG & 0x100)
/// - Supplementary alignments (FLAG & 0x800)
/// - Duplicate read IDs (only counts each read once)
pub fn count_reads_per_reference(sam_path: &Path) -> Result<HashMap<String, usize>> {
    let file = File::open(sam_path)
        .with_context(|| format!("Cannot open SAM: {}", sam_path.display()))?;
    let reader = BufReader::new(file);

    let mut counts: HashMap<String, usize> = HashMap::new();
    let mut seen_reads: HashSet<String> = HashSet::new();

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('@') {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 3 {
            continue;
        }

        let qname = fields[0];
        let flag: u16 = fields[1].parse().unwrap_or(0);
        let rname = fields[2];

        // Skip unmapped
        if rname == "*" {
            continue;
        }

        // Skip secondary (0x100 = 256)
        if flag & 0x100 != 0 {
            continue;
        }

        // Skip supplementary (0x800 = 2048)
        if flag & 0x800 != 0 {
            continue;
        }

        // Deduplicate by read ID
        if seen_reads.insert(qname.to_string()) {
            *counts.entry(rname.to_string()).or_insert(0) += 1;
        }
    }

    Ok(counts)
}

/// Parse a SAM file and return a vec of (read_name, mapped_reference) pairs.
///
/// Only primary alignments are included (skips secondary 0x100,
/// supplementary 0x800, and unmapped). Each read appears at most once.
pub fn get_read_mappings(sam_path: &Path) -> Result<Vec<(String, String)>> {
    let file = File::open(sam_path)
        .with_context(|| format!("Cannot open SAM: {}", sam_path.display()))?;
    let reader = BufReader::new(file);

    let mut mappings = Vec::new();
    let mut seen_reads: HashSet<String> = HashSet::new();

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('@') {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 3 {
            continue;
        }

        let qname = fields[0];
        let flag: u16 = fields[1].parse().unwrap_or(0);
        let rname = fields[2];

        if rname == "*" {
            continue;
        }
        if flag & 0x100 != 0 || flag & 0x800 != 0 {
            continue;
        }

        if seen_reads.insert(qname.to_string()) {
            mappings.push((qname.to_string(), rname.to_string()));
        }
    }

    Ok(mappings)
}

/// Parse SAM file and return the set of read IDs that map (RNAME != *).
/// Used for host filtering.
pub fn get_mapped_read_ids(sam_path: &Path) -> Result<HashSet<String>> {
    let file = File::open(sam_path)
        .with_context(|| format!("Cannot open SAM: {}", sam_path.display()))?;
    let reader = BufReader::new(file);

    let mut mapped = HashSet::new();

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('@') {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 3 {
            continue;
        }

        let rname = fields[2];
        if rname != "*" {
            mapped.insert(fields[0].to_string());
        }
    }

    Ok(mapped)
}
