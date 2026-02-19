use anyhow::{Context, Result};
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

/// Parse a plain-text ID list file into a HashSet.
/// Skips empty lines and lines starting with '#'.
/// Takes the first whitespace-delimited token from each line.
pub fn parse_id_set(path: &Path) -> Result<HashSet<String>> {
    let file = File::open(path)
        .with_context(|| format!("Cannot open ID file: {}", path.display()))?;
    let reader = BufReader::new(file);
    let mut ids = HashSet::new();

    for line in reader.lines() {
        let line = line?;
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        if let Some(id) = line.split_whitespace().next() {
            ids.insert(id.to_string());
        }
    }

    Ok(ids)
}

/// Extract the source sequence ID from a fragment read name.
///
/// Read names follow the pattern `{seq_id}_fragment_{n}` (with optional
/// trailing fields like `start=... length=...`). We find the last occurrence
/// of `_fragment_` and take everything before it as the source ID.
pub fn extract_source_id(read_name: &str) -> Option<&str> {
    // Take the first whitespace-delimited token (the ID portion)
    let id_part = read_name.split_whitespace().next().unwrap_or(read_name);
    // Find the last occurrence of "_fragment_" to handle IDs that contain that substring
    id_part.rfind("_fragment_").map(|pos| &id_part[..pos])
}

/// Parse a sample manifest TSV file into a HashMap of id -> weight.
///
/// Format: `id<tab>weight` (weight is optional, defaults to 1.0).
/// Skips empty lines and lines starting with '#'.
pub fn parse_sample_manifest(path: &Path) -> Result<HashMap<String, f64>> {
    let file = File::open(path)
        .with_context(|| format!("Cannot open sample manifest: {}", path.display()))?;
    let reader = BufReader::new(file);
    let mut samples = HashMap::new();

    for line in reader.lines() {
        let line = line?;
        let line = line.trim().to_string();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let parts: Vec<&str> = line.split('\t').collect();
        let id = parts[0].trim().to_string();
        let weight = if parts.len() >= 2 {
            parts[1].trim().parse::<f64>()
                .with_context(|| format!("Invalid weight for '{}': {}", id, parts[1].trim()))?
        } else {
            1.0
        };
        samples.insert(id, weight);
    }

    Ok(samples)
}

/// Write a list of IDs to a file, one per line.
pub fn write_id_list(ids: &[String], path: &Path) -> Result<()> {
    let file = File::create(path)
        .with_context(|| format!("Cannot create ID file: {}", path.display()))?;
    let mut writer = BufWriter::new(file);
    for id in ids {
        writeln!(writer, "{}", id)?;
    }
    writer.flush()?;
    Ok(())
}
