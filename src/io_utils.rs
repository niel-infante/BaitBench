use anyhow::{Context, Result};
use std::collections::HashSet;
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
