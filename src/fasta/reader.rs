use anyhow::{Context, Result};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

/// Parse a FASTA file into a map of sequence ID -> sequence string.
/// Multi-line sequences are concatenated and uppercased.
pub fn parse_fasta(path: &Path) -> Result<HashMap<String, String>> {
    let file = File::open(path).with_context(|| format!("Cannot open FASTA: {}", path.display()))?;
    let reader = BufReader::new(file);

    let mut sequences = HashMap::new();
    let mut current_id: Option<String> = None;
    let mut current_seq = Vec::new();

    for line in reader.lines() {
        let line = line?;
        let line = line.trim_end();
        if line.starts_with('>') {
            if let Some(id) = current_id.take() {
                sequences.insert(id, current_seq.join(""));
                current_seq.clear();
            }
            let id = line[1..].split_whitespace().next().unwrap_or("").to_string();
            current_id = Some(id);
        } else if !line.is_empty() {
            current_seq.push(line.to_uppercase());
        }
    }

    if let Some(id) = current_id {
        sequences.insert(id, current_seq.join(""));
    }

    Ok(sequences)
}

/// Parse FASTA file and return just the sequence IDs in order.
pub fn parse_fasta_ids(path: &Path) -> Result<Vec<String>> {
    let file = File::open(path).with_context(|| format!("Cannot open FASTA: {}", path.display()))?;
    let reader = BufReader::new(file);

    let mut ids = Vec::new();
    for line in reader.lines() {
        let line = line?;
        if line.starts_with('>') {
            let id = line[1..].split_whitespace().next().unwrap_or("").to_string();
            ids.push(id);
        }
    }

    Ok(ids)
}

/// Count the number of sequences in a FASTA file (counts '>' lines).
pub fn count_sequences(path: &Path) -> Result<usize> {
    let file = File::open(path).with_context(|| format!("Cannot open FASTA: {}", path.display()))?;
    let reader = BufReader::new(file);

    let mut count = 0;
    for line in reader.lines() {
        let line = line?;
        if line.starts_with('>') {
            count += 1;
        }
    }

    Ok(count)
}
