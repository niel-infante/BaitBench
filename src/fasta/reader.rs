use anyhow::{Context, Result};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

/// Parse a FASTA file into a map of sequence ID -> sequence string.
/// Multi-line sequences are concatenated and uppercased.
pub fn parse_fasta(path: &Path) -> Result<HashMap<String, String>> {
    let file = File::open(path).with_context(|| format!("Cannot open FASTA: {}", path.display()))?;
    let mut reader = BufReader::new(file);

    let mut sequences = HashMap::new();
    let mut current_id: Option<String> = None;
    // One growing buffer, not a Vec<String> of lines joined at the end. A 2 GB genome
    // at 60 bp per line is ~35 million lines, and the old form allocated a String for
    // every one of them, then held both the fragments and the joined copy at once.
    // That is the difference between a few seconds and minutes of allocator churn, and
    // it is what pushed a large run into swap.
    let mut current_seq = String::new();
    let mut line = String::new();

    loop {
        line.clear();
        if reader.read_line(&mut line)? == 0 {
            break;
        }
        let trimmed = line.trim_end();
        if let Some(header) = trimmed.strip_prefix('>') {
            if let Some(id) = current_id.take() {
                // Growth doubles capacity, so an unshrunk buffer can retain twice the
                // sequence. Across hundreds of contigs that is gigabytes held for
                // nothing.
                current_seq.shrink_to_fit();
                sequences.insert(id, std::mem::take(&mut current_seq));
            }
            current_id = Some(
                header
                    .split_whitespace()
                    .next()
                    .unwrap_or("")
                    .to_string(),
            );
        } else if !trimmed.is_empty() {
            // Uppercase in place rather than allocating a per-line String. FASTA is
            // ASCII, so this matches the previous to_uppercase() behaviour.
            current_seq.extend(trimmed.chars().map(|c| c.to_ascii_uppercase()));
        }
    }

    if let Some(id) = current_id {
        current_seq.shrink_to_fit();
        sequences.insert(id, current_seq);
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
            let id = line.strip_prefix('>').unwrap_or("").split_whitespace().next().unwrap_or("").to_string();
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

// ── Format-aware reads helpers ────────────────────────────────────────────────

/// Detect whether a file is FASTQ by peeking at its first byte.
/// Returns true if the file starts with '@' (FASTQ), false if '>' (FASTA).
pub fn is_fastq(path: &Path) -> Result<bool> {
    let mut file = File::open(path)
        .with_context(|| format!("Cannot open: {}", path.display()))?;
    let mut buf = [0u8; 1];
    match file.read(&mut buf)? {
        0 => Ok(false), // empty file → treat as FASTA
        _ => Ok(buf[0] == b'@'),
    }
}

/// Parse a FASTQ file into a map of sequence ID → sequence string.
/// Quality scores are discarded (only sequences are needed for alignment).
fn parse_fastq(path: &Path) -> Result<HashMap<String, String>> {
    let file = File::open(path).with_context(|| format!("Cannot open FASTQ: {}", path.display()))?;
    let reader = BufReader::new(file);
    let mut sequences = HashMap::new();
    let mut lines = reader.lines();
    while let Some(header) = lines.next() {
        let header = header?;
        let seq    = lines.next().ok_or_else(|| anyhow::anyhow!("Truncated FASTQ: {}", path.display()))??;
        let _plus  = lines.next().ok_or_else(|| anyhow::anyhow!("Truncated FASTQ: {}", path.display()))??;
        let _qual  = lines.next().ok_or_else(|| anyhow::anyhow!("Truncated FASTQ: {}", path.display()))??;
        let id = header.strip_prefix('@').unwrap_or(&header)
            .split_whitespace().next().unwrap_or("").to_string();
        if !id.is_empty() {
            sequences.insert(id, seq.to_uppercase());
        }
    }
    Ok(sequences)
}

/// Parse a FASTQ file and return sequence IDs in order.
fn parse_fastq_ids(path: &Path) -> Result<Vec<String>> {
    let file = File::open(path).with_context(|| format!("Cannot open FASTQ: {}", path.display()))?;
    let reader = BufReader::new(file);
    let mut ids = Vec::new();
    let mut lines = reader.lines();
    while let Some(header) = lines.next() {
        let header = header?;
        let _seq  = lines.next().ok_or_else(|| anyhow::anyhow!("Truncated FASTQ: {}", path.display()))??;
        let _plus = lines.next().ok_or_else(|| anyhow::anyhow!("Truncated FASTQ: {}", path.display()))??;
        let _qual = lines.next().ok_or_else(|| anyhow::anyhow!("Truncated FASTQ: {}", path.display()))??;
        let id = header.strip_prefix('@').unwrap_or(&header)
            .split_whitespace().next().unwrap_or("").to_string();
        if !id.is_empty() {
            ids.push(id);
        }
    }
    Ok(ids)
}

/// Count sequences in a FASTQ file (each record is 4 lines).
fn count_fastq_sequences(path: &Path) -> Result<usize> {
    let file = File::open(path).with_context(|| format!("Cannot open FASTQ: {}", path.display()))?;
    let reader = BufReader::new(file);
    let mut count = 0usize;
    let mut lines = reader.lines();
    while let Some(header) = lines.next() {
        let _h    = header?;
        let _seq  = lines.next().ok_or_else(|| anyhow::anyhow!("Truncated FASTQ: {}", path.display()))??;
        let _plus = lines.next().ok_or_else(|| anyhow::anyhow!("Truncated FASTQ: {}", path.display()))??;
        let _qual = lines.next().ok_or_else(|| anyhow::anyhow!("Truncated FASTQ: {}", path.display()))??;
        count += 1;
    }
    Ok(count)
}

/// Parse a file as FASTA or FASTQ, returning sequence ID → sequence.
/// Format is detected by peeking at the first byte.
/// Quality scores are discarded when the input is FASTQ.
pub fn parse_reads(path: &Path) -> Result<HashMap<String, String>> {
    if is_fastq(path)? { parse_fastq(path) } else { parse_fasta(path) }
}

/// Return sequence IDs in order from a FASTA or FASTQ file.
pub fn parse_reads_ids(path: &Path) -> Result<Vec<String>> {
    if is_fastq(path)? { parse_fastq_ids(path) } else { parse_fasta_ids(path) }
}

/// Count sequences in a FASTA or FASTQ file.
pub fn count_reads(path: &Path) -> Result<usize> {
    if is_fastq(path)? { count_fastq_sequences(path) } else { count_sequences(path) }
}
