use anyhow::{Context, Result};
use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use super::reader::is_fastq;

/// Extract sequences from a FASTA file by ID set (replaces seqtk subseq).
/// Single-pass streaming: never loads full file into memory.
/// Returns the number of sequences extracted.
pub fn extract_by_ids(fasta_path: &Path, ids: &HashSet<String>, output_path: &Path) -> Result<usize> {
    let file = File::open(fasta_path)
        .with_context(|| format!("Cannot open FASTA: {}", fasta_path.display()))?;
    let reader = BufReader::new(file);

    let out_file = File::create(output_path)
        .with_context(|| format!("Cannot create output: {}", output_path.display()))?;
    let mut writer = BufWriter::new(out_file);

    let mut writing = false;
    let mut count = 0;

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('>') {
            let id = line.strip_prefix('>').unwrap_or("").split_whitespace().next().unwrap_or("");
            writing = ids.contains(id);
            if writing {
                writeln!(writer, "{}", line)?;
                count += 1;
            }
        } else if writing {
            writeln!(writer, "{}", line)?;
        }
    }

    writer.flush()?;
    Ok(count)
}

/// Extract records from a FASTA or FASTQ file by ID set.
/// Format is detected automatically. FASTQ records are preserved intact (all 4 lines).
/// Returns the number of records extracted.
pub fn extract_reads_by_ids(path: &Path, ids: &HashSet<String>, output: &Path) -> Result<usize> {
    if is_fastq(path)? {
        extract_fastq_by_ids(path, ids, output)
    } else {
        extract_by_ids(path, ids, output)
    }
}

fn extract_fastq_by_ids(fastq_path: &Path, ids: &HashSet<String>, output_path: &Path) -> Result<usize> {
    let file = File::open(fastq_path)
        .with_context(|| format!("Cannot open FASTQ: {}", fastq_path.display()))?;
    let reader = BufReader::new(file);

    let out_file = File::create(output_path)
        .with_context(|| format!("Cannot create output: {}", output_path.display()))?;
    let mut writer = BufWriter::new(out_file);

    let mut count = 0usize;
    let mut lines = reader.lines();
    while let Some(header) = lines.next() {
        let header = header?;
        let seq    = lines.next().ok_or_else(|| anyhow::anyhow!("Truncated FASTQ: {}", fastq_path.display()))??;
        let plus   = lines.next().ok_or_else(|| anyhow::anyhow!("Truncated FASTQ: {}", fastq_path.display()))??;
        let qual   = lines.next().ok_or_else(|| anyhow::anyhow!("Truncated FASTQ: {}", fastq_path.display()))??;
        let id = header.strip_prefix('@').unwrap_or(&header)
            .split_whitespace().next().unwrap_or("");
        if ids.contains(id) {
            writeln!(writer, "{}", header)?;
            writeln!(writer, "{}", seq)?;
            writeln!(writer, "{}", plus)?;
            writeln!(writer, "{}", qual)?;
            count += 1;
        }
    }

    writer.flush()?;
    Ok(count)
}

/// Concatenate multiple FASTA files into one output file.
pub fn concatenate_fastas(inputs: &[&Path], output: &Path) -> Result<()> {
    let out_file = File::create(output)
        .with_context(|| format!("Cannot create output: {}", output.display()))?;
    let mut writer = BufWriter::new(out_file);

    for input in inputs {
        let f = File::open(input).with_context(|| format!("Cannot open: {}", input.display()))?;
        let reader = BufReader::new(f);
        for line in reader.lines() {
            writeln!(writer, "{}", line?)?;
        }
    }

    writer.flush()?;
    Ok(())
}
