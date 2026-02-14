use anyhow::{Context, Result};
use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

/// Write a single FASTA record to a writer.
#[allow(dead_code)]
pub fn write_fasta_record(writer: &mut impl Write, id: &str, sequence: &str) -> Result<()> {
    writeln!(writer, ">{}", id)?;
    // Write sequence in 80-character lines
    for chunk in sequence.as_bytes().chunks(80) {
        writer.write_all(chunk)?;
        writeln!(writer)?;
    }
    Ok(())
}

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
            let id = line[1..].split_whitespace().next().unwrap_or("");
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

/// Concatenate two FASTA files into one output file.
pub fn concatenate_fastas(fasta1: &Path, fasta2: &Path, output: &Path) -> Result<()> {
    let out_file = File::create(output)
        .with_context(|| format!("Cannot create output: {}", output.display()))?;
    let mut writer = BufWriter::new(out_file);

    // Write first file
    let f1 = File::open(fasta1).with_context(|| format!("Cannot open: {}", fasta1.display()))?;
    let reader1 = BufReader::new(f1);
    for line in reader1.lines() {
        writeln!(writer, "{}", line?)?;
    }

    // Write second file
    let f2 = File::open(fasta2).with_context(|| format!("Cannot open: {}", fasta2.display()))?;
    let reader2 = BufReader::new(f2);
    for line in reader2.lines() {
        writeln!(writer, "{}", line?)?;
    }

    writer.flush()?;
    Ok(())
}
