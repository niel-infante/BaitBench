use anyhow::{Context, Result};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use crate::sdust;

/// Count sequences and total bases in a FASTA file (streaming).
pub fn count_fasta_stats(path: &Path) -> Result<(usize, usize)> {
    let file =
        File::open(path).with_context(|| format!("Cannot open FASTA: {}", path.display()))?;
    let reader = BufReader::new(file);

    let mut num_seqs = 0;
    let mut total_bases = 0;

    for line in reader.lines() {
        let line = line?;
        let trimmed = line.trim_end();
        if trimmed.starts_with('>') {
            num_seqs += 1;
        } else if !trimmed.is_empty() {
            total_bases += trimmed.len();
        }
    }

    Ok((num_seqs, total_bases))
}

fn compute_n_frac(seq: &str) -> f64 {
    if seq.is_empty() {
        return 0.0;
    }
    let n_count = seq
        .chars()
        .filter(|c| !matches!(c, 'A' | 'C' | 'G' | 'T' | 'a' | 'c' | 'g' | 't'))
        .count();
    n_count as f64 / seq.len() as f64
}

/// Filter FASTA sequences by N (ambiguous base) content (streaming).
///
/// Keeps sequences whose N fraction is at most `max_n_frac`.
pub fn filter_n_content(input: &Path, output: &Path, max_n_frac: f64) -> Result<()> {
    let file =
        File::open(input).with_context(|| format!("Cannot open FASTA: {}", input.display()))?;
    let reader = BufReader::new(file);

    let out_file =
        File::create(output).with_context(|| format!("Cannot create: {}", output.display()))?;
    let mut writer = BufWriter::new(out_file);

    let mut current_id: Option<String> = None;
    let mut current_seq = String::new();

    for line in reader.lines() {
        let line = line?;
        let trimmed = line.trim_end();
        if trimmed.starts_with('>') {
            if let Some(ref id) = current_id {
                if compute_n_frac(&current_seq) <= max_n_frac {
                    writeln!(writer, ">{}", id)?;
                    writeln!(writer, "{}", current_seq)?;
                }
            }
            current_id = Some(trimmed.strip_prefix('>').unwrap_or("").to_string());
            current_seq.clear();
        } else if !trimmed.is_empty() {
            current_seq.push_str(&trimmed.to_uppercase());
        }
    }

    if let Some(ref id) = current_id {
        if compute_n_frac(&current_seq) <= max_n_frac {
            writeln!(writer, ">{}", id)?;
            writeln!(writer, "{}", current_seq)?;
        }
    }

    writer.flush()?;
    Ok(())
}

/// Filter FASTA sequences shorter than `min_length` (streaming).
pub fn filter_short_sequences(input: &Path, output: &Path, min_length: usize) -> Result<()> {
    let file =
        File::open(input).with_context(|| format!("Cannot open FASTA: {}", input.display()))?;
    let reader = BufReader::new(file);

    let out_file =
        File::create(output).with_context(|| format!("Cannot create: {}", output.display()))?;
    let mut writer = BufWriter::new(out_file);

    let mut current_id: Option<String> = None;
    let mut current_seq = String::new();

    for line in reader.lines() {
        let line = line?;
        let trimmed = line.trim_end();
        if trimmed.starts_with('>') {
            if let Some(ref id) = current_id {
                if current_seq.len() >= min_length {
                    writeln!(writer, ">{}", id)?;
                    writeln!(writer, "{}", current_seq)?;
                }
            }
            current_id = Some(trimmed.strip_prefix('>').unwrap_or("").to_string());
            current_seq.clear();
        } else if !trimmed.is_empty() {
            current_seq.push_str(&trimmed.to_uppercase());
        }
    }

    if let Some(ref id) = current_id {
        if current_seq.len() >= min_length {
            writeln!(writer, ">{}", id)?;
            writeln!(writer, "{}", current_seq)?;
        }
    }

    writer.flush()?;
    Ok(())
}

fn compute_gc(seq: &str) -> f64 {
    if seq.is_empty() {
        return 0.0;
    }
    let gc_count = seq
        .chars()
        .filter(|c| matches!(c, 'G' | 'C' | 'g' | 'c'))
        .count();
    gc_count as f64 / seq.len() as f64
}

/// Filter FASTA sequences by GC content (streaming).
pub fn filter_gc(input: &Path, output: &Path, min_gc: f64, max_gc: f64) -> Result<()> {
    let file =
        File::open(input).with_context(|| format!("Cannot open FASTA: {}", input.display()))?;
    let reader = BufReader::new(file);

    let out_file =
        File::create(output).with_context(|| format!("Cannot create: {}", output.display()))?;
    let mut writer = BufWriter::new(out_file);

    let mut current_id: Option<String> = None;
    let mut current_seq = String::new();

    for line in reader.lines() {
        let line = line?;
        let trimmed = line.trim_end();
        if trimmed.starts_with('>') {
            if let Some(ref id) = current_id {
                let gc = compute_gc(&current_seq);
                if gc >= min_gc && gc <= max_gc {
                    writeln!(writer, ">{}", id)?;
                    writeln!(writer, "{}", current_seq)?;
                }
            }
            current_id = Some(trimmed.strip_prefix('>').unwrap_or("").to_string());
            current_seq.clear();
        } else if !trimmed.is_empty() {
            current_seq.push_str(&trimmed.to_uppercase());
        }
    }

    if let Some(ref id) = current_id {
        let gc = compute_gc(&current_seq);
        if gc >= min_gc && gc <= max_gc {
            writeln!(writer, ">{}", id)?;
            writeln!(writer, "{}", current_seq)?;
        }
    }

    writer.flush()?;
    Ok(())
}

/// Filter FASTA sequences by low-complexity content using sDUST (streaming).
///
/// Keeps sequences whose masked fraction is at most `max_masked_frac`.
pub fn filter_complexity(
    input: &Path,
    output: &Path,
    dust_threshold: f64,
    dust_window: usize,
    max_masked_frac: f64,
) -> Result<()> {
    let file =
        File::open(input).with_context(|| format!("Cannot open FASTA: {}", input.display()))?;
    let reader = BufReader::new(file);

    let out_file =
        File::create(output).with_context(|| format!("Cannot create: {}", output.display()))?;
    let mut writer = BufWriter::new(out_file);

    let mut current_id: Option<String> = None;
    let mut current_seq = String::new();

    for line in reader.lines() {
        let line = line?;
        let trimmed = line.trim_end();
        if trimmed.starts_with('>') {
            if let Some(ref id) = current_id {
                let frac = sdust::masked_fraction(current_seq.as_bytes(), dust_threshold, dust_window);
                if frac <= max_masked_frac {
                    writeln!(writer, ">{}", id)?;
                    writeln!(writer, "{}", current_seq)?;
                }
            }
            current_id = Some(trimmed.strip_prefix('>').unwrap_or("").to_string());
            current_seq.clear();
        } else if !trimmed.is_empty() {
            current_seq.push_str(&trimmed.to_uppercase());
        }
    }

    if let Some(ref id) = current_id {
        let frac = sdust::masked_fraction(current_seq.as_bytes(), dust_threshold, dust_window);
        if frac <= max_masked_frac {
            writeln!(writer, ">{}", id)?;
            writeln!(writer, "{}", current_seq)?;
        }
    }

    writer.flush()?;
    Ok(())
}
