use anyhow::{Context, Result, bail};
use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;
use std::process::Command;

/// Check that blastn is available on PATH.
pub fn check_available() -> Result<()> {
    let status = Command::new("blastn")
        .arg("-version")
        .stdout(std::process::Stdio::null())
        .stderr(std::process::Stdio::null())
        .status();

    match status {
        Ok(s) if s.success() => Ok(()),
        _ => bail!("blastn not found on PATH. Please install BLAST+."),
    }
}

/// Run blastn for probe capture.
///
/// Uses blastn-short task with tabular output format.
pub fn capture_align(
    blast_db: &str,
    reads: &Path,
    output_tsv: &Path,
    log_file: &Path,
    threads: usize,
) -> Result<()> {
    let log = File::create(log_file)
        .with_context(|| format!("Cannot create log: {}", log_file.display()))?;

    let status = Command::new("blastn")
        .args(["-task", "blastn-short"])
        .args(["-db", blast_db])
        .arg("-query")
        .arg(reads)
        .arg("-out")
        .arg(output_tsv)
        .args(["-outfmt", "6 qseqid sseqid nident length mismatch gapopen"])
        .args(["-max_hsps", "1"])
        .args(["-max_target_seqs", "1"])
        .args(["-evalue", "1000"])
        .args(["-word_size", "7"])
        .args(["-dust", "no"])
        .args(["-soft_masking", "false"])
        .args(["-num_threads", &threads.to_string()])
        .stderr(log)
        .status()
        .context("Failed to execute blastn")?;

    if !status.success() {
        bail!("blastn failed (exit code {:?})", status.code());
    }

    Ok(())
}

/// Filter BLAST tabular output and return passing read IDs.
///
/// Filters: first hit per read, no gaps (gapopen == 0), nident >= min_match_bases.
pub fn filter_blast_results(
    blast_tsv: &Path,
    min_match_bases: u32,
) -> Result<HashSet<String>> {
    let file = File::open(blast_tsv)
        .with_context(|| format!("Cannot open BLAST results: {}", blast_tsv.display()))?;
    let reader = BufReader::new(file);

    let mut passing = HashSet::new();
    let mut total_hits = 0u64;

    for line in reader.lines() {
        let line = line?;
        let line = line.trim();
        if line.is_empty() {
            continue;
        }

        total_hits += 1;
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 6 {
            continue;
        }

        let qseqid = fields[0];
        let nident: u32 = fields[2].parse().unwrap_or(0);
        let gapopen: u32 = fields[5].parse().unwrap_or(1);

        // First hit per read only, no gaps, sufficient identity
        if !passing.contains(qseqid) && gapopen == 0 && nident >= min_match_bases {
            passing.insert(qseqid.to_string());
        }
    }

    log::info!(
        "BLAST filtering: {} total hits, {} reads passing",
        total_hits,
        passing.len()
    );

    Ok(passing)
}
