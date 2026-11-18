use anyhow::{bail, Context, Result};
use std::fs::{self, File};
use std::path::Path;
use std::process::Command;

use crate::fasta;

/// Check that blastn and makeblastdb are available on PATH.
pub fn check_available() -> Result<()> {
    let blastn_ok = Command::new("blastn")
        .arg("-version")
        .stdout(std::process::Stdio::null())
        .stderr(std::process::Stdio::null())
        .status()
        .map(|s| s.success())
        .unwrap_or(false);
    if !blastn_ok {
        bail!("blastn not found on PATH. Install BLAST+ (conda install -c bioconda blast).");
    }

    let makeblastdb_ok = Command::new("makeblastdb")
        .arg("-version")
        .stdout(std::process::Stdio::null())
        .stderr(std::process::Stdio::null())
        .status()
        .map(|s| s.success())
        .unwrap_or(false);
    if !makeblastdb_ok {
        bail!("makeblastdb not found on PATH. Install BLAST+ (conda install -c bioconda blast).");
    }

    Ok(())
}

/// Align `query` against `reference` with blastn-short, reporting every HSP
/// (not just the best hit per query) so cross-reactivity analysis sees all
/// homologous regions, mirroring `minimap2::xreact_align`'s use of secondary
/// alignments.
///
/// Builds a temporary nucleotide BLAST database from `reference` next to
/// `output_tsv`, runs the search, then removes the database files.
///
/// Output columns (tab-separated, no header): qseqid, sseqid, qlen, qstart,
/// qend, nident, length — parsed by `alignment::blast_tab::parse_blast_hits`.
/// Align probes to a reference for capture simulation, emitting the aligned sequences.
///
/// Output columns: `qseqid sseqid sstart send sstrand qseq sseq`.
///
/// `qseq`/`sseq` are the gapped alignment strings, so the probe/reference base pairs
/// that `delta_g()` needs can be read straight off them — no CIGAR or MD
/// reconstruction, and no orientation handling, because BLAST reports both strands in
/// alignment orientation.
///
/// This exists because minimizer seeding loses divergent probes. On the MTEC bacterial
/// panel against *B. burgdorferi* B31, blastn-short finds 14 of the nifS/clpA probes
/// where the embedded aligner finds 6 — and the 6 are a strict subset. Probe density
/// drives per-locus capture efficiency, so dropping 57% of the probes at a locus
/// distorts the simulation in exactly the regime divergent panels care about.
pub fn probe_align(
    reference: &Path,
    probes: &Path,
    output_tsv: &Path,
    log_file: &Path,
    threads: usize,
) -> Result<()> {
    if fasta::count_sequences(reference)? == 0 || fasta::count_sequences(probes)? == 0 {
        File::create(output_tsv)
            .with_context(|| format!("Cannot create: {}", output_tsv.display()))?;
        File::create(log_file)
            .with_context(|| format!("Cannot create log: {}", log_file.display()))?;
        return Ok(());
    }

    let db_prefix = output_tsv.with_extension("blastdb");
    let log = File::create(log_file)
        .with_context(|| format!("Cannot create log: {}", log_file.display()))?;

    let makeblastdb_status = Command::new("makeblastdb")
        .arg("-in")
        .arg(reference)
        .args(["-dbtype", "nucl"])
        .arg("-out")
        .arg(&db_prefix)
        .stdout(log.try_clone().context("Cannot duplicate log handle")?)
        .stderr(log.try_clone().context("Cannot duplicate log handle")?)
        .status()
        .context("Failed to execute makeblastdb")?;
    if !makeblastdb_status.success() {
        cleanup_blast_db(&db_prefix);
        bail!(
            "makeblastdb failed (exit code {:?}); see {}",
            makeblastdb_status.code(),
            log_file.display()
        );
    }

    let blastn_status = Command::new("blastn")
        .args(["-task", "blastn-short"])
        .arg("-db")
        .arg(&db_prefix)
        .arg("-query")
        .arg(probes)
        .arg("-out")
        .arg(output_tsv)
        .args(["-outfmt", "6 qseqid sseqid sstart send sstrand qseq sseq"])
        .args(["-evalue", "1000"])
        .args(["-dust", "no"])
        // Without these, blastn-short (word size 7) returns an enormous number of tiny
        // spurious HSPs — 882k for a 228-probe panel against a single 0.9 Mb
        // chromosome. Such fragments cannot nucleate a duplex, but they would still
        // enter level-2 sampling and dilute the genuine hits.
        //
        // The thresholds match those used to derive the reference probe intervals for
        // this panel (`probes_merged.bed`): at least half the probe aligned, at no
        // worse than 80% identity. 80% over 120 bp is 24 mismatches, far more
        // permissive than minimizer seeding, so the divergent probes are retained.
        .args(["-perc_identity", "80"])
        .args(["-qcov_hsp_perc", "50"])
        .args(["-num_threads", &threads.to_string()])
        .stdout(std::process::Stdio::null())
        .stderr(log)
        .status()
        .context("Failed to execute blastn")?;

    cleanup_blast_db(&db_prefix);
    if !blastn_status.success() {
        bail!(
            "blastn failed (exit code {:?}); see {}",
            blastn_status.code(),
            log_file.display()
        );
    }
    Ok(())
}

pub fn xreact_align(
    reference: &Path,
    query: &Path,
    output_tsv: &Path,
    log_file: &Path,
    threads: usize,
) -> Result<()> {
    // makeblastdb/blastn error out on empty FASTA input (e.g. build-probes
    // filtered every probe away). minimap2 handles this gracefully (0 PAF
    // records), so match that behavior: skip the search and write an empty
    // results file.
    if fasta::count_sequences(reference)? == 0 || fasta::count_sequences(query)? == 0 {
        File::create(output_tsv)
            .with_context(|| format!("Cannot create: {}", output_tsv.display()))?;
        File::create(log_file)
            .with_context(|| format!("Cannot create log: {}", log_file.display()))?;
        return Ok(());
    }

    let db_prefix = output_tsv.with_extension("blastdb");
    let log = File::create(log_file)
        .with_context(|| format!("Cannot create log: {}", log_file.display()))?;

    let makeblastdb_status = Command::new("makeblastdb")
        .arg("-in")
        .arg(reference)
        .args(["-dbtype", "nucl"])
        .arg("-out")
        .arg(&db_prefix)
        .stdout(log.try_clone().context("Cannot duplicate log handle")?)
        .stderr(log.try_clone().context("Cannot duplicate log handle")?)
        .status()
        .context("Failed to execute makeblastdb")?;

    if !makeblastdb_status.success() {
        cleanup_blast_db(&db_prefix);
        bail!(
            "makeblastdb failed (exit code {:?}); see {}",
            makeblastdb_status.code(),
            log_file.display()
        );
    }

    let blastn_status = Command::new("blastn")
        .args(["-task", "blastn-short"])
        .arg("-db")
        .arg(&db_prefix)
        .arg("-query")
        .arg(query)
        .arg("-out")
        .arg(output_tsv)
        .args(["-outfmt", "6 qseqid sseqid qlen qstart qend nident length"])
        .args(["-evalue", "1000"])
        .args(["-dust", "no"])
        .args(["-num_threads", &threads.to_string()])
        .stdout(std::process::Stdio::null())
        .stderr(log)
        .status()
        .context("Failed to execute blastn")?;

    cleanup_blast_db(&db_prefix);

    if !blastn_status.success() {
        bail!(
            "blastn failed (exit code {:?}); see {}",
            blastn_status.code(),
            log_file.display()
        );
    }

    Ok(())
}

/// Remove all BLAST database files sharing `db_prefix` (extension varies by
/// BLAST+ version: .ndb/.nhr/.nin/.njs/.not/.nsq/.ntf/.nto or the older
/// .nhr/.nin/.nsq trio).
fn cleanup_blast_db(db_prefix: &Path) {
    let dir = db_prefix.parent().unwrap_or_else(|| Path::new("."));
    let prefix_name = match db_prefix.file_name().and_then(|n| n.to_str()) {
        Some(n) => n.to_string(),
        None => return,
    };

    let entries = match fs::read_dir(dir) {
        Ok(e) => e,
        Err(_) => return,
    };

    for entry in entries.flatten() {
        if let Some(name) = entry.file_name().to_str() {
            if name.starts_with(&prefix_name) && name != prefix_name {
                let _ = fs::remove_file(entry.path());
            }
        }
    }
}
