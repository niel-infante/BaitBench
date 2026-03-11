use anyhow::{Context, Result, bail};
use std::fs::File;
use std::path::Path;
use std::process::Command;

/// Check that minimap2 is available on PATH.
pub fn check_available() -> Result<()> {
    let status = Command::new("minimap2")
        .arg("--version")
        .stdout(std::process::Stdio::null())
        .stderr(std::process::Stdio::null())
        .status();

    match status {
        Ok(s) if s.success() => Ok(()),
        _ => bail!("minimap2 not found on PATH. Please install minimap2."),
    }
}

/// Run minimap2 for probe capture (PAF output with CIGAR).
///
/// `minimap2 -x sr -c --cs -A 4 -B 2 -O 12,32 --secondary=yes <probes> <reads> > <output>`
pub fn capture_align(
    probes: &Path,
    reads: &Path,
    output_paf: &Path,
    log_file: &Path,
) -> Result<()> {
    let log = File::create(log_file)
        .with_context(|| format!("Cannot create log: {}", log_file.display()))?;
    let out = File::create(output_paf)
        .with_context(|| format!("Cannot create PAF: {}", output_paf.display()))?;

    let status = Command::new("minimap2")
        .args(["-x", "sr", "-c", "--cs"])
        .args(["-A", "4", "-B", "2", "-O", "12,32"])
        .arg("--secondary=yes")
        .arg(probes)
        .arg(reads)
        .stdout(out)
        .stderr(log)
        .status()
        .context("Failed to execute minimap2")?;

    if !status.success() {
        bail!("minimap2 capture alignment failed (exit code {:?})", status.code());
    }

    Ok(())
}

/// Run minimap2 for read mapping (SAM output).
///
/// `minimap2 -ax <preset> --secondary=no <reference> <reads> > <output>`
pub fn map_reads(
    preset: &str,
    reference: &Path,
    reads: &Path,
    output_sam: &Path,
    log_file: &Path,
) -> Result<()> {
    let out = File::create(output_sam)
        .with_context(|| format!("Cannot create SAM: {}", output_sam.display()))?;
    let log = File::create(log_file)
        .with_context(|| format!("Cannot create log: {}", log_file.display()))?;

    let status = Command::new("minimap2")
        .args(["-ax", preset])
        .arg("--secondary=no")
        .arg(reference)
        .arg(reads)
        .stdout(out)
        .stderr(log)
        .status()
        .context("Failed to execute minimap2")?;

    if !status.success() {
        bail!("minimap2 mapping failed (exit code {:?})", status.code());
    }

    Ok(())
}

/// Run minimap2 for probe-to-target mapping (SAM output with secondary alignments).
///
/// `minimap2 -ax <preset> --secondary=yes -N 1000 <targets> <probes> > <output>`
///
/// Secondary alignments are kept because a single probe can legitimately
/// tile conserved regions across multiple target sequences.
pub fn probe_align(
    preset: &str,
    targets: &Path,
    probes: &Path,
    output_sam: &Path,
    log_file: &Path,
) -> Result<()> {
    let out = File::create(output_sam)
        .with_context(|| format!("Cannot create SAM: {}", output_sam.display()))?;
    let log = File::create(log_file)
        .with_context(|| format!("Cannot create log: {}", log_file.display()))?;

    let status = Command::new("minimap2")
        .args(["-ax", preset])
        .arg("--secondary=yes")
        .args(["-N", "1000"])
        .arg(targets)
        .arg(probes)
        .stdout(out)
        .stderr(log)
        .status()
        .context("Failed to execute minimap2")?;

    if !status.success() {
        bail!(
            "minimap2 probe alignment failed (exit code {:?})",
            status.code()
        );
    }

    Ok(())
}

/// Run minimap2 for cross-reactivity analysis (PAF output with secondary alignments).
///
/// `minimap2 -x <preset> -c --secondary=yes -N 10000 -p 0.5 <reference> <probes> > <output>`
///
/// Uses PAF output to get lightweight tabular format with matching base counts.
/// Secondary alignments are kept to capture all cross-reactive hits.
/// -p 0.5 lowers the secondary score threshold to catch weaker homology.
pub fn xreact_align(
    preset: &str,
    reference: &Path,
    probes: &Path,
    output_paf: &Path,
    log_file: &Path,
) -> Result<()> {
    let log = File::create(log_file)
        .with_context(|| format!("Cannot create log: {}", log_file.display()))?;
    let out = File::create(output_paf)
        .with_context(|| format!("Cannot create PAF: {}", output_paf.display()))?;

    let status = Command::new("minimap2")
        .args(["-x", preset, "-c"])
        .arg("--secondary=yes")
        .args(["-N", "10000"])
        .args(["-p", "0.5"])
        .arg(reference)
        .arg(probes)
        .stdout(out)
        .stderr(log)
        .status()
        .context("Failed to execute minimap2")?;

    if !status.success() {
        bail!(
            "minimap2 xreact alignment failed (exit code {:?})",
            status.code()
        );
    }

    Ok(())
}

/// Run minimap2 for host filtering (SAM output).
///
/// `minimap2 -ax <preset> <host> <reads> > <output>`
pub fn host_align(
    preset: &str,
    host: &Path,
    reads: &Path,
    output_sam: &Path,
    log_file: &Path,
) -> Result<()> {
    let out = File::create(output_sam)
        .with_context(|| format!("Cannot create SAM: {}", output_sam.display()))?;
    let log = File::create(log_file)
        .with_context(|| format!("Cannot create log: {}", log_file.display()))?;

    let status = Command::new("minimap2")
        .args(["-ax", preset])
        .arg(host)
        .arg(reads)
        .stdout(out)
        .stderr(log)
        .status()
        .context("Failed to execute minimap2")?;

    if !status.success() {
        bail!("minimap2 host alignment failed (exit code {:?})", status.code());
    }

    Ok(())
}
