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
) -> Result<()> {
    let out = File::create(output_sam)
        .with_context(|| format!("Cannot create SAM: {}", output_sam.display()))?;

    let status = Command::new("minimap2")
        .args(["-ax", preset])
        .arg("--secondary=no")
        .arg(reference)
        .arg(reads)
        .stdout(out)
        .stderr(std::process::Stdio::null())
        .status()
        .context("Failed to execute minimap2")?;

    if !status.success() {
        bail!("minimap2 mapping failed (exit code {:?})", status.code());
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
