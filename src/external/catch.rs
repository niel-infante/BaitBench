// Wrapper for the external CATCH probe design tool (Metsky et al. 2019).
// Requires the `catch` conda package (bioconda): https://github.com/broadinstitute/catch
//
// The main command is `design_probes.py`, installed to PATH by the catch package.

use anyhow::{Context, Result, bail};
use std::fs::{self, File};
use std::io::{BufRead, BufReader};
use std::path::Path;
use std::process::Command;

/// Check that design_probes.py (external CATCH) is available on PATH.
pub fn check_available() -> Result<()> {
    let status = Command::new("design_probes.py")
        .arg("--help")
        .stdout(std::process::Stdio::null())
        .stderr(std::process::Stdio::null())
        .status();

    match status {
        Ok(_) => Ok(()),
        _ => bail!(
            "design_probes.py not found on PATH. \
             Install the external CATCH tool: conda install -c bioconda catch"
        ),
    }
}

/// Run external CATCH to design probes from an input FASTA.
///
/// Calls: `design_probes.py <input> -pl <probe_len> -s <stride> -m <mismatches> -o <output>`
///
/// Returns the number of probes written to the output file.
pub fn design_probes(
    input: &Path,
    output: &Path,
    probe_len: usize,
    stride: usize,
    mismatches: usize,
) -> Result<usize> {
    let log_path = output.with_extension("catch.log");
    let log = File::create(&log_path)
        .with_context(|| format!("Cannot create log: {}", log_path.display()))?;

    let status = Command::new("design_probes.py")
        .arg(input)
        .arg("-pl")
        .arg(probe_len.to_string())
        .arg("-s")
        .arg(stride.to_string())
        .arg("-m")
        .arg(mismatches.to_string())
        .arg("-o")
        .arg(output)
        .stdout(log.try_clone()?)
        .stderr(log)
        .status()
        .context("Failed to execute design_probes.py")?;

    if !status.success() {
        bail!("design_probes.py failed (exit code {:?})", status.code());
    }

    let count = count_fasta_sequences(output)?;
    log::info!("  External CATCH: {} probes written", count);
    let _ = fs::remove_file(&log_path);
    Ok(count)
}

fn count_fasta_sequences(path: &Path) -> Result<usize> {
    let file = File::open(path).with_context(|| format!("Cannot open: {}", path.display()))?;
    let reader = BufReader::new(file);
    let mut count = 0;
    for line in reader.lines() {
        if line?.starts_with('>') {
            count += 1;
        }
    }
    Ok(count)
}
