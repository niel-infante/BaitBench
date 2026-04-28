use anyhow::{Context, Result, bail};
use std::fs::File;
use std::path::Path;
use std::process::Command;

/// Check that cd-hit-est is available on PATH.
pub fn check_available() -> Result<()> {
    if is_available() {
        Ok(())
    } else {
        bail!("cd-hit-est not found on PATH. Please install cd-hit.")
    }
}

/// Return true if cd-hit-est is available on PATH.
pub fn is_available() -> bool {
    Command::new("cd-hit-est")
        .arg("-h")
        .stdout(std::process::Stdio::null())
        .stderr(std::process::Stdio::null())
        .status()
        .is_ok()
}

/// Run cd-hit-est to cluster sequences by identity.
///
/// `cd-hit-est -i <input> -o <output> -c <threshold> -T <threads>`
pub fn cluster(
    input: &Path,
    output: &Path,
    threshold: f64,
    threads: usize,
    log_file: &Path,
) -> Result<()> {
    let log = File::create(log_file)
        .with_context(|| format!("Cannot create log: {}", log_file.display()))?;

    let status = Command::new("cd-hit-est")
        .arg("-i")
        .arg(input)
        .arg("-o")
        .arg(output)
        .arg("-c")
        .arg(format!("{:.2}", threshold))
        .arg("-T")
        .arg(threads.to_string())
        .stdout(log.try_clone()?)
        .stderr(log)
        .status()
        .context("Failed to execute cd-hit-est")?;

    if !status.success() {
        bail!(
            "cd-hit-est failed (exit code {:?})",
            status.code()
        );
    }

    Ok(())
}
