use anyhow::{Context, Result, bail};
use std::path::Path;
use std::process::Command;

/// Check that Rscript is available on PATH.
pub fn check_available() -> bool {
    Command::new("Rscript")
        .arg("--version")
        .stdout(std::process::Stdio::null())
        .stderr(std::process::Stdio::null())
        .status()
        .is_ok_and(|s| s.success())
}

/// Find the R scripts directory. Search order:
/// 1. $BAITBENCH_R_DIR environment variable
/// 2. ../share/baitbench/R/ relative to binary
/// 3. ./R/ in current working directory
pub fn find_r_dir() -> Option<std::path::PathBuf> {
    // Check environment variable
    if let Ok(dir) = std::env::var("BAITBENCH_R_DIR") {
        let path = std::path::PathBuf::from(dir);
        if path.is_dir() {
            return Some(path);
        }
    }

    // Check relative to binary
    if let Ok(exe) = std::env::current_exe() {
        let share_dir = exe
            .parent()
            .and_then(|p| p.parent())
            .map(|p| p.join("share/baitbench/R"));
        if let Some(dir) = share_dir {
            if dir.is_dir() {
                return Some(dir);
            }
        }
    }

    // Check current working directory
    let cwd_r = std::path::PathBuf::from("R");
    if cwd_r.is_dir() {
        return Some(cwd_r);
    }

    None
}

/// Run an R script with the given arguments.
pub fn run_rscript(script: &Path, args: &[&str]) -> Result<()> {
    let status = Command::new("Rscript")
        .arg(script)
        .args(args)
        .status()
        .context("Failed to execute Rscript")?;

    if !status.success() {
        bail!("Rscript failed (exit code {:?})", status.code());
    }

    Ok(())
}
