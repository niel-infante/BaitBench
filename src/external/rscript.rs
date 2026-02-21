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
/// 2. Walk up from the resolved binary location looking for R/report.R
///    (resolves symlinks via canonicalize, works on both macOS and Linux)
/// 3. ./R/ in current working directory
pub fn find_r_dir() -> Option<std::path::PathBuf> {
    // Check environment variable
    if let Ok(dir) = std::env::var("BAITBENCH_R_DIR") {
        let path = std::path::PathBuf::from(dir);
        if path.is_dir() {
            return Some(path);
        }
    }

    // Resolve symlinks (current_exe doesn't resolve on macOS) then walk up
    if let Ok(exe) = std::env::current_exe() {
        let resolved = std::fs::canonicalize(&exe).unwrap_or(exe);
        let mut dir = resolved.parent();
        for _ in 0..4 {
            match dir {
                Some(d) => {
                    let r_dir = d.join("R");
                    if r_dir.join("report.R").is_file() {
                        return Some(r_dir);
                    }
                    let share_dir = d.join("share/baitbench/R");
                    if share_dir.join("report.R").is_file() {
                        return Some(share_dir);
                    }
                    dir = d.parent();
                }
                None => break,
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
