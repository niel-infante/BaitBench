use anyhow::Result;
use std::path::Path;

/// No-op: cd-hit-est is now implemented natively; no external binary required.
pub fn check_available() -> Result<()> {
    Ok(())
}

/// Cluster nucleotide sequences by identity.
///
/// Delegates to the native Rust cd-hit-est implementation in `crate::cdhit_est`.
/// The `threads` argument is accepted for API compatibility but ignored
/// (the native implementation is single-threaded).
pub fn cluster(
    input: &Path,
    output: &Path,
    threshold: f64,
    _threads: usize,
    log_file: &Path,
) -> Result<()> {
    crate::cdhit_est::cluster(input, output, threshold, log_file)
}
