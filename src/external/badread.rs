use anyhow::{bail, Context, Result};
use std::fs::File;
use std::path::Path;
use std::process::Command;

/// Check that badread is available on PATH.
pub fn check_available() -> Result<()> {
    let status = Command::new("badread")
        .arg("--version")
        .stdout(std::process::Stdio::null())
        .stderr(std::process::Stdio::null())
        .status();
    match status {
        Ok(_) => Ok(()),
        Err(_) => bail!(
            "badread not found on PATH. Install it with: conda install -c conda-forge badread"
        ),
    }
}

/// Map a user-facing profile name to badread's (error_model, qscore_model, len_mean, len_sd).
///
/// | Profile     | Chemistry                  | Error model    |
/// |-------------|---------------------------|----------------|
/// | `ont`       | ONT R10.4.1 / Kit14        | nanopore2023   |
/// | `ont-2020`  | ONT R9.4.1                 | nanopore2020   |
/// | `pacbio`    | PacBio CLR                 | pacbio2016     |
pub fn profile_params(profile: &str) -> Result<(&'static str, &'static str, u32, u32)> {
    match profile {
        "ont"      => Ok(("nanopore2023", "nanopore2023", 9000,  7000)),
        "ont-2020" => Ok(("nanopore2020", "nanopore2020", 9000,  7000)),
        "pacbio"   => Ok(("pacbio2016",   "pacbio2016",   15000, 13000)),
        other => bail!(
            "Unknown badread profile '{}'. Valid values: ont, ont-2020, pacbio",
            other
        ),
    }
}

/// Run badread to simulate long reads from a multi-FASTA of captured fragments.
///
/// Command:
/// ```text
/// badread simulate --reference <input> --quantity <depth>x
///                  --error_model <model> --qscore_model <model>
///                  --length <mean>,<sd> [--seed <seed>]
/// ```
///
/// Reads are written in FASTQ format to `output_fastq` (uncompressed).
/// Read names embed the source sequence ID: `@{uuid} {ref_name},{start}-{end} ...`
pub fn run_simulation(
    input: &Path,
    output_fastq: &Path,
    profile: &str,
    coverage_depth: f64,
    seed: Option<u64>,
    log_file: &Path,
) -> Result<()> {
    let (error_model, qscore_model, len_mean, len_sd) = profile_params(profile)?;
    let log = File::create(log_file)
        .with_context(|| format!("Cannot create log: {}", log_file.display()))?;
    let output = File::create(output_fastq)
        .with_context(|| format!("Cannot create output: {}", output_fastq.display()))?;

    let mut cmd = Command::new("badread");
    cmd.arg("simulate")
        .arg("--reference").arg(input)
        .arg("--quantity").arg(format!("{:.1}x", coverage_depth))
        .arg("--error_model").arg(error_model)
        .arg("--qscore_model").arg(qscore_model)
        .arg("--length").arg(format!("{},{}", len_mean, len_sd));

    if let Some(s) = seed {
        cmd.arg("--seed").arg(s.to_string());
    }

    cmd.stdout(output).stderr(log);

    let status = cmd.status().context("Failed to execute badread")?;
    if !status.success() {
        bail!("badread exited with status: {}", status);
    }
    Ok(())
}
