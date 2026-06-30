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
/// `len_mean` / `len_sd`: override the profile's default read length distribution.
/// Realism controls (`glitches`, `junk_reads`, `random_reads`, `chimeras`) are only
/// passed to badread when `Some`; when `None`, badread's own defaults apply.
pub fn run_simulation(
    input: &Path,
    output_fastq: &Path,
    profile: &str,
    coverage_depth: f64,
    seed: Option<u64>,
    log_file: &Path,
    len_mean: Option<usize>,
    len_sd: Option<usize>,
    glitches: Option<&str>,
    junk_reads: Option<f64>,
    random_reads: Option<f64>,
    chimeras: Option<f64>,
) -> Result<()> {
    let (error_model, qscore_model, profile_mean, profile_sd) = profile_params(profile)?;
    let final_mean = len_mean.unwrap_or(profile_mean as usize);
    let final_sd   = len_sd.unwrap_or(profile_sd as usize);

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
        .arg("--length").arg(format!("{},{}", final_mean, final_sd));

    if let Some(s) = seed {
        cmd.arg("--seed").arg(s.to_string());
    }
    if let Some(g) = glitches {
        cmd.arg("--glitches").arg(g);
    }
    if let Some(j) = junk_reads {
        cmd.arg("--junk_reads").arg(j.to_string());
    }
    if let Some(r) = random_reads {
        cmd.arg("--random_reads").arg(r.to_string());
    }
    if let Some(c) = chimeras {
        cmd.arg("--chimeras").arg(c.to_string());
    }

    cmd.stdout(output).stderr(log);

    let status = cmd.status().context("Failed to execute badread")?;
    if !status.success() {
        bail!("badread exited with status: {}", status);
    }
    Ok(())
}
