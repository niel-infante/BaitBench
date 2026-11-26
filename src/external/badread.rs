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

/// Map a user-facing profile name to badread's
/// (error_model, qscore_model, len_mean, len_sd, identity).
///
/// | Profile        | Chemistry                | Error model  | Identity      |
/// |----------------|--------------------------|--------------|---------------|
/// | `ont`          | ONT R10.4.1 / Kit14      | nanopore2023 | badread deflt |
/// | `ont-2020`     | ONT R9.4.1               | nanopore2020 | badread deflt |
/// | `pacbio`       | PacBio CLR               | pacbio2016   | badread deflt |
/// | `pacbio-hifi`  | PacBio HiFi / CCS        | pacbio2021   | 99.5,100,1.5  |
///
/// **Identity is a separate axis from the error model.** The error model shapes *what
/// kind* of errors appear; `--identity` sets *how many*. Badread's default is
/// `95,99,2.5`, so before `pacbio-hifi` existed every PacBio simulation ran at ~95%
/// mean identity — CLR-grade — no matter which error model was chosen. That is why
/// `--read-simulator perfect` had been the closer match to HiFi data, and it was a
/// limitation of this wrapper rather than of badread.
///
/// For reference, measured on the D6331 mock community's real CCS reads (342k
/// alignments, identity from NM over aligned length): median 99.83%, mean 98.53%,
/// sd 2.82, p5 91.53%. The distribution is strongly left-skewed — most reads near
/// perfect with a thin poor tail — which is what badread's beta form (mean,max,stdev)
/// reproduces. Pass `--read-identity` to match a specific run rather than the generic
/// HiFi default.
pub fn profile_params(
    profile: &str,
) -> Result<(&'static str, &'static str, u32, u32, Option<&'static str>)> {
    match profile {
        "ont" => Ok(("nanopore2023", "nanopore2023", 9000, 7000, None)),
        "ont-2020" => Ok(("nanopore2020", "nanopore2020", 9000, 7000, None)),
        "pacbio" => Ok(("pacbio2016", "pacbio2016", 15000, 13000, None)),
        // Q20+ is the HiFi specification, with most reads well above it.
        "pacbio-hifi" => Ok((
            "pacbio2021",
            "pacbio2021",
            15000,
            13000,
            Some("99.5,100,1.5"),
        )),
        other => bail!(
            "Unknown badread profile '{}'. Valid values: ont, ont-2020, pacbio, pacbio-hifi",
            other
        ),
    }
}

/// Run badread to simulate long reads from a multi-FASTA of captured fragments.
///
/// `len_mean` / `len_sd`: override the profile's default read length distribution.
/// `identity`: override the profile's identity distribution, as badread's
/// `mean,max,stdev` (beta) or `mean,stdev` (normal qscore). `None` uses the profile's
/// default, and a profile with no default leaves badread's own (`95,99,2.5`) in place.
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
    identity: Option<&str>,
    glitches: Option<&str>,
    junk_reads: Option<f64>,
    random_reads: Option<f64>,
    chimeras: Option<f64>,
) -> Result<()> {
    let (error_model, qscore_model, profile_mean, profile_sd, profile_identity) =
        profile_params(profile)?;
    let final_mean = len_mean.unwrap_or(profile_mean as usize);
    let final_sd = len_sd.unwrap_or(profile_sd as usize);
    let final_identity = identity.or(profile_identity);

    let log = File::create(log_file)
        .with_context(|| format!("Cannot create log: {}", log_file.display()))?;
    let output = File::create(output_fastq)
        .with_context(|| format!("Cannot create output: {}", output_fastq.display()))?;

    let mut cmd = Command::new("badread");
    cmd.arg("simulate")
        .arg("--reference")
        .arg(input)
        .arg("--quantity")
        .arg(format!("{:.1}x", coverage_depth))
        .arg("--error_model")
        .arg(error_model)
        .arg("--qscore_model")
        .arg(qscore_model)
        .arg("--length")
        .arg(format!("{},{}", final_mean, final_sd));

    if let Some(i) = final_identity {
        cmd.arg("--identity").arg(i);
    }
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
