use anyhow::{Context, Result, bail};
use std::fs::File;
use std::path::Path;
use std::process::Command;

/// Check that art_modern is available on PATH.
pub fn check_available() -> Result<()> {
    let status = Command::new("art_modern")
        .arg("--help")
        .stdout(std::process::Stdio::null())
        .stderr(std::process::Stdio::null())
        .status();

    match status {
        Ok(_) => Ok(()),
        Err(_) => bail!(
            "art_modern not found on PATH. Install it with: conda install -c bioconda art_modern"
        ),
    }
}

/// Run art_modern to simulate Illumina reads from captured fragments.
///
/// Single-end:
/// ```
/// art_modern --i-file <input> --o-fastq <out_r1> --o-sam <out_sam>
///            --read_len <n> --builtin_qual_file <profile>
///            --i-fcov <depth> --lc se
/// ```
///
/// Paired-end:
/// ```
/// art_modern --i-file <input> --o-fastq <out_r1> --o-fastq2 <out_r2> --o-sam <out_sam>
///            --read_len <n> --builtin_qual_file <profile>
///            --i-fcov <depth> --lc pe
///            --pe_frag_dist_mean <pe_mean> --pe_frag_dist_std_dev <pe_sd>
/// ```
///
/// Returns paths to the generated FASTQ file(s): `(r1_path, r2_path_if_pe)`.
#[allow(clippy::too_many_arguments)]
pub fn run_simulation(
    input: &Path,
    output_r1: &Path,
    output_r2: Option<&Path>,
    output_sam: &Path,
    log_file: &Path,
    read_len: usize,
    profile: &str,
    coverage_depth: f64,
    paired_end: bool,
    pe_frag_len_mean: usize,
    pe_frag_len_sd: usize,
) -> Result<()> {
    let log = File::create(log_file)
        .with_context(|| format!("Cannot create log: {}", log_file.display()))?;

    let mut cmd = Command::new("art_modern");
    cmd.arg("--i-file").arg(input)
        .arg("--o-fastq").arg(output_r1)
        .arg("--o-sam").arg(output_sam)
        .arg("--read_len").arg(read_len.to_string())
        .arg("--builtin_qual_file").arg(profile)
        .arg("--i-fcov").arg(coverage_depth.to_string());

    if paired_end {
        let r2 = output_r2.ok_or_else(|| anyhow::anyhow!("output_r2 required for paired-end"))?;
        cmd.arg("--o-fastq2").arg(r2)
            .arg("--lc").arg("pe")
            .arg("--pe_frag_dist_mean").arg(pe_frag_len_mean.to_string())
            .arg("--pe_frag_dist_std_dev").arg(pe_frag_len_sd.to_string());
    } else {
        cmd.arg("--lc").arg("se");
    }

    cmd.stderr(log);

    let status = cmd.status().context("Failed to execute art_modern")?;

    if !status.success() {
        bail!("art_modern failed (exit code {:?}). Check {}", status.code(), log_file.display());
    }

    Ok(())
}
