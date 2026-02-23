use anyhow::{Context, Result, bail};
use std::path::Path;

use crate::external::rscript;

pub struct ReportArgs<'a> {
    pub summary: &'a Path,
    pub detail: &'a Path,
    pub params: &'a Path,
    pub run_name: &'a str,
    pub output: &'a Path,
}

pub fn execute(args: &ReportArgs) -> Result<()> {
    if !rscript::check_available() {
        bail!(
            "Rscript not found on PATH. Install R to generate reports. \
             All other outputs (TSV, JSON) are still available."
        );
    }

    let r_dir = rscript::find_r_dir().ok_or_else(|| {
        anyhow::anyhow!(
            "Cannot find R scripts directory. Set BAITBENCH_R_DIR or ensure ./R/ exists."
        )
    })?;

    let report_script = r_dir.join("report.R");
    if !report_script.exists() {
        bail!("R report script not found: {}", report_script.display());
    }

    log::info!("Generating HTML report...");
    let summary_abs = std::fs::canonicalize(args.summary)
        .with_context(|| format!("Cannot find summary file: {}", args.summary.display()))?;
    let detail_abs = std::fs::canonicalize(args.detail)
        .with_context(|| format!("Cannot find detail file: {}", args.detail.display()))?;
    let params_abs = std::fs::canonicalize(args.params)
        .with_context(|| format!("Cannot find params file: {}", args.params.display()))?;
    let output_abs = if args.output.is_absolute() {
        args.output.to_path_buf()
    } else {
        std::env::current_dir()?.join(args.output)
    };
    let summary_str = summary_abs.to_str().unwrap_or("");
    let detail_str = detail_abs.to_str().unwrap_or("");
    let params_str = params_abs.to_str().unwrap_or("");
    let output_str = output_abs.to_str().unwrap_or("");

    rscript::run_rscript(
        &report_script,
        &[
            "--summary", summary_str,
            "--detail", detail_str,
            "--params", params_str,
            "--run-name", args.run_name,
            "--output", output_str,
        ],
    )?;

    log::info!("Report generated: {}", args.output.display());
    Ok(())
}
