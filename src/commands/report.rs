use anyhow::{Result, bail};
use std::path::Path;

use crate::external::rscript;

pub struct ReportArgs<'a> {
    pub summary: &'a Path,
    pub detail: &'a Path,
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
    let summary_str = args.summary.to_str().unwrap_or("");
    let detail_str = args.detail.to_str().unwrap_or("");
    let output_str = args.output.to_str().unwrap_or("");

    rscript::run_rscript(
        &report_script,
        &[
            "--summary", summary_str,
            "--detail", detail_str,
            "--run-name", args.run_name,
            "--output", output_str,
        ],
    )?;

    log::info!("Report generated: {}", args.output.display());
    Ok(())
}
