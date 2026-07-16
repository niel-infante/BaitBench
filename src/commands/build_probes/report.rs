use anyhow::{Context, Result};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use crate::commands::report::substitute_rmd_params;
use crate::external::rscript;
use crate::io_utils::abs_path_str;

use super::BuildProbesArgs;
use super::StepStats;

pub fn write_stats_tsv(path: &Path, stats: &[StepStats]) -> Result<()> {
    let file = File::create(path).with_context(|| format!("Cannot create: {}", path.display()))?;
    let mut w = BufWriter::new(file);

    writeln!(w, "step\tsequences\tbases")?;
    for s in stats {
        writeln!(w, "{}\t{}\t{}", s.step, s.sequences, s.bases)?;
    }

    w.flush()?;
    Ok(())
}

pub fn write_run_params(path: &Path, args: &BuildProbesArgs, cdhit_available: bool) -> Result<()> {
    let file = File::create(path).with_context(|| format!("Cannot create: {}", path.display()))?;
    let mut w = BufWriter::new(file);

    writeln!(w, "parameter\tflag\tvalue")?;
    writeln!(w, "targets\t--targets\t{}", args.targets.display())?;
    writeln!(w, "method\t--method\t{:?}", args.method)?;
    writeln!(w, "probe_length\t--probe-length\t{}", args.probe_length)?;
    writeln!(w, "step\t--step\t{}", args.step)?;
    writeln!(w, "catch_probe_stride\t--catch-stride\t{}", args.catch_probe_stride)?;
    writeln!(w, "catch_mismatches\t--catch-mismatches\t{}", args.catch_mismatches)?;
    writeln!(w, "catch_extension\t--catch-extension\t{}", args.catch_extension)?;
    writeln!(w, "catch_coverage\t--catch-coverage\t{:.2}", args.catch_coverage)?;
    writeln!(w, "catch_minhash_threshold\t--catch-minhash-threshold\t{:.2}", args.catch_minhash_threshold)?;
    writeln!(w, "min_gc\t--min-gc\t{:.2}", args.min_gc)?;
    writeln!(w, "max_gc\t--max-gc\t{:.2}", args.max_gc)?;
    writeln!(w, "max_n_frac\t--max-n-frac\t{:.2}", args.max_n_frac)?;
    writeln!(w, "dust_threshold\t--dust-threshold\t{}", args.dust_threshold)?;
    writeln!(w, "dust_window\t--dust-window\t{}", args.dust_window)?;
    writeln!(w, "max_masked_frac\t--max-masked-frac\t{:.2}", args.max_masked_frac)?;
    writeln!(w, "collapse_threshold\t--collapse-threshold\t{:.2}", args.collapse_threshold)?;
    writeln!(w, "dedup_threshold\t--dedup-threshold\t{:.2}", args.dedup_threshold)?;
    writeln!(w, "threads\t--threads\t{}", args.threads)?;
    writeln!(w, "outdir\t--outdir\t{}", args.outdir.display())?;
    writeln!(w, "cdhit_available\t--n/a\t{}", cdhit_available)?;
    writeln!(w, "pt_step\t--pt-step\t{}", args.pt_step)?;
    writeln!(w, "pt_identity\t--pt-identity\t{:.2}", args.pt_identity)?;
    writeln!(w, "pt_coverage\t--pt-coverage\t{:.2}", args.pt_coverage)?;
    writeln!(w, "pt_batch_size\t--pt-batch-size\t{}", args.pt_batch_size)?;
    writeln!(w, "pt_max_panel_size\t--pt-max-panel-size\t{}", args.pt_max_panel_size.map(|v| v.to_string()).unwrap_or_default())?;
    writeln!(w, "pt_min_depth\t--pt-min-depth\t{}", args.pt_min_depth)?;
    writeln!(w, "pt_max_iterations\t--pt-max-iterations\t{}", args.pt_max_iterations)?;
    writeln!(w, "pt_min_coverage_gain\t--pt-min-coverage-gain\t{:.4}", args.pt_min_coverage_gain)?;

    w.flush()?;
    Ok(())
}

pub fn generate_report(stats_path: &Path, params_path: &Path, output_path: &Path) -> Result<()> {
    let r_dir =
        rscript::find_r_dir().ok_or_else(|| anyhow::anyhow!("Cannot find R scripts directory."))?;

    let script = r_dir.join("build_probes.R");
    if !script.exists() {
        anyhow::bail!("Build probes R script not found: {}", script.display());
    }

    rscript::run_rscript(
        &script,
        &[
            "--stats",
            &abs_path_str(stats_path)?,
            "--params",
            &abs_path_str(params_path)?,
            "--output",
            &abs_path_str(output_path)?,
        ],
    )
}

pub fn write_rmd(stats_path: &Path, params_path: &Path, rmd_path: &Path) -> Result<()> {
    let r_dir =
        rscript::find_r_dir().ok_or_else(|| anyhow::anyhow!("Cannot find R scripts directory."))?;

    let template_path = r_dir.join("build_probes.Rmd");
    if !template_path.exists() {
        anyhow::bail!(
            "Build probes Rmd template not found: {}",
            template_path.display()
        );
    }

    let template = std::fs::read_to_string(&template_path)
        .with_context(|| format!("Failed to read template: {}", template_path.display()))?;

    let stats_abs = abs_path_str(stats_path)?;
    let params_abs = abs_path_str(params_path)?;

    let substituted = substitute_rmd_params(
        &template,
        &[
            ("stats_file", stats_abs.as_str()),
            ("params_file", params_abs.as_str()),
        ],
    );

    std::fs::write(rmd_path, substituted)?;
    Ok(())
}
