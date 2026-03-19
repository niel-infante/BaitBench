use anyhow::{Context, Result, bail};
use std::path::Path;

use crate::cli::ReportMode;
use crate::external::rscript;

pub struct ReportArgs<'a> {
    pub summary: &'a Path,
    pub detail: &'a Path,
    pub params: &'a Path,
    pub coverage: Option<&'a Path>,
    pub species_calls: Option<&'a Path>,
    pub run_name: &'a str,
    pub output: &'a Path,
    pub report: ReportMode,
}

pub fn execute(args: &ReportArgs) -> Result<()> {
    match args.report {
        ReportMode::None => {
            log::info!("Skipping report generation (--report none)");
            Ok(())
        }
        ReportMode::Full => execute_full(args),
        ReportMode::Rmd => execute_rmd(args),
    }
}

fn execute_full(args: &ReportArgs) -> Result<()> {
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

    let mut r_args: Vec<&str> = vec![
        "--summary", summary_str,
        "--detail", detail_str,
        "--params", params_str,
        "--run-name", args.run_name,
        "--output", output_str,
    ];

    let coverage_str;
    if let Some(cov_path) = args.coverage {
        if cov_path.exists() {
            let coverage_abs = std::fs::canonicalize(cov_path)
                .with_context(|| format!("Cannot find coverage file: {}", cov_path.display()))?;
            coverage_str = coverage_abs.to_str().unwrap_or("").to_string();
            r_args.push("--coverage");
            r_args.push(&coverage_str);
        }
    }

    let species_calls_str;
    if let Some(sc_path) = args.species_calls {
        if sc_path.exists() {
            let sc_abs = std::fs::canonicalize(sc_path)
                .with_context(|| format!("Cannot find species calls file: {}", sc_path.display()))?;
            species_calls_str = sc_abs.to_str().unwrap_or("").to_string();
            r_args.push("--species-calls");
            r_args.push(&species_calls_str);
        }
    }

    rscript::run_rscript(&report_script, &r_args)?;

    log::info!("Report generated: {}", args.output.display());
    Ok(())
}

fn execute_rmd(args: &ReportArgs) -> Result<()> {
    let r_dir = rscript::find_r_dir().ok_or_else(|| {
        anyhow::anyhow!(
            "Cannot find R scripts directory. Set BAITBENCH_R_DIR or ensure ./R/ exists."
        )
    })?;

    let rmd_template = r_dir.join("report.Rmd");
    if !rmd_template.exists() {
        bail!("RMarkdown template not found: {}", rmd_template.display());
    }

    let summary_abs = std::fs::canonicalize(args.summary)
        .with_context(|| format!("Cannot find summary file: {}", args.summary.display()))?;
    let detail_abs = std::fs::canonicalize(args.detail)
        .with_context(|| format!("Cannot find detail file: {}", args.detail.display()))?;
    let params_abs = std::fs::canonicalize(args.params)
        .with_context(|| format!("Cannot find params file: {}", args.params.display()))?;

    let coverage_val = if let Some(cov_path) = args.coverage {
        if cov_path.exists() {
            std::fs::canonicalize(cov_path)
                .with_context(|| format!("Cannot find coverage file: {}", cov_path.display()))?
                .to_str()
                .unwrap_or("")
                .to_string()
        } else {
            String::new()
        }
    } else {
        String::new()
    };

    let species_calls_val = if let Some(sc_path) = args.species_calls {
        if sc_path.exists() {
            std::fs::canonicalize(sc_path)
                .with_context(|| format!("Cannot find species calls file: {}", sc_path.display()))?
                .to_str()
                .unwrap_or("")
                .to_string()
        } else {
            String::new()
        }
    } else {
        String::new()
    };

    let params = vec![
        ("summary_file", summary_abs.to_str().unwrap_or("")),
        ("detail_file", detail_abs.to_str().unwrap_or("")),
        ("params_file", params_abs.to_str().unwrap_or("")),
        ("coverage_file", &coverage_val),
        ("species_calls_file", &species_calls_val),
        ("run_name", args.run_name),
    ];

    let template_content = std::fs::read_to_string(&rmd_template)
        .with_context(|| format!("Failed to read template: {}", rmd_template.display()))?;

    let output_content = substitute_rmd_params(&template_content, &params);

    let output_path = rmd_output_path(args.output);
    std::fs::write(&output_path, output_content)
        .with_context(|| format!("Failed to write RMarkdown file: {}", output_path.display()))?;

    log::info!("RMarkdown file written: {}", output_path.display());
    log::info!(
        "Edit and render with: Rscript -e 'rmarkdown::render(\"{}\")'",
        output_path.display()
    );
    Ok(())
}

/// Replace the `params:` block in an RMarkdown YAML front matter with actual values.
///
/// Handles both simple string params (`key: "value"` or `key: value`) and
/// R expression params (`key: !r c()`). All values are quoted as strings
/// except for integer/numeric values which are left unquoted.
pub fn substitute_rmd_params(template: &str, params: &[(&str, &str)]) -> String {
    let mut result = String::with_capacity(template.len());
    let mut in_yaml = false;
    let mut in_params = false;
    let mut yaml_count = 0;

    for line in template.lines() {
        let trimmed = line.trim();

        // Track YAML front matter boundaries
        if trimmed == "---" {
            yaml_count += 1;
            if yaml_count == 1 {
                in_yaml = true;
            } else {
                in_yaml = false;
                in_params = false;
            }
            result.push_str(line);
            result.push('\n');
            continue;
        }

        if in_yaml && trimmed == "params:" {
            in_params = true;
            result.push_str(line);
            result.push('\n');
            continue;
        }

        // Inside the params block: lines starting with spaces followed by a key
        if in_params {
            // Check if this line is still a param (indented key: value)
            if line.starts_with(' ') || line.starts_with('\t') {
                // Try to match a param key
                let key_part = trimmed.split(':').next().unwrap_or("");
                if let Some((_, value)) = params.iter().find(|(k, _)| *k == key_part) {
                    // Get the indentation from the original line
                    let indent: String = line.chars().take_while(|c| c.is_whitespace()).collect();
                    // Check if value is a pure integer or float
                    if value.parse::<i64>().is_ok() || value.parse::<f64>().is_ok() {
                        result.push_str(&format!("{}{}: {}\n", indent, key_part, value));
                    } else {
                        result.push_str(&format!("{}{}: \"{}\"\n", indent, key_part, value));
                    }
                    continue;
                }
            } else {
                // No longer indented — we've left the params block
                in_params = false;
            }
        }

        result.push_str(line);
        result.push('\n');
    }

    result
}

/// Convert an output path from .html to .Rmd (or append .Rmd if no .html extension).
pub fn rmd_output_path(html_path: &Path) -> std::path::PathBuf {
    let mut path = html_path.to_path_buf();
    if path.extension().map_or(false, |ext| ext == "html") {
        path.set_extension("Rmd");
    } else {
        let mut name = path
            .file_name()
            .unwrap_or_default()
            .to_string_lossy()
            .to_string();
        name.push_str(".Rmd");
        path.set_file_name(name);
    }
    path
}
