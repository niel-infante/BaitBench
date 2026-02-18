use anyhow::{Context, Result, bail};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

/// Parse a weights file (tab-separated: id<tab>weight, ignoring # comments).
pub fn parse_weights(path: &Path) -> Result<HashMap<String, f64>> {
    let file = File::open(path).with_context(|| format!("Cannot open weights: {}", path.display()))?;
    let reader = BufReader::new(file);

    let mut weights = HashMap::new();

    for line in reader.lines() {
        let line = line?;
        let line = line.trim().to_string();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() >= 2 {
            let id = parts[0].to_string();
            let weight: f64 = parts[1].parse()
                .with_context(|| format!("Invalid weight value: {}", parts[1]))?;
            weights.insert(id, weight);
        }
    }

    Ok(weights)
}

/// Generate weights for targets and distractors.
///
/// Targets get weight 1.0 each.
/// Distractors are weighted to produce distractor_fraction of total reads.
/// Formula: distractor_weight = (fraction * n_targets) / (n_distractors * (1 - fraction))
pub fn generate_weights(
    target_ids: &[String],
    distractor_ids: &[String],
    distractor_fraction: f64,
    output: &Path,
) -> Result<(f64, f64)> {
    let n_targets = target_ids.len();
    let n_distractors = distractor_ids.len();

    if n_targets == 0 {
        bail!("No target sequences found!");
    }

    let target_weight = 1.0;

    let distractor_weight = if n_distractors > 0 && distractor_fraction > 0.0 {
        if distractor_fraction >= 1.0 {
            bail!("distractor_fraction must be less than 1.0");
        }
        (distractor_fraction * n_targets as f64 * target_weight)
            / (n_distractors as f64 * (1.0 - distractor_fraction))
    } else {
        0.0
    };

    let file = File::create(output)
        .with_context(|| format!("Cannot create weights: {}", output.display()))?;
    let mut writer = BufWriter::new(file);

    writeln!(writer, "# BaitBench weights file")?;
    writeln!(writer, "# Targets: {}, Distractors: {}", n_targets, n_distractors)?;
    writeln!(writer, "# Distractor fraction: {}", distractor_fraction)?;
    writeln!(writer, "# Target weight: {}, Distractor weight: {:.6}", target_weight, distractor_weight)?;
    writeln!(writer, "#")?;

    for id in target_ids {
        writeln!(writer, "{}\t{}", id, target_weight)?;
    }

    for id in distractor_ids {
        writeln!(writer, "{}\t{:.6}", id, distractor_weight)?;
    }

    writer.flush()?;
    Ok((target_weight, distractor_weight))
}
