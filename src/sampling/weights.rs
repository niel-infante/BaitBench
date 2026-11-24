use anyhow::{bail, Context, Result};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

/// Parse a weights file (tab-separated: id<tab>weight, ignoring # comments).
pub fn parse_weights(path: &Path) -> Result<HashMap<String, f64>> {
    let file =
        File::open(path).with_context(|| format!("Cannot open weights: {}", path.display()))?;
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
            let weight: f64 = parts[1]
                .parse()
                .with_context(|| format!("Invalid weight value: {}", parts[1]))?;
            weights.insert(id, weight);
        }
    }

    Ok(weights)
}

/// Generate weights for targets and distractors.
///
/// Sample targets use weights from `sample_weights` map.
/// Non-sample targets get weight 0.0 (no reads generated).
/// Distractors are weighted to produce `distractor_fraction` of total reads.
///
/// Formula: distractor_weight = (fraction * total_sample_weight) / (n_distractors * (1 - fraction))
///
/// When no `--sample` is provided, all targets are in `sample_weights` with weight 1.0,
/// so this reduces to the original formula.
pub fn generate_weights(
    target_ids: &[String],
    distractor_ids: &[String],
    sample_weights: &HashMap<String, f64>,
    distractor_fraction: f64,
    output: &Path,
) -> Result<()> {
    let n_targets = target_ids.len();
    let n_sample = sample_weights.len();
    let n_distractors = distractor_ids.len();

    if n_sample == 0 {
        bail!("No sample sequences found!");
    }

    let total_sample_weight: f64 = sample_weights.values().sum();

    let distractor_weight = if n_distractors > 0 && distractor_fraction > 0.0 {
        if distractor_fraction >= 1.0 {
            bail!("distractor_fraction must be less than 1.0");
        }
        (distractor_fraction * total_sample_weight)
            / (n_distractors as f64 * (1.0 - distractor_fraction))
    } else {
        0.0
    };

    let file = File::create(output)
        .with_context(|| format!("Cannot create weights: {}", output.display()))?;
    let mut writer = BufWriter::new(file);

    writeln!(writer, "# BaitBench weights file")?;
    writeln!(
        writer,
        "# Targets: {} ({} in sample), Distractors: {}",
        n_targets, n_sample, n_distractors
    )?;
    writeln!(writer, "# Distractor fraction: {}", distractor_fraction)?;
    writeln!(
        writer,
        "# Total sample weight: {}, Distractor weight: {:.6}",
        total_sample_weight, distractor_weight
    )?;
    writeln!(writer, "#")?;

    for id in target_ids {
        let w = sample_weights.get(id).copied().unwrap_or(0.0);
        // `{:e}`, not `{:.6}`: a fixed six decimals silently rounds any weight below
        // 5e-7 to zero. That is not a corner case — weighting a genome by
        // `abundance / genome_length` puts weights around 1e-8, and the whole sample
        // then lands at zero weight and the run dies with "no reference sequence can
        // yield a fragment", pointing at fragment length rather than at the real cause.
        writeln!(writer, "{}\t{:e}", id, w)?;
    }

    for id in distractor_ids {
        writeln!(writer, "{}\t{:e}", id, distractor_weight)?;
    }

    writer.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn strs(ids: &[&str]) -> Vec<String> {
        ids.iter().map(|s| s.to_string()).collect()
    }

    fn sample(pairs: &[(&str, f64)]) -> HashMap<String, f64> {
        pairs.iter().map(|(k, v)| (k.to_string(), *v)).collect()
    }

    fn written_weights(path: &std::path::Path) -> HashMap<String, f64> {
        parse_weights(path).unwrap()
    }

    #[test]
    fn weight_sample_target_uses_given_weight() {
        let tmp = tempfile::NamedTempFile::new().unwrap();
        generate_weights(
            &strs(&["t1"]),
            &[],
            &sample(&[("t1", 2.5)]),
            0.0,
            tmp.path(),
        )
        .unwrap();
        let w = written_weights(tmp.path());
        assert!((w["t1"] - 2.5).abs() < 1e-5);
    }

    #[test]
    fn weight_nonsample_target_is_zero() {
        let tmp = tempfile::NamedTempFile::new().unwrap();
        generate_weights(
            &strs(&["t1", "t2"]),
            &[],
            &sample(&[("t1", 1.0)]),
            0.0,
            tmp.path(),
        )
        .unwrap();
        let w = written_weights(tmp.path());
        assert!((w["t2"] - 0.0).abs() < 1e-10);
    }

    #[test]
    fn weight_distractor_formula() {
        // 1 sample (w=1.0), 1 distractor, df=0.5 → (0.5×1.0)/(1×0.5) = 1.0
        let tmp = tempfile::NamedTempFile::new().unwrap();
        generate_weights(
            &strs(&["t1"]),
            &strs(&["d1"]),
            &sample(&[("t1", 1.0)]),
            0.5,
            tmp.path(),
        )
        .unwrap();
        let w = written_weights(tmp.path());
        assert!((w["d1"] - 1.0).abs() < 1e-5);
    }

    #[test]
    fn weight_two_distractors_halved() {
        // 1 sample (w=1.0), 2 distractors, df=0.5 → (0.5×1.0)/(2×0.5) = 0.5 each
        let tmp = tempfile::NamedTempFile::new().unwrap();
        generate_weights(
            &strs(&["t1"]),
            &strs(&["d1", "d2"]),
            &sample(&[("t1", 1.0)]),
            0.5,
            tmp.path(),
        )
        .unwrap();
        let w = written_weights(tmp.path());
        assert!((w["d1"] - 0.5).abs() < 1e-5);
        assert!((w["d2"] - 0.5).abs() < 1e-5);
    }

    #[test]
    fn weight_distractor_zero_when_fraction_zero() {
        let tmp = tempfile::NamedTempFile::new().unwrap();
        generate_weights(
            &strs(&["t1"]),
            &strs(&["d1"]),
            &sample(&[("t1", 1.0)]),
            0.0,
            tmp.path(),
        )
        .unwrap();
        let w = written_weights(tmp.path());
        assert!((w["d1"] - 0.0).abs() < 1e-10);
    }

    #[test]
    fn weight_empty_sample_errors() {
        let tmp = tempfile::NamedTempFile::new().unwrap();
        assert!(generate_weights(&strs(&["t1"]), &[], &HashMap::new(), 0.0, tmp.path()).is_err());
    }

    #[test]
    fn weight_distractor_fraction_ge_1_errors() {
        let tmp = tempfile::NamedTempFile::new().unwrap();
        assert!(generate_weights(
            &strs(&["t1"]),
            &strs(&["d1"]),
            &sample(&[("t1", 1.0)]),
            1.0,
            tmp.path()
        )
        .is_err());
    }

    /// Small weights must survive being written and read back.
    ///
    /// Weighting a genome by `abundance / genome_length` lands around 1e-8. Written
    /// with a fixed six decimals those all became 0.000000, every sequence was dropped
    /// as zero-weight, and the run failed with a message about fragment length.
    #[test]
    fn small_weights_survive_the_round_trip() {
        use std::io::Write as _;
        let dir = std::env::temp_dir();
        let path = dir.join(format!("bb_small_weights_{}.txt", std::process::id()));
        {
            let mut f = std::fs::File::create(&path).unwrap();
            writeln!(f, "# BaitBench weights file").unwrap();
            for (id, w) in [("a", 3.44e-8_f64), ("b", 6.7e-9), ("c", 1.0)] {
                writeln!(f, "{}\t{:e}", id, w).unwrap();
            }
        }
        let got = parse_weights(&path).unwrap();
        assert!((got["a"] - 3.44e-8).abs() < 1e-18, "got {}", got["a"]);
        assert!((got["b"] - 6.7e-9).abs() < 1e-18, "got {}", got["b"]);
        assert!(
            got["a"] > 0.0 && got["b"] > 0.0,
            "small weights must not become zero"
        );
        // And the ratio, which is what actually drives composition, is preserved.
        assert!(((got["a"] / got["b"]) - (3.44e-8 / 6.7e-9)).abs() < 1e-6);
        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn parse_weights_round_trip() {
        let mut tmp = tempfile::NamedTempFile::new().unwrap();
        writeln!(tmp, "# comment").unwrap();
        writeln!(tmp, "seq1\t0.500000").unwrap();
        writeln!(tmp, "seq2\t2.000000").unwrap();
        let w = written_weights(tmp.path());
        assert!((w["seq1"] - 0.5).abs() < 1e-5);
        assert!((w["seq2"] - 2.0).abs() < 1e-5);
        assert!(!w.contains_key("# comment"));
    }
}
