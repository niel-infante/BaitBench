use anyhow::{bail, Context, Result};
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

/// Parse a plain-text ID list file into a HashSet.
/// Skips empty lines and lines starting with '#'.
/// Takes the first whitespace-delimited token from each line.
pub fn parse_id_set(path: &Path) -> Result<HashSet<String>> {
    let file = File::open(path)
        .with_context(|| format!("Cannot open ID file: {}", path.display()))?;
    let reader = BufReader::new(file);
    let mut ids = HashSet::new();

    for line in reader.lines() {
        let line = line?;
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        if let Some(id) = line.split_whitespace().next() {
            ids.insert(id.to_string());
        }
    }

    Ok(ids)
}

/// Extract the source sequence ID from a fragment read name.
///
/// Read names follow the pattern `{seq_id}_fragment_{n}` (with optional
/// trailing fields like `start=... length=...`). We find the last occurrence
/// of `_fragment_` and take everything before it as the source ID.
pub fn extract_source_id(read_name: &str) -> Option<&str> {
    // Take the first whitespace-delimited token (the ID portion)
    let id_part = read_name.split_whitespace().next().unwrap_or(read_name);
    // Find the last occurrence of "_fragment_" to handle IDs that contain that substring
    id_part.rfind("_fragment_").map(|pos| &id_part[..pos])
}

/// Parse a sample manifest TSV file into a HashMap of id -> weight.
///
/// Format: `id<tab>weight` (weight is optional, defaults to 1.0).
/// Skips empty lines and lines starting with '#'.
pub fn parse_sample_manifest(path: &Path) -> Result<HashMap<String, f64>> {
    let file = File::open(path)
        .with_context(|| format!("Cannot open sample manifest: {}", path.display()))?;
    let reader = BufReader::new(file);
    let mut samples = HashMap::new();

    for line in reader.lines() {
        let line = line?;
        let line = line.trim().to_string();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let parts: Vec<&str> = line.split('\t').collect();
        let id = parts[0].trim().to_string();
        let weight = if parts.len() >= 2 {
            parts[1].trim().parse::<f64>()
                .with_context(|| format!("Invalid weight for '{}': {}", id, parts[1].trim()))?
        } else {
            1.0
        };
        samples.insert(id, weight);
    }

    Ok(samples)
}

/// Resolve a `--sample` CLI argument into a sample manifest HashMap.
///
/// If exactly one token is provided and it is a path to an existing file,
/// the file is parsed as a TSV sample manifest (`id<tab>weight`).
/// Otherwise, the tokens are parsed as an inline list: each token is treated
/// as a sample ID (weight 1.0) unless it parses as a float, in which case it
/// sets the weight of the preceding ID.
///
/// Examples:
///   `--sample manifest.tsv`             → parse file
///   `--sample t1 t2 t3`                 → {t1: 1.0, t2: 1.0, t3: 1.0}
///   `--sample t1 t2 t3 5 t4`           → {t1: 1.0, t2: 1.0, t3: 5.0, t4: 1.0}
pub fn resolve_sample_arg(tokens: &[String]) -> Result<HashMap<String, f64>> {
    if tokens.len() == 1 && Path::new(&tokens[0]).is_file() {
        parse_sample_manifest(Path::new(&tokens[0]))
    } else {
        parse_sample_inline(tokens)
    }
}

/// Parse an inline sample list into a HashMap of id -> weight.
///
/// Tokens are interpreted sequentially: each token is either a sample ID
/// (if it cannot be parsed as f64) or a weight for the preceding ID
/// (if it can be parsed as f64). IDs default to weight 1.0.
fn parse_sample_inline(tokens: &[String]) -> Result<HashMap<String, f64>> {
    if tokens.is_empty() {
        bail!("No sample IDs provided");
    }

    let mut samples = HashMap::new();
    let mut last_id: Option<String> = None;

    for token in tokens {
        if let Ok(weight) = token.parse::<f64>() {
            match &last_id {
                Some(id) => {
                    samples.insert(id.clone(), weight);
                }
                None => {
                    bail!(
                        "Weight value '{}' appears before any sample ID. \
                         The first token must be a sample ID.",
                        token
                    );
                }
            }
        } else {
            // It's a new ID — insert previous ID with default weight if not already set
            if let Some(ref id) = last_id {
                samples.entry(id.clone()).or_insert(1.0);
            }
            last_id = Some(token.clone());
        }
    }
    // Insert the last ID if not yet inserted
    if let Some(ref id) = last_id {
        samples.entry(id.clone()).or_insert(1.0);
    }

    Ok(samples)
}

/// Format a sample manifest HashMap as a human-readable string for display.
///
/// Output format: `t1 t2 t3(5.0) t4` — IDs with weight 1.0 are shown plain,
/// others show their weight in parentheses.
pub fn format_sample_display(samples: &HashMap<String, f64>) -> String {
    let mut entries: Vec<(&String, &f64)> = samples.iter().collect();
    entries.sort_by(|a, b| a.0.cmp(b.0));
    entries
        .iter()
        .map(|(id, &w)| {
            if (w - 1.0).abs() < f64::EPSILON {
                id.to_string()
            } else {
                format!("{}({:.1})", id, w)
            }
        })
        .collect::<Vec<_>>()
        .join(" ")
}

/// Write a list of IDs to a file, one per line.
pub fn write_id_list(ids: &[String], path: &Path) -> Result<()> {
    let file = File::create(path)
        .with_context(|| format!("Cannot create ID file: {}", path.display()))?;
    let mut writer = BufWriter::new(file);
    for id in ids {
        writeln!(writer, "{}", id)?;
    }
    writer.flush()?;
    Ok(())
}

/// Parse a sample-target-map TSV file into a HashMap of genome_id -> Vec<target_id>.
///
/// Format: `genome_id<tab>target_id` (one mapping per line, multiple lines per genome).
/// Skips empty lines and lines starting with '#'.
pub fn parse_sample_target_map(path: &Path) -> Result<HashMap<String, Vec<String>>> {
    let file = File::open(path)
        .with_context(|| format!("Cannot open sample-target-map: {}", path.display()))?;
    let reader = BufReader::new(file);
    let mut map: HashMap<String, Vec<String>> = HashMap::new();

    for line in reader.lines() {
        let line = line?;
        let line = line.trim().to_string();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 2 || parts[1].trim().is_empty() {
            bail!(
                "Invalid sample-target-map line (expected genome_id<tab>target_id): '{}'",
                line
            );
        }
        let genome_id = parts[0].trim().trim_start_matches('>').to_string();
        let target_id = parts[1].trim().trim_start_matches('>').to_string();
        map.entry(genome_id).or_default().push(target_id);
    }

    Ok(map)
}

/// Write a sample-target-map to a TSV file.
///
/// Format: `genome_id<tab>target_id` (one line per mapping).
pub fn write_sample_target_map(
    map: &HashMap<String, Vec<String>>,
    path: &Path,
) -> Result<()> {
    let file = File::create(path)
        .with_context(|| format!("Cannot create sample-target-map: {}", path.display()))?;
    let mut writer = BufWriter::new(file);
    writeln!(writer, "# BaitBench sample-target-map")?;
    writeln!(writer, "# genome_id\ttarget_id")?;

    let mut genome_ids: Vec<&String> = map.keys().collect();
    genome_ids.sort();
    for genome_id in genome_ids {
        let mut targets = map[genome_id].clone();
        targets.sort();
        for target_id in &targets {
            writeln!(writer, "{}\t{}", genome_id, target_id)?;
        }
    }

    writer.flush()?;
    Ok(())
}
