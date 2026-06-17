use anyhow::{Context, Result};
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use crate::alignment::paf;
use crate::external::minimap2;

/// A single pairwise similarity record between two targets.
pub struct TargetSimilarity {
    pub target_a: String,
    pub target_b: String,
    pub identity_pct: f64,
    pub matching_bases: u32,
    pub len_a: u32,
    pub len_b: u32,
}

/// Per-species discriminability summary.
pub struct SpeciesDiscriminability {
    pub species_id: String,
    pub total_targets: usize,
    pub unique_targets: usize,
    pub shared_targets: usize,
    pub discriminability_score: f64,
    pub confusable_with: Vec<(String, usize)>,
}

/// Computed similarity context for species-level analysis.
pub struct SimilarityContext {
    /// target_id -> set of similar target_ids from OTHER species (bidirectional).
    pub cross_species_similar: HashMap<String, HashSet<String>>,
    /// target_id -> true if unique (no cross-species similarity).
    pub target_is_unique: HashMap<String, bool>,
    /// species_id -> list of its target_ids.
    pub species_targets: HashMap<String, Vec<String>>,
    /// target_id -> species that own it.
    pub target_to_species: HashMap<String, Vec<String>>,
}

/// Compute target-vs-target similarity using sequence alignment.
///
/// Aligns all targets against themselves and returns pairwise similarities
/// above the identity threshold. Self-hits (same target) are excluded.
pub fn compute_target_similarity(
    targets_fasta: &Path,
    minimap_preset: &str,
    identity_threshold: f64,
    work_dir: &Path,
) -> Result<Vec<TargetSimilarity>> {
    minimap2::check_available()?;

    let paf_path = work_dir.join("target_vs_target.paf");
    let log_path = work_dir.join("target_vs_target.log");

    log::info!("Aligning targets against themselves for similarity analysis...");
    minimap2::xreact_align(
        minimap_preset,
        targets_fasta,
        targets_fasta,
        &paf_path,
        &log_path,
    )?;

    let records = paf::parse_paf_records(&paf_path)?;
    log::info!("Parsed {} PAF records from target-vs-target alignment", records.len());

    let mut similarities = Vec::new();
    for rec in &records {
        // Skip self-hits
        if rec.query_name == rec.target_name {
            continue;
        }

        let min_len = rec.query_length.min(rec.target_length);
        if min_len == 0 {
            continue;
        }
        let identity = rec.matching_bases as f64 / min_len as f64 * 100.0;

        if identity >= identity_threshold {
            similarities.push(TargetSimilarity {
                target_a: rec.query_name.clone(),
                target_b: rec.target_name.clone(),
                identity_pct: identity,
                matching_bases: rec.matching_bases,
                len_a: rec.query_length,
                len_b: rec.target_length,
            });
        }
    }

    log::info!(
        "Found {} target pairs above {:.1}% identity threshold",
        similarities.len(),
        identity_threshold
    );

    Ok(similarities)
}

/// Write target similarity records to a TSV file.
pub fn write_target_similarity(
    similarities: &[TargetSimilarity],
    path: &Path,
) -> Result<()> {
    let file = File::create(path)
        .with_context(|| format!("Cannot create target similarity file: {}", path.display()))?;
    let mut w = BufWriter::new(file);

    writeln!(w, "target_a\ttarget_b\tidentity_pct\tmatching_bases\tlen_a\tlen_b")?;
    for s in similarities {
        writeln!(
            w,
            "{}\t{}\t{:.2}\t{}\t{}\t{}",
            s.target_a, s.target_b, s.identity_pct, s.matching_bases, s.len_a, s.len_b
        )?;
    }
    w.flush()?;
    Ok(())
}

/// Parse a target_similarity.tsv file.
pub fn parse_target_similarity(path: &Path) -> Result<Vec<TargetSimilarity>> {
    let file = File::open(path)
        .with_context(|| format!("Cannot open target similarity file: {}", path.display()))?;
    let reader = BufReader::new(file);
    let mut sims = Vec::new();

    for (i, line) in reader.lines().enumerate() {
        let line = line?;
        let line = line.trim();
        if line.is_empty() || i == 0 {
            // Skip header
            continue;
        }
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 6 {
            continue;
        }
        sims.push(TargetSimilarity {
            target_a: parts[0].to_string(),
            target_b: parts[1].to_string(),
            identity_pct: parts[2].parse().unwrap_or(0.0),
            matching_bases: parts[3].parse().unwrap_or(0),
            len_a: parts[4].parse().unwrap_or(0),
            len_b: parts[5].parse().unwrap_or(0),
        });
    }

    Ok(sims)
}

/// Build the full similarity context from raw similarities and a genome-to-targets map.
///
/// The `genome_to_targets` map uses species/genome IDs as keys and target ID lists as values.
/// This is the same format as the sample-target-map.
pub fn build_similarity_context(
    similarities: &[TargetSimilarity],
    genome_to_targets: &HashMap<String, Vec<String>>,
) -> SimilarityContext {
    // Build target -> species reverse mapping
    let mut target_to_species: HashMap<String, Vec<String>> = HashMap::new();
    let mut species_targets: HashMap<String, Vec<String>> = HashMap::new();

    for (species, targets) in genome_to_targets {
        species_targets.insert(species.clone(), targets.clone());
        for t in targets {
            target_to_species
                .entry(t.clone())
                .or_default()
                .push(species.clone());
        }
    }

    // Build cross-species similarity map
    let mut cross_species_similar: HashMap<String, HashSet<String>> = HashMap::new();

    for sim in similarities {
        let species_a = target_to_species.get(&sim.target_a);
        let species_b = target_to_species.get(&sim.target_b);

        if let (Some(sp_a), Some(sp_b)) = (species_a, species_b) {
            // Check if these targets belong to different species
            let different_species = sp_a.iter().any(|sa| sp_b.iter().all(|sb| sa != sb));
            if different_species {
                cross_species_similar
                    .entry(sim.target_a.clone())
                    .or_default()
                    .insert(sim.target_b.clone());
                cross_species_similar
                    .entry(sim.target_b.clone())
                    .or_default()
                    .insert(sim.target_a.clone());
            }
        }
    }

    // Determine uniqueness for each target
    let all_targets: HashSet<&String> = target_to_species.keys().collect();
    let mut target_is_unique: HashMap<String, bool> = HashMap::new();
    for t in &all_targets {
        let is_unique = !cross_species_similar.contains_key(*t);
        target_is_unique.insert((*t).clone(), is_unique);
    }

    SimilarityContext {
        cross_species_similar,
        target_is_unique,
        species_targets,
        target_to_species,
    }
}

/// Compute per-species discriminability from a similarity context.
pub fn compute_discriminability(
    ctx: &SimilarityContext,
) -> Vec<SpeciesDiscriminability> {
    let mut results = Vec::new();

    for (species, targets) in &ctx.species_targets {
        let total = targets.len();
        let unique = targets
            .iter()
            .filter(|t| *ctx.target_is_unique.get(*t).unwrap_or(&true))
            .count();
        let shared = total - unique;

        // Which other species share targets?
        let mut confusable: HashMap<String, usize> = HashMap::new();
        for t in targets {
            if let Some(similar_targets) = ctx.cross_species_similar.get(t) {
                for st in similar_targets {
                    if let Some(other_species) = ctx.target_to_species.get(st) {
                        for os in other_species {
                            if os != species {
                                *confusable.entry(os.clone()).or_default() += 1;
                            }
                        }
                    }
                }
            }
        }
        let mut confusable_with: Vec<(String, usize)> = confusable.into_iter().collect();
        confusable_with.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));

        let score = if total > 0 {
            unique as f64 / total as f64
        } else {
            0.0
        };

        results.push(SpeciesDiscriminability {
            species_id: species.clone(),
            total_targets: total,
            unique_targets: unique,
            shared_targets: shared,
            discriminability_score: score,
            confusable_with,
        });
    }

    results.sort_by(|a, b| {
        a.discriminability_score
            .partial_cmp(&b.discriminability_score)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then(a.species_id.cmp(&b.species_id))
    });

    results
}

/// Build a species-by-species confusion matrix.
///
/// Returns (species_ids_sorted, matrix) where matrix[i][j] = number of targets
/// in species i that are similar to targets in species j.
pub fn build_confusion_matrix(
    ctx: &SimilarityContext,
) -> (Vec<String>, Vec<Vec<usize>>) {
    let mut species_ids: Vec<String> = ctx.species_targets.keys().cloned().collect();
    species_ids.sort();

    let n = species_ids.len();
    let species_idx: HashMap<&str, usize> = species_ids
        .iter()
        .enumerate()
        .map(|(i, s)| (s.as_str(), i))
        .collect();

    let mut matrix = vec![vec![0usize; n]; n];

    for (species, targets) in &ctx.species_targets {
        let si = match species_idx.get(species.as_str()) {
            Some(&i) => i,
            None => continue,
        };
        for t in targets {
            if let Some(similar) = ctx.cross_species_similar.get(t) {
                for st in similar {
                    if let Some(other_species_list) = ctx.target_to_species.get(st) {
                        for os in other_species_list {
                            if os != species {
                                if let Some(&sj) = species_idx.get(os.as_str()) {
                                    matrix[si][sj] += 1;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    (species_ids, matrix)
}

/// Write species discriminability results to a TSV file.
pub fn write_discriminability(
    results: &[SpeciesDiscriminability],
    path: &Path,
) -> Result<()> {
    let file = File::create(path)
        .with_context(|| format!("Cannot create discriminability file: {}", path.display()))?;
    let mut w = BufWriter::new(file);

    writeln!(
        w,
        "species_id\ttotal_targets\tunique_targets\tshared_targets\tdiscriminability_score\tconfusable_species"
    )?;
    for r in results {
        let confusable_str = if r.confusable_with.is_empty() {
            "-".to_string()
        } else {
            r.confusable_with
                .iter()
                .map(|(s, n)| format!("{}({})", s, n))
                .collect::<Vec<_>>()
                .join(",")
        };
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{:.3}\t{}",
            r.species_id,
            r.total_targets,
            r.unique_targets,
            r.shared_targets,
            r.discriminability_score,
            confusable_str
        )?;
    }
    w.flush()?;
    Ok(())
}

/// Write species confusion matrix to a TSV file.
pub fn write_confusion_matrix(
    species_ids: &[String],
    matrix: &[Vec<usize>],
    path: &Path,
) -> Result<()> {
    let file = File::create(path)
        .with_context(|| format!("Cannot create confusion matrix file: {}", path.display()))?;
    let mut w = BufWriter::new(file);

    // Header
    write!(w, "species")?;
    for s in species_ids {
        write!(w, "\t{}", s)?;
    }
    writeln!(w)?;

    // Rows
    for (i, s) in species_ids.iter().enumerate() {
        write!(w, "{}", s)?;
        for val in &matrix[i] {
            write!(w, "\t{}", val)?;
        }
        writeln!(w)?;
    }

    w.flush()?;
    Ok(())
}
