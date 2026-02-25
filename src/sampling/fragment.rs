use anyhow::{Context, Result, bail};
use rand::distributions::WeightedIndex;
use rand::prelude::*;
use rand_distr::Normal;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

/// Generate random fragments from sequences with weighted sampling.
///
/// Sequences are weighted by both their explicit weight AND sequence length.
/// Fragment lengths follow a normal distribution clamped to [min, max].
pub fn generate_fragments(
    sequences: &HashMap<String, String>,
    weights: &HashMap<String, f64>,
    num_fragments: usize,
    output: &Path,
    seed: Option<u64>,
    fragment_length_mean: f64,
    fragment_length_min: usize,
    fragment_length_max: usize,
) -> Result<usize> {
    // Filter to sequences with positive weights
    let weighted_seqs: Vec<(&String, &String)> = sequences
        .iter()
        .filter(|(id, _)| weights.get(*id).copied().unwrap_or(0.0) > 0.0)
        .collect();

    if weighted_seqs.is_empty() {
        bail!("No sequences with positive weights found!");
    }

    // Calculate sampling weights: explicit_weight * sequence_length
    let seq_ids: Vec<&String> = weighted_seqs.iter().map(|(id, _)| *id).collect();
    let seq_weights: Vec<f64> = weighted_seqs
        .iter()
        .map(|(id, seq)| weights[*id] * seq.len() as f64)
        .collect();

    log::info!(
        "Loaded {} sequences, generating {} fragments from {} weighted sequences",
        sequences.len(),
        num_fragments,
        weighted_seqs.len()
    );
    log::info!(
        "Fragment length: {}-{} (mean {})",
        fragment_length_min,
        fragment_length_max,
        fragment_length_mean
    );

    let dist = WeightedIndex::new(&seq_weights)?;
    let sd = (fragment_length_max as f64 - fragment_length_min as f64) / 6.0;
    let normal = Normal::new(fragment_length_mean, sd)?;

    let mut rng: StdRng = match seed {
        Some(s) => StdRng::seed_from_u64(s),
        None => StdRng::from_entropy(),
    };

    let mut fragments_generated = 0usize;
    let mut source_counts: HashMap<&String, usize> = HashMap::new();

    let file = File::create(output)
        .with_context(|| format!("Cannot create fragments output: {}", output.display()))?;
    let mut writer = BufWriter::new(file);

    for i in 0..num_fragments {
        let idx = dist.sample(&mut rng);
        let seq_id = seq_ids[idx];
        let seq = &weighted_seqs[idx].1;

        *source_counts.entry(seq_id).or_insert(0) += 1;

        // Fragment length from normal distribution, clamped
        let frag_len_f: f64 = normal.sample(&mut rng);
        let frag_len = (frag_len_f.round() as usize)
            .max(fragment_length_min)
            .min(fragment_length_max);

        // Skip if sequence is shorter than fragment
        if seq.len() < frag_len {
            continue;
        }

        let max_start = seq.len() - frag_len;
        let start = rng.gen_range(0..=max_start);
        let fragment = &seq[start..start + frag_len];

        writeln!(
            writer,
            ">{}_fragment_{} start={} length={}",
            seq_id,
            i + 1,
            start,
            frag_len
        )?;
        writeln!(writer, "{}", fragment)?;

        fragments_generated += 1;

        if (i + 1) % 10000 == 0 {
            log::info!("Generated {} fragments...", i + 1);
        }
    }

    writer.flush()?;
    log::info!("Successfully generated {} fragments", fragments_generated);

    // Log source distribution
    log::info!("Source distribution:");
    let mut sorted_counts: Vec<_> = source_counts.iter().collect();
    sorted_counts.sort_by_key(|(id, _)| id.as_str());
    for (id, count) in sorted_counts {
        log::info!("  {}: {} fragments", id, count);
    }

    Ok(fragments_generated)
}
