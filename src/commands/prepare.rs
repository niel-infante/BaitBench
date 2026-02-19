use anyhow::{bail, Result};
use std::collections::{HashMap, HashSet};
use std::fs;
use std::path::{Path, PathBuf};

use crate::fasta;
use crate::io_utils;
use crate::sampling::weights;

pub struct PrepareArgs<'a> {
    pub targets: &'a Path,
    pub distractors: &'a [PathBuf],
    pub sample: Option<&'a Path>,
    pub distractor_fraction: f64,
    pub outdir: &'a Path,
}

pub fn execute(args: &PrepareArgs) -> Result<()> {
    fs::create_dir_all(args.outdir)?;

    let output_fasta = args.outdir.join("combined_reference.fa");
    let output_weights = args.outdir.join("weights.txt");
    let targets_list = args.outdir.join("targets.txt");
    let distractors_list = args.outdir.join("distractors.txt");
    let sample_list = args.outdir.join("sample.txt");

    // Parse target FASTA for IDs
    log::info!("Reading targets from {}...", args.targets.display());
    let target_ids = fasta::parse_fasta_ids(args.targets)?;
    log::info!("  Found {} target sequences", target_ids.len());

    // Parse all distractor FASTAs for IDs
    let mut distractor_ids: Vec<String> = Vec::new();
    for distractor_path in args.distractors {
        log::info!("Reading distractors from {}...", distractor_path.display());
        let ids = fasta::parse_fasta_ids(distractor_path)?;
        log::info!("  Found {} distractor sequences", ids.len());
        distractor_ids.extend(ids);
    }
    log::info!("  Total distractor sequences: {}", distractor_ids.len());

    // Check for overlapping IDs
    let target_set: HashSet<&str> = target_ids.iter().map(|s| s.as_str()).collect();
    let distractor_set: HashSet<&str> = distractor_ids.iter().map(|s| s.as_str()).collect();
    let overlap: Vec<&&str> = target_set.intersection(&distractor_set).collect();
    if !overlap.is_empty() {
        log::warn!(
            "{} sequence IDs appear in both targets and distractors",
            overlap.len()
        );
        for id in overlap.iter().take(5) {
            log::warn!("  {}", id);
        }
        if overlap.len() > 5 {
            log::warn!("  ... and {} more", overlap.len() - 5);
        }
    }

    // Parse sample manifest or default to all targets
    let sample_weights: HashMap<String, f64> = if let Some(sample_path) = args.sample {
        log::info!("Reading sample manifest from {}...", sample_path.display());
        let manifest = io_utils::parse_sample_manifest(sample_path)?;
        log::info!("  Found {} sample entries", manifest.len());

        // Validate all sample IDs exist in target_ids
        for id in manifest.keys() {
            if !target_set.contains(id.as_str()) {
                bail!(
                    "Sample ID '{}' not found in targets FASTA. All sample IDs must be target sequences.",
                    id
                );
            }
        }

        manifest
    } else {
        // Default: all targets in sample with weight 1.0
        target_ids.iter().map(|id| (id.clone(), 1.0)).collect()
    };

    let sample_ids: Vec<String> = sample_weights.keys().cloned().collect();
    log::info!("  Sample contains {} of {} targets", sample_ids.len(), target_ids.len());

    // Concatenate FASTAs (targets + all distractor files)
    log::info!("Combining FASTAs to {}...", output_fasta.display());
    let mut input_paths: Vec<&Path> = vec![args.targets];
    for d in args.distractors {
        input_paths.push(d.as_path());
    }
    fasta::concatenate_fastas(&input_paths, &output_fasta)?;

    // Generate weights
    log::info!("Generating weights to {}...", output_weights.display());
    weights::generate_weights(
        &target_ids,
        &distractor_ids,
        &sample_weights,
        args.distractor_fraction,
        &output_weights,
    )?;

    // Write ID lists
    log::info!("Writing ID lists...");
    io_utils::write_id_list(&target_ids, &targets_list)?;
    io_utils::write_id_list(&distractor_ids, &distractors_list)?;
    io_utils::write_id_list(&sample_ids, &sample_list)?;

    log::info!("Prepare step complete.");
    Ok(())
}
