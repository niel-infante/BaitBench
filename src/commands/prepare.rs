use anyhow::Result;
use std::collections::HashSet;
use std::fs;
use std::path::Path;

use crate::fasta;
use crate::io_utils;
use crate::sampling::weights;

pub struct PrepareArgs<'a> {
    pub targets: &'a Path,
    pub distractors: &'a Path,
    pub distractor_fraction: f64,
    pub outdir: &'a Path,
}

pub fn execute(args: &PrepareArgs) -> Result<()> {
    fs::create_dir_all(args.outdir)?;

    let output_fasta = args.outdir.join("combined_reference.fa");
    let output_weights = args.outdir.join("weights.txt");
    let targets_list = args.outdir.join("targets.txt");
    let distractors_list = args.outdir.join("distractors.txt");

    // Parse input FASTAs for IDs
    log::info!("Reading targets from {}...", args.targets.display());
    let target_ids = fasta::parse_fasta_ids(args.targets)?;
    log::info!("  Found {} target sequences", target_ids.len());

    log::info!("Reading distractors from {}...", args.distractors.display());
    let distractor_ids = fasta::parse_fasta_ids(args.distractors)?;
    log::info!("  Found {} distractor sequences", distractor_ids.len());

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

    // Concatenate FASTAs
    log::info!("Combining FASTAs to {}...", output_fasta.display());
    fasta::concatenate_fastas(args.targets, args.distractors, &output_fasta)?;

    // Generate weights
    log::info!("Generating weights to {}...", output_weights.display());
    let (target_w, distractor_w) = weights::generate_weights(
        &target_ids,
        &distractor_ids,
        args.distractor_fraction,
        &output_weights,
    )?;
    log::info!("  Target weight: {}", target_w);
    log::info!("  Distractor weight: {:.6}", distractor_w);

    // Write ID lists
    log::info!("Writing ID lists...");
    io_utils::write_id_list(&target_ids, &targets_list)?;
    io_utils::write_id_list(&distractor_ids, &distractors_list)?;

    log::info!("Prepare step complete.");
    Ok(())
}
