use anyhow::{bail, Result};
use std::collections::{HashMap, HashSet};
use std::fs;
use std::path::{Path, PathBuf};

use crate::fasta;
use crate::io_utils;
use crate::io_utils::prefixed_join;
use crate::sampling::weights;

pub struct PrepareArgs<'a> {
    pub targets: &'a Path,
    pub genomes: Option<&'a Path>,
    pub distractors: &'a [PathBuf],
    pub sample: Option<&'a HashMap<String, f64>>,
    pub sample_target_map: Option<&'a HashMap<String, Vec<String>>>,
    /// Optional target groups file (seq_id → group_name TSV).
    pub groups: Option<&'a Path>,
    /// Optional distractor groups file (contig_id → group_name TSV).
    /// When absent, contigs are grouped by their source FASTA file stem.
    pub distractor_groups: Option<&'a Path>,
    pub distractor_fraction: f64,
    pub outdir: &'a Path,
    pub output_prefix: &'a str,
}

pub fn execute(args: &PrepareArgs) -> Result<()> {
    fs::create_dir_all(args.outdir)?;

    let pfx = args.output_prefix;
    let output_fasta = prefixed_join(args.outdir, pfx, "combined_reference.fa");
    let output_weights = prefixed_join(args.outdir, pfx, "weights.txt");
    let targets_list = prefixed_join(args.outdir, pfx, "targets.txt");
    let distractors_list = prefixed_join(args.outdir, pfx, "distractors.txt");
    let sample_list = prefixed_join(args.outdir, pfx, "sample.txt");

    // Parse target FASTA for IDs
    log::info!("Reading targets from {}...", args.targets.display());
    let target_ids = fasta::parse_fasta_ids(args.targets)?;
    log::info!("  Found {} target sequences", target_ids.len());
    let target_set: HashSet<&str> = target_ids.iter().map(|s| s.as_str()).collect();

    // Parse genome FASTA for IDs (if provided)
    let genome_ids: Option<Vec<String>> = if let Some(genomes_path) = args.genomes {
        log::info!("Reading genomes from {}...", genomes_path.display());
        let ids = fasta::parse_fasta_ids(genomes_path)?;
        log::info!("  Found {} genome sequences", ids.len());
        Some(ids)
    } else {
        None
    };

    // Parse all distractor FASTAs for IDs, recording which file each contig came from
    let mut distractor_ids: Vec<String> = Vec::new();
    // Maps contig_id -> file stem (for default distractor grouping)
    let mut distractor_file_group: HashMap<String, String> = HashMap::new();
    for distractor_path in args.distractors {
        log::info!("Reading distractors from {}...", distractor_path.display());
        let ids = fasta::parse_fasta_ids(distractor_path)?;
        log::info!("  Found {} distractor sequences", ids.len());
        let stem = distractor_path
            .file_stem()
            .map(|s| s.to_string_lossy().into_owned())
            .unwrap_or_else(|| distractor_path.display().to_string());
        for id in &ids {
            distractor_file_group.insert(id.clone(), stem.clone());
        }
        distractor_ids.extend(ids);
    }
    log::info!("  Total distractor sequences: {}", distractor_ids.len());

    // Determine the fragment source IDs (genomes if provided, otherwise targets)
    let source_ids: &[String] = genome_ids.as_deref().unwrap_or(&target_ids);
    let source_set: HashSet<&str> = source_ids.iter().map(|s| s.as_str()).collect();

    // Check for overlapping IDs between sources and distractors
    let distractor_set: HashSet<&str> = distractor_ids.iter().map(|s| s.as_str()).collect();
    let overlap: Vec<&&str> = source_set.intersection(&distractor_set).collect();
    if !overlap.is_empty() {
        log::warn!(
            "{} sequence IDs appear in both {} and distractors",
            overlap.len(),
            if args.genomes.is_some() { "genomes" } else { "targets" }
        );
        for id in overlap.iter().take(5) {
            log::warn!("  {}", id);
        }
        if overlap.len() > 5 {
            log::warn!("  ... and {} more", overlap.len() - 5);
        }
    }

    // Use provided sample manifest or default to all sources
    let sample_weights: HashMap<String, f64> = if let Some(sample_map) = args.sample {
        log::info!("Using sample with {} entries", sample_map.len());

        // Validate all sample IDs exist in source sequences (genomes or targets)
        for id in sample_map.keys() {
            if !source_set.contains(id.as_str()) {
                bail!(
                    "Sample ID '{}' not found in {} FASTA. All sample IDs must be {} sequences.",
                    id,
                    if args.genomes.is_some() { "genomes" } else { "targets" },
                    if args.genomes.is_some() { "genome" } else { "target" },
                );
            }
        }

        sample_map.clone()
    } else {
        // Default: all sources in sample with weight 1.0
        source_ids.iter().map(|id| (id.clone(), 1.0)).collect()
    };

    let sample_ids: Vec<String> = sample_weights.keys().cloned().collect();
    log::info!(
        "  Sample contains {} of {} {}",
        sample_ids.len(),
        source_ids.len(),
        if args.genomes.is_some() { "genomes" } else { "targets" }
    );

    if args.genomes.is_some() {
        let genomes_path = args.genomes.unwrap();
        let genome_id_list = genome_ids.as_ref().unwrap();

        // Build combined_reference.fa = genomes + distractors (for fragment generation)
        log::info!("Combining genomes + distractors to {}...", output_fasta.display());
        let mut input_paths: Vec<&Path> = vec![genomes_path];
        for d in args.distractors {
            input_paths.push(d.as_path());
        }
        fasta::concatenate_fastas(&input_paths, &output_fasta)?;

        // Build mapping_reference.fa = targets + distractors (for read mapping)
        let mapping_ref = prefixed_join(args.outdir, pfx, "mapping_reference.fa");
        log::info!("Combining targets + distractors to {}...", mapping_ref.display());
        let mut mapping_paths: Vec<&Path> = vec![args.targets];
        for d in args.distractors {
            mapping_paths.push(d.as_path());
        }
        fasta::concatenate_fastas(&mapping_paths, &mapping_ref)?;

        // Generate weights for genomes + distractors
        log::info!("Generating weights to {}...", output_weights.display());
        weights::generate_weights(
            genome_id_list,
            &distractor_ids,
            &sample_weights,
            args.distractor_fraction,
            &output_weights,
        )?;

        // Resolve sample-target-map: explicit mappings + auto-inferred identity matches
        let mut resolved_map: HashMap<String, Vec<String>> = HashMap::new();

        if let Some(explicit_map) = args.sample_target_map {
            // Validate and copy explicit mappings
            for (genome_id, target_list) in explicit_map {
                if !source_set.contains(genome_id.as_str()) {
                    bail!(
                        "Sample-target-map genome ID '{}' not found in genomes FASTA.",
                        genome_id
                    );
                }
                for target_id in target_list {
                    if !target_set.contains(target_id.as_str()) {
                        bail!(
                            "Sample-target-map target ID '{}' not found in targets FASTA.",
                            target_id
                        );
                    }
                }
                resolved_map.insert(genome_id.clone(), target_list.clone());
            }
        }

        // Auto-infer mappings for genomes not in the explicit map
        // Match by: (1) exact name match, or (2) target starts with "{genome_id}|"
        let mut auto_linked = 0usize;
        for genome_id in genome_id_list {
            if resolved_map.contains_key(genome_id) {
                continue;
            }
            let mut matched_targets: Vec<String> = Vec::new();

            // Exact match
            if target_set.contains(genome_id.as_str()) {
                matched_targets.push(genome_id.clone());
            }

            // Prefix match: target IDs like "{genome_id}|gene_name"
            let prefix = format!("{}|", genome_id);
            for target_id in &target_ids {
                if target_id.starts_with(&prefix) {
                    matched_targets.push(target_id.clone());
                }
            }

            if !matched_targets.is_empty() {
                log::info!(
                    "  Auto-linked genome '{}' to {} target(s): {}",
                    genome_id,
                    matched_targets.len(),
                    matched_targets.join(", ")
                );
                resolved_map.insert(genome_id.clone(), matched_targets);
                auto_linked += 1;
            }
        }

        let untargeted: Vec<&String> = genome_id_list
            .iter()
            .filter(|id| !resolved_map.contains_key(*id))
            .collect();

        log::info!(
            "Sample-target-map: {} explicit, {} auto-linked, {} untargeted",
            resolved_map.len() - auto_linked,
            auto_linked,
            untargeted.len()
        );
        if !untargeted.is_empty() {
            for id in untargeted.iter().take(5) {
                log::info!("  Untargeted genome: {}", id);
            }
            if untargeted.len() > 5 {
                log::info!("  ... and {} more untargeted", untargeted.len() - 5);
            }
        }

        // Write sample-target-map
        let map_path = prefixed_join(args.outdir, pfx, "sample_target_map.txt");
        log::info!("Writing sample-target-map to {}...", map_path.display());
        io_utils::write_sample_target_map(&resolved_map, &map_path)?;

        // Write ID lists
        log::info!("Writing ID lists...");
        io_utils::write_id_list(genome_id_list, &prefixed_join(args.outdir, pfx, "genomes.txt"))?;
        io_utils::write_id_list(&target_ids, &targets_list)?;
        io_utils::write_id_list(&distractor_ids, &distractors_list)?;
        io_utils::write_id_list(&sample_ids, &sample_list)?;
    } else {
        // Standard mode: targets are the fragment source

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
    }

    // Write target groups file (if --groups provided)
    if let Some(groups_path) = args.groups {
        log::info!("Validating and writing target groups...");
        let groups_map = io_utils::parse_groups_file(groups_path)?;
        // Validate all seq IDs in groups file exist in targets
        for seq_id in groups_map.keys() {
            if !target_set.contains(seq_id.as_str()) {
                bail!(
                    "Groups file seq ID '{}' not found in targets FASTA.",
                    seq_id
                );
            }
        }
        let out_groups = prefixed_join(args.outdir, pfx, "target_groups.tsv");
        io_utils::write_groups_file(&groups_map, &out_groups)?;
        log::info!(
            "  Wrote target groups: {} mappings ({} unique groups) to {}",
            groups_map.len(),
            groups_map.values().collect::<HashSet<_>>().len(),
            out_groups.display()
        );
    }

    // Write distractor groups file (always)
    let distractor_groups_map: HashMap<String, String> = if let Some(dg_path) = args.distractor_groups {
        log::info!("Loading distractor groups from {}...", dg_path.display());
        let map = io_utils::parse_groups_file(dg_path)?;
        // Validate all IDs exist in distractor set
        let dist_set: HashSet<&str> = distractor_ids.iter().map(|s| s.as_str()).collect();
        for id in map.keys() {
            if !dist_set.contains(id.as_str()) {
                bail!(
                    "Distractor groups file ID '{}' not found in any distractor FASTA.",
                    id
                );
            }
        }
        map
    } else {
        // Default: group by source FASTA file stem
        distractor_file_group
    };

    let out_distractor_groups = prefixed_join(args.outdir, pfx, "distractor_groups.tsv");
    io_utils::write_groups_file(&distractor_groups_map, &out_distractor_groups)?;
    {
        let n_groups = distractor_groups_map.values().collect::<HashSet<_>>().len();
        log::info!(
            "  Wrote distractor groups: {} contigs in {} group(s) to {}",
            distractor_groups_map.len(),
            n_groups,
            out_distractor_groups.display()
        );
    }

    log::info!("Prepare step complete.");
    Ok(())
}
