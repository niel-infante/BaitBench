use anyhow::{Result, bail};
use rand::prelude::*;
use rand::seq::SliceRandom;
use std::collections::HashSet;
use std::path::Path;

use crate::fasta;
use crate::io_utils;

pub struct EnrichArgs<'a> {
    pub captured: &'a Path,
    pub fragments: &'a Path,
    pub targets: &'a Path,
    pub distractors: &'a Path,
    pub fold_enrichment: f64,
    pub seed: Option<u64>,
    pub output: &'a Path,
}

pub fn execute(args: &EnrichArgs) -> Result<()> {
    if args.fold_enrichment < 1.0 {
        bail!("Fold enrichment must be >= 1.0 (1.0 = no enrichment)");
    }

    // Load target and distractor ID sets
    let target_ids = io_utils::parse_id_set(args.targets)?;
    let distractor_ids = io_utils::parse_id_set(args.distractors)?;

    // Get all fragment IDs and classify them
    let all_fragment_ids = fasta::parse_fasta_ids(args.fragments)?;
    let captured_ids: HashSet<String> = fasta::parse_fasta_ids(args.captured)?
        .into_iter()
        .collect();

    let mut n_target: usize = 0;
    let mut n_distractor: usize = 0;
    let mut captured_target_ids: Vec<String> = Vec::new();
    let mut captured_distractor_ids: Vec<String> = Vec::new();
    let mut uncaptured_distractor_ids: Vec<String> = Vec::new();

    for frag_id in &all_fragment_ids {
        let source = io_utils::extract_source_id(frag_id);
        let is_target = source.map_or(false, |s| target_ids.contains(s));
        let is_distractor = source.map_or(false, |s| distractor_ids.contains(s));

        if is_target {
            n_target += 1;
            if captured_ids.contains(frag_id) {
                captured_target_ids.push(frag_id.clone());
            }
        } else if is_distractor {
            n_distractor += 1;
            if captured_ids.contains(frag_id) {
                captured_distractor_ids.push(frag_id.clone());
            } else {
                uncaptured_distractor_ids.push(frag_id.clone());
            }
        }
    }

    let c_target = captured_target_ids.len();
    let c_distractor = captured_distractor_ids.len();

    log::info!("Pre-capture: {} target fragments, {} distractor fragments", n_target, n_distractor);
    log::info!("Binary capture: {} target captured, {} distractor captured", c_target, c_distractor);

    if n_target == 0 || n_distractor == 0 {
        log::warn!("Cannot apply fold enrichment: no target or no distractor fragments. Outputting binary capture as-is.");
        std::fs::copy(args.captured, args.output)?;
        return Ok(());
    }

    if c_target == 0 {
        log::warn!("No target fragments were captured. Outputting binary capture as-is.");
        std::fs::copy(args.captured, args.output)?;
        return Ok(());
    }

    // Calculate current fold enrichment from binary capture
    let r_pre = n_target as f64 / n_distractor as f64;
    let fe_current = if c_distractor == 0 {
        f64::INFINITY
    } else {
        (c_target as f64 / c_distractor as f64) / r_pre
    };

    log::info!("Binary capture fold enrichment: {:.1}", fe_current);
    log::info!("Requested fold enrichment: {:.1}", args.fold_enrichment);

    if fe_current <= args.fold_enrichment {
        log::info!(
            "Binary capture enrichment ({:.1}) is already at or below requested ({:.1}). Outputting binary capture as-is.",
            fe_current,
            args.fold_enrichment
        );
        std::fs::copy(args.captured, args.output)?;
        return Ok(());
    }

    // Calculate desired distractor count: D_post = C_target / (FE * R_pre)
    let d_post = (c_target as f64 / (args.fold_enrichment * r_pre)).round() as usize;
    let d_post = d_post.max(1); // at least 1 distractor to avoid division issues

    log::info!("Target distractor count in output: {}", d_post);

    // Build the final ID set
    let mut rng: StdRng = match args.seed {
        Some(s) => StdRng::seed_from_u64(s),
        None => StdRng::from_entropy(),
    };

    let mut final_ids: HashSet<String> = HashSet::new();

    // All captured targets are kept
    for id in &captured_target_ids {
        final_ids.insert(id.clone());
    }

    if d_post <= c_distractor {
        // Subsample captured distractors
        let mut captured_dist = captured_distractor_ids.clone();
        captured_dist.shuffle(&mut rng);
        for id in captured_dist.into_iter().take(d_post) {
            final_ids.insert(id);
        }
    } else {
        // Keep all captured distractors + add uncaptured ones
        for id in &captured_distractor_ids {
            final_ids.insert(id.clone());
        }
        let extra_needed = d_post - c_distractor;
        if extra_needed <= uncaptured_distractor_ids.len() {
            uncaptured_distractor_ids.shuffle(&mut rng);
            for id in uncaptured_distractor_ids.into_iter().take(extra_needed) {
                final_ids.insert(id);
            }
        } else {
            // Use all available uncaptured distractors
            for id in uncaptured_distractor_ids {
                final_ids.insert(id);
            }
            let actual_d = c_distractor + final_ids.len() - c_target;
            let actual_fe = if actual_d == 0 {
                f64::INFINITY
            } else {
                (c_target as f64 / actual_d as f64) / r_pre
            };
            log::warn!(
                "Not enough distractor fragments to achieve {:.1}x enrichment. Using all available. Actual enrichment: {:.1}x",
                args.fold_enrichment,
                actual_fe
            );
        }
    }

    // Extract final set from fragments.fa
    let count = fasta::extract_by_ids(args.fragments, &final_ids, args.output)?;

    let final_target = c_target;
    let final_distractor = count - final_target;
    let actual_fe = if final_distractor == 0 {
        f64::INFINITY
    } else {
        (final_target as f64 / final_distractor as f64) / r_pre
    };

    log::info!(
        "Enriched output: {} target + {} distractor = {} total fragments (actual FE: {:.1}x)",
        final_target,
        final_distractor,
        count,
        actual_fe
    );

    Ok(())
}
