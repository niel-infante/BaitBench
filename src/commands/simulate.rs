use anyhow::Result;
use std::path::Path;

use crate::fasta;
use crate::sampling::{self, weights};

pub struct SimulateArgs<'a> {
    pub reference: &'a Path,
    pub weights: &'a Path,
    pub num_fragments: usize,
    pub seed: Option<u64>,
    pub output: &'a Path,
    pub fragment_length_mean: f64,
    pub fragment_length_min: usize,
    pub fragment_length_max: usize,
}

pub fn execute(args: &SimulateArgs) -> Result<()> {
    if let Some(seed) = args.seed {
        log::info!("Using random seed: {}", seed);
    }

    log::info!("Parsing FASTA file...");
    let sequences = fasta::parse_fasta(args.reference)?;

    log::info!("Parsing weights file...");
    let wts = weights::parse_weights(args.weights)?;

    sampling::generate_fragments(
        &sequences,
        &wts,
        args.num_fragments,
        args.output,
        args.seed,
        args.fragment_length_mean,
        args.fragment_length_min,
        args.fragment_length_max,
    )?;

    Ok(())
}
