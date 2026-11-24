use anyhow::Result;
use std::path::Path;

use crate::external::minimap2;

pub struct MapArgs<'a> {
    pub reference: &'a Path,
    pub reads: &'a Path,
    /// R2 reads for paired-end mode (optional).
    pub reads_r2: Option<&'a Path>,
    pub minimap_preset: &'a str,
    pub output: &'a Path,
    pub log_file: &'a Path,
    /// Rayon workers for the mapping step. Mapping is the dominant cost of a run.
    pub threads: usize,
}

pub fn execute(args: &MapArgs) -> Result<()> {
    minimap2::check_available()?;

    log::info!("Mapping reads to reference...");
    minimap2::map_reads(
        args.minimap_preset,
        args.reference,
        args.reads,
        args.reads_r2,
        args.output,
        args.log_file,
        args.threads,
    )?;

    log::info!("Mapping complete: {}", args.output.display());
    Ok(())
}
