use anyhow::Result;
use std::path::Path;

use crate::external::minimap2;

pub struct MapArgs<'a> {
    pub reference: &'a Path,
    pub reads: &'a Path,
    pub minimap_preset: &'a str,
    pub output: &'a Path,
    pub log_file: &'a Path,
}

pub fn execute(args: &MapArgs) -> Result<()> {
    minimap2::check_available()?;

    log::info!("Mapping reads to reference...");
    minimap2::map_reads(
        args.minimap_preset,
        args.reference,
        args.reads,
        args.output,
        args.log_file,
    )?;

    log::info!("Mapping complete: {}", args.output.display());
    Ok(())
}
