use anyhow::Result;
use std::path::Path;

use crate::alignment::sam;
use crate::external::minimap2;
use crate::fasta;

pub struct FilterArgs<'a> {
    pub host: &'a Path,
    pub reads: &'a Path,
    pub minimap_preset: &'a str,
    pub output: &'a Path,
    pub log_file: &'a Path,
}

pub fn execute(args: &FilterArgs) -> Result<()> {
    minimap2::check_available()?;

    let host_sam = args.output.with_extension("host.sam");

    // Map captured reads to host genome
    log::info!("Mapping reads to host genome...");
    minimap2::host_align(args.minimap_preset, args.host, args.reads, &host_sam, args.log_file)?;

    // Get IDs of reads that map to host
    let host_ids = sam::get_mapped_read_ids(&host_sam)?;
    log::info!("Found {} reads mapping to host", host_ids.len());

    // Get all read IDs from captured FASTA
    let all_ids: std::collections::HashSet<String> =
        fasta::parse_fasta_ids(args.reads)?.into_iter().collect();

    // Passing = all - host
    let passing: std::collections::HashSet<String> =
        all_ids.difference(&host_ids).cloned().collect();
    log::info!("{} reads passing host filter", passing.len());

    // Extract passing reads
    let count = fasta::extract_by_ids(args.reads, &passing, args.output)?;
    log::info!("Wrote {} filtered sequences", count);

    // Clean up temp SAM
    let _ = std::fs::remove_file(&host_sam);

    Ok(())
}
