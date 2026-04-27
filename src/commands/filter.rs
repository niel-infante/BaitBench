use anyhow::Result;
use std::path::Path;

use crate::alignment::sam;
use crate::external::minimap2;
use crate::fasta;

pub struct FilterArgs<'a> {
    pub host: &'a Path,
    pub reads: &'a Path,
    /// R2 reads FASTA for paired-end mode (optional).
    pub reads_r2: Option<&'a Path>,
    pub minimap_preset: &'a str,
    pub output: &'a Path,
    /// R2 filtered output path — required when `reads_r2` is `Some`.
    pub output_r2: Option<&'a Path>,
    pub log_file: &'a Path,
}

pub fn execute(args: &FilterArgs) -> Result<()> {
    minimap2::check_available()?;

    let host_sam = args.output.with_extension("host.sam");

    // Map R1 (and optionally R2) to host genome; use R1 IDs to drive filtering for both.
    log::info!("Mapping reads to host genome...");
    minimap2::host_align(args.minimap_preset, args.host, args.reads, args.reads_r2, &host_sam, args.log_file)?;

    // Get IDs of reads that map to host
    let host_ids = sam::get_mapped_read_ids(&host_sam)?;
    log::info!("Found {} reads mapping to host", host_ids.len());

    // Passing = all R1 IDs minus host-mapped
    let all_ids: std::collections::HashSet<String> =
        fasta::parse_fasta_ids(args.reads)?.into_iter().collect();
    let passing: std::collections::HashSet<String> =
        all_ids.difference(&host_ids).cloned().collect();
    log::info!("{} reads passing host filter", passing.len());

    // Extract passing reads (R1)
    let count = fasta::extract_by_ids(args.reads, &passing, args.output)?;
    log::info!("Wrote {} filtered sequences (R1)", count);

    // Extract passing R2 reads using the same ID set (paired reads share names)
    if let (Some(reads_r2), Some(output_r2)) = (args.reads_r2, args.output_r2) {
        let count_r2 = fasta::extract_by_ids(reads_r2, &passing, output_r2)?;
        log::info!("Wrote {} filtered sequences (R2)", count_r2);
    }

    // Clean up temp SAM
    let _ = std::fs::remove_file(&host_sam);

    Ok(())
}
