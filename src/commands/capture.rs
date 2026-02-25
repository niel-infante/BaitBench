use anyhow::Result;
use std::path::Path;

use crate::alignment::paf;
use crate::external::{blastn, minimap2};
use crate::fasta;

#[derive(Clone, Copy, PartialEq)]
pub enum CaptureMethod {
    Minimap2,
    Blast,
}

pub struct CaptureArgs<'a> {
    pub method: CaptureMethod,
    pub probes: &'a Path,
    pub fragments: &'a Path,
    pub max_mismatches: u32,
    pub min_match_bases: u32,
    pub blast_db: Option<&'a str>,
    pub output: &'a Path,
    pub log_file: &'a Path,
    pub threads: usize,
}

pub fn execute(args: &CaptureArgs) -> Result<()> {
    let passing_ids = match args.method {
        CaptureMethod::Minimap2 => capture_minimap2(args)?,
        CaptureMethod::Blast => capture_blast(args)?,
    };

    // Extract passing fragments from FASTA
    let count = fasta::extract_by_ids(args.fragments, &passing_ids, args.output)?;
    log::info!("Captured {} sequences", count);

    Ok(())
}

fn capture_minimap2(args: &CaptureArgs) -> Result<std::collections::HashSet<String>> {
    minimap2::check_available()?;

    let paf_path = args.output.with_extension("paf");

    log::info!("Running minimap2 capture alignment...");
    minimap2::capture_align(args.probes, args.fragments, &paf_path, args.log_file)?;

    log::info!("Filtering PAF results...");
    let passing = paf::filter_paf(&paf_path, args.max_mismatches, args.min_match_bases)?;

    // Clean up temp PAF
    let _ = std::fs::remove_file(&paf_path);

    Ok(passing)
}

fn capture_blast(args: &CaptureArgs) -> Result<std::collections::HashSet<String>> {
    blastn::check_available()?;

    let blast_db = args
        .blast_db
        .ok_or_else(|| anyhow::anyhow!("--blast-db is required when using BLAST capture method"))?;

    let blast_tsv = args.output.with_extension("blast.tsv");

    log::info!("Running BLAST capture alignment...");
    blastn::capture_align(blast_db, args.fragments, &blast_tsv, args.log_file, args.threads)?;

    log::info!("Filtering BLAST results...");
    let passing = blastn::filter_blast_results(&blast_tsv, args.min_match_bases)?;

    // Clean up temp file
    let _ = std::fs::remove_file(&blast_tsv);

    Ok(passing)
}
