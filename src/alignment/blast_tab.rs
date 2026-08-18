use anyhow::{Context, Result};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

/// A parsed blastn tabular (outfmt 6) hit, restricted to the fields needed
/// for cross-reactivity homology/coverage calculations. `query_start`/`query_end`
/// are normalized to 0-based half-open coordinates to match `alignment::paf::PafRecord`,
/// so downstream code can treat minimap2 and BLAST hits identically.
pub struct BlastHit {
    pub query_name: String,
    pub target_name: String,
    pub query_length: u32,
    pub query_start: u32,
    pub query_end: u32,
    pub matching_bases: u32,
    pub alignment_length: u32,
}

/// Parse blastn `-outfmt "6 qseqid sseqid qlen qstart qend nident length"` output.
///
/// Returns all HSPs without filtering. Skips malformed lines (< 7 columns).
pub fn parse_blast_hits(tsv_path: &Path) -> Result<Vec<BlastHit>> {
    let file = File::open(tsv_path)
        .with_context(|| format!("Cannot open BLAST results: {}", tsv_path.display()))?;
    let reader = BufReader::new(file);

    let mut hits = Vec::new();

    for line in reader.lines() {
        let line = line?;
        let line = line.trim();
        if line.is_empty() {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 7 {
            log::debug!("Skipping malformed BLAST line ({} columns)", fields.len());
            continue;
        }

        // BLAST qstart/qend are 1-based inclusive with qstart <= qend (query
        // coordinates never flip for blastn, unlike sstart/send). Convert to
        // 0-based half-open to match PafRecord's convention.
        let qstart: u32 = fields[3].parse().unwrap_or(0);
        let qend: u32 = fields[4].parse().unwrap_or(0);

        hits.push(BlastHit {
            query_name: fields[0].to_string(),
            target_name: fields[1].to_string(),
            query_length: fields[2].parse().unwrap_or(0),
            query_start: qstart.saturating_sub(1),
            query_end: qend,
            matching_bases: fields[5].parse().unwrap_or(0),
            alignment_length: fields[6].parse().unwrap_or(0),
        });
    }

    Ok(hits)
}
