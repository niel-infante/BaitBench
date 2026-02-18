use anyhow::Result;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use crate::alignment::sam;

pub struct ListArgs<'a> {
    pub sam: &'a Path,
    pub output: &'a Path,
}

pub fn execute(args: &ListArgs) -> Result<()> {
    log::info!("Parsing SAM file...");
    let counts = sam::count_reads_per_reference(args.sam)?;

    // Sort by count ascending (matching original pipeline behavior)
    let mut sorted: Vec<_> = counts.into_iter().collect();
    sorted.sort_by_key(|(_, count)| *count);

    let file = File::create(args.output)?;
    let mut writer = BufWriter::new(file);

    for (ref_id, count) in &sorted {
        writeln!(writer, "{}\t{}", ref_id, count)?;
    }

    writer.flush()?;
    log::info!(
        "Detection list: {} references detected",
        sorted.len()
    );

    Ok(())
}
