use anyhow::{Context, Result};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};

use crate::io_utils::extract_source_id;

pub struct SequenceArgs<'a> {
    pub input: &'a Path,
    pub output: &'a Path,
    pub read_length: usize,
    pub num_sequences: Option<usize>,
    pub seed: Option<u64>,
}

pub fn execute(args: &SequenceArgs) -> Result<()> {
    match args.num_sequences {
        Some(n) => execute_sampled(args, n),
        None => execute_passthrough(args),
    }
}

/// Original 1:1 behavior: trim every fragment to read_length.
fn execute_passthrough(args: &SequenceArgs) -> Result<()> {
    log::info!(
        "Sequencing captured fragments (trim to {}bp)...",
        args.read_length
    );

    let file = File::open(args.input)
        .with_context(|| format!("Cannot open input: {}", args.input.display()))?;
    let reader = BufReader::new(file);

    let out_file = File::create(args.output)
        .with_context(|| format!("Cannot create output: {}", args.output.display()))?;
    let mut writer = BufWriter::new(out_file);

    let mut current_header: Option<String> = None;
    let mut current_seq = String::new();
    let mut count = 0usize;

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('>') {
            if let Some(ref header) = current_header {
                write_trimmed_record(&mut writer, header, &current_seq, args.read_length)?;
                count += 1;
                current_seq.clear();
            }
            current_header = Some(line);
        } else if !line.is_empty() {
            current_seq.push_str(&line);
        }
    }
    // Write last record
    if let Some(ref header) = current_header {
        write_trimmed_record(&mut writer, header, &current_seq, args.read_length)?;
        count += 1;
    }

    writer.flush()?;
    log::info!("Sequenced {} fragments into reads", count);
    Ok(())
}

/// Sample num_sequences fragments with replacement, renumber, trim.
fn execute_sampled(args: &SequenceArgs, num_sequences: usize) -> Result<()> {
    log::info!(
        "Sequencing {} reads (sampled with replacement, trim to {}bp)...",
        num_sequences,
        args.read_length
    );

    // Load all fragments into memory
    let fragments = load_fragments(args.input)?;
    if fragments.is_empty() {
        log::warn!("No fragments found in input");
        File::create(args.output)?;
        return Ok(());
    }
    log::info!("Loaded {} fragments for sampling", fragments.len());

    let mut rng: StdRng = match args.seed {
        Some(s) => StdRng::seed_from_u64(s),
        None => StdRng::from_entropy(),
    };

    let out_file = File::create(args.output)
        .with_context(|| format!("Cannot create output: {}", args.output.display()))?;
    let mut writer = BufWriter::new(out_file);

    for i in 0..num_sequences {
        let idx = rng.gen_range(0..fragments.len());
        let header: &str = &fragments[idx].0;
        let seq: &str = &fragments[idx].1;

        // Extract source ID and metadata from original header
        let id_part = header
            .strip_prefix('>')
            .unwrap_or(header)
            .split_whitespace()
            .next()
            .unwrap_or("");
        let source_id = extract_source_id(id_part).unwrap_or(id_part);
        let metadata = parse_header_metadata(header);

        // Write with new sequential ID
        write!(writer, ">{}_fragment_{}", source_id, i + 1)?;
        if !metadata.is_empty() {
            write!(writer, " {}", metadata)?;
        }
        writeln!(writer)?;

        let trimmed = if seq.len() > args.read_length {
            &seq[..args.read_length]
        } else {
            seq
        };
        writeln!(writer, "{}", trimmed)?;
    }

    writer.flush()?;
    log::info!(
        "Sampled {} reads from {} fragments",
        num_sequences,
        fragments.len()
    );
    Ok(())
}

/// Load all FASTA records into memory as (header_line, sequence) pairs.
fn load_fragments(path: &Path) -> Result<Vec<(String, String)>> {
    let file = File::open(path)
        .with_context(|| format!("Cannot open input: {}", path.display()))?;
    let reader = BufReader::new(file);

    let mut fragments = Vec::new();
    let mut current_header: Option<String> = None;
    let mut current_seq = String::new();

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('>') {
            if let Some(header) = current_header.take() {
                fragments.push((header, std::mem::take(&mut current_seq)));
            }
            current_header = Some(line);
        } else if !line.is_empty() {
            current_seq.push_str(&line);
        }
    }
    if let Some(header) = current_header {
        fragments.push((header, current_seq));
    }

    Ok(fragments)
}

/// Extract the metadata portion (everything after the first whitespace token) from a FASTA header.
fn parse_header_metadata(header: &str) -> &str {
    let h = header.strip_prefix('>').unwrap_or(header);
    match h.find(char::is_whitespace) {
        Some(pos) => h[pos..].trim_start(),
        None => "",
    }
}

fn write_trimmed_record(
    writer: &mut impl Write,
    header: &str,
    sequence: &str,
    read_length: usize,
) -> Result<()> {
    let trimmed = if sequence.len() > read_length {
        &sequence[..read_length]
    } else {
        sequence
    };
    writeln!(writer, "{}", header)?;
    writeln!(writer, "{}", trimmed)?;
    Ok(())
}
