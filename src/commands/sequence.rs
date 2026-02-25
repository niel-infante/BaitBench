use anyhow::{Context, Result};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

pub struct SequenceArgs<'a> {
    pub input: &'a Path,
    pub output: &'a Path,
    pub read_length: usize,
}

pub fn execute(args: &SequenceArgs) -> Result<()> {
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
