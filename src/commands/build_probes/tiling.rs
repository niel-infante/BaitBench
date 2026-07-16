use anyhow::{Context, Result};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

/// Tile probes across each sequence using a sliding window.
///
/// For each sequence, probes of `probe_length` bp are generated starting from
/// position 0 and advancing by `stride = probe_length + step` each time.
/// - `step < 0`: probes overlap (e.g., -60 with length 120 → stride 60, 50% overlap)
/// - `step = 0`: probes are perfectly tiled (no overlap, no gap)
/// - `step > 0`: gaps between probes
///
/// A final probe is always placed at the end of the sequence (last `probe_length` bp),
/// regardless of overlap with the previous probe. Sequences shorter than `probe_length`
/// produce a single probe of whatever length is available.
pub fn tile_probes(input: &Path, output: &Path, probe_length: usize, step: i64) -> Result<()> {
    let stride = (probe_length as i64 + step) as usize;

    let file =
        File::open(input).with_context(|| format!("Cannot open FASTA: {}", input.display()))?;
    let reader = BufReader::new(file);

    let out_file =
        File::create(output).with_context(|| format!("Cannot create: {}", output.display()))?;
    let mut writer = BufWriter::new(out_file);

    let mut current_id: Option<String> = None;
    let mut current_seq = String::new();

    let emit_tiles =
        |id: &str, seq: &str, writer: &mut BufWriter<File>| -> Result<()> {
            let seq_len = seq.len();
            if seq_len == 0 {
                return Ok(());
            }

            if seq_len <= probe_length {
                writeln!(writer, ">probe_{}|tile_1", id)?;
                writeln!(writer, "{}", seq)?;
                return Ok(());
            }

            let mut tile_num = 0usize;
            let mut start = 0usize;

            while start + probe_length <= seq_len {
                tile_num += 1;
                writeln!(writer, ">probe_{}|tile_{}", id, tile_num)?;
                writeln!(writer, "{}", &seq[start..start + probe_length])?;
                start += stride;
            }

            // Final probe anchored to end; only emit if it doesn't duplicate the last one.
            let final_start = seq_len - probe_length;
            let last_emitted_start = start - stride;
            if final_start != last_emitted_start {
                tile_num += 1;
                writeln!(writer, ">probe_{}|tile_{}", id, tile_num)?;
                writeln!(writer, "{}", &seq[final_start..seq_len])?;
            }

            Ok(())
        };

    for line in reader.lines() {
        let line = line?;
        let trimmed = line.trim_end();
        if trimmed.starts_with('>') {
            if let Some(ref id) = current_id {
                emit_tiles(id, &current_seq, &mut writer)?;
            }
            current_id = Some(
                trimmed.strip_prefix('>').unwrap_or("").split_whitespace().next().unwrap_or("").to_string()
            );
            current_seq.clear();
        } else if !trimmed.is_empty() {
            current_seq.push_str(&trimmed.to_uppercase());
        }
    }

    if let Some(ref id) = current_id {
        emit_tiles(id, &current_seq, &mut writer)?;
    }

    writer.flush()?;
    Ok(())
}
