use anyhow::{Context, Result};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

/// Replace every N/n in each probe with a non-N base (streaming).
///
/// Preferred replacement is T; falls back through A → C → G if T would be
/// immediately adjacent to an existing T on either side.
/// Returns `(probes_with_n, total_n_replaced)`.
pub fn fix_n_bases(input: &Path, output: &Path) -> Result<(usize, usize)> {
    let file =
        File::open(input).with_context(|| format!("Cannot open FASTA: {}", input.display()))?;
    let reader = BufReader::new(file);

    let out_file =
        File::create(output).with_context(|| format!("Cannot create: {}", output.display()))?;
    let mut writer = BufWriter::new(out_file);

    let mut current_id: Option<String> = None;
    let mut current_seq = String::new();
    let mut probes_with_n = 0usize;
    let mut total_n_replaced = 0usize;

    for line in reader.lines() {
        let line = line?;
        let trimmed = line.trim_end();
        if trimmed.starts_with('>') {
            if let Some(ref id) = current_id {
                let (fixed, count) = replace_ns(&current_seq);
                if count > 0 {
                    probes_with_n += 1;
                    total_n_replaced += count;
                }
                writeln!(writer, ">{}", id)?;
                writeln!(writer, "{}", fixed)?;
            }
            current_id = Some(trimmed.strip_prefix('>').unwrap_or("").to_string());
            current_seq.clear();
        } else if !trimmed.is_empty() {
            current_seq.push_str(&trimmed.to_uppercase());
        }
    }

    if let Some(ref id) = current_id {
        let (fixed, count) = replace_ns(&current_seq);
        if count > 0 {
            probes_with_n += 1;
            total_n_replaced += count;
        }
        writeln!(writer, ">{}", id)?;
        writeln!(writer, "{}", fixed)?;
    }

    writer.flush()?;
    Ok((probes_with_n, total_n_replaced))
}

/// Replace each N/n in `seq` with a non-N base.
///
/// Scans left-to-right. For each N, tries T first; if T would be immediately
/// adjacent to the left or right neighbor, tries A, then C, then G. The left
/// neighbor is taken from the already-modified output (so consecutive N's are
/// handled correctly); the right neighbor is taken from the original sequence.
fn replace_ns(seq: &str) -> (String, usize) {
    let mut result: Vec<u8> = seq.bytes().collect();
    let mut count = 0usize;
    for i in 0..result.len() {
        if result[i] == b'N' {
            let left = if i > 0 { result[i - 1] } else { b'X' };
            let right = if i + 1 < result.len() {
                result[i + 1]
            } else {
                b'X'
            };
            let replacement = b"TACG"
                .iter()
                .copied()
                .find(|&b| b != left && b != right)
                .unwrap_or(b'T');
            result[i] = replacement;
            count += 1;
        }
    }
    (String::from_utf8(result).unwrap_or_default(), count)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// With neither neighbour a T, the first candidate wins.
    #[test]
    fn n_becomes_t_when_unambiguous() {
        assert_eq!(replace_ns("ANA"), ("ATA".to_string(), 1));
    }

    /// T is skipped when it would duplicate a neighbour, falling through to A.
    #[test]
    fn replacement_avoids_both_neighbours() {
        assert_eq!(replace_ns("TNG"), ("TAG".to_string(), 1));
        assert_eq!(replace_ns("GNT"), ("GAT".to_string(), 1));
        // T and A both excluded by the flanks, so C is taken.
        assert_eq!(replace_ns("TNA"), ("TCA".to_string(), 1));
    }

    /// The left neighbour comes from the already-modified output, so a run of Ns does
    /// not collapse into a homopolymer.
    #[test]
    fn consecutive_ns_do_not_repeat() {
        let (fixed, n) = replace_ns("NN");
        assert_eq!(n, 2);
        assert_eq!(fixed, "TA");
        let (fixed, n) = replace_ns("ANNA");
        assert_eq!(n, 2);
        assert_eq!(fixed.len(), 4);
        let b = fixed.as_bytes();
        assert_ne!(b[1], b[2], "adjacent replacements must differ: {fixed}");
        assert!(!fixed.contains('N'));
    }

    /// Ends have no neighbour on one side; the sentinel must not block a candidate.
    #[test]
    fn ns_at_the_ends_are_replaced() {
        assert_eq!(replace_ns("NA"), ("TA".to_string(), 1));
        assert_eq!(replace_ns("AN"), ("AT".to_string(), 1));
        assert_eq!(replace_ns("N"), ("T".to_string(), 1));
    }

    /// A sequence without N is returned untouched and reports no replacements.
    #[test]
    fn sequence_without_n_is_unchanged() {
        assert_eq!(replace_ns("ACGT"), ("ACGT".to_string(), 0));
    }

    /// `replace_ns` matches uppercase only; `fix_n_bases` uppercases first (see its
    /// read loop), so lowercase input is handled there rather than here.
    #[test]
    fn lowercase_n_is_not_matched_here() {
        assert_eq!(replace_ns("ana"), ("ana".to_string(), 0));
    }
}
