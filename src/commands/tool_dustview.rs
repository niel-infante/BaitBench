use anyhow::Result;
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::Path;

use crate::sdust::sdust;

pub fn execute(input: Option<&Path>, dust_threshold: f64, dust_window: usize) -> Result<()> {
    let reader: Box<dyn BufRead> = match input {
        Some(path) => Box::new(BufReader::new(File::open(path)?)),
        None => Box::new(BufReader::new(io::stdin())),
    };

    let records = parse_fasta(reader);

    for (i, (header, seq)) in records.iter().enumerate() {
        if i > 0 {
            println!();
        }

        let bytes = seq.as_bytes();
        let regions = sdust(bytes, dust_threshold, dust_window);

        let mut masked = bytes.to_vec();
        for &(start, end) in &regions {
            for b in &mut masked[start..end] {
                *b = b'X';
            }
        }

        let scores: Vec<f64> = (0..bytes.len())
            .map(|i| position_score(bytes, i, dust_window))
            .collect();

        let masked_bases: usize = regions.iter().map(|&(s, e)| e - s).sum();
        let frac_masked = masked_bases as f64 / bytes.len().max(1) as f64;

        let avg_score = if scores.is_empty() {
            0.0
        } else {
            scores.iter().sum::<f64>() / scores.len() as f64
        };

        let below: Vec<f64> = scores
            .iter()
            .copied()
            .filter(|&s| s < dust_threshold)
            .collect();
        let avg_below = if below.is_empty() {
            0.0
        } else {
            below.iter().sum::<f64>() / below.len() as f64
        };

        println!("{}", header);
        println!("{}", seq);
        println!("{}", String::from_utf8(masked).unwrap());
        println!(
            "masked: {:.1}%  avg_score: {:.2}  avg_score_below_threshold: {:.2}",
            frac_masked * 100.0,
            avg_score,
            avg_below
        );
    }

    Ok(())
}

fn parse_fasta(reader: impl BufRead) -> Vec<(String, String)> {
    let mut records = Vec::new();
    let mut current_header: Option<String> = None;
    let mut current_seq = String::new();

    for line in reader.lines() {
        let line = line.unwrap_or_default();
        let line = line.trim_end();
        if line.starts_with('>') {
            if let Some(header) = current_header.take() {
                records.push((header, current_seq.clone()));
                current_seq.clear();
            }
            current_header = Some(line.to_string());
        } else if current_header.is_some() {
            current_seq.push_str(line);
        }
    }
    if let Some(header) = current_header {
        records.push((header, current_seq));
    }
    records
}

fn position_score(seq: &[u8], pos: usize, window: usize) -> f64 {
    let len = seq.len();
    if len < 3 {
        return 0.0;
    }
    let half = window / 2;
    let start = pos.saturating_sub(half);
    let end = (start + window).min(len);
    let start = end.saturating_sub(window);
    let sub = &seq[start..end];

    if sub.len() < 3 {
        return 0.0;
    }

    let mut counts = [0u32; 64];
    let mut valid_triplets: usize = 0;

    for i in 0..sub.len() - 2 {
        let b0 = encode_base(sub[i]);
        let b1 = encode_base(sub[i + 1]);
        let b2 = encode_base(sub[i + 2]);
        if let (Some(b0), Some(b1), Some(b2)) = (b0, b1, b2) {
            let t = ((b0 << 4) | (b1 << 2) | b2) as usize;
            counts[t] += 1;
            valid_triplets += 1;
        }
    }

    if valid_triplets < 2 {
        return 0.0;
    }

    let r: u32 = counts.iter().map(|&c| c.saturating_sub(1) * c / 2).sum();
    r as f64 / (valid_triplets - 1) as f64
}

#[inline]
fn encode_base(b: u8) -> Option<u8> {
    match b {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}
