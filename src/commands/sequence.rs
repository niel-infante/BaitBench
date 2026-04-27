use anyhow::{bail, Context, Result};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};

use crate::external::{art_modern, badread};
use crate::io_utils::extract_source_id;

/// Which tool generates reads from captured fragments.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum ReadSimulator {
    /// Trim each fragment to `read_length` bp from its start. No errors. One read per fragment.
    Perfect,
    /// Illumina error simulation via ART-modern. Requires `art_modern` on PATH.
    Art,
    /// Long-read simulation via badread (ONT / PacBio CLR). Requires `badread` on PATH.
    Badread,
}

impl ReadSimulator {
    pub fn from_str(s: &str) -> Result<Self> {
        match s {
            "perfect" => Ok(Self::Perfect),
            "art"     => Ok(Self::Art),
            "badread" => Ok(Self::Badread),
            other => bail!(
                "Unknown read simulator '{}'. Valid values: perfect, art, badread",
                other
            ),
        }
    }
}

pub struct SequenceArgs<'a> {
    pub input: &'a Path,
    /// Primary output (single-end reads, or R1 for paired-end).
    pub output: &'a Path,
    /// R2 output path — only used when `paired_end = true`.
    pub output_r2: Option<&'a Path>,
    pub read_length: usize,
    pub num_sequences: Option<usize>,
    pub seed: Option<u64>,
    pub simulator: ReadSimulator,
    /// Error/quality profile. Ignored for `Perfect`.
    /// For `Art`: ART-modern built-in profile (e.g. `HiSeq2500_150bp`).
    /// For `Pbsim`: abstract preset (`ont`, `ont-hq`, `pacbio-rsii`, `pacbio-sequel`, `hifi`).
    pub sequencer_profile: String,
    /// Reads generated per fragment relative to fragment length (Art/Pbsim only).
    pub coverage_depth: f64,
    /// Enable paired-end output (Art only).
    pub paired_end: bool,
    /// Mean insert size for paired-end (Art + paired_end only).
    pub pe_frag_len_mean: usize,
    /// Insert size std-dev for paired-end (Art + paired_end only).
    pub pe_frag_len_sd: usize,
}

pub fn execute(args: &SequenceArgs) -> Result<()> {
    match args.simulator {
        ReadSimulator::Perfect => execute_perfect(args),
        ReadSimulator::Art     => execute_art(args),
        ReadSimulator::Badread => execute_badread(args),
    }
}

// ── Perfect (original behavior) ───────────────────────────────────────────────

fn execute_perfect(args: &SequenceArgs) -> Result<()> {
    match args.num_sequences {
        Some(n) => execute_perfect_sampled(args, n),
        None    => execute_perfect_passthrough(args),
    }
}

fn execute_perfect_passthrough(args: &SequenceArgs) -> Result<()> {
    log::info!("Sequencing fragments (perfect trim to {}bp)...", args.read_length);
    let reader = BufReader::new(
        File::open(args.input)
            .with_context(|| format!("Cannot open input: {}", args.input.display()))?,
    );
    let mut writer = BufWriter::new(
        File::create(args.output)
            .with_context(|| format!("Cannot create output: {}", args.output.display()))?,
    );
    let mut header: Option<String> = None;
    let mut seq = String::new();
    let mut count = 0usize;
    for line in reader.lines() {
        let line = line?;
        if line.starts_with('>') {
            if let Some(ref h) = header {
                write_trimmed(&mut writer, h, &seq, args.read_length)?;
                count += 1;
                seq.clear();
            }
            header = Some(line);
        } else if !line.is_empty() {
            seq.push_str(&line);
        }
    }
    if let Some(ref h) = header {
        write_trimmed(&mut writer, h, &seq, args.read_length)?;
        count += 1;
    }
    writer.flush()?;
    log::info!("Sequenced {} fragments into reads", count);
    Ok(())
}

fn execute_perfect_sampled(args: &SequenceArgs, n: usize) -> Result<()> {
    log::info!("Sequencing {} reads (sampled, trim to {}bp)...", n, args.read_length);
    let frags = load_fasta_records(args.input)?;
    let mut rng: StdRng = match args.seed {
        Some(s) => StdRng::seed_from_u64(s),
        None    => StdRng::from_entropy(),
    };
    let mut writer = BufWriter::new(
        File::create(args.output)
            .with_context(|| format!("Cannot create output: {}", args.output.display()))?,
    );
    for i in 0..n {
        let idx = rng.gen_range(0..frags.len());
        let (header, seq) = &frags[idx];
        let id_part = header.strip_prefix('>').unwrap_or(header)
            .split_whitespace().next().unwrap_or("");
        let source = extract_source_id(id_part).unwrap_or(id_part);
        let meta = header_metadata(header);
        write!(writer, ">{}_fragment_{}", source, i + 1)?;
        if !meta.is_empty() { write!(writer, " {}", meta)?; }
        writeln!(writer)?;
        write_trimmed_seq(&mut writer, seq, args.read_length)?;
    }
    writer.flush()?;
    log::info!("Sampled {} reads", n);
    Ok(())
}

// ── ART-modern ────────────────────────────────────────────────────────────────

fn execute_art(args: &SequenceArgs) -> Result<()> {
    if args.paired_end {
        log::info!(
            "Sequencing with ART-modern (PE, profile={}, read_len={}, depth={})...",
            args.sequencer_profile, args.read_length, args.coverage_depth
        );
    } else {
        log::info!(
            "Sequencing with ART-modern (SE, profile={}, read_len={}, depth={})...",
            args.sequencer_profile, args.read_length, args.coverage_depth
        );
    }
    art_modern::check_available()?;

    let tmp_r1  = tmp_path(args.output, "art_r1.fq");
    let tmp_r2  = tmp_path(args.output, "art_r2.fq");
    let tmp_sam = tmp_path(args.output, "art.sam");
    let log     = tmp_path(args.output, "art.log");

    art_modern::run_simulation(
        args.input,
        &tmp_r1,
        if args.paired_end { Some(&tmp_r2) } else { None },
        &tmp_sam,
        &log,
        args.read_length,
        &args.sequencer_profile,
        args.coverage_depth,
        args.paired_end,
        args.pe_frag_len_mean,
        args.pe_frag_len_sd,
    )?;

    let qname_to_frag = parse_sam_qname_to_rname(&tmp_sam)?;

    // Write full FASTA (before optional sampling)
    let full_r1 = match args.num_sequences {
        Some(_) => tmp_path(args.output, "art_r1_full.fa"),
        None    => args.output.to_path_buf(),
    };
    fastq_to_renamed_fasta(&tmp_r1, &full_r1, &qname_to_frag)?;

    let full_r2 = if args.paired_end {
        let r2_out = args.output_r2.ok_or_else(|| {
            anyhow::anyhow!("output_r2 required when paired_end = true")
        })?;
        let path = match args.num_sequences {
            Some(_) => tmp_path(r2_out, "art_r2_full.fa"),
            None    => r2_out.to_path_buf(),
        };
        fastq_to_renamed_fasta(&tmp_r2, &path, &qname_to_frag)?;
        Some(path)
    } else {
        None
    };

    // Optional sampling pass
    if let Some(n) = args.num_sequences {
        sample_fasta_reads(&full_r1, args.output, n, args.seed)?;
        let _ = std::fs::remove_file(&full_r1);
        if args.paired_end {
            let r2_full = full_r2.as_ref().unwrap();
            let r2_out  = args.output_r2.unwrap();
            sample_fasta_reads(r2_full, r2_out, n, args.seed)?;
            let _ = std::fs::remove_file(r2_full);
        }
    }

    for p in [&tmp_r1, &tmp_r2, &tmp_sam, &log] {
        let _ = std::fs::remove_file(p);
    }
    Ok(())
}

// ── badread (ONT / PacBio CLR) ───────────────────────────────────────────────

fn execute_badread(args: &SequenceArgs) -> Result<()> {
    log::info!(
        "Sequencing with badread (profile={}, depth={})...",
        args.sequencer_profile, args.coverage_depth
    );
    badread::check_available()?;

    let tmp_fastq = tmp_path(args.output, "badread.fq");
    let log_path  = tmp_path(args.output, "badread.log");

    badread::run_simulation(
        args.input,
        &tmp_fastq,
        &args.sequencer_profile,
        args.coverage_depth,
        args.seed,
        &log_path,
    )?;

    let full_out = match args.num_sequences {
        Some(_) => tmp_path(args.output, "badread_full.fa"),
        None    => args.output.to_path_buf(),
    };
    badread_fastq_to_fasta(&tmp_fastq, &full_out)?;

    for p in [&tmp_fastq, &log_path] {
        let _ = std::fs::remove_file(p);
    }

    if let Some(n) = args.num_sequences {
        sample_fasta_reads(&full_out, args.output, n, args.seed)?;
        let _ = std::fs::remove_file(&full_out);
    }
    Ok(())
}

/// Convert a badread FASTQ output to renamed FASTA.
///
/// Badread read descriptions contain the source reference ID:
///   `@{uuid} {ref_name},{start}-{end} strand=...`
/// We extract `ref_name` (the original fragment header ID), apply
/// `extract_source_id` to strip the `_fragment_N` suffix, and emit
/// `>{source_id}_fragment_{counter}` entries.
fn badread_fastq_to_fasta(fastq: &Path, fasta_out: &Path) -> Result<()> {
    let reader = BufReader::new(
        File::open(fastq)
            .with_context(|| format!("Cannot open badread FASTQ: {}", fastq.display()))?,
    );
    let mut writer = BufWriter::new(
        File::create(fasta_out)
            .with_context(|| format!("Cannot create FASTA: {}", fasta_out.display()))?,
    );

    let mut counter = 0usize;
    let mut lines = reader.lines();
    while let Some(header) = lines.next() {
        let header   = header?;
        let seq      = lines.next().ok_or_else(|| anyhow::anyhow!("Truncated FASTQ"))??;
        let _plus    = lines.next().ok_or_else(|| anyhow::anyhow!("Truncated FASTQ"))??;
        let _qual    = lines.next().ok_or_else(|| anyhow::anyhow!("Truncated FASTQ"))??;

        // "@{uuid} {ref_name},{start}-{end} strand=..."
        let desc      = header.strip_prefix('@').unwrap_or(&header);
        let ref_field = desc.split_whitespace().nth(1).unwrap_or("");
        let frag_id   = ref_field.split(',').next().unwrap_or(ref_field);
        let source_id = extract_source_id(frag_id).unwrap_or(frag_id);

        counter += 1;
        writeln!(writer, ">{}_fragment_{}", source_id, counter)?;
        writeln!(writer, "{}", seq)?;
    }
    writer.flush()?;
    log::debug!("badread: converted {} reads → {}", counter, fasta_out.display());
    Ok(())
}

// ── Shared helpers ─────────────────────────────────────────────────────────────

fn write_trimmed(writer: &mut impl Write, header: &str, seq: &str, read_length: usize) -> Result<()> {
    writeln!(writer, "{}", header)?;
    write_trimmed_seq(writer, seq, read_length)
}

fn write_trimmed_seq(writer: &mut impl Write, seq: &str, read_length: usize) -> Result<()> {
    let trimmed = if seq.len() > read_length { &seq[..read_length] } else { seq };
    writeln!(writer, "{}", trimmed)?;
    Ok(())
}

fn load_fasta_records(path: &Path) -> Result<Vec<(String, String)>> {
    let reader = BufReader::new(
        File::open(path).with_context(|| format!("Cannot open: {}", path.display()))?,
    );
    let mut records = Vec::new();
    let mut cur_header: Option<String> = None;
    let mut cur_seq = String::new();
    for line in reader.lines() {
        let line = line?;
        if line.starts_with('>') {
            if let Some(h) = cur_header.take() {
                records.push((h, std::mem::take(&mut cur_seq)));
            }
            cur_header = Some(line);
        } else if !line.is_empty() {
            cur_seq.push_str(&line);
        }
    }
    if let Some(h) = cur_header {
        records.push((h, cur_seq));
    }
    Ok(records)
}


fn header_metadata(header: &str) -> &str {
    let h = header.strip_prefix('>').unwrap_or(header);
    match h.find(char::is_whitespace) {
        Some(pos) => h[pos..].trim_start(),
        None      => "",
    }
}

/// Parse a SAM file and return a map of read QNAME → source fragment RNAME.
/// Only primary alignments to mapped references are included.
fn parse_sam_qname_to_rname(sam_path: &Path) -> Result<HashMap<String, String>> {
    let reader = BufReader::new(
        File::open(sam_path)
            .with_context(|| format!("Cannot open SAM: {}", sam_path.display()))?,
    );
    let mut map = HashMap::new();
    for line in reader.lines() {
        let line = line?;
        if line.starts_with('@') { continue; }
        let fields: Vec<&str> = line.splitn(6, '\t').collect();
        if fields.len() < 4 { continue; }
        let flag: u16 = fields[1].parse().unwrap_or(0);
        if flag & 0x4 != 0 || flag & 0x100 != 0 || flag & 0x800 != 0 { continue; }
        let qname = fields[0].to_string();
        let rname = fields[2].to_string();
        if rname == "*" { continue; }
        map.entry(qname).or_insert(rname);
    }
    Ok(map)
}

/// Read a plain FASTQ and write a renamed FASTA using the qname→fragment map.
/// Each output read is named `{source_id}_fragment_{counter}`.
fn fastq_to_renamed_fasta(
    fastq_path: &Path,
    fasta_out: &Path,
    qname_to_frag: &HashMap<String, String>,
) -> Result<()> {
    let reader = BufReader::new(
        File::open(fastq_path)
            .with_context(|| format!("Cannot open FASTQ: {}", fastq_path.display()))?,
    );
    let mut writer = BufWriter::new(
        File::create(fasta_out)
            .with_context(|| format!("Cannot create FASTA: {}", fasta_out.display()))?,
    );

    let mut counter = 0usize;
    let mut lines = reader.lines();
    while let Some(header_line) = lines.next() {
        let header_line = header_line?;
        let seq_line    = lines.next().ok_or_else(|| anyhow::anyhow!("Truncated FASTQ"))??;
        let _plus       = lines.next().ok_or_else(|| anyhow::anyhow!("Truncated FASTQ"))??;
        let _qual       = lines.next().ok_or_else(|| anyhow::anyhow!("Truncated FASTQ"))??;

        let qname = header_line.strip_prefix('@').unwrap_or(&header_line)
            .split_whitespace().next().unwrap_or("");
        // Strip /1 /2 paired-end suffix if present
        let qname_base = qname.trim_end_matches("/1").trim_end_matches("/2");

        let source_id = if let Some(frag_id) = qname_to_frag.get(qname_base) {
            extract_source_id(frag_id).unwrap_or(frag_id.as_str()).to_string()
        } else {
            extract_source_id(qname_base).unwrap_or(qname_base).to_string()
        };

        counter += 1;
        writeln!(writer, ">{}_fragment_{}", source_id, counter)?;
        writeln!(writer, "{}", seq_line)?;
    }
    writer.flush()?;
    log::debug!("Converted {} reads → {}", counter, fasta_out.display());
    Ok(())
}


/// Sample `n` reads with replacement from a FASTA file into another FASTA file.
fn sample_fasta_reads(input: &Path, output: &Path, n: usize, seed: Option<u64>) -> Result<()> {
    let records = load_fasta_records(input)?;
    if records.is_empty() {
        File::create(output)?;
        return Ok(());
    }
    let mut rng: StdRng = match seed {
        Some(s) => StdRng::seed_from_u64(s),
        None    => StdRng::from_entropy(),
    };
    let mut writer = BufWriter::new(
        File::create(output)
            .with_context(|| format!("Cannot create: {}", output.display()))?,
    );
    for i in 0..n {
        let idx = rng.gen_range(0..records.len());
        let (header, seq) = &records[idx];
        let id_part = header.strip_prefix('>').unwrap_or(header)
            .split_whitespace().next().unwrap_or("");
        let source = extract_source_id(id_part).unwrap_or(id_part);
        write!(writer, ">{}_fragment_{}", source, i + 1)?;
        let meta = header_metadata(header);
        if !meta.is_empty() { write!(writer, " {}", meta)?; }
        writeln!(writer)?;
        writeln!(writer, "{}", seq)?;
    }
    writer.flush()?;
    Ok(())
}

/// Construct a sibling temp file path with a descriptive suffix.
fn tmp_path(base: &Path, suffix: &str) -> PathBuf {
    let stem = base.file_stem().unwrap_or_default().to_string_lossy();
    let dir  = base.parent().unwrap_or(Path::new("."));
    dir.join(format!("{}_{}", stem, suffix))
}
