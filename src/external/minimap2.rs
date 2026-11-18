use anyhow::{Context, Result, bail};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use rammap::{Aligner, MapOpts, Preset, Strand};
use rammap::align::map::AlignFlags;

use crate::fasta;

/// (ref name → sequence length, ref name → per-position depth vector), as returned by
/// `probe_depth_in_memory`.
type RefLengthsAndCoverage = (
    std::collections::HashMap<String, usize>,
    std::collections::HashMap<String, Vec<u32>>,
);

// ---- helpers ----

fn to_preset(preset: &str) -> Result<Preset> {
    match preset {
        "sr"        => Ok(Preset::Sr),
        "map-ont"   => Ok(Preset::MapOnt),
        "map-pb"    => Ok(Preset::MapPb),
        "map-hifi"  => Ok(Preset::MapHifi),
        "map-iclr"  => Ok(Preset::MapIclr),
        "lr:hq"     => Ok(Preset::LrHq),
        "lr:hqae"   => Ok(Preset::LrHqae),
        "splice"    => Ok(Preset::Splice),
        "splice:hq" => Ok(Preset::SpliceHq),
        "splice:sr" => Ok(Preset::SpliceSr),
        "cdna"      => Ok(Preset::Cdna),
        "asm5"      => Ok(Preset::Asm5),
        "asm10"     => Ok(Preset::Asm10),
        "asm20"     => Ok(Preset::Asm20),
        "ava-ont"   => Ok(Preset::AvaOnt),
        "ava-pb"    => Ok(Preset::AvaPb),
        other => bail!("Unknown alignment preset: '{}'", other),
    }
}

fn path_str(path: &Path) -> Result<&str> {
    path.to_str()
        .ok_or_else(|| anyhow::anyhow!("Path has non-UTF-8 characters: {}", path.display()))
}

fn write_stub_log(log_file: &Path, label: &str) -> Result<()> {
    let mut f = File::create(log_file)
        .with_context(|| format!("Cannot create log: {}", log_file.display()))?;
    writeln!(f, "[rammap] {}", label)?;
    Ok(())
}

fn write_sam_header<W: Write>(w: &mut W, aligner: &Aligner) -> Result<()> {
    writeln!(w, "@HD\tVN:1.6\tSO:unsorted")?;
    for seq in &aligner.index().seqs {
        writeln!(w, "@SQ\tSN:{}\tLN:{}", seq.name, seq.len)?;
    }
    Ok(())
}

/// Write one mapped SAM record.
///
/// SEQ is emitted in **reference orientation**: for a reverse-strand mapping the query
/// is reverse-complemented, as the SAM specification requires. rammap returns CIGAR and
/// MD in reference orientation, so writing the forward query alongside them would
/// describe an alignment that does not exist — and `sampling::thermo_sim` would then
/// compute ΔG by pairing bases that were never opposite each other.
///
/// Soft clips are emitted so the CIGAR's query length matches SEQ. rammap's `cigar_str`
/// covers only the aligned block; without the clips a record whose alignment starts
/// past the first base is invalid SAM, and every downstream query index is shifted.
fn write_sam_mapped<W: Write>(
    w: &mut W,
    qname: &str,
    qseq: &str,
    flag: u16,
    m: &rammap::Mapping,
) -> Result<()> {
    let is_reverse = m.strand == Strand::Reverse;
    let oriented: String = if is_reverse {
        String::from_utf8_lossy(&crate::seq_utils::reverse_complement(qseq.as_bytes())).into_owned()
    } else {
        qseq.to_string()
    };

    let cigar = match m.cigar.as_deref() {
        Some(core) => {
            let qlen = qseq.len();
            // Clip lengths are expressed in reference orientation.
            let (head, tail) = if is_reverse {
                (qlen.saturating_sub(m.query_end), m.query_start)
            } else {
                (m.query_start, qlen.saturating_sub(m.query_end))
            };
            let mut c = String::with_capacity(core.len() + 12);
            if head > 0 {
                c.push_str(&format!("{}S", head));
            }
            c.push_str(core);
            if tail > 0 {
                c.push_str(&format!("{}S", tail));
            }
            c
        }
        None => "*".to_string(),
    };

    write!(
        w,
        "{}\t{}\t{}\t{}\t{}\t{}\t*\t0\t0\t{}\t*\tNM:i:{}",
        qname,
        flag,
        m.target_name,
        m.target_start + 1,
        m.mapq,
        cigar,
        oriented,
        m.edit_distance,
    )?;
    if let Some(ref md) = m.md {
        write!(w, "\tMD:Z:{}", md)?;
    }
    writeln!(w)?;
    Ok(())
}

fn write_sam_unmapped<W: Write>(w: &mut W, qname: &str, qseq: &str) -> Result<()> {
    writeln!(w, "{}\t4\t*\t0\t0\t*\t*\t0\t0\t{}\t*", qname, qseq)?;
    Ok(())
}

fn write_paf_record<W: Write>(
    w: &mut W,
    qname: &str,
    qlen: usize,
    m: &rammap::Mapping,
) -> Result<()> {
    let strand = if m.strand == Strand::Forward { '+' } else { '-' };
    write!(
        w,
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        qname, qlen,
        m.query_start, m.query_end,
        strand,
        m.target_name, m.target_len,
        m.target_start, m.target_end,
        m.matches, m.block_len, m.mapq,
    )?;
    write!(w, "\tNM:i:{}", m.edit_distance)?;
    if let Some(ref cigar) = m.cigar {
        write!(w, "\tcg:Z:{}", cigar)?;
    }
    writeln!(w)?;
    Ok(())
}

// SAM FLAG constants
const FLAG_REVERSE: u16 = 0x10;
const FLAG_SECONDARY: u16 = 0x100;

// ---- public API (signatures unchanged from subprocess version) ----

/// rammap is compiled into the binary — always available.
pub fn check_available() -> Result<()> {
    Ok(())
}

/// Map reads to a reference (SAM output, primary alignments only).
///
/// Replicates `minimap2 -ax <preset> --secondary=no <reference> <reads> [<reads_r2>]`
pub fn map_reads(
    preset: &str,
    reference: &Path,
    reads: &Path,
    reads_r2: Option<&Path>,
    output_sam: &Path,
    log_file: &Path,
) -> Result<()> {
    write_stub_log(log_file, "map_reads")?;

    let aligner = Aligner::from_fasta(path_str(reference)?, to_preset(preset)?)
        .with_context(|| format!("Failed to build alignment index for {}", reference.display()))?;

    let out = File::create(output_sam)
        .with_context(|| format!("Cannot create SAM: {}", output_sam.display()))?;
    let mut w = BufWriter::new(out);
    write_sam_header(&mut w, &aligner)?;

    for (name, seq) in fasta::parse_reads(reads)? {
        map_and_write_primary(&aligner, &mut w, &name, &seq)?;
    }
    if let Some(r2) = reads_r2 {
        for (name, seq) in fasta::parse_reads(r2)? {
            map_and_write_primary(&aligner, &mut w, &name, &seq)?;
        }
    }

    Ok(())
}

fn map_and_write_primary<W: Write>(
    aligner: &Aligner,
    w: &mut W,
    name: &str,
    seq: &str,
) -> Result<()> {
    let result = aligner.map_seq(name, seq.as_bytes());
    match result.mappings.iter().find(|m| m.is_primary && !m.is_supplementary) {
        Some(m) => {
            let flag = if m.strand == Strand::Reverse { FLAG_REVERSE } else { 0 };
            write_sam_mapped(w, name, seq, flag, m)
        }
        None => write_sam_unmapped(w, name, seq),
    }
}

/// Align probes to targets (SAM output with secondary alignments and MD tags).
///
/// Replicates `minimap2 -ax <preset> --secondary=yes -N <max_secondary> -t <threads> <targets> <probes>`
///
/// MD tags are enabled for thermodynamic scoring in `sampling/thermo_sim.rs`.
/// Secondary alignments are kept because a single probe can legitimately
/// tile conserved regions across multiple target sequences.
///
/// `max_secondary`: maximum secondary alignments per query. Use 1000 for normal
/// pangenome coverage; use a larger value (e.g. 10 000) when computing per-target
/// coverage without competition (--all-individual-targets).
pub fn probe_align(
    preset: &str,
    targets: &Path,
    probes: &Path,
    output_sam: &Path,
    log_file: &Path,
    threads: usize,
    max_secondary: usize,
) -> Result<()> {
    write_stub_log(log_file, "probe_align")?;

    let mut aligner = Aligner::from_fasta(path_str(targets)?, to_preset(preset)?)
        .with_context(|| format!("Failed to build alignment index for {}", targets.display()))?;

    // Enable secondary alignments (the SR preset suppresses them by default).
    aligner.options_mut().flags.remove(AlignFlags::NO_PRINT_2ND);
    aligner.options_mut().filtering.best_n = max_secondary as i32;

    // Retune for divergent probes. The `sr` preset is built for near-identical short
    // reads: mismatch_penalty = 8 means a 120-mer at 80% identity scores
    // 96*2 - 24*8 = 0 and is discarded before alignment is attempted. Probes designed
    // against relatives of the sampled organism routinely sit in that range, and
    // because probe density drives per-locus capture efficiency, silently dropping
    // them distorts the simulation.
    //
    // k/w are index-build parameters and not reachable through rammap's public API,
    // but the scoring and chaining thresholds are. Lowering them recovers divergent
    // hits without changing the aligner or paying BLAST's runtime.
    {
        let o = aligner.options_mut();
        o.scoring.mismatch_penalty = 4; // rammap default; `sr` doubles it
        o.chaining.min_cnt = 1; // one seed is enough to anchor a 120-mer
        o.chaining.min_chain_score = 10;
        o.alignment.min_dp_max = 20;
        o.filtering.pri_ratio = 0.1; // keep weaker secondaries for competition
    }
    // Enable MD tags for thermo_sim.rs.
    aligner.output_config_mut().do_md = true;

    let opts = MapOpts { cs: None, md: Some(true) };

    let out = File::create(output_sam)
        .with_context(|| format!("Cannot create SAM: {}", output_sam.display()))?;
    let mut w = BufWriter::new(out);
    write_sam_header(&mut w, &aligner)?;

    let probe_seqs = fasta::parse_fasta(probes)?;

    if threads > 1 {
        use std::sync::Arc;
        use rayon::prelude::*;

        let _ = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global();

        let shared = Arc::new(aligner);
        let results: Vec<(String, String, rammap::MapResult)> = probe_seqs
            .into_par_iter()
            .map(|(name, seq)| {
                let result = shared.map_seq_with(&name, seq.as_bytes(), opts);
                (name, seq, result)
            })
            .collect();

        for (name, seq, result) in results {
            write_probe_alignments(&mut w, &name, &seq, &result.mappings, max_secondary)?;
        }
    } else {
        for (name, seq) in probe_seqs {
            let result = aligner.map_seq_with(&name, seq.as_bytes(), opts);
            write_probe_alignments(&mut w, &name, &seq, &result.mappings, max_secondary)?;
        }
    }

    Ok(())
}

fn write_probe_alignments<W: Write>(
    w: &mut W,
    name: &str,
    seq: &str,
    mappings: &[rammap::Mapping],
    max_secondary: usize,
) -> Result<()> {
    if mappings.is_empty() {
        return write_sam_unmapped(w, name, seq);
    }

    // Rank secondaries before truncating at `max_secondary`. rammap returns mappings in
    // its own order, so taking the first N kept an arbitrary subset rather than the
    // best ones whenever a probe had more hits than the cap — which is exactly the
    // case for conserved regions tiled across many references.
    let mut ordered: Vec<&rammap::Mapping> = mappings.iter().filter(|m| !m.is_supplementary).collect();
    ordered.sort_by(|a, b| {
        // Primary first, then by identity (fewest edits per aligned base), then by
        // alignment length. Ties broken on target name/position so the output is
        // deterministic regardless of the order rammap produced.
        let ident = |m: &rammap::Mapping| {
            if m.block_len > 0 {
                m.edit_distance as f64 / m.block_len as f64
            } else {
                f64::MAX
            }
        };
        b.is_primary
            .cmp(&a.is_primary)
            .then(ident(a).partial_cmp(&ident(b)).unwrap_or(std::cmp::Ordering::Equal))
            .then(b.block_len.cmp(&a.block_len))
            .then(a.target_name.cmp(&b.target_name))
            .then(a.target_start.cmp(&b.target_start))
    });

    let mut secondary_count = 0usize;
    for m in ordered {
        if m.is_primary {
            let flag = if m.strand == Strand::Reverse { FLAG_REVERSE } else { 0 };
            write_sam_mapped(w, name, seq, flag, m)?;
        } else if secondary_count < max_secondary {
            let flag = FLAG_SECONDARY | if m.strand == Strand::Reverse { FLAG_REVERSE } else { 0 };
            write_sam_mapped(w, name, seq, flag, m)?;
            secondary_count += 1;
        }
    }
    Ok(())
}

/// Align pre-loaded probe sequences against a single-target FASTA and return per-position depth.
///
/// No SAM file is written. Probe sequences are passed pre-loaded to avoid re-reading the
/// probe FASTA on every call (used by individual target coverage computation).
///
/// Probe mapping is parallelised across `threads` rayon workers; coverage accumulation
/// is serial (no atomic overhead). Pass `threads = 1` to disable parallelism.
pub fn probe_depth_in_memory(
    preset: &str,
    target_path: &Path,
    probe_seqs: &std::collections::HashMap<String, String>,
    max_secondary: usize,
    threads: usize,
) -> anyhow::Result<RefLengthsAndCoverage> {
    use crate::alignment::cigar::{parse_cigar, CigarOp};
    use std::sync::Arc;
    use rayon::prelude::*;

    let mut aligner = Aligner::from_fasta(path_str(target_path)?, to_preset(preset)?)
        .with_context(|| format!("Failed to build index for {}", target_path.display()))?;
    aligner.options_mut().flags.remove(AlignFlags::NO_PRINT_2ND);
    aligner.options_mut().filtering.best_n = max_secondary as i32;
    aligner.output_config_mut().do_md = true;
    let opts = MapOpts { cs: None, md: Some(true) };

    let ref_lengths: std::collections::HashMap<String, usize> = aligner
        .index()
        .seqs
        .iter()
        .map(|s| (s.name.to_string(), s.len))
        .collect();
    let mut coverage: std::collections::HashMap<String, Vec<u32>> = ref_lengths
        .iter()
        .map(|(n, &l)| (n.clone(), vec![0u32; l]))
        .collect();

    // Parallel phase: map all probes, collecting raw mappings.
    // Aligner is Sync so Arc<Aligner> can be shared across rayon workers.
    let shared = Arc::new(aligner);
    let probe_vec: Vec<(&String, &String)> = probe_seqs.iter().collect();

    let all_mappings: Vec<Vec<rammap::Mapping>> = if threads > 1 {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap_or_else(|_| rayon::ThreadPoolBuilder::new().build().unwrap());
        pool.install(|| {
            probe_vec
                .par_iter()
                .map(|(name, seq)| shared.map_seq_with(name, seq.as_bytes(), opts).mappings)
                .collect()
        })
    } else {
        probe_vec
            .iter()
            .map(|(name, seq)| shared.map_seq_with(name, seq.as_bytes(), opts).mappings)
            .collect()
    };

    // Serial phase: accumulate coverage from collected mappings.
    for mappings in &all_mappings {
        let mut secondary_count = 0usize;
        for m in mappings {
            if m.is_supplementary {
                continue;
            }
            if !m.is_primary {
                if secondary_count >= max_secondary {
                    continue;
                }
                secondary_count += 1;
            }
            let cov = match coverage.get_mut(m.target_name.as_ref()) {
                Some(c) => c,
                None => continue,
            };
            let cigar = m.cigar.as_deref().unwrap_or("*");
            let mut ref_pos = m.target_start;
            for (len, op) in parse_cigar(cigar) {
                match op {
                    CigarOp::Match | CigarOp::Equal | CigarOp::Mismatch => {
                        for _ in 0..len {
                            if ref_pos < cov.len() {
                                cov[ref_pos] += 1;
                            }
                            ref_pos += 1;
                        }
                    }
                    CigarOp::Del | CigarOp::Skip => ref_pos += len as usize,
                    _ => {}
                }
            }
        }
    }

    Ok((ref_lengths, coverage))
}

/// Cross-reactivity alignment (PAF output with all secondary alignments).
///
/// Replicates `minimap2 -x <preset> -c --secondary=yes -N 10000 -p 0.5 <reference> <probes>`
///
/// Uses PAF format for lightweight tabular output with matching base counts.
/// `-p 0.5` is replicated via `filtering.pri_ratio = 0.5` to catch weaker homology.
pub fn xreact_align(
    preset: &str,
    reference: &Path,
    probes: &Path,
    output_paf: &Path,
    log_file: &Path,
) -> Result<()> {
    write_stub_log(log_file, "xreact_align")?;

    let mut aligner = Aligner::from_fasta(path_str(reference)?, to_preset(preset)?)
        .with_context(|| format!("Failed to build alignment index for {}", reference.display()))?;

    // Replicate minimap2 --secondary=yes -N 10000 -p 0.5.
    aligner.options_mut().flags.remove(AlignFlags::NO_PRINT_2ND);
    aligner.options_mut().filtering.best_n = 10000;
    aligner.options_mut().filtering.pri_ratio = 0.5;

    let out = File::create(output_paf)
        .with_context(|| format!("Cannot create PAF: {}", output_paf.display()))?;
    let mut w = BufWriter::new(out);

    for (name, seq) in fasta::parse_fasta(probes)? {
        let qlen = seq.len();
        let result = aligner.map_seq(&name, seq.as_bytes());
        for m in &result.mappings {
            if !m.is_supplementary {
                write_paf_record(&mut w, &name, qlen, m)?;
            }
        }
    }

    Ok(())
}

/// Host read filtering alignment (SAM output).
///
/// Replicates `minimap2 -ax <preset> <host> <reads> [<reads_r2>]`
pub fn host_align(
    preset: &str,
    host: &Path,
    reads: &Path,
    reads_r2: Option<&Path>,
    output_sam: &Path,
    log_file: &Path,
) -> Result<()> {
    write_stub_log(log_file, "host_align")?;

    let aligner = Aligner::from_fasta(path_str(host)?, to_preset(preset)?)
        .with_context(|| format!("Failed to build alignment index for {}", host.display()))?;

    let out = File::create(output_sam)
        .with_context(|| format!("Cannot create SAM: {}", output_sam.display()))?;
    let mut w = BufWriter::new(out);
    write_sam_header(&mut w, &aligner)?;

    for (name, seq) in fasta::parse_reads(reads)? {
        map_and_write_primary(&aligner, &mut w, &name, &seq)?;
    }
    if let Some(r2) = reads_r2 {
        for (name, seq) in fasta::parse_reads(r2)? {
            map_and_write_primary(&aligner, &mut w, &name, &seq)?;
        }
    }

    Ok(())
}
