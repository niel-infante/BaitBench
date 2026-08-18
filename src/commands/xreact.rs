use anyhow::{bail, Context, Result};
use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use crate::alignment::{blast_tab, paf};
use crate::cleanup;
use crate::cli::{Aligner, ReportMode};
use crate::commands::report::{rmd_output_path, substitute_rmd_params};
use crate::external::{blastn, minimap2, rscript};
use crate::fasta;
use crate::io_utils::{abs_path_str, prefixed_join};

pub struct XreactArgs<'a> {
    pub probes: &'a Path,
    pub against: &'a [PathBuf],
    pub self_mode: bool,
    pub threshold: f64,
    pub aligner: Aligner,
    pub outdir: &'a Path,
    pub output_prefix: &'a str,
    pub minimap_preset: &'a str,
    pub threads: usize,
    pub report: ReportMode,
    pub cleanup: bool,
}

/// Alignment hit reduced to the fields cross-reactivity scoring needs,
/// common to both the minimap2 (PAF) and BLAST (tabular) backends.
struct AlignHit {
    query_name: String,
    target_name: String,
    query_length: u32,
    query_start: u32,
    query_end: u32,
    matching_bases: u32,
    block_length: u32,
}

impl From<paf::PafRecord> for AlignHit {
    fn from(r: paf::PafRecord) -> Self {
        AlignHit {
            query_name: r.query_name,
            target_name: r.target_name,
            query_length: r.query_length,
            query_start: r.query_start,
            query_end: r.query_end,
            matching_bases: r.matching_bases,
            block_length: r.block_length,
        }
    }
}

impl From<blast_tab::BlastHit> for AlignHit {
    fn from(h: blast_tab::BlastHit) -> Self {
        AlignHit {
            query_name: h.query_name,
            target_name: h.target_name,
            query_length: h.query_length,
            query_start: h.query_start,
            query_end: h.query_end,
            matching_bases: h.matching_bases,
            block_length: h.alignment_length,
        }
    }
}

/// Run the selected aligner and return hits reduced to `AlignHit`, cleaning
/// up the intermediate alignment file afterward.
#[allow(clippy::too_many_arguments)]
fn align(
    aligner: Aligner,
    minimap_preset: &str,
    threads: usize,
    reference: &Path,
    query: &Path,
    outdir: &Path,
    pfx: &str,
    label: &str,
) -> Result<Vec<AlignHit>> {
    let log_path = prefixed_join(outdir, pfx, &format!("{label}.log"));

    let hits = match aligner {
        Aligner::Minimap2 => {
            let paf_path = prefixed_join(outdir, pfx, &format!("{label}.paf"));
            minimap2::xreact_align(minimap_preset, reference, query, &paf_path, &log_path)?;
            let records = paf::parse_paf_records(&paf_path)?;
            log::info!("PAF records ({label}): {}", records.len());
            let hits: Vec<AlignHit> = records.into_iter().map(AlignHit::from).collect();
            if let Err(e) = fs::remove_file(&paf_path) {
                log::debug!("Could not remove {}: {}", paf_path.display(), e);
            }
            hits
        }
        Aligner::Blast => {
            let tsv_path = prefixed_join(outdir, pfx, &format!("{label}.blast.tsv"));
            blastn::xreact_align(reference, query, &tsv_path, &log_path, threads)?;
            let records = blast_tab::parse_blast_hits(&tsv_path)?;
            log::info!("BLAST hits ({label}): {}", records.len());
            let hits: Vec<AlignHit> = records.into_iter().map(AlignHit::from).collect();
            if let Err(e) = fs::remove_file(&tsv_path) {
                log::debug!("Could not remove {}: {}", tsv_path.display(), e);
            }
            hits
        }
    };

    Ok(hits)
}

struct HitRecord {
    probe_id: String,
    target_id: String,
    homology_pct: f64,
    identity_pct: f64,
    query_coverage_pct: f64,
    matching_bases: u32,
    alignment_length: u32,
    probe_length: u32,
    mode: &'static str,
}

struct SummaryRecord {
    probe_id: String,
    mode: &'static str,
    max_homology_pct: f64,
    best_hit: String,
    num_hits: usize,
}

pub fn execute(args: &XreactArgs) -> Result<()> {
    if args.against.is_empty() && !args.self_mode {
        bail!("At least one of --against or --self must be specified");
    }
    if !args.probes.exists() {
        bail!("Probes file not found: {}", args.probes.display());
    }
    for path in args.against {
        if !path.exists() {
            bail!("Reference file not found: {}", path.display());
        }
    }

    fs::create_dir_all(args.outdir)?;
    match args.aligner {
        Aligner::Minimap2 => minimap2::check_available()?,
        Aligner::Blast => blastn::check_available()?,
    }

    log::info!("=============================================");
    log::info!("BaitBench - Cross-Reactivity Analysis");
    log::info!("=============================================");
    log::info!("Probes    : {}", args.probes.display());
    if !args.against.is_empty() {
        for p in args.against {
            log::info!("Against   : {}", p.display());
        }
    }
    log::info!("Self mode : {}", args.self_mode);
    log::info!("Threshold : {:.1}%", args.threshold);
    match args.aligner {
        Aligner::Minimap2 => log::info!("Aligner   : minimap2 (preset {})", args.minimap_preset),
        Aligner::Blast => log::info!("Aligner   : blast ({} threads)", args.threads),
    }
    log::info!("Output    : {}", args.outdir.display());

    let pfx = args.output_prefix;
    let n_probes = fasta::count_sequences(args.probes)?;
    let probe_ids = fasta::parse_fasta_ids(args.probes)?;
    log::info!("Probe sequences: {}", n_probes);

    let mut all_hits: Vec<HitRecord> = Vec::new();

    // Probe-to-genome mode
    if !args.against.is_empty() {
        let reference_path = if args.against.len() == 1 {
            args.against[0].clone()
        } else {
            let combined = prefixed_join(args.outdir, pfx, "against_combined.fa");
            let refs: Vec<&Path> = args.against.iter().map(|p| p.as_path()).collect();
            fasta::concatenate_fastas(&refs, &combined)?;
            combined
        };

        let n_refs = fasta::count_sequences(&reference_path)?;
        log::info!("Reference sequences: {}", n_refs);

        log::info!("Aligning probes against reference...");
        let records = align(
            args.aligner,
            args.minimap_preset,
            args.threads,
            &reference_path,
            args.probes,
            args.outdir,
            pfx,
            "against",
        )?;

        for rec in &records {
            if rec.query_length == 0 {
                continue;
            }
            let homology = rec.matching_bases as f64 / rec.query_length as f64 * 100.0;
            if homology >= args.threshold {
                let identity = if rec.block_length > 0 {
                    rec.matching_bases as f64 / rec.block_length as f64 * 100.0
                } else {
                    0.0
                };
                let query_cov =
                    (rec.query_end - rec.query_start) as f64 / rec.query_length as f64 * 100.0;
                all_hits.push(HitRecord {
                    probe_id: rec.query_name.clone(),
                    target_id: rec.target_name.clone(),
                    homology_pct: homology,
                    identity_pct: identity,
                    query_coverage_pct: query_cov,
                    matching_bases: rec.matching_bases,
                    alignment_length: rec.block_length,
                    probe_length: rec.query_length,
                    mode: "against",
                });
            }
        }

        if args.against.len() > 1 {
            let combined = prefixed_join(args.outdir, pfx, "against_combined.fa");
            if let Err(e) = fs::remove_file(&combined) {
                log::debug!("Could not remove {}: {}", combined.display(), e);
            }
        }
    }

    // Self mode
    if args.self_mode {
        log::info!("Aligning probes against themselves...");
        let records = align(
            args.aligner,
            args.minimap_preset,
            args.threads,
            args.probes,
            args.probes,
            args.outdir,
            pfx,
            "self",
        )?;
        log::info!("Alignment hits (self, including self-hits): {}", records.len());

        for rec in &records {
            // Skip self-hits
            if rec.query_name == rec.target_name {
                continue;
            }
            if rec.query_length == 0 {
                continue;
            }
            let homology = rec.matching_bases as f64 / rec.query_length as f64 * 100.0;
            if homology >= args.threshold {
                let identity = if rec.block_length > 0 {
                    rec.matching_bases as f64 / rec.block_length as f64 * 100.0
                } else {
                    0.0
                };
                let query_cov =
                    (rec.query_end - rec.query_start) as f64 / rec.query_length as f64 * 100.0;
                all_hits.push(HitRecord {
                    probe_id: rec.query_name.clone(),
                    target_id: rec.target_name.clone(),
                    homology_pct: homology,
                    identity_pct: identity,
                    query_coverage_pct: query_cov,
                    matching_bases: rec.matching_bases,
                    alignment_length: rec.block_length,
                    probe_length: rec.query_length,
                    mode: "self",
                });
            }
        }
    }

    // Write run_params.tsv
    let params_path = prefixed_join(args.outdir, pfx, "run_params.tsv");
    write_run_params(&params_path, args)?;

    // Write hits.tsv
    let hits_path = prefixed_join(args.outdir, pfx, "hits.tsv");
    write_hits_tsv(&hits_path, &all_hits)?;
    log::info!(
        "Hits above {:.1}% threshold: {} (written to {})",
        args.threshold,
        all_hits.len(),
        hits_path.display()
    );

    // Build and write summary.tsv
    let summaries = build_summaries(
        &probe_ids,
        &all_hits,
        args.self_mode,
        !args.against.is_empty(),
    );
    let summary_path = prefixed_join(args.outdir, pfx, "summary.tsv");
    write_summary_tsv(&summary_path, &summaries)?;

    // Console summary
    if !args.against.is_empty() {
        let against_flagged = summaries
            .iter()
            .filter(|s| s.mode == "against" && s.num_hits > 0)
            .count();
        log::info!(
            "Probes with cross-reactive genome hits: {}/{}",
            against_flagged,
            n_probes
        );
    }
    if args.self_mode {
        let self_flagged = summaries
            .iter()
            .filter(|s| s.mode == "self" && s.num_hits > 0)
            .count();
        log::info!(
            "Probes with cross-reactive probe hits: {}/{}",
            self_flagged,
            n_probes
        );
    }

    // Report generation
    match args.report {
        ReportMode::None => {
            log::info!("Skipping report generation (--report none)");
        }
        ReportMode::Full => {
            if rscript::check_available() {
                let report_path = prefixed_join(args.outdir, pfx, "xreact_report.html");
                log::info!("Generating cross-reactivity report...");
                match generate_xreact_report(
                    &hits_path,
                    &summary_path,
                    &params_path,
                    args.threshold,
                    &report_path,
                ) {
                    Ok(()) => log::info!("Report generated: {}", report_path.display()),
                    Err(e) => log::warn!("Report generation failed (non-fatal): {}", e),
                }
            } else {
                log::warn!("Rscript not found -- skipping HTML report.");
            }
        }
        ReportMode::Rmd => {
            let report_path = prefixed_join(args.outdir, pfx, "xreact_report.html");
            log::info!("Generating cross-reactivity RMarkdown file...");
            match write_xreact_rmd(
                &hits_path,
                &summary_path,
                &params_path,
                args.threshold,
                &report_path,
            ) {
                Ok(()) => {}
                Err(e) => log::warn!("RMarkdown generation failed (non-fatal): {}", e),
            }
        }
        ReportMode::BothR => {
            let report_path = prefixed_join(args.outdir, pfx, "xreact_report.html");
            log::info!("Generating cross-reactivity RMarkdown file...");
            match write_xreact_rmd(
                &hits_path,
                &summary_path,
                &params_path,
                args.threshold,
                &report_path,
            ) {
                Ok(()) => {}
                Err(e) => log::warn!("RMarkdown generation failed (non-fatal): {}", e),
            }
            if rscript::check_available() {
                log::info!("Generating cross-reactivity HTML report...");
                match generate_xreact_report(
                    &hits_path,
                    &summary_path,
                    &params_path,
                    args.threshold,
                    &report_path,
                ) {
                    Ok(()) => log::info!("Report generated: {}", report_path.display()),
                    Err(e) => log::warn!("Report generation failed (non-fatal): {}", e),
                }
            } else {
                log::warn!("Rscript not found — skipping HTML report (Rmd still written).");
            }
        }
    }

    // Cleanup intermediate files if requested
    if args.cleanup {
        log::info!("Cleaning up intermediate files...");
        let cleanup_names: Vec<String> = ["against.log", "self.log"]
            .iter()
            .map(|f| format!("{}{}", pfx, f))
            .collect();
        let cleanup_refs: Vec<&str> = cleanup_names.iter().map(|s| s.as_str()).collect();
        cleanup::cleanup_files(args.outdir, &cleanup_refs);
    }

    log::info!("=============================================");
    log::info!("Cross-reactivity analysis complete!");
    log::info!("Results in {}", args.outdir.display());
    log::info!("=============================================");

    Ok(())
}

fn write_run_params(path: &Path, args: &XreactArgs) -> Result<()> {
    let file =
        File::create(path).with_context(|| format!("Cannot create params file: {}", path.display()))?;
    let mut w = BufWriter::new(file);

    writeln!(w, "parameter\tflag\tvalue")?;
    writeln!(w, "probes\t--probes\t{}", args.probes.display())?;
    for p in args.against {
        writeln!(w, "against\t--against\t{}", p.display())?;
    }
    writeln!(w, "self\t--self\t{}", args.self_mode)?;
    writeln!(w, "threshold\t--threshold\t{:.1}", args.threshold)?;
    let aligner_str = match args.aligner {
        Aligner::Minimap2 => "minimap2",
        Aligner::Blast => "blast",
    };
    writeln!(w, "aligner\t--aligner\t{}", aligner_str)?;
    match args.aligner {
        Aligner::Minimap2 => {
            writeln!(w, "minimap_preset\t--minimap-preset\t{}", args.minimap_preset)?;
        }
        Aligner::Blast => {
            writeln!(w, "threads\t--threads\t{}", args.threads)?;
        }
    }
    writeln!(w, "outdir\t-o\t{}", args.outdir.display())?;

    w.flush()?;
    Ok(())
}

fn generate_xreact_report(
    hits_path: &Path,
    summary_path: &Path,
    params_path: &Path,
    threshold: f64,
    output_path: &Path,
) -> Result<()> {
    let r_dir = rscript::find_r_dir()
        .ok_or_else(|| anyhow::anyhow!("Cannot find R scripts directory."))?;

    let script = r_dir.join("xreact.R");
    if !script.exists() {
        bail!("Xreact R script not found: {}", script.display());
    }

    let hits_str = abs_path_str(hits_path)?;
    let summary_str = abs_path_str(summary_path)?;
    let params_str = abs_path_str(params_path)?;
    let output_str = abs_path_str(output_path)?;
    let threshold_str = format!("{:.1}", threshold);

    rscript::run_rscript(
        &script,
        &[
            "--hits",
            &hits_str,
            "--summary",
            &summary_str,
            "--params",
            &params_str,
            "--threshold",
            &threshold_str,
            "--output",
            &output_str,
        ],
    )
}

fn write_xreact_rmd(
    hits_path: &Path,
    summary_path: &Path,
    params_path: &Path,
    threshold: f64,
    output_path: &Path,
) -> Result<()> {
    let r_dir = rscript::find_r_dir()
        .ok_or_else(|| anyhow::anyhow!("Cannot find R scripts directory."))?;

    let rmd_template = r_dir.join("xreact.Rmd");
    if !rmd_template.exists() {
        bail!("RMarkdown template not found: {}", rmd_template.display());
    }

    let hits_abs = abs_path_str(hits_path)?;
    let summary_abs = abs_path_str(summary_path)?;
    let params_abs = abs_path_str(params_path)?;
    let threshold_str = format!("{:.1}", threshold);

    let params = vec![
        ("hits_file", hits_abs.as_str()),
        ("summary_file", summary_abs.as_str()),
        ("params_file", params_abs.as_str()),
        ("threshold", threshold_str.as_str()),
    ];

    let template_content = std::fs::read_to_string(&rmd_template)
        .with_context(|| format!("Failed to read template: {}", rmd_template.display()))?;

    let output_content = substitute_rmd_params(&template_content, &params);

    let rmd_path = rmd_output_path(output_path);
    std::fs::write(&rmd_path, output_content)
        .with_context(|| format!("Failed to write RMarkdown file: {}", rmd_path.display()))?;

    log::info!("RMarkdown file written: {}", rmd_path.display());
    log::info!(
        "Edit and render with: Rscript -e 'rmarkdown::render(\"{}\")'",
        rmd_path.display()
    );
    Ok(())
}

fn write_hits_tsv(path: &Path, hits: &[HitRecord]) -> Result<()> {
    let file =
        File::create(path).with_context(|| format!("Cannot create hits file: {}", path.display()))?;
    let mut w = BufWriter::new(file);

    writeln!(
        w,
        "probe_id\ttarget_id\thomology_pct\tidentity_pct\tquery_coverage_pct\tmatching_bases\talignment_length\tprobe_length\tmode"
    )?;

    // Sort by probe_id, then mode, then descending homology
    let mut sorted: Vec<usize> = (0..hits.len()).collect();
    sorted.sort_by(|&a, &b| {
        hits[a]
            .probe_id
            .cmp(&hits[b].probe_id)
            .then(hits[a].mode.cmp(hits[b].mode))
            .then(
                hits[b]
                    .homology_pct
                    .partial_cmp(&hits[a].homology_pct)
                    .unwrap_or(std::cmp::Ordering::Equal),
            )
    });

    for &i in &sorted {
        let h = &hits[i];
        writeln!(
            w,
            "{}\t{}\t{:.1}\t{:.1}\t{:.1}\t{}\t{}\t{}\t{}",
            h.probe_id,
            h.target_id,
            h.homology_pct,
            h.identity_pct,
            h.query_coverage_pct,
            h.matching_bases,
            h.alignment_length,
            h.probe_length,
            h.mode,
        )?;
    }

    w.flush()?;
    Ok(())
}

fn build_summaries(
    probe_ids: &[String],
    hits: &[HitRecord],
    self_mode: bool,
    against_mode: bool,
) -> Vec<SummaryRecord> {
    // Group hits by (probe_id, mode) → track max homology + count
    let mut map: HashMap<(&str, &str), (f64, String, usize)> = HashMap::new();

    for h in hits {
        let entry = map
            .entry((&h.probe_id, h.mode))
            .or_insert((0.0, "NA".to_string(), 0));
        entry.2 += 1;
        if h.homology_pct > entry.0 {
            entry.0 = h.homology_pct;
            entry.1 = h.target_id.clone();
        }
    }

    let mut summaries = Vec::new();

    let modes: Vec<&str> = {
        let mut m = Vec::new();
        if against_mode {
            m.push("against");
        }
        if self_mode {
            m.push("self");
        }
        m
    };

    for probe_id in probe_ids {
        for &mode in &modes {
            let (max_hom, best_hit, num_hits) = map
                .get(&(probe_id.as_str(), mode))
                .cloned()
                .unwrap_or((0.0, "NA".to_string(), 0));
            summaries.push(SummaryRecord {
                probe_id: probe_id.clone(),
                mode,
                max_homology_pct: max_hom,
                best_hit,
                num_hits,
            });
        }
    }

    summaries
}

fn write_summary_tsv(path: &Path, summaries: &[SummaryRecord]) -> Result<()> {
    let file = File::create(path)
        .with_context(|| format!("Cannot create summary file: {}", path.display()))?;
    let mut w = BufWriter::new(file);

    writeln!(
        w,
        "probe_id\tmode\tmax_homology_pct\tbest_hit\tnum_hits_above_threshold"
    )?;

    for s in summaries {
        writeln!(
            w,
            "{}\t{}\t{:.1}\t{}\t{}",
            s.probe_id, s.mode, s.max_homology_pct, s.best_hit, s.num_hits,
        )?;
    }

    w.flush()?;
    Ok(())
}
