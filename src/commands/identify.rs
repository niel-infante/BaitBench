use anyhow::{bail, Context, Result};
use std::collections::{HashMap, HashSet};
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use crate::io_utils;
use crate::io_utils::prefixed_join;
use crate::target_similarity::{self, SimilarityContext};

pub struct IdentifyArgs<'a> {
    pub detected_detail: &'a Path,
    pub sample_target_map: &'a Path,
    pub target_similarity_file: Option<&'a Path>,
    pub targets_fasta: Option<&'a Path>,
    pub identity_threshold: f64,
    pub minimap_preset: &'a str,
    pub min_unique_targets: usize,
    pub outdir: &'a Path,
    pub output_prefix: &'a str,
}

/// Species call result.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum SpeciesCall {
    Present,
    Absent,
    Ambiguous,
}

impl std::fmt::Display for SpeciesCall {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            SpeciesCall::Present => write!(f, "PRESENT"),
            SpeciesCall::Absent => write!(f, "ABSENT"),
            SpeciesCall::Ambiguous => write!(f, "AMBIGUOUS"),
        }
    }
}

/// Full evidence record for a species call.
pub struct SpeciesCallResult {
    pub species_id: String,
    pub call: SpeciesCall,
    pub total_targets: usize,
    pub unique_targets: usize,
    pub shared_targets: usize,
    pub unique_detected: usize,
    pub shared_detected: usize,
    pub total_detected: usize,
    pub total_reads: usize,
    pub explained_by: Vec<String>,
    pub reason: String,
}

/// Parsed row from detected_detail.tsv.
pub struct DetailRow {
    reference_id: String,
    detected: bool,
    reads_assigned: usize,
}

pub fn execute(args: &IdentifyArgs) -> Result<()> {
    if !args.detected_detail.exists() {
        bail!(
            "Detected detail file not found: {}",
            args.detected_detail.display()
        );
    }
    if !args.sample_target_map.exists() {
        bail!(
            "Sample-target-map not found: {}",
            args.sample_target_map.display()
        );
    }

    if !(0.0..=100.0).contains(&args.identity_threshold) {
        bail!(
            "--identity-threshold must be in [0, 100], got {}",
            args.identity_threshold
        );
    }

    fs::create_dir_all(args.outdir)?;

    log::info!("=============================================");
    log::info!("BaitBench - Species Identification");
    log::info!("=============================================");

    // Parse genome-to-target mapping
    let genome_to_targets = io_utils::parse_sample_target_map(args.sample_target_map)?;

    // Get or compute target similarity
    let similarities = if let Some(sim_path) = args.target_similarity_file {
        if !sim_path.exists() {
            bail!(
                "Target similarity file not found: {}",
                sim_path.display()
            );
        }
        log::info!(
            "Loading pre-computed target similarity from {}",
            sim_path.display()
        );
        target_similarity::parse_target_similarity(sim_path)?
    } else if let Some(targets) = args.targets_fasta {
        if !targets.exists() {
            bail!("Targets FASTA not found: {}", targets.display());
        }
        log::info!("Computing target similarity on-the-fly...");
        target_similarity::compute_target_similarity(
            targets,
            args.minimap_preset,
            args.identity_threshold,
            args.outdir,
        )?
    } else {
        bail!(
            "Either --target-similarity or --targets must be provided \
             for species identification."
        );
    };

    // Build similarity context
    let ctx = target_similarity::build_similarity_context(&similarities, &genome_to_targets);

    // Parse detected_detail.tsv
    let detail_rows = parse_detected_detail(args.detected_detail)?;

    // Run species calling algorithm
    let calls = call_species(&ctx, &detail_rows, args.min_unique_targets);

    // Write results
    let calls_tsv_path = prefixed_join(args.outdir, args.output_prefix, "species_calls.tsv");
    write_species_calls_tsv(&calls, &calls_tsv_path)?;

    let calls_json_path = prefixed_join(args.outdir, args.output_prefix, "species_calls.json");
    write_species_calls_json(&calls, &calls_json_path)?;

    // Console summary
    let n_present = calls.iter().filter(|c| c.call == SpeciesCall::Present).count();
    let n_absent = calls.iter().filter(|c| c.call == SpeciesCall::Absent).count();
    let n_ambiguous = calls
        .iter()
        .filter(|c| c.call == SpeciesCall::Ambiguous)
        .count();

    log::info!("Species calls: {} PRESENT, {} ABSENT, {} AMBIGUOUS", n_present, n_absent, n_ambiguous);

    for c in &calls {
        let marker = match c.call {
            SpeciesCall::Present => "+",
            SpeciesCall::Absent => "-",
            SpeciesCall::Ambiguous => "?",
        };
        log::info!(
            "  [{}] {} — {} ({} unique/{} shared detected, {} reads)",
            marker,
            c.species_id,
            c.reason,
            c.unique_detected,
            c.shared_detected,
            c.total_reads,
        );
    }

    log::info!("=============================================");
    log::info!("Species identification complete!");
    log::info!("Results in {}", args.outdir.display());
    log::info!("=============================================");

    Ok(())
}

/// Parse detected_detail.tsv to extract reference_id, detected, reads_assigned.
fn parse_detected_detail(path: &Path) -> Result<Vec<DetailRow>> {
    let file = File::open(path)
        .with_context(|| format!("Cannot open detected detail: {}", path.display()))?;
    let reader = BufReader::new(file);
    let mut rows = Vec::new();
    let mut header_idx: HashMap<String, usize> = HashMap::new();

    for (i, line) in reader.lines().enumerate() {
        let line = line?;
        let line = line.trim().to_string();
        if line.is_empty() {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if i == 0 {
            for (j, &f) in fields.iter().enumerate() {
                header_idx.insert(f.to_string(), j);
            }
            continue;
        }

        let ref_id_col = header_idx.get("reference_id").copied().unwrap_or(0);
        let detected_col = header_idx.get("detected").copied().unwrap_or(3);
        let reads_col = header_idx.get("reads_assigned").copied().unwrap_or(6);

        let reference_id = fields.get(ref_id_col).unwrap_or(&"").to_string();
        let detected = fields.get(detected_col).unwrap_or(&"false") == &"true";
        let reads_assigned: usize = fields
            .get(reads_col)
            .unwrap_or(&"0")
            .parse()
            .unwrap_or(0);

        rows.push(DetailRow {
            reference_id,
            detected,
            reads_assigned,
        });
    }

    Ok(rows)
}

/// Core species calling algorithm.
///
/// Uses a weighted unique marker approach with cross-reactivity explanation:
/// 1. Classify each target as unique or shared (from SimilarityContext)
/// 2. Collect detection evidence from detail rows
/// 3. Sort species by evidence strength (unique detections, then total reads)
/// 4. Process in order: call PRESENT/ABSENT/AMBIGUOUS with explanation
pub fn call_species(
    ctx: &SimilarityContext,
    detail_rows: &[DetailRow],
    min_unique_targets: usize,
) -> Vec<SpeciesCallResult> {
    // Build detection lookup: target_id -> (detected, reads)
    let detection: HashMap<&str, (bool, usize)> = detail_rows
        .iter()
        .map(|r| (r.reference_id.as_str(), (r.detected, r.reads_assigned)))
        .collect();

    // Build per-species evidence
    let mut evidence: Vec<SpeciesCallResult> = Vec::new();

    for (species, targets) in &ctx.species_targets {
        let mut unique_detected = 0usize;
        let mut shared_detected = 0usize;
        let mut unique_reads = 0usize;
        let mut shared_reads = 0usize;
        let mut unique_count = 0usize;
        let mut shared_count = 0usize;

        for t in targets {
            let is_unique = *ctx.target_is_unique.get(t.as_str()).unwrap_or(&true);
            let (det, reads) = detection
                .get(t.as_str())
                .copied()
                .unwrap_or((false, 0));

            if is_unique {
                unique_count += 1;
                if det {
                    unique_detected += 1;
                    unique_reads += reads;
                }
            } else {
                shared_count += 1;
                if det {
                    shared_detected += 1;
                    shared_reads += reads;
                }
            }
        }

        evidence.push(SpeciesCallResult {
            species_id: species.clone(),
            call: SpeciesCall::Ambiguous, // placeholder
            total_targets: targets.len(),
            unique_targets: unique_count,
            shared_targets: shared_count,
            unique_detected,
            shared_detected,
            total_detected: unique_detected + shared_detected,
            total_reads: unique_reads + shared_reads,
            explained_by: Vec::new(),
            reason: String::new(),
        });
    }

    // Sort by evidence strength: unique_detected DESC, then total_reads DESC
    evidence.sort_by(|a, b| {
        b.unique_detected
            .cmp(&a.unique_detected)
            .then(b.total_reads.cmp(&a.total_reads))
            .then(a.species_id.cmp(&b.species_id))
    });

    // Track which targets are "explained" by species already called present
    let mut explained_targets: HashSet<String> = HashSet::new();
    let mut called_present: Vec<String> = Vec::new();

    for ev in &mut evidence {
        let targets = match ctx.species_targets.get(&ev.species_id) {
            Some(t) => t,
            None => continue,
        };

        if ev.total_detected == 0 {
            // No targets detected at all
            ev.call = SpeciesCall::Absent;
            ev.reason = "no_detections".to_string();
            continue;
        }

        // Count unexplained detections for this species
        let unexplained_unique: Vec<&String> = targets
            .iter()
            .filter(|t| {
                *ctx.target_is_unique.get(t.as_str()).unwrap_or(&true)
                    && detection
                        .get(t.as_str())
                        .map_or(false, |(det, _)| *det)
                    && !explained_targets.contains(t.as_str())
            })
            .collect();

        let unexplained_shared: Vec<&String> = targets
            .iter()
            .filter(|t| {
                !*ctx.target_is_unique.get(t.as_str()).unwrap_or(&true)
                    && detection
                        .get(t.as_str())
                        .map_or(false, |(det, _)| *det)
                    && !explained_targets.contains(t.as_str())
            })
            .collect();

        if unexplained_unique.len() >= min_unique_targets {
            // Strong evidence: unique targets detected
            ev.call = SpeciesCall::Present;
            ev.reason = "unique_markers_detected".to_string();

            // Mark all this species' targets as explained
            for t in targets {
                explained_targets.insert(t.clone());
                // Also mark similar targets from other species as explained
                if let Some(similar) = ctx.cross_species_similar.get(t) {
                    for st in similar {
                        explained_targets.insert(st.clone());
                    }
                }
            }
            called_present.push(ev.species_id.clone());
        } else if ev.unique_targets == 0
            && unexplained_unique.is_empty()
            && !unexplained_shared.is_empty()
        {
            // No unique markers at all — can't confirm
            ev.call = SpeciesCall::Ambiguous;
            ev.reason = "no_unique_markers".to_string();
        } else if unexplained_unique.is_empty() && unexplained_shared.is_empty() {
            // All detected targets are explained by species already called present
            ev.call = SpeciesCall::Absent;
            ev.reason = "cross_reactivity_explained".to_string();

            // Find which present species explain the shared hits
            let mut explainers: Vec<String> = Vec::new();
            for t in targets {
                let (det, _) = detection.get(t.as_str()).copied().unwrap_or((false, 0));
                if det {
                    if let Some(similar) = ctx.cross_species_similar.get(t) {
                        for st in similar {
                            if let Some(species_list) = ctx.target_to_species.get(st) {
                                for s in species_list {
                                    if called_present.contains(s) && !explainers.contains(s) {
                                        explainers.push(s.clone());
                                    }
                                }
                            }
                        }
                    }
                }
            }
            ev.explained_by = explainers;
        } else {
            // Some unique targets exist but not enough detected
            ev.call = SpeciesCall::Ambiguous;
            ev.reason = "insufficient_unique_evidence".to_string();
        }
    }

    // Re-sort alphabetically for output
    evidence.sort_by(|a, b| a.species_id.cmp(&b.species_id));

    evidence
}

fn write_species_calls_tsv(calls: &[SpeciesCallResult], path: &Path) -> Result<()> {
    let file = File::create(path)
        .with_context(|| format!("Cannot create species calls file: {}", path.display()))?;
    let mut w = BufWriter::new(file);

    writeln!(
        w,
        "species_id\tcall\ttotal_targets\tunique_targets\tshared_targets\tunique_detected\tshared_detected\ttotal_detected\ttotal_reads\texplained_by\treason"
    )?;

    for c in calls {
        let explained = if c.explained_by.is_empty() {
            "-".to_string()
        } else {
            c.explained_by.join(",")
        };
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            c.species_id,
            c.call,
            c.total_targets,
            c.unique_targets,
            c.shared_targets,
            c.unique_detected,
            c.shared_detected,
            c.total_detected,
            c.total_reads,
            explained,
            c.reason,
        )?;
    }

    w.flush()?;
    Ok(())
}

fn write_species_calls_json(calls: &[SpeciesCallResult], path: &Path) -> Result<()> {
    let file = File::create(path)
        .with_context(|| format!("Cannot create species calls JSON: {}", path.display()))?;
    let mut w = BufWriter::new(file);

    writeln!(w, "{{")?;
    writeln!(w, "  \"species_calls\": [")?;

    for (i, c) in calls.iter().enumerate() {
        let comma = if i < calls.len() - 1 { "," } else { "" };
        let explained_json: Vec<String> = c
            .explained_by
            .iter()
            .map(|s| format!("\"{}\"", s))
            .collect();
        writeln!(w, "    {{")?;
        writeln!(w, "      \"species_id\": \"{}\",", c.species_id)?;
        writeln!(w, "      \"call\": \"{}\",", c.call)?;
        writeln!(w, "      \"total_targets\": {},", c.total_targets)?;
        writeln!(w, "      \"unique_targets\": {},", c.unique_targets)?;
        writeln!(w, "      \"shared_targets\": {},", c.shared_targets)?;
        writeln!(w, "      \"unique_detected\": {},", c.unique_detected)?;
        writeln!(w, "      \"shared_detected\": {},", c.shared_detected)?;
        writeln!(w, "      \"total_detected\": {},", c.total_detected)?;
        writeln!(w, "      \"total_reads\": {},", c.total_reads)?;
        writeln!(
            w,
            "      \"explained_by\": [{}],",
            explained_json.join(", ")
        )?;
        writeln!(w, "      \"reason\": \"{}\"", c.reason)?;
        writeln!(w, "    }}{}", comma)?;
    }

    writeln!(w, "  ]")?;
    writeln!(w, "}}")?;

    w.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_row(id: &str, detected: bool, reads: usize) -> DetailRow {
        DetailRow { reference_id: id.to_string(), detected, reads_assigned: reads }
    }

    fn unique_ctx(species: &str, target: &str) -> SimilarityContext {
        SimilarityContext {
            cross_species_similar: HashMap::new(),
            target_is_unique: [(target.to_string(), true)].into_iter().collect(),
            species_targets: [(species.to_string(), vec![target.to_string()])].into_iter().collect(),
            target_to_species: [(target.to_string(), vec![species.to_string()])].into_iter().collect(),
        }
    }

    #[test]
    fn call_absent_no_detections() {
        let ctx = unique_ctx("SpA", "t1");
        let calls = call_species(&ctx, &[make_row("t1", false, 0)], 1);
        assert_eq!(calls.len(), 1);
        assert!(matches!(calls[0].call, SpeciesCall::Absent));
        assert_eq!(calls[0].reason, "no_detections");
    }

    #[test]
    fn call_present_unique_detected() {
        let ctx = unique_ctx("SpA", "t1");
        let calls = call_species(&ctx, &[make_row("t1", true, 50)], 1);
        assert_eq!(calls.len(), 1);
        assert!(matches!(calls[0].call, SpeciesCall::Present));
        assert_eq!(calls[0].reason, "unique_markers_detected");
    }

    #[test]
    fn call_ambiguous_no_unique_markers() {
        // Species with only shared (non-unique) targets, one detected → AMBIGUOUS
        let mut similar: HashMap<String, HashSet<String>> = HashMap::new();
        similar.insert("t_shared".to_string(), ["t_other".to_string()].into_iter().collect());
        let ctx = SimilarityContext {
            cross_species_similar: similar,
            target_is_unique: [("t_shared".to_string(), false)].into_iter().collect(),
            species_targets: [("SpA".to_string(), vec!["t_shared".to_string()])].into_iter().collect(),
            target_to_species: [("t_shared".to_string(), vec!["SpA".to_string()])].into_iter().collect(),
        };
        let calls = call_species(&ctx, &[make_row("t_shared", true, 10)], 1);
        assert_eq!(calls.len(), 1);
        assert!(matches!(calls[0].call, SpeciesCall::Ambiguous));
        assert_eq!(calls[0].reason, "no_unique_markers");
    }

    #[test]
    fn call_cross_reactivity_explained() {
        // SpA: unique t1 + shared t2 — detected first, called PRESENT.
        // SpB: only shared t2 — all detections explained by SpA → ABSENT.
        let mut similar: HashMap<String, HashSet<String>> = HashMap::new();
        similar.insert("t1".to_string(), ["t2".to_string()].into_iter().collect());
        similar.insert("t2".to_string(), ["t1".to_string()].into_iter().collect());
        let ctx = SimilarityContext {
            cross_species_similar: similar,
            target_is_unique: [("t1".to_string(), true), ("t2".to_string(), false)].into_iter().collect(),
            species_targets: [
                ("SpA".to_string(), vec!["t1".to_string(), "t2".to_string()]),
                ("SpB".to_string(), vec!["t2".to_string()]),
            ].into_iter().collect(),
            target_to_species: [
                ("t1".to_string(), vec!["SpA".to_string()]),
                ("t2".to_string(), vec!["SpA".to_string(), "SpB".to_string()]),
            ].into_iter().collect(),
        };
        let rows = vec![make_row("t1", true, 50), make_row("t2", true, 30)];
        let calls = call_species(&ctx, &rows, 1);
        let spa = calls.iter().find(|c| c.species_id == "SpA").unwrap();
        let spb = calls.iter().find(|c| c.species_id == "SpB").unwrap();
        assert!(matches!(spa.call, SpeciesCall::Present), "SpA should be PRESENT");
        assert!(matches!(spb.call, SpeciesCall::Absent), "SpB should be ABSENT (explained)");
        assert_eq!(spb.reason, "cross_reactivity_explained");
        assert!(spb.explained_by.contains(&"SpA".to_string()));
    }

    #[test]
    fn call_ambiguous_insufficient_unique() {
        // 2 unique targets, only 1 detected, min_unique=2 → AMBIGUOUS
        let ctx = SimilarityContext {
            cross_species_similar: HashMap::new(),
            target_is_unique: [("t1".to_string(), true), ("t2".to_string(), true)].into_iter().collect(),
            species_targets: [("SpA".to_string(), vec!["t1".to_string(), "t2".to_string()])].into_iter().collect(),
            target_to_species: [
                ("t1".to_string(), vec!["SpA".to_string()]),
                ("t2".to_string(), vec!["SpA".to_string()]),
            ].into_iter().collect(),
        };
        let rows = vec![make_row("t1", true, 20), make_row("t2", false, 0)];
        let calls = call_species(&ctx, &rows, 2);
        assert_eq!(calls.len(), 1);
        assert!(matches!(calls[0].call, SpeciesCall::Ambiguous));
        assert_eq!(calls[0].reason, "insufficient_unique_evidence");
    }
}
