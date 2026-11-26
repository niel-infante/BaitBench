use anyhow::Result;
use std::path::Path;

use crate::alignment::probe_cache;
use crate::external::{blastn, minimap2};
use crate::fasta;
use crate::sampling::thermo_sim::SimulateMode;
use crate::sampling::{thermo_sim, weights, FragmentLengthParams};
use crate::thermodynamics::ThermoModel;

pub struct SimulateArgs<'a> {
    pub reference: &'a Path,
    pub weights: &'a Path,
    pub probes: &'a Path,
    pub num_fragments: usize,
    /// Capture efficiency (0.0–1.0): of the fragments carrying enough of a probe's
    /// site, the fraction the probes pull down. A preset — the probe concentration
    /// achieving it is back-solved, unless `probe_concentration` overrides it. Not an
    /// output composition target; see `solve_probe_concentration`.
    pub capture_fraction: f64,
    /// Probability a fragment overlapping no probe is retained anyway.
    pub background_retention: f64,
    /// Physical probe concentration, bypassing the `capture_fraction` solver.
    pub probe_concentration: Option<f64>,
    /// Fraction of a probe's site a fragment must carry to be captured by it.
    pub min_overlap_frac: f64,
    /// Optional path to a probe-placement cache. Reused when it matches this probes
    /// file, reference and aligner; written after alignment when it does not exist.
    pub probe_hits: Option<&'a Path>,
    /// Pairs in the cooperative melting unit; 0 uses the whole duplex.
    pub melting_window: usize,
    pub simulate_mode: SimulateMode,
    /// Hybridization temperature in °C (thermodynamic mode only)
    pub hybridization_temperature: f64,
    /// Na+ concentration in molar (e.g. 0.05 for 50 mM); thermodynamic mode only.
    /// SantaLucia (1998) tables assume 1 M; other values apply an entropy correction.
    pub na_conc_m: f64,
    pub seed: Option<u64>,
    pub output: &'a Path,
    pub fragment_length_mean: f64,
    pub fragment_length_min: usize,
    pub fragment_length_max: usize,
    /// Standard deviation of the fragment length distribution. `None` uses
    /// `(max - min) / 6`.
    pub fragment_length_sd: Option<f64>,
    /// Shape of the fragment length distribution.
    pub fragment_length_shape: crate::sampling::FragmentLengthShape,
    pub threads: usize,
    /// Aligner used to place probes on the reference for scoring.
    pub aligner: crate::cli::ProbeAligner,
}

pub fn execute(args: &SimulateArgs) -> Result<()> {
    if args.fragment_length_mean < args.fragment_length_min as f64
        || args.fragment_length_mean > args.fragment_length_max as f64
    {
        anyhow::bail!(
            "--fragment-length-mean ({}) must be within [--fragment-length-min ({}), --fragment-length-max ({})]",
            args.fragment_length_mean,
            args.fragment_length_min,
            args.fragment_length_max,
        );
    }

    if let Some(seed) = args.seed {
        log::info!("Using random seed: {}", seed);
    }

    log::info!("Parsing FASTA reference...");
    let sequences = fasta::parse_fasta(args.reference)?;

    log::info!("Parsing weights file...");
    let wts = weights::parse_weights(args.weights)?;

    // --- Step 1: Place probes on the reference ---
    //
    // The aligner choice matters more here than it looks. Minimizer seeding loses
    // divergent probes, and probe density drives per-locus capture, so a panel built
    // against relatives of the sampled organism can be substantially truncated before
    // the thermodynamic model sees it. `--probe-aligner blast` recovers those.
    //
    // This is also the most expensive step by a wide margin, and it depends only on the
    // probes and the reference. `--probe-hits` persists the result so a titration
    // series, a replicate set or a parameter sweep pays for it once.
    use crate::cli::ProbeAligner;
    let aligner_name = match args.aligner {
        ProbeAligner::Minimap2 => "minimap2",
        ProbeAligner::SeedExtend => "seed-extend",
        ProbeAligner::Blast => "blast",
    };
    let thermo_model = ThermoModel::new(args.hybridization_temperature, args.na_conc_m);

    // Parsed once and reused: for the cache key below and for the alignment-rate
    // report further down.
    let probe_seqs = fasta::parse_fasta(args.probes)?;
    let cache_key = args
        .probe_hits
        .map(|_| probe_cache::CacheKey::new(&probe_seqs, &sequences, &wts, aligner_name));

    let mut placements: Option<Vec<probe_cache::Placement>> = None;
    let mut replacing_stale = false;
    if let (Some(path), Some(key)) = (args.probe_hits, cache_key.as_ref()) {
        match probe_cache::read(path, key) {
            Ok(p) => {
                log::info!(
                    "Reusing {} cached probe placements from {}",
                    p.len(),
                    path.display()
                );
                placements = Some(p);
            }
            Err(miss) => {
                replacing_stale = matches!(miss, probe_cache::Miss::Stale { .. });
                log::info!("Computing probe placements ({}): {}", path.display(), miss);
            }
        }
    }

    let placements = match placements {
        Some(p) => p,
        None => {
            let temp_hits = match args.aligner {
                ProbeAligner::Blast => args.output.with_extension("probe_hits.tsv"),
                _ => args.output.with_extension("probe_hits.sam"),
            };
            let temp_log = args.output.with_extension("probe_align.log");

            log::info!(
                "Aligning probes to reference ({})...",
                match args.aligner {
                    ProbeAligner::Minimap2 => "rammap, sr preset retuned for divergent probes",
                    ProbeAligner::SeedExtend => "built-in ungapped local seed-and-extend",
                    ProbeAligner::Blast => "blastn-short",
                }
            );
            let computed = match args.aligner {
                ProbeAligner::Blast => {
                    blastn::check_available()?;
                    blastn::probe_align(
                        args.reference,
                        args.probes,
                        &temp_hits,
                        &temp_log,
                        args.threads,
                    )?;
                    thermo_sim::load_probe_hits_blast(&temp_hits, &wts)?
                }
                ProbeAligner::Minimap2 => {
                    minimap2::probe_align(
                        "sr",
                        args.reference,
                        args.probes,
                        &temp_hits,
                        &temp_log,
                        args.threads,
                        1000,
                    )?;
                    thermo_sim::load_probe_hits(&temp_hits, &sequences, &wts)?
                }
                ProbeAligner::SeedExtend => {
                    thermo_sim::probe_hits_seed_extend(&sequences, args.probes, &wts, args.threads)?
                }
            };

            let _ = std::fs::remove_file(&temp_hits);
            let _ = std::fs::remove_file(&temp_log);

            if let (Some(path), Some(key)) = (args.probe_hits, cache_key.as_ref()) {
                if replacing_stale {
                    // One file holds one (probes, reference, aligner) combination.
                    // Alternating two panels through a single path is not wrong, but
                    // every run then rebuilds the alignment it just discarded, which
                    // looks like the cache simply not working.
                    log::warn!(
                        "Replacing the cache at {} — it was built for different inputs. \
                         Give each probes/reference combination its own --probe-hits \
                         path, or they will take turns overwriting each other.",
                        path.display()
                    );
                }
                probe_cache::write(path, key, &computed)?;
                log::info!(
                    "Wrote {} probe placements to {} for reuse",
                    computed.len(),
                    path.display()
                );
            }
            computed
        }
    };

    // --- Step 2: Score the placements at this temperature and salt ---
    //
    // Cheap, and deliberately not cached: it is what makes one stored alignment valid
    // across a --hybridization-temperature or --simulate-mode sweep.
    let hits_by_probe = thermo_sim::score_placements(
        &placements,
        args.simulate_mode,
        &thermo_model,
        args.melting_window,
    );

    let n_probes_with_hits = hits_by_probe.len();
    let n_total_hits: usize = hits_by_probe.values().map(|v| v.len()).sum();

    // Report the probe alignment RATE, not just the count. A panel designed against
    // relatives of the sampled organism can lose a large share of its probes to the
    // aligner's seeding before the thermodynamic model ever sees them, and probe
    // density is the dominant driver of per-locus capture. Silently simulating with a
    // fraction of the panel produces plausible-looking but wrong coverage.
    let n_probes_total = probe_seqs.len();
    let align_rate = if n_probes_total > 0 {
        n_probes_with_hits as f64 / n_probes_total as f64
    } else {
        0.0
    };
    log::info!(
        "Found {} of {} probes with hits ({:.1}%), {} total hits",
        n_probes_with_hits,
        n_probes_total,
        100.0 * align_rate,
        n_total_hits
    );
    if align_rate < 0.5 && n_probes_total > 0 {
        if matches!(args.aligner, ProbeAligner::Blast | ProbeAligner::SeedExtend) {
            // Already on the sensitive path, so a low rate is informative rather than
            // actionable: a cross-species panel legitimately has probes with no
            // detectable homology to this reference.
            log::info!(
                "{:.1}% of probes aligned. For a panel designed against relatives of \
                 these sequences that is expected; the remainder have no homology at \
                 80% identity over half the probe.",
                100.0 * align_rate
            );
        } else {
            log::warn!(
                "Only {:.1}% of probes aligned to the reference. If the panel targets \
                 relatives of these sequences rather than the sequences themselves, \
                 much of it may be too divergent for minimizer seeding. Try \
                 --probe-aligner blast, which is markedly more sensitive here.",
                100.0 * align_rate
            );
        }
    }

    // --- Step 3: Build the retention model ---
    //
    // One probability per fragment, p_retain = 1 - (1 - p_ns) * prod_i (1 - theta_i),
    // integrated exactly over each reference because it is piecewise constant in
    // position. Composition emerges from abundance x efficiency instead of being
    // imposed as a quota, which is what makes a dilution series titrate and what keeps
    // --distractor-fraction meaningful above 0.5.
    let capture_fraction = args.capture_fraction.clamp(0.0, 1.0);
    let by_seq = thermo_sim::group_hits_by_sequence(&hits_by_probe);

    let mut params = thermo_sim::RetentionParams {
        probe_concentration: 0.0,
        background_retention: args.background_retention,
        min_overlap_frac: args.min_overlap_frac,
        capture_efficiency: capture_fraction,
        frag_len: (args.fragment_length_mean.round() as usize).max(1),
    };

    params.probe_concentration = args
        .probe_concentration
        .unwrap_or(thermo_sim::DEFAULT_PROBE_CONCENTRATION);
    log::info!(
        "Capture: efficiency {:.3} (--capture-fraction), [probe]={:.4e}, \
         background_retention={:.2e}",
        capture_fraction,
        params.probe_concentration,
        params.background_retention
    );

    let (short_frac, n_short, n_refs) =
        thermo_sim::short_contig_mass_fraction(&sequences, &wts, params.frag_len);
    // Triggered on either the mass share or the reference count, because they catch
    // different shapes of the same problem. The D6331 panel's Candida albicans
    // reference is 13,062 of 13,145 references but only ~0.9% of weighted mass, so a
    // mass-only threshold of 1% stayed silent on the assembly that motivated this.
    let short_ref_frac = if n_refs > 0 {
        n_short as f64 / n_refs as f64
    } else {
        0.0
    };
    if short_frac > 0.001 || short_ref_frac > 0.05 {
        log::warn!(
            "{} of {} weighted references ({:.1}%) are shorter than the fragment \
             length ({} bp), carrying {:.2}% of weighted input mass. Those yield \
             truncated fragments, which shifts the library's length distribution. \
             Check --fragment-length-mean against the assembly's contig sizes.",
            n_short,
            n_refs,
            100.0 * short_ref_frac,
            params.frag_len,
            100.0 * short_frac
        );
    }

    let profiles = thermo_sim::build_profiles(&sequences, &wts, &by_seq, &params, &thermo_model);
    let stats = thermo_sim::retention_stats(&profiles, &wts);
    log::info!(
        "Capture model: [probe]={:.4e}, background_retention={:.2e} \u{2192} \
         on-target {:.2}%, yield {:.4}%, enrichment {:.1}x over {} references",
        params.probe_concentration,
        params.background_retention,
        100.0 * stats.on_target_fraction,
        100.0 * stats.yield_fraction,
        stats.enrichment,
        stats.n_references
    );

    let len_params = FragmentLengthParams {
        mean: args.fragment_length_mean,
        min: args.fragment_length_min,
        max: args.fragment_length_max,
        sd: args.fragment_length_sd,
        shape: args.fragment_length_shape,
    };

    // The model's own account of the capture reaction, for metrics to report. These
    // are analytic expectations, not counts: under per-fragment retention every
    // fragment written was retained, so counting the output cannot recover the yield.
    let stats_path = args.output.with_extension("model_stats.tsv");
    {
        use std::io::Write as _;
        let mut f = std::fs::File::create(&stats_path)?;
        writeln!(f, "parameter\tvalue")?;
        writeln!(f, "probe_concentration\t{:.6e}", params.probe_concentration)?;
        writeln!(
            f,
            "background_retention\t{:.6e}",
            params.background_retention
        )?;
        writeln!(f, "capture_efficiency_target\t{:.6}", capture_fraction)?;
        writeln!(f, "yield_fraction\t{:.8}", stats.yield_fraction)?;
        writeln!(f, "on_target_fraction\t{:.6}", stats.on_target_fraction)?;
        writeln!(f, "enrichment\t{:.4}", stats.enrichment)?;
    }

    // The model's own expectation, per reference, before sampling noise. Per-locus
    // efficiency is what the Borrelia LOD data constrains, and it is invisible in the
    // aggregate on-target fraction.
    let retention_report = args.output.with_extension("retention.tsv");
    thermo_sim::write_retention_report(&profiles, &wts, &retention_report)?;
    log::info!(
        "Wrote per-reference retention report to {}",
        retention_report.display()
    );

    log::info!("Sampling {} fragments by retention...", args.num_fragments);
    let all_fragments = thermo_sim::sample_retention_fragments(
        &profiles,
        &sequences,
        &wts,
        args.num_fragments,
        len_params,
        args.seed,
        0,
    )?;

    // --- Step 4: Write output ---
    log::info!(
        "Writing {} fragments to {}",
        all_fragments.len(),
        args.output.display()
    );
    thermo_sim::write_fragments(&all_fragments, args.output)?;

    log::info!(
        "Simulate complete: {} fragments written",
        all_fragments.len()
    );
    Ok(())
}
