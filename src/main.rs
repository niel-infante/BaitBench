mod alignment;
mod cleanup;
mod cli;
mod commands;
mod external;
mod fasta;
mod io_utils;
mod probes;
mod sampling;
mod sdust;
mod seq_utils;
mod target_similarity;
mod thermodynamics;

use anyhow::{bail, Result};
use clap::Parser;

use cli::{Cli, Commands, ToolCommands};
use commands::{assess_probes, build_probes, coverage_curve, filter, generate_list, identify, map_reads, metrics, panel_qc, prepare, probe_coverage, report, run, sequence, simulate, xreact};
use probes::{catch, syotti};
use sequence::ReadSimulator;
use sampling::thermo_sim::SimulateMode as ThSimMode;
use io_utils::resolve_sample_arg;

/// Default distractor fraction when neither --distractor-fraction nor --ct is specified.
const DEFAULT_DISTRACTOR_FRACTION: f64 = 0.9;

/// Pick the best minimap2 preset for read mapping/filtering based on the simulator and profile.
fn auto_minimap_preset(simulator: &ReadSimulator, profile: &str) -> &'static str {
    match simulator {
        ReadSimulator::Badread => match profile {
            "pacbio" => "map-pb",
            _        => "map-ont",  // ont, ont-2020
        },
        _ => "sr",
    }
}

/// Convert a CT (cycle threshold) score to a distractor fraction.
///
/// Formula: target_fraction = baseline_fraction * (1 + efficiency)^(baseline_ct - ct)
///          distractor_fraction = 1 - target_fraction
pub(crate) fn ct_to_distractor_fraction(ct: f64, baseline_ct: f64, baseline_fraction: f64, efficiency: f64) -> Result<f64> {
    let target_fraction = baseline_fraction * (1.0 + efficiency).powf(baseline_ct - ct);
    if target_fraction >= 1.0 {
        bail!(
            "CT {:.1} with baseline CT {:.1} / fraction {:.4} / efficiency {:.4} yields target_fraction {:.4} >= 1.0. \
             Use a higher CT value or adjust --ct-baseline / --ct-baseline-fraction.",
            ct, baseline_ct, baseline_fraction, efficiency, target_fraction
        );
    }
    if target_fraction <= 0.0 {
        bail!(
            "CT {:.1} yields target_fraction <= 0. Check --ct-baseline and --ct-baseline-fraction.",
            ct
        );
    }
    let distractor_fraction = 1.0 - target_fraction;
    log::info!(
        "CT {:.1} (efficiency {:.4}) → target fraction {:.6} → distractor fraction {:.6}",
        ct, efficiency, target_fraction, distractor_fraction
    );
    Ok(distractor_fraction)
}

/// Parse a "CT,fraction" calibration point string.
fn parse_calibration_point(s: &str) -> Result<(f64, f64)> {
    let parts: Vec<&str> = s.splitn(2, ',').collect();
    if parts.len() != 2 {
        bail!("Calibration point must be 'CT,fraction' (e.g. '25.0,0.001'), got: {}", s);
    }
    let ct = parts[0].trim().parse::<f64>()
        .map_err(|_| anyhow::anyhow!("Invalid CT value '{}' in calibration point '{}'", parts[0].trim(), s))?;
    let fraction = parts[1].trim().parse::<f64>()
        .map_err(|_| anyhow::anyhow!("Invalid fraction '{}' in calibration point '{}'", parts[1].trim(), s))?;
    Ok((ct, fraction))
}

/// Derive (baseline_ct, baseline_fraction, efficiency) from two calibration points.
///
/// Formula: 1 + E = (f1 / f2)^(1 / (ct2 - ct1))
fn derive_calibration_params(ct1: f64, f1: f64, ct2: f64, f2: f64) -> Result<(f64, f64, f64)> {
    if (ct1 - ct2).abs() < 1e-10 {
        bail!("CT values in --ct-calibration must be different (got {} and {})", ct1, ct2);
    }
    if f1 <= 0.0 || f1 >= 1.0 {
        bail!("Calibration fraction must be in (0, 1), got {}", f1);
    }
    if f2 <= 0.0 || f2 >= 1.0 {
        bail!("Calibration fraction must be in (0, 1), got {}", f2);
    }
    let efficiency = (f1 / f2).powf(1.0 / (ct2 - ct1)) - 1.0;
    if efficiency <= 0.0 {
        bail!(
            "Derived PCR efficiency {:.4} ≤ 0. Check calibration points — fractions may be inverted \
             (lower CT should correspond to higher target fraction).",
            efficiency
        );
    }
    if efficiency > 1.0 {
        bail!(
            "Derived PCR efficiency {:.4} > 1.0, which is physically impossible. \
             Check calibration points for errors.",
            efficiency
        );
    }
    log::info!(
        "CT calibration: ({:.1}, {:.6}) and ({:.1}, {:.6}) → derived efficiency {:.4}",
        ct1, f1, ct2, f2, efficiency
    );
    Ok((ct1, f1, efficiency))
}

/// Resolve (baseline_ct, baseline_fraction, efficiency) from CLI flags.
/// If --ct-calibration is provided, derives all three from the two calibration points.
fn resolve_ct_params(
    ct_baseline: f64,
    ct_baseline_fraction: f64,
    ct_efficiency: f64,
    ct_calibration: Option<&[String]>,
) -> Result<(f64, f64, f64)> {
    if let Some(cal) = ct_calibration {
        let (ct1, f1) = parse_calibration_point(&cal[0])?;
        let (ct2, f2) = parse_calibration_point(&cal[1])?;
        derive_calibration_params(ct1, f1, ct2, f2)
    } else {
        if ct_efficiency <= 0.0 || ct_efficiency > 1.0 {
            bail!("--ct-efficiency must be in (0, 1], got {}", ct_efficiency);
        }
        Ok((ct_baseline, ct_baseline_fraction, ct_efficiency))
    }
}

/// Resolve distractor fraction from --distractor-fraction and --ct flags.
/// Clap enforces mutual exclusivity; this resolves the final value.
fn resolve_distractor_fraction(
    distractor_fraction: Option<f64>,
    ct: Option<f64>,
    baseline_ct: f64,
    baseline_fraction: f64,
    efficiency: f64,
) -> Result<f64> {
    match ct {
        Some(ct_val) => ct_to_distractor_fraction(ct_val, baseline_ct, baseline_fraction, efficiency),
        None => {
            let df = distractor_fraction.unwrap_or(DEFAULT_DISTRACTOR_FRACTION);
            if !(0.0..=1.0).contains(&df) {
                bail!("--distractor-fraction must be in [0, 1], got {}", df);
            }
            Ok(df)
        }
    }
}

struct ResolvedCommonArgs {
    distractor_fraction: f64,
    sample: Option<std::collections::HashMap<String, f64>>,
    sample_target_map: Option<std::collections::HashMap<String, Vec<String>>>,
}

fn resolve_common_pipeline_args(
    ct_baseline: f64,
    ct_baseline_fraction: f64,
    ct_efficiency: f64,
    ct_calibration: Option<&[String]>,
    distractor_fraction: Option<f64>,
    ct: Option<f64>,
    sample: Option<&[String]>,
    sample_target_map: Option<&std::path::Path>,
) -> Result<ResolvedCommonArgs> {
    let (baseline_ct, baseline_frac, efficiency) =
        resolve_ct_params(ct_baseline, ct_baseline_fraction, ct_efficiency, ct_calibration)?;
    let resolved_df = resolve_distractor_fraction(
        distractor_fraction, ct, baseline_ct, baseline_frac, efficiency,
    )?;
    let resolved_sample = sample.map(resolve_sample_arg).transpose()?;
    let resolved_stm = sample_target_map
        .map(io_utils::parse_sample_target_map)
        .transpose()?;
    Ok(ResolvedCommonArgs {
        distractor_fraction: resolved_df,
        sample: resolved_sample,
        sample_target_map: resolved_stm,
    })
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    // Initialize logging
    let log_level = if cli.verbose { "debug" } else { "info" };
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or(log_level))
        .format_timestamp(None)
        .format_target(false)
        .init();

    match cli.command {
        Commands::Run {
            targets,
            genomes,
            distractors,
            probes,
            sample,
            sample_target_map,
            groups,
            distractor_groups,
            host_fasta,
            run_name,
            num_fragments,
            distractor_fraction,
            ct,
            ct_baseline,
            ct_baseline_fraction,
            ct_efficiency,
            ct_calibration,
            seed,
            simulate_mode,
            hybridization_temperature,
            salt_concentration,
            capture_fraction,
            minimap_preset,
            host_minimap_preset,
            fragment_length_mean,
            fragment_length_min,
            fragment_length_max,
            read_length,
            num_sequences,
            read_simulator,
            sequencer_profile,
            coverage_depth,
            paired_end,
            pe_frag_len_mean,
            pe_frag_len_sd,
            output_format,
            long_read_length_mean,
            long_read_length_sd,
            badread_glitches,
            badread_junk_reads,
            badread_random_reads,
            badread_chimeras,
            outdir,
            output_prefix,
            threads,
            identify,
            identity_threshold,
            min_unique_targets,
            report,
            cleanup,
        } => {
            let resolved = resolve_common_pipeline_args(
                ct_baseline, ct_baseline_fraction, ct_efficiency, ct_calibration.as_deref(),
                distractor_fraction, ct,
                sample.as_deref(), sample_target_map.as_deref(),
            )?;

            let run_name = run_name.unwrap_or_else(|| {
                format!("run_{}", chrono::Local::now().format("%Y%m%d_%H%M%S"))
            });
            let full_outdir = outdir.join(&run_name);

            let sim_mode = match simulate_mode {
                cli::SimulateMode::Thermodynamic => ThSimMode::Thermodynamic,
                cli::SimulateMode::Simple => ThSimMode::Simple,
            };

            let simulator = ReadSimulator::from_str(&read_simulator)?;
            let sequencer_profile = sequencer_profile.unwrap_or_else(|| match simulator {
                ReadSimulator::Art    => "HiSeq2500_150bp".to_string(),
                ReadSimulator::Badread => "ont".to_string(),
                ReadSimulator::Perfect => String::new(),
            });
            let preset = auto_minimap_preset(&simulator, &sequencer_profile);
            let minimap_preset = minimap_preset.unwrap_or_else(|| preset.to_string());
            let host_minimap_preset = host_minimap_preset.unwrap_or_else(|| preset.to_string());

            run::execute(&run::RunArgs {
                targets: &targets,
                genomes: genomes.as_deref(),
                distractors: &distractors,
                probes: &probes,
                sample: resolved.sample.as_ref(),
                sample_target_map: resolved.sample_target_map.as_ref(),
                groups: groups.as_deref(),
                distractor_groups: distractor_groups.as_deref(),
                host_fasta: host_fasta.as_deref(),
                run_name,
                num_fragments,
                distractor_fraction: resolved.distractor_fraction,
                ct,
                ct_baseline,
                ct_baseline_fraction,
                ct_efficiency,
                ct_calibration,
                seed,
                simulate_mode: sim_mode,
                hybridization_temperature,
                na_conc_m: salt_concentration / 1000.0,
                capture_fraction,
                minimap_preset,
                host_minimap_preset,
                fragment_length_mean,
                fragment_length_min,
                fragment_length_max,
                read_length,
                num_sequences,
                simulator,
                sequencer_profile,
                coverage_depth,
                paired_end,
                pe_frag_len_mean,
                pe_frag_len_sd,
                output_format,
                long_read_length_mean,
                long_read_length_sd,
                badread_glitches,
                badread_junk_reads,
                badread_random_reads,
                badread_chimeras,
                outdir: full_outdir,
                threads,
                identify,
                identity_threshold,
                min_unique_targets,
                report,
                cleanup,
                output_prefix,
            })?;
        }

        Commands::Prepare {
            targets,
            genomes,
            distractors,
            sample,
            sample_target_map,
            groups,
            distractor_groups,
            distractor_fraction,
            ct,
            ct_baseline,
            ct_baseline_fraction,
            ct_efficiency,
            ct_calibration,
            outdir,
            output_prefix,
        } => {
            let resolved = resolve_common_pipeline_args(
                ct_baseline, ct_baseline_fraction, ct_efficiency, ct_calibration.as_deref(),
                distractor_fraction, ct,
                sample.as_deref(), sample_target_map.as_deref(),
            )?;

            prepare::execute(&prepare::PrepareArgs {
                targets: &targets,
                genomes: genomes.as_deref(),
                distractors: &distractors,
                sample: resolved.sample.as_ref(),
                sample_target_map: resolved.sample_target_map.as_ref(),
                groups: groups.as_deref(),
                distractor_groups: distractor_groups.as_deref(),
                distractor_fraction: resolved.distractor_fraction,
                outdir: &outdir,
                output_prefix: &output_prefix,
            })?;
        }

        Commands::Simulate {
            reference,
            weights,
            probes,
            num_fragments,
            capture_fraction,
            simulate_mode,
            hybridization_temperature,
            salt_concentration,
            seed,
            output,
            fragment_length_mean,
            fragment_length_min,
            fragment_length_max,
            threads,
        } => {
            let sim_mode = match simulate_mode {
                cli::SimulateMode::Thermodynamic => ThSimMode::Thermodynamic,
                cli::SimulateMode::Simple => ThSimMode::Simple,
            };
            simulate::execute(&simulate::SimulateArgs {
                reference: &reference,
                weights: &weights,
                probes: &probes,
                num_fragments,
                capture_fraction,
                simulate_mode: sim_mode,
                hybridization_temperature,
                na_conc_m: salt_concentration / 1000.0,
                seed,
                output: &output,
                fragment_length_mean,
                fragment_length_min,
                fragment_length_max,
                threads,
            })?;
        }

        Commands::Sequence {
            input,
            output,
            output_r2,
            read_length,
            num_sequences,
            seed,
            read_simulator,
            sequencer_profile,
            coverage_depth,
            paired_end,
            pe_frag_len_mean,
            pe_frag_len_sd,
            output_format,
            long_read_length_mean,
            long_read_length_sd,
            badread_glitches,
            badread_junk_reads,
            badread_random_reads,
            badread_chimeras,
        } => {
            let simulator = sequence::ReadSimulator::from_str(&read_simulator)?;
            let profile = sequencer_profile.unwrap_or_else(|| match simulator {
                sequence::ReadSimulator::Art     => "HiSeq2500_150bp".to_string(),
                sequence::ReadSimulator::Badread => "ont".to_string(),
                sequence::ReadSimulator::Perfect => String::new(),
            });
            sequence::execute(&sequence::SequenceArgs {
                input: &input,
                output: &output,
                output_r2: output_r2.as_deref(),
                read_length,
                num_sequences,
                seed,
                simulator,
                sequencer_profile: profile,
                coverage_depth,
                paired_end,
                pe_frag_len_mean,
                pe_frag_len_sd,
                output_format,
                long_read_length_mean,
                long_read_length_sd,
                badread_glitches,
                badread_junk_reads,
                badread_random_reads,
                badread_chimeras,
            })?;
        }

        Commands::Filter {
            host,
            reads,
            minimap_preset,
            output,
            log_file,
        } => {
            filter::execute(&filter::FilterArgs {
                host: &host,
                reads: &reads,
                reads_r2: None,
                minimap_preset: &minimap_preset,
                output: &output,
                output_r2: None,
                log_file: &log_file,
            })?;
        }

        Commands::Map {
            reference,
            reads,
            minimap_preset,
            output,
            log_file,
        } => {
            map_reads::execute(&map_reads::MapArgs {
                reference: &reference,
                reads: &reads,
                reads_r2: None,
                minimap_preset: &minimap_preset,
                output: &output,
                log_file: &log_file,
            })?;
        }

        Commands::List { sam, output } => {
            generate_list::execute(&generate_list::ListArgs {
                sam: &sam,
                output: &output,
            })?;
        }

        Commands::Metrics {
            targets,
            distractors,
            sample,
            detected,
            fragments,
            captured,
            sam,
            run_name,
            num_fragments,
            seed,
            output_summary,
            output_detail,
            output_json,
            output_coverage,
            reads_sequenced,
            reads_after_filter,
            target_groups,
            distractor_groups,
        } => {
            metrics::execute(&metrics::MetricsArgs {
                targets: &targets,
                distractors: &distractors,
                sample: &sample,
                sample_target_map: None,
                target_groups: target_groups.as_deref(),
                distractor_groups: distractor_groups.as_deref(),
                detected: &detected,
                fragments: &fragments,
                captured: &captured,
                sam: &sam,
                run_name: &run_name,
                num_fragments,
                seed: &seed,
                output_summary: &output_summary,
                output_detail: &output_detail,
                output_group_detail: None,
                output_json: output_json.as_deref(),
                output_coverage: output_coverage.as_deref(),
                reads_sequenced,
                reads_after_filter,
            })?;
        }

        Commands::ProbeCoverage {
            targets,
            probes,
            outdir,
            minimap_preset,
            proximity,
            output_prefix,
            report,
            cleanup,
        } => {
            probe_coverage::execute(&probe_coverage::ProbeCoverageArgs {
                targets: &targets,
                probes: &probes,
                outdir: &outdir,
                output_prefix: &output_prefix,
                minimap_preset: &minimap_preset,
                proximity,
                report,
                cleanup,
            })?;
        }

        Commands::Report {
            summary,
            detail,
            params,
            coverage,
            run_name,
            output,
            report,
        } => {
            report::execute(&report::ReportArgs {
                summary: &summary,
                detail: &detail,
                params: &params,
                coverage: coverage.as_deref(),
                species_calls: None,
                group_detail: None,
                run_name: &run_name,
                output: &output,
                report,
            })?;
        }

        Commands::Xreact {
            probes,
            against,
            self_mode,
            threshold,
            minimap_preset,
            outdir,
            output_prefix,
            report,
            cleanup,
        } => {
            xreact::execute(&xreact::XreactArgs {
                probes: &probes,
                against: &against.unwrap_or_default(),
                self_mode,
                threshold,
                minimap_preset: &minimap_preset,
                outdir: &outdir,
                output_prefix: &output_prefix,
                report,
                cleanup,
            })?;
        }

        Commands::PanelQc {
            targets,
            sample_target_map,
            identity_threshold,
            minimap_preset,
            outdir,
            output_prefix,
            report,
            cleanup,
        } => {
            panel_qc::execute(&panel_qc::PanelQcArgs {
                targets: &targets,
                sample_target_map: &sample_target_map,
                identity_threshold,
                minimap_preset: &minimap_preset,
                outdir: &outdir,
                output_prefix: &output_prefix,
                report,
                cleanup,
            })?;
        }

        Commands::Identify {
            detected_detail,
            sample_target_map,
            target_similarity,
            targets,
            identity_threshold,
            minimap_preset,
            min_unique_targets,
            outdir,
            output_prefix,
        } => {
            identify::execute(&identify::IdentifyArgs {
                detected_detail: &detected_detail,
                sample_target_map: &sample_target_map,
                target_similarity_file: target_similarity.as_deref(),
                targets_fasta: targets.as_deref(),
                identity_threshold,
                minimap_preset: &minimap_preset,
                min_unique_targets,
                outdir: &outdir,
                output_prefix: &output_prefix,
            })?;
        }

        Commands::BuildProbes {
            targets,
            method,
            probe_length,
            step,
            catch_probe_stride,
            catch_mismatches,
            catch_extension,
            catch_coverage,
            catch_minhash_threshold,
            min_gc,
            max_gc,
            max_n_frac,
            dust_threshold,
            dust_window,
            max_masked_frac,
            collapse_threshold,
            dedup_threshold,
            threads,
            genomes,
            threshold,
            minimap_preset,
            proximity,
            skip_assess,
            outdir,
            output_prefix,
            report,
            cleanup,
            refine_threshold,
            refine_iterations,
            refine_until_stable,
            syotti_mismatches,
            syotti_seed_len,
            pt_step,
            pt_identity,
            pt_coverage,
            pt_batch_size,
            pt_max_panel_size,
            pt_min_depth,
            pt_max_iterations,
            pt_min_coverage_gain,
            no_n_in_probes,
        } => {
            build_probes::execute(&build_probes::BuildProbesArgs {
                targets: &targets,
                method,
                probe_length,
                step,
                catch_probe_stride,
                catch_mismatches,
                catch_extension,
                catch_coverage,
                catch_minhash_threshold,
                min_gc,
                max_gc,
                max_n_frac,
                dust_threshold,
                dust_window,
                max_masked_frac,
                collapse_threshold,
                dedup_threshold,
                threads,
                genomes: &genomes.unwrap_or_default(),
                threshold,
                minimap_preset: &minimap_preset,
                proximity,
                skip_assess,
                outdir: &outdir,
                output_prefix: &output_prefix,
                report,
                cleanup,
                refine_threshold,
                refine_iterations,
                refine_until_stable,
                syotti_mismatches,
                syotti_seed_len,
                pt_step,
                pt_identity,
                pt_coverage,
                pt_batch_size,
                pt_max_panel_size,
                pt_min_depth,
                pt_max_iterations,
                pt_min_coverage_gain,
                no_n_in_probes,
            })?;
        }

        Commands::AssessProbes {
            targets,
            probes,
            genomes,
            threshold,
            minimap_preset,
            proximity,
            outdir,
            output_prefix,
            report,
            cleanup,
            refine_threshold,
            refine_iterations,
            refine_until_stable,
            gap_min_length,
            no_individual_targets,
            threads,
        } => {
            let threads = threads.unwrap_or_else(|| {
                std::thread::available_parallelism().map(|n| n.get()).unwrap_or(4)
            });
            assess_probes::execute(&assess_probes::AssessProbesArgs {
                targets: &targets,
                probes: &probes,
                genomes: &genomes.unwrap_or_default(),
                threshold,
                minimap_preset: &minimap_preset,
                proximity,
                outdir: &outdir,
                output_prefix: &output_prefix,
                report,
                cleanup,
                build_stats_file: None,
                build_params_file: None,
                refine_threshold,
                refine_iterations,
                refine_until_stable,
                no_individual_targets,
                gap_min_length,
                threads,
            })?;
        }

        Commands::CoverageCurve {
            targets,
            genomes,
            distractors,
            probes,
            sample,
            sample_target_map,
            ct_values,
            ct,
            distractor_fraction,
            distractor_fraction_values,
            ct_baseline,
            ct_baseline_fraction,
            ct_efficiency,
            ct_calibration,
            capture_fraction_values,
            capture_fraction,
            simulate_mode,
            hybridization_temperature_values,
            hybridization_temperature,
            salt_concentration,
            num_sequences_values,
            num_sequences,
            num_fragments,
            read_length,
            seed,
            fragment_length_mean,
            fragment_length_min,
            fragment_length_max,
            host_fasta,
            minimap_preset,
            host_minimap_preset,
            threads,
            outdir,
            output_prefix,
            report,
            cleanup,
        } => {
            let resolved_sample = resolve_sample_arg(&sample)?;
            let resolved_stm = sample_target_map
                .as_ref()
                .map(|p| io_utils::parse_sample_target_map(p))
                .transpose()?;

            let (baseline_ct, baseline_frac, efficiency) = resolve_ct_params(
                ct_baseline, ct_baseline_fraction, ct_efficiency, ct_calibration.as_deref(),
            )?;

            // Resolve CT/DF dimension: --distractor-fraction-values (sweep DF) or
            // --ct-values (sweep CT) or --ct/--distractor-fraction (fixed single value)
            let mut swept_params = Vec::new();
            let (ct_display_values, ct_distractor_fractions): (Vec<f64>, Vec<f64>) =
                if let Some(df_vals) = distractor_fraction_values {
                    swept_params.push("distractor_fraction".to_string());
                    (df_vals.clone(), df_vals)
                } else if let Some(ct_vals) = ct_values {
                    swept_params.push("ct".to_string());
                    let mut dfs = Vec::new();
                    for ct_val in &ct_vals {
                        dfs.push(ct_to_distractor_fraction(
                            *ct_val, baseline_ct, baseline_frac, efficiency,
                        )?);
                    }
                    (ct_vals, dfs)
                } else if let Some(ct_val) = ct {
                    let df = ct_to_distractor_fraction(ct_val, baseline_ct, baseline_frac, efficiency)?;
                    (vec![ct_val], vec![df])
                } else {
                    let df = distractor_fraction.unwrap_or(DEFAULT_DISTRACTOR_FRACTION);
                    (vec![df], vec![df])
                };

            // Resolve capture-fraction dimension: --capture-fraction-values (sweep) or --capture-fraction (fixed)
            let resolved_cf_values: Vec<f64> = if let Some(cf_vals) = capture_fraction_values {
                swept_params.push("capture_fraction".to_string());
                cf_vals
            } else {
                vec![capture_fraction]
            };

            // Resolve NS dimension: --num-sequences-values (sweep) or --num-sequences (fixed)
            let resolved_ns_values: Vec<Option<usize>> = if let Some(ns_vals) = num_sequences_values {
                swept_params.push("num_sequences".to_string());
                ns_vals.into_iter().map(Some).collect()
            } else {
                vec![num_sequences]
            };

            // Resolve hybridization-temperature dimension
            let resolved_temp_values: Vec<f64> = if let Some(temp_vals) = hybridization_temperature_values {
                swept_params.push("hybridization_temperature".to_string());
                temp_vals
            } else {
                vec![hybridization_temperature]
            };

            let sim_mode = match simulate_mode {
                cli::SimulateMode::Thermodynamic => ThSimMode::Thermodynamic,
                cli::SimulateMode::Simple => ThSimMode::Simple,
            };

            coverage_curve::execute(&coverage_curve::CoverageCurveArgs {
                targets: &targets,
                genomes: genomes.as_deref(),
                distractors: &distractors,
                probes: &probes,
                sample: &resolved_sample,
                sample_target_map: resolved_stm.as_ref(),
                ct_display_values,
                ct_distractor_fractions,
                cf_values: resolved_cf_values,
                ns_values: resolved_ns_values,
                swept_params,
                simulate_mode: sim_mode,
                temp_values: resolved_temp_values,
                na_conc_m: salt_concentration / 1000.0,
                num_fragments,
                read_length,
                seed,
                fragment_length_mean,
                fragment_length_min,
                fragment_length_max,
                host_fasta: host_fasta.as_deref(),
                minimap_preset: &minimap_preset,
                host_minimap_preset: &host_minimap_preset,
                threads,
                outdir,
                output_prefix,
                report,
                cleanup,
            })?;
        }

        Commands::Tool { command } => match command {
            ToolCommands::Syotti {
                targets,
                output,
                probe_length,
                mismatches,
                seed_len,
            } => {
                syotti::design_probes(&targets, &output, probe_length, mismatches, seed_len)?;
            }

            ToolCommands::Catch {
                targets,
                output,
                probe_length,
                stride,
                mismatches,
                extension,
                coverage,
                minhash_threshold,
            } => {
                catch::design_probes(
                    &targets,
                    &output,
                    probe_length,
                    stride,
                    mismatches,
                    extension,
                    coverage,
                    minhash_threshold,
                )?;
            }

            ToolCommands::Dustview {
                input,
                dust_threshold,
                dust_window,
            } => {
                commands::tool_dustview::execute(input.as_deref(), dust_threshold, dust_window)?;
            }

            ToolCommands::Collapse {
                input,
                output,
                threshold,
                threads,
                log_file,
            } => {
                external::cdhit::check_available()?;
                external::cdhit::cluster(&input, &output, threshold, threads, &log_file)?;
            }
        },
    }

    Ok(())
}
