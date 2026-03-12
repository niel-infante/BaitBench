mod alignment;
mod cli;
mod commands;
mod external;
mod fasta;
mod io_utils;
mod sampling;

use anyhow::{bail, Result};
use clap::Parser;

use cli::{Cli, Commands};
use commands::{capture, coverage_curve, enrich, filter, generate_list, map_reads, metrics, prepare, probe_coverage, report, run, sequence, simulate, xreact};
use io_utils::resolve_sample_arg;

/// Default distractor fraction when neither --distractor-fraction nor --ct is specified.
const DEFAULT_DISTRACTOR_FRACTION: f64 = 0.9;

/// Convert a CT (cycle threshold) score to a distractor fraction.
///
/// Formula: target_fraction = baseline_fraction * 2^(baseline_ct - ct)
///          distractor_fraction = 1 - target_fraction
pub(crate) fn ct_to_distractor_fraction(ct: f64, baseline_ct: f64, baseline_fraction: f64) -> Result<f64> {
    let target_fraction = baseline_fraction * 2.0_f64.powf(baseline_ct - ct);
    if target_fraction >= 1.0 {
        bail!(
            "CT {:.1} with baseline CT {:.1} / fraction {:.4} yields target_fraction {:.4} >= 1.0. \
             Use a higher CT value or adjust --ct-baseline / --ct-baseline-fraction.",
            ct, baseline_ct, baseline_fraction, target_fraction
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
        "CT {:.1} → target fraction {:.6} → distractor fraction {:.6}",
        ct, target_fraction, distractor_fraction
    );
    Ok(distractor_fraction)
}

/// Resolve distractor fraction from --distractor-fraction and --ct flags.
/// Clap enforces mutual exclusivity; this resolves the final value.
fn resolve_distractor_fraction(
    distractor_fraction: Option<f64>,
    ct: Option<f64>,
    ct_baseline: f64,
    ct_baseline_fraction: f64,
) -> Result<f64> {
    match ct {
        Some(ct_val) => ct_to_distractor_fraction(ct_val, ct_baseline, ct_baseline_fraction),
        None => Ok(distractor_fraction.unwrap_or(DEFAULT_DISTRACTOR_FRACTION)),
    }
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
            host_fasta,
            run_name,
            num_fragments,
            distractor_fraction,
            ct,
            ct_baseline,
            ct_baseline_fraction,
            seed,
            capture_method,
            max_mismatches,
            min_match_bases,
            blast_db,
            minimap_preset,
            host_minimap_preset,
            fragment_length_mean,
            fragment_length_min,
            fragment_length_max,
            read_length,
            num_sequences,
            outdir,
            threads,
            fold_enrichment,
            report,
        } => {
            let resolved_df = resolve_distractor_fraction(
                distractor_fraction, ct, ct_baseline, ct_baseline_fraction,
            )?;

            let resolved_sample = sample
                .as_ref()
                .map(|s| resolve_sample_arg(s))
                .transpose()?;

            let resolved_stm = sample_target_map
                .as_ref()
                .map(|p| io_utils::parse_sample_target_map(p))
                .transpose()?;

            let run_name = run_name.unwrap_or_else(|| {
                format!("run_{}", chrono::Local::now().format("%Y%m%d_%H%M%S"))
            });
            let full_outdir = outdir.join(&run_name);

            run::execute(&run::RunArgs {
                targets: &targets,
                genomes: genomes.as_deref(),
                distractors: &distractors,
                probes: &probes,
                sample: resolved_sample.as_ref(),
                sample_target_map: resolved_stm.as_ref(),
                host_fasta: host_fasta.as_deref(),
                run_name,
                num_fragments,
                distractor_fraction: resolved_df,
                ct,
                ct_baseline,
                ct_baseline_fraction,
                seed,
                capture_method: capture_method.into(),
                max_mismatches,
                min_match_bases,
                blast_db,
                minimap_preset,
                host_minimap_preset,
                fragment_length_mean,
                fragment_length_min,
                fragment_length_max,
                read_length,
                num_sequences,
                outdir: full_outdir,
                threads,
                fold_enrichment,
                report,
            })?;
        }

        Commands::Prepare {
            targets,
            genomes,
            distractors,
            sample,
            sample_target_map,
            distractor_fraction,
            ct,
            ct_baseline,
            ct_baseline_fraction,
            outdir,
        } => {
            let resolved_df = resolve_distractor_fraction(
                distractor_fraction, ct, ct_baseline, ct_baseline_fraction,
            )?;

            let resolved_sample = sample
                .as_ref()
                .map(|s| resolve_sample_arg(s))
                .transpose()?;

            let resolved_stm = sample_target_map
                .as_ref()
                .map(|p| io_utils::parse_sample_target_map(p))
                .transpose()?;

            prepare::execute(&prepare::PrepareArgs {
                targets: &targets,
                genomes: genomes.as_deref(),
                distractors: &distractors,
                sample: resolved_sample.as_ref(),
                sample_target_map: resolved_stm.as_ref(),
                distractor_fraction: resolved_df,
                outdir: &outdir,
            })?;
        }

        Commands::Simulate {
            reference,
            weights,
            num_fragments,
            seed,
            output,
            fragment_length_mean,
            fragment_length_min,
            fragment_length_max,
        } => {
            simulate::execute(&simulate::SimulateArgs {
                reference: &reference,
                weights: &weights,
                num_fragments,
                seed,
                output: &output,
                fragment_length_mean,
                fragment_length_min,
                fragment_length_max,
            })?;
        }

        Commands::Capture {
            probes,
            fragments,
            method,
            max_mismatches,
            min_match_bases,
            blast_db,
            output,
            log_file,
            threads,
        } => {
            capture::execute(&capture::CaptureArgs {
                method: method.into(),
                probes: &probes,
                fragments: &fragments,
                max_mismatches,
                min_match_bases,
                blast_db: blast_db.as_deref(),
                output: &output,
                log_file: &log_file,
                threads,
            })?;
        }

        Commands::Enrich {
            captured,
            fragments,
            targets,
            distractors,
            fold_enrichment,
            seed,
            output,
        } => {
            enrich::execute(&enrich::EnrichArgs {
                captured: &captured,
                fragments: &fragments,
                targets: &targets,
                distractors: &distractors,
                fold_enrichment,
                seed,
                output: &output,
            })?;
        }

        Commands::Sequence {
            input,
            output,
            read_length,
            num_sequences,
            seed,
        } => {
            sequence::execute(&sequence::SequenceArgs {
                input: &input,
                output: &output,
                read_length,
                num_sequences,
                seed,
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
                minimap_preset: &minimap_preset,
                output: &output,
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
        } => {
            metrics::execute(&metrics::MetricsArgs {
                targets: &targets,
                distractors: &distractors,
                sample: &sample,
                sample_target_map: None,
                detected: &detected,
                fragments: &fragments,
                captured: &captured,
                sam: &sam,
                run_name: &run_name,
                num_fragments,
                seed: &seed,
                output_summary: &output_summary,
                output_detail: &output_detail,
                output_json: output_json.as_deref(),
                output_coverage: output_coverage.as_deref(),
            })?;
        }

        Commands::ProbeCoverage {
            targets,
            probes,
            outdir,
            minimap_preset,
            proximity,
            report,
        } => {
            probe_coverage::execute(&probe_coverage::ProbeCoverageArgs {
                targets: &targets,
                probes: &probes,
                outdir: &outdir,
                minimap_preset: &minimap_preset,
                proximity,
                report,
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
            report,
        } => {
            xreact::execute(&xreact::XreactArgs {
                probes: &probes,
                against: &against.unwrap_or_default(),
                self_mode,
                threshold,
                minimap_preset: &minimap_preset,
                outdir: &outdir,
                report,
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
            ct_baseline,
            ct_baseline_fraction,
            fold_enrichment_values,
            fold_enrichment,
            num_sequences_values,
            num_sequences,
            num_fragments,
            read_length,
            seed,
            fragment_length_mean,
            fragment_length_min,
            fragment_length_max,
            capture_method,
            max_mismatches,
            min_match_bases,
            blast_db,
            host_fasta,
            minimap_preset,
            host_minimap_preset,
            threads,
            outdir,
            report,
        } => {
            let resolved_sample = resolve_sample_arg(&sample)?;
            let resolved_stm = sample_target_map
                .as_ref()
                .map(|p| io_utils::parse_sample_target_map(p))
                .transpose()?;

            // Resolve CT dimension: --ct-values (sweep) or --ct/--distractor-fraction (fixed)
            let mut swept_params = Vec::new();
            let (ct_display_values, ct_distractor_fractions): (Vec<f64>, Vec<f64>) =
                if let Some(ct_vals) = ct_values {
                    swept_params.push("ct".to_string());
                    let mut dfs = Vec::new();
                    for ct_val in &ct_vals {
                        dfs.push(ct_to_distractor_fraction(
                            *ct_val, ct_baseline, ct_baseline_fraction,
                        )?);
                    }
                    (ct_vals, dfs)
                } else if let Some(ct_val) = ct {
                    let df = ct_to_distractor_fraction(ct_val, ct_baseline, ct_baseline_fraction)?;
                    (vec![ct_val], vec![df])
                } else {
                    let df = distractor_fraction.unwrap_or(DEFAULT_DISTRACTOR_FRACTION);
                    (vec![df], vec![df])
                };

            // Resolve FE dimension: --fold-enrichment-values (sweep) or --fold-enrichment (fixed)
            let resolved_fe_values: Vec<Option<f64>> = if let Some(fe_vals) = fold_enrichment_values {
                swept_params.push("fold_enrichment".to_string());
                fe_vals.into_iter().map(Some).collect()
            } else {
                vec![fold_enrichment]
            };

            // Resolve NS dimension: --num-sequences-values (sweep) or --num-sequences (fixed)
            let resolved_ns_values: Vec<Option<usize>> = if let Some(ns_vals) = num_sequences_values {
                swept_params.push("num_sequences".to_string());
                ns_vals.into_iter().map(Some).collect()
            } else {
                vec![num_sequences]
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
                fe_values: resolved_fe_values,
                ns_values: resolved_ns_values,
                swept_params,
                num_fragments,
                read_length,
                seed,
                fragment_length_mean,
                fragment_length_min,
                fragment_length_max,
                capture_method: capture_method.into(),
                max_mismatches,
                min_match_bases,
                blast_db: blast_db.as_deref(),
                host_fasta: host_fasta.as_deref(),
                minimap_preset: &minimap_preset,
                host_minimap_preset: &host_minimap_preset,
                threads,
                outdir,
                report,
            })?;
        }
    }

    Ok(())
}
