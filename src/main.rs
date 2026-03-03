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
use commands::{capture, enrich, filter, generate_list, map_reads, metrics, prepare, probe_coverage, report, run, sequence, simulate};

/// Default distractor fraction when neither --distractor-fraction nor --ct is specified.
const DEFAULT_DISTRACTOR_FRACTION: f64 = 0.9;

/// Convert a CT (cycle threshold) score to a distractor fraction.
///
/// Formula: target_fraction = baseline_fraction * 2^(baseline_ct - ct)
///          distractor_fraction = 1 - target_fraction
fn ct_to_distractor_fraction(ct: f64, baseline_ct: f64, baseline_fraction: f64) -> Result<f64> {
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
            distractors,
            probes,
            sample,
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
            outdir,
            threads,
            fold_enrichment,
            no_report,
        } => {
            let resolved_df = resolve_distractor_fraction(
                distractor_fraction, ct, ct_baseline, ct_baseline_fraction,
            )?;

            let run_name = run_name.unwrap_or_else(|| {
                format!("run_{}", chrono::Local::now().format("%Y%m%d_%H%M%S"))
            });
            let full_outdir = outdir.join(&run_name);

            run::execute(&run::RunArgs {
                targets: &targets,
                distractors: &distractors,
                probes: &probes,
                sample: sample.as_deref(),
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
                outdir: full_outdir,
                threads,
                fold_enrichment,
                no_report,
            })?;
        }

        Commands::Prepare {
            targets,
            distractors,
            sample,
            distractor_fraction,
            ct,
            ct_baseline,
            ct_baseline_fraction,
            outdir,
        } => {
            let resolved_df = resolve_distractor_fraction(
                distractor_fraction, ct, ct_baseline, ct_baseline_fraction,
            )?;

            prepare::execute(&prepare::PrepareArgs {
                targets: &targets,
                distractors: &distractors,
                sample: sample.as_deref(),
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
        } => {
            sequence::execute(&sequence::SequenceArgs {
                input: &input,
                output: &output,
                read_length,
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
            no_report,
        } => {
            probe_coverage::execute(&probe_coverage::ProbeCoverageArgs {
                targets: &targets,
                probes: &probes,
                outdir: &outdir,
                minimap_preset: &minimap_preset,
                proximity,
                no_report,
            })?;
        }

        Commands::Report {
            summary,
            detail,
            params,
            coverage,
            run_name,
            output,
        } => {
            report::execute(&report::ReportArgs {
                summary: &summary,
                detail: &detail,
                params: &params,
                coverage: coverage.as_deref(),
                run_name: &run_name,
                output: &output,
            })?;
        }
    }

    Ok(())
}
