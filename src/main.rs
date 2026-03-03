mod alignment;
mod cli;
mod commands;
mod external;
mod fasta;
mod io_utils;
mod sampling;

use anyhow::Result;
use clap::Parser;

use cli::{Cli, Commands};
use commands::{capture, enrich, filter, generate_list, map_reads, metrics, prepare, probe_coverage, report, run, sequence, simulate};

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
                distractor_fraction,
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
            outdir,
        } => {
            prepare::execute(&prepare::PrepareArgs {
                targets: &targets,
                distractors: &distractors,
                sample: sample.as_deref(),
                distractor_fraction,
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
