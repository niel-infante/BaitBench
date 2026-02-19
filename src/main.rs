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
use commands::{capture, filter, generate_list, map_reads, metrics, prepare, report, run, simulate};

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
            host_fasta,
            run_name,
            num_reads,
            distractor_fraction,
            seed,
            capture_method,
            max_mismatches,
            min_match_bases,
            blast_db,
            minimap_preset,
            host_minimap_preset,
            read_length_mean,
            read_length_min,
            read_length_max,
            outdir,
            threads,
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
                host_fasta: host_fasta.as_deref(),
                run_name,
                num_reads,
                distractor_fraction,
                seed,
                capture_method: capture_method.into(),
                max_mismatches,
                min_match_bases,
                blast_db,
                minimap_preset,
                host_minimap_preset,
                read_length_mean,
                read_length_min,
                read_length_max,
                outdir: full_outdir,
                threads,
                no_report,
            })?;
        }

        Commands::Prepare {
            targets,
            distractors,
            distractor_fraction,
            outdir,
        } => {
            prepare::execute(&prepare::PrepareArgs {
                targets: &targets,
                distractors: &distractors,
                distractor_fraction,
                outdir: &outdir,
            })?;
        }

        Commands::Simulate {
            reference,
            weights,
            num_reads,
            seed,
            output,
            read_length_mean,
            read_length_min,
            read_length_max,
        } => {
            simulate::execute(&simulate::SimulateArgs {
                reference: &reference,
                weights: &weights,
                num_reads,
                seed,
                output: &output,
                read_length_mean,
                read_length_min,
                read_length_max,
            })?;
        }

        Commands::Capture {
            probes,
            reads,
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
                reads: &reads,
                max_mismatches,
                min_match_bases,
                blast_db: blast_db.as_deref(),
                output: &output,
                log_file: &log_file,
                threads,
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
            detected,
            reads,
            captured,
            run_name,
            num_reads,
            seed,
            output_summary,
            output_detail,
            output_json,
        } => {
            metrics::execute(&metrics::MetricsArgs {
                targets: &targets,
                distractors: &distractors,
                detected: &detected,
                reads: &reads,
                captured: &captured,
                run_name: &run_name,
                num_reads,
                seed: &seed,
                output_summary: &output_summary,
                output_detail: &output_detail,
                output_json: output_json.as_deref(),
            })?;
        }

        Commands::Report {
            summary,
            detail,
            run_name,
            output,
        } => {
            report::execute(&report::ReportArgs {
                summary: &summary,
                detail: &detail,
                run_name: &run_name,
                output: &output,
            })?;
        }
    }

    Ok(())
}
