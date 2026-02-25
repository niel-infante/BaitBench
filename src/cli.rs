use clap::{Parser, Subcommand, ValueEnum};
use std::path::PathBuf;

#[derive(Parser)]
#[command(name = "baitbench")]
#[command(about = "Probe capture efficiency testing via in-silico simulation")]
#[command(version)]
pub struct Cli {
    /// Enable verbose logging
    #[arg(short, long, global = true)]
    pub verbose: bool,

    #[command(subcommand)]
    pub command: Commands,
}

#[derive(Subcommand)]
pub enum Commands {
    /// Run the full pipeline
    Run {
        /// Target sequences FASTA
        #[arg(short, long)]
        targets: PathBuf,

        /// Distractor sequences FASTA (can be specified multiple times)
        #[arg(short, long, num_args = 1..)]
        distractors: Vec<PathBuf>,

        /// Probe sequences FASTA
        #[arg(short, long)]
        probes: PathBuf,

        /// Sample manifest TSV (id<tab>weight; subset of targets present in sample)
        #[arg(long)]
        sample: Option<PathBuf>,

        /// Host genome FASTA for filtering (optional)
        #[arg(long)]
        host_fasta: Option<PathBuf>,

        /// Run name (auto-generated if not specified)
        #[arg(long)]
        run_name: Option<String>,

        /// Number of fragments to simulate
        #[arg(short, long, default_value = "10000")]
        num_fragments: usize,

        /// Fraction of reads from distractors (0-1)
        #[arg(long, default_value = "0.9")]
        distractor_fraction: f64,

        /// Random seed for reproducibility
        #[arg(short, long)]
        seed: Option<u64>,

        /// Capture method
        #[arg(long, default_value = "minimap2")]
        capture_method: CaptureMethodArg,

        /// Max mismatches for minimap2 capture
        #[arg(long, default_value = "10")]
        max_mismatches: u32,

        /// Min matching bases required
        #[arg(long, default_value = "60")]
        min_match_bases: u32,

        /// BLAST database path (required if capture-method=blast)
        #[arg(long)]
        blast_db: Option<String>,

        /// Minimap2 preset for read mapping
        #[arg(long, default_value = "sr")]
        minimap_preset: String,

        /// Minimap2 preset for host filtering
        #[arg(long, default_value = "sr")]
        host_minimap_preset: String,

        /// Mean fragment length
        #[arg(long, default_value = "175")]
        fragment_length_mean: f64,

        /// Minimum fragment length
        #[arg(long, default_value = "150")]
        fragment_length_min: usize,

        /// Maximum fragment length
        #[arg(long, default_value = "200")]
        fragment_length_max: usize,

        /// Sequencing read length (trim captured fragments to this length)
        #[arg(long, default_value = "120")]
        read_length: usize,

        /// Output directory
        #[arg(short, long, default_value = "./results")]
        outdir: PathBuf,

        /// Number of threads for external tools
        #[arg(long, default_value = "1")]
        threads: usize,

        /// Skip HTML report generation
        #[arg(long)]
        no_report: bool,
    },

    /// Combine target and distractor FASTAs, generate weights
    Prepare {
        /// Target sequences FASTA
        #[arg(short, long)]
        targets: PathBuf,

        /// Distractor sequences FASTA (can be specified multiple times)
        #[arg(short, long, num_args = 1..)]
        distractors: Vec<PathBuf>,

        /// Sample manifest TSV (id<tab>weight; subset of targets present in sample)
        #[arg(long)]
        sample: Option<PathBuf>,

        /// Fraction of reads from distractors (0-1)
        #[arg(short = 'f', long, default_value = "0.9")]
        distractor_fraction: f64,

        /// Output directory
        #[arg(short, long, default_value = ".")]
        outdir: PathBuf,
    },

    /// Generate weighted random fragments from FASTA
    Simulate {
        /// Combined reference FASTA
        #[arg(short, long)]
        reference: PathBuf,

        /// Weights file
        #[arg(short, long)]
        weights: PathBuf,

        /// Number of fragments to generate
        #[arg(short, long)]
        num_fragments: usize,

        /// Random seed
        #[arg(short, long)]
        seed: Option<u64>,

        /// Output FASTA file
        #[arg(short, long)]
        output: PathBuf,

        /// Mean fragment length
        #[arg(long, default_value = "175")]
        fragment_length_mean: f64,

        /// Minimum fragment length
        #[arg(long, default_value = "150")]
        fragment_length_min: usize,

        /// Maximum fragment length
        #[arg(long, default_value = "200")]
        fragment_length_max: usize,
    },

    /// Simulate probe capture using minimap2 or BLAST
    Capture {
        /// Probe sequences FASTA
        #[arg(short, long)]
        probes: PathBuf,

        /// Fragments FASTA to capture
        #[arg(short, long)]
        fragments: PathBuf,

        /// Capture method
        #[arg(long, default_value = "minimap2")]
        method: CaptureMethodArg,

        /// Max mismatches (minimap2)
        #[arg(long, default_value = "10")]
        max_mismatches: u32,

        /// Min matching bases
        #[arg(long, default_value = "60")]
        min_match_bases: u32,

        /// BLAST database path
        #[arg(long)]
        blast_db: Option<String>,

        /// Output captured fragments FASTA
        #[arg(short, long)]
        output: PathBuf,

        /// Log file
        #[arg(long, default_value = "capture.log")]
        log_file: PathBuf,

        /// Number of threads
        #[arg(long, default_value = "1")]
        threads: usize,
    },

    /// Simulate sequencing of captured fragments (trim to read length)
    Sequence {
        /// Input captured fragments FASTA
        #[arg(short, long)]
        input: PathBuf,

        /// Output reads FASTA
        #[arg(short, long)]
        output: PathBuf,

        /// Read length to trim to
        #[arg(long, default_value = "120")]
        read_length: usize,
    },

    /// Filter out host reads using minimap2
    Filter {
        /// Host genome FASTA
        #[arg(long)]
        host: PathBuf,

        /// Reads FASTA
        #[arg(short, long)]
        reads: PathBuf,

        /// Minimap2 preset
        #[arg(long, default_value = "sr")]
        minimap_preset: String,

        /// Output filtered reads FASTA
        #[arg(short, long)]
        output: PathBuf,

        /// Log file
        #[arg(long, default_value = "host_filter.log")]
        log_file: PathBuf,
    },

    /// Map reads back to reference
    Map {
        /// Reference FASTA
        #[arg(short, long)]
        reference: PathBuf,

        /// Reads FASTA
        #[arg(long)]
        reads: PathBuf,

        /// Minimap2 preset
        #[arg(long, default_value = "sr")]
        minimap_preset: String,

        /// Output SAM file
        #[arg(short, long)]
        output: PathBuf,

        /// Log file
        #[arg(long, default_value = "mapping.log")]
        log_file: PathBuf,
    },

    /// Generate detection list from SAM file
    List {
        /// Input SAM file
        #[arg(short, long)]
        sam: PathBuf,

        /// Output detection list
        #[arg(short, long)]
        output: PathBuf,
    },

    /// Calculate TP/FP/FN/TN metrics
    Metrics {
        /// Targets ID file
        #[arg(long)]
        targets: PathBuf,

        /// Distractors ID file
        #[arg(long)]
        distractors: PathBuf,

        /// Sample ID file (subset of targets present in sample)
        #[arg(long)]
        sample: PathBuf,

        /// Detection list file
        #[arg(long)]
        detected: PathBuf,

        /// Generated fragments FASTA
        #[arg(long)]
        fragments: PathBuf,

        /// Captured fragments FASTA
        #[arg(long)]
        captured: PathBuf,

        /// Mapped reads SAM file
        #[arg(long)]
        sam: PathBuf,

        /// Run name
        #[arg(long)]
        run_name: String,

        /// Number of fragments requested
        #[arg(long)]
        num_fragments: usize,

        /// Random seed used
        #[arg(long, default_value = "NA")]
        seed: String,

        /// Output summary TSV
        #[arg(long)]
        output_summary: PathBuf,

        /// Output detail TSV
        #[arg(long)]
        output_detail: PathBuf,

        /// Output JSON (optional)
        #[arg(long)]
        output_json: Option<PathBuf>,
    },

    /// Generate HTML report with ggplot2 figures
    Report {
        /// Summary TSV file (results.tsv)
        #[arg(long)]
        summary: PathBuf,

        /// Detail TSV file (detected_detail.tsv)
        #[arg(long)]
        detail: PathBuf,

        /// Run parameters file (run_params.tsv)
        #[arg(long)]
        params: PathBuf,

        /// Run name
        #[arg(long, default_value = "BaitBench Run")]
        run_name: String,

        /// Output HTML file
        #[arg(short, long)]
        output: PathBuf,
    },
}

#[derive(Clone, Copy, ValueEnum)]
pub enum CaptureMethodArg {
    Minimap2,
    Blast,
}

impl From<CaptureMethodArg> for crate::commands::capture::CaptureMethod {
    fn from(arg: CaptureMethodArg) -> Self {
        match arg {
            CaptureMethodArg::Minimap2 => crate::commands::capture::CaptureMethod::Minimap2,
            CaptureMethodArg::Blast => crate::commands::capture::CaptureMethod::Blast,
        }
    }
}
