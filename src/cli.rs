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

        /// Genome sequences FASTA for fragment generation (optional).
        /// When provided, fragments are generated from genomes instead of targets.
        /// Use for organisms where the genome is much larger than the target regions
        /// (e.g., bacteria with specific gene targets).
        #[arg(short, long)]
        genomes: Option<PathBuf>,

        /// Distractor sequences FASTA (can be specified multiple times)
        #[arg(short, long, num_args = 1..)]
        distractors: Vec<PathBuf>,

        /// Probe sequences FASTA
        #[arg(short, long)]
        probes: PathBuf,

        /// Sample targets: either a manifest TSV file path, or inline IDs with optional weights.
        /// Examples: --sample manifest.tsv | --sample t1 t2 t3 | --sample t1 t2 t3 5 t4
        /// Inline: IDs default to weight 1.0; a number after an ID sets that ID's weight.
        /// When --genomes is used, IDs refer to genome sequences instead of targets.
        #[arg(long, num_args = 1..)]
        sample: Option<Vec<String>>,

        /// Sample-to-target mapping TSV (optional, used with --genomes).
        /// Format: genome_id<tab>target_id (one mapping per line).
        /// Maps genome sequences to their target regions for metrics classification.
        /// Genomes with matching target IDs are auto-linked when this file is absent.
        #[arg(long)]
        sample_target_map: Option<PathBuf>,

        /// Host genome FASTA for filtering (optional)
        #[arg(long)]
        host_fasta: Option<PathBuf>,

        /// Run name (auto-generated if not specified)
        #[arg(long)]
        run_name: Option<String>,

        /// Number of fragments to simulate
        #[arg(short, long, default_value = "10000")]
        num_fragments: usize,

        /// Fraction of reads from distractors (0-1). Mutually exclusive with --ct.
        #[arg(long, conflicts_with = "ct")]
        distractor_fraction: Option<f64>,

        /// CT (cycle threshold) score to determine distractor fraction.
        /// Mutually exclusive with --distractor-fraction.
        /// Converts via: target_fraction = ct_baseline_fraction * 2^(ct_baseline - ct)
        #[arg(long, conflicts_with = "distractor_fraction")]
        ct: Option<f64>,

        /// CT baseline value (CT at which target fraction equals --ct-baseline-fraction)
        #[arg(long, default_value = "20.0")]
        ct_baseline: f64,

        /// Target fraction at the baseline CT value
        #[arg(long, default_value = "0.01")]
        ct_baseline_fraction: f64,

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

        /// Number of sequences to sample in sequencing step (with replacement). If not specified, all captured fragments become reads.
        #[arg(long)]
        num_sequences: Option<usize>,

        /// Output directory
        #[arg(short, long, default_value = "./results")]
        outdir: PathBuf,

        /// Number of threads for external tools
        #[arg(long, default_value = "1")]
        threads: usize,

        /// Fold enrichment for capture (e.g. 100 = 100x more target relative to distractor post-capture vs pre-capture). Omit for binary capture (default behavior).
        #[arg(long)]
        fold_enrichment: Option<f64>,

        /// Report output mode: full (HTML report), none (skip), rmd (editable RMarkdown file)
        #[arg(long, default_value = "full")]
        report: ReportMode,

        /// Delete intermediate files after pipeline completes, keeping only report inputs and final outputs
        #[arg(long)]
        cleanup: bool,
    },

    /// Combine target and distractor FASTAs, generate weights
    Prepare {
        /// Target sequences FASTA
        #[arg(short, long)]
        targets: PathBuf,

        /// Genome sequences FASTA for fragment generation (optional)
        #[arg(short, long)]
        genomes: Option<PathBuf>,

        /// Distractor sequences FASTA (can be specified multiple times)
        #[arg(short, long, num_args = 1..)]
        distractors: Vec<PathBuf>,

        /// Sample targets: either a manifest TSV file path, or inline IDs with optional weights.
        /// Examples: --sample manifest.tsv | --sample t1 t2 t3 | --sample t1 t2 t3 5 t4
        /// When --genomes is used, IDs refer to genome sequences instead of targets.
        #[arg(long, num_args = 1..)]
        sample: Option<Vec<String>>,

        /// Sample-to-target mapping TSV (optional, used with --genomes)
        #[arg(long)]
        sample_target_map: Option<PathBuf>,

        /// Fraction of reads from distractors (0-1). Mutually exclusive with --ct.
        #[arg(short = 'f', long, conflicts_with = "ct")]
        distractor_fraction: Option<f64>,

        /// CT (cycle threshold) score to determine distractor fraction.
        /// Mutually exclusive with --distractor-fraction.
        #[arg(long, conflicts_with = "distractor_fraction")]
        ct: Option<f64>,

        /// CT baseline value (CT at which target fraction equals --ct-baseline-fraction)
        #[arg(long, default_value = "20.0")]
        ct_baseline: f64,

        /// Target fraction at the baseline CT value
        #[arg(long, default_value = "0.01")]
        ct_baseline_fraction: f64,

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

    /// Adjust captured fragment pool to a target fold enrichment
    Enrich {
        /// Captured fragments FASTA (from capture step)
        #[arg(short, long)]
        captured: PathBuf,

        /// All fragments FASTA (from simulate step)
        #[arg(short, long)]
        fragments: PathBuf,

        /// Targets ID file
        #[arg(short, long)]
        targets: PathBuf,

        /// Distractors ID file
        #[arg(short, long)]
        distractors: PathBuf,

        /// Fold enrichment (>= 1.0; 1.0 = no enrichment)
        #[arg(long)]
        fold_enrichment: f64,

        /// Random seed
        #[arg(short, long)]
        seed: Option<u64>,

        /// Output enriched fragments FASTA
        #[arg(short, long)]
        output: PathBuf,
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

        /// Number of sequences to sample (with replacement). If not specified, all fragments pass through.
        #[arg(long)]
        num_sequences: Option<usize>,

        /// Random seed for sampling reproducibility
        #[arg(short, long)]
        seed: Option<u64>,
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

        /// Output per-position coverage TSV (optional)
        #[arg(long)]
        output_coverage: Option<PathBuf>,

        /// Number of reads after sequencing step (for pipeline flow tracking)
        #[arg(long)]
        reads_sequenced: Option<usize>,

        /// Number of reads after host filtering (for pipeline flow tracking)
        #[arg(long)]
        reads_after_filter: Option<usize>,
    },

    /// Analyze probe tiling and coverage across target sequences (probe design QC)
    ProbeCoverage {
        /// Target sequences FASTA
        #[arg(short, long)]
        targets: PathBuf,

        /// Probe sequences FASTA
        #[arg(short, long)]
        probes: PathBuf,

        /// Output directory
        #[arg(short, long, default_value = "./probe_coverage")]
        outdir: PathBuf,

        /// Minimap2 alignment preset
        #[arg(long, default_value = "sr")]
        minimap_preset: String,

        /// Proximity distance (bp) for pull-down zone metric
        #[arg(long, default_value = "50")]
        proximity: usize,

        /// Report output mode: full (HTML report), none (skip), rmd (editable RMarkdown file)
        #[arg(long, default_value = "full")]
        report: ReportMode,

        /// Delete intermediate files after completion, keeping only report inputs and final outputs
        #[arg(long)]
        cleanup: bool,
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

        /// Coverage profile TSV (optional, for coverage plots)
        #[arg(long)]
        coverage: Option<PathBuf>,

        /// Run name
        #[arg(long, default_value = "BaitBench Run")]
        run_name: String,

        /// Output file (HTML for full mode, RMarkdown for rmd mode)
        #[arg(short, long)]
        output: PathBuf,

        /// Report output mode: full (HTML report), rmd (editable RMarkdown file)
        #[arg(long, default_value = "full")]
        report: ReportMode,
    },

    /// Analyze probe cross-reactivity against genomes or other probes
    Xreact {
        /// Probe sequences FASTA
        #[arg(short, long)]
        probes: PathBuf,

        /// Reference genome FASTA(s) to check cross-reactivity against (can be specified multiple times)
        #[arg(long, num_args = 1..)]
        against: Option<Vec<PathBuf>>,

        /// Check probe-vs-probe cross-reactivity (self-hits excluded)
        #[arg(long = "self")]
        self_mode: bool,

        /// Minimum homology percentage to report (matching_bases / probe_length * 100)
        #[arg(long, default_value = "80.0")]
        threshold: f64,

        /// Minimap2 alignment preset
        #[arg(long, default_value = "sr")]
        minimap_preset: String,

        /// Output directory
        #[arg(short, long, default_value = "./xreact_results")]
        outdir: PathBuf,

        /// Report output mode: full (HTML report), none (skip), rmd (editable RMarkdown file)
        #[arg(long, default_value = "full")]
        report: ReportMode,

        /// Delete intermediate files after completion, keeping only report inputs and final outputs
        #[arg(long)]
        cleanup: bool,
    },

    /// Generate coverage depth curves, optionally sweeping CT, fold-enrichment, and/or num-sequences
    CoverageCurve {
        /// Target sequences FASTA
        #[arg(short, long)]
        targets: PathBuf,

        /// Genome sequences FASTA for fragment generation (optional)
        #[arg(short, long)]
        genomes: Option<PathBuf>,

        /// Distractor sequences FASTA (can be specified multiple times)
        #[arg(short, long, num_args = 1..)]
        distractors: Vec<PathBuf>,

        /// Probe sequences FASTA
        #[arg(short, long)]
        probes: PathBuf,

        /// Sample targets (required): either a manifest TSV file path, or inline IDs with optional weights.
        /// Examples: --sample manifest.tsv | --sample t1 t2 t3 | --sample t1 t2 t3 5 t4
        /// When --genomes is used, IDs refer to genome sequences instead of targets.
        #[arg(long, num_args = 1.., required = true)]
        sample: Vec<String>,

        /// Sample-to-target mapping TSV (optional, used with --genomes)
        #[arg(long)]
        sample_target_map: Option<PathBuf>,

        /// CT values to sweep (space-separated). Conflicts with --ct and --distractor-fraction.
        #[arg(long, num_args = 1.., conflicts_with_all = ["ct", "distractor_fraction"])]
        ct_values: Option<Vec<f64>>,

        /// Fixed CT value (when not sweeping CT). Conflicts with --ct-values and --distractor-fraction.
        #[arg(long, conflicts_with_all = ["ct_values", "distractor_fraction"])]
        ct: Option<f64>,

        /// Fixed distractor fraction (when not sweeping CT). Conflicts with --ct-values and --ct.
        #[arg(long, conflicts_with_all = ["ct_values", "ct"])]
        distractor_fraction: Option<f64>,

        /// CT baseline value
        #[arg(long, default_value = "20.0")]
        ct_baseline: f64,

        /// Target fraction at baseline CT
        #[arg(long, default_value = "0.01")]
        ct_baseline_fraction: f64,

        /// Fold-enrichment values to sweep (space-separated). Conflicts with --fold-enrichment.
        #[arg(long, num_args = 1.., conflicts_with = "fold_enrichment")]
        fold_enrichment_values: Option<Vec<f64>>,

        /// Fixed fold enrichment (when not sweeping). Conflicts with --fold-enrichment-values.
        #[arg(long, conflicts_with = "fold_enrichment_values")]
        fold_enrichment: Option<f64>,

        /// Num-sequences values to sweep (space-separated). Conflicts with --num-sequences.
        #[arg(long, num_args = 1.., conflicts_with = "num_sequences")]
        num_sequences_values: Option<Vec<usize>>,

        /// Fixed num-sequences (when not sweeping). Conflicts with --num-sequences-values.
        #[arg(long, conflicts_with = "num_sequences_values")]
        num_sequences: Option<usize>,

        /// Number of fragments to simulate
        #[arg(short, long, default_value = "10000")]
        num_fragments: usize,

        /// Sequencing read length
        #[arg(long, default_value = "120")]
        read_length: usize,

        /// Random seed
        #[arg(short, long)]
        seed: Option<u64>,

        /// Mean fragment length
        #[arg(long, default_value = "175")]
        fragment_length_mean: f64,

        /// Minimum fragment length
        #[arg(long, default_value = "150")]
        fragment_length_min: usize,

        /// Maximum fragment length
        #[arg(long, default_value = "200")]
        fragment_length_max: usize,

        /// Capture method
        #[arg(long, default_value = "minimap2")]
        capture_method: CaptureMethodArg,

        /// Max mismatches for minimap2 capture
        #[arg(long, default_value = "10")]
        max_mismatches: u32,

        /// Min matching bases required
        #[arg(long, default_value = "60")]
        min_match_bases: u32,

        /// BLAST database path
        #[arg(long)]
        blast_db: Option<String>,

        /// Host genome FASTA for filtering
        #[arg(long)]
        host_fasta: Option<PathBuf>,

        /// Minimap2 preset for read mapping
        #[arg(long, default_value = "sr")]
        minimap_preset: String,

        /// Minimap2 preset for host filtering
        #[arg(long, default_value = "sr")]
        host_minimap_preset: String,

        /// Number of threads for external tools
        #[arg(long, default_value = "1")]
        threads: usize,

        /// Output directory
        #[arg(short, long, default_value = "./coverage_curve_results")]
        outdir: PathBuf,

        /// Report output mode: full (HTML report), none (skip), rmd (editable RMarkdown file)
        #[arg(long, default_value = "full")]
        report: ReportMode,

        /// Delete intermediate files after completion, keeping only report inputs and final outputs
        #[arg(long)]
        cleanup: bool,
    },
}

#[derive(Clone, Copy, ValueEnum, Debug, PartialEq, Eq)]
pub enum ReportMode {
    /// Generate full HTML report (default)
    Full,
    /// Skip report generation
    None,
    /// Output parameterized RMarkdown file for manual editing and rendering
    Rmd,
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
