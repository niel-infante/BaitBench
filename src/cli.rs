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

        /// String to prepend to every output filename
        #[arg(long, default_value = "")]
        output_prefix: String,

        /// Number of threads for external tools
        #[arg(long, default_value = "1")]
        threads: usize,

        /// Fold enrichment for capture (e.g. 100 = 100x more target relative to distractor post-capture vs pre-capture). Omit for binary capture (default behavior).
        #[arg(long)]
        fold_enrichment: Option<f64>,

        /// Enable species-level identification after metrics (genome mode only, requires --sample-target-map).
        /// Calls species PRESENT/ABSENT/AMBIGUOUS based on multi-target detection patterns.
        #[arg(long)]
        identify: bool,

        /// Minimum sequence identity % to consider targets "similar" for species identification
        #[arg(long, default_value = "90.0")]
        identity_threshold: f64,

        /// Minimum unique targets detected to call a species PRESENT
        #[arg(long, default_value = "1")]
        min_unique_targets: usize,

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

        /// String to prepend to every output filename
        #[arg(long, default_value = "")]
        output_prefix: String,
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

        /// String to prepend to every output filename
        #[arg(long, default_value = "")]
        output_prefix: String,

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

        /// String to prepend to every output filename
        #[arg(long, default_value = "")]
        output_prefix: String,

        /// Report output mode: full (HTML report), none (skip), rmd (editable RMarkdown file)
        #[arg(long, default_value = "full")]
        report: ReportMode,

        /// Delete intermediate files after completion, keeping only report inputs and final outputs
        #[arg(long)]
        cleanup: bool,
    },

    /// Assess target panel discriminability for species-level identification
    PanelQc {
        /// Target sequences FASTA
        #[arg(short, long)]
        targets: PathBuf,

        /// Sample-to-target mapping TSV (genome_id<tab>target_id).
        /// Maps species/genome IDs to their target regions.
        #[arg(long)]
        sample_target_map: PathBuf,

        /// Minimum sequence identity percentage to consider two targets "similar"
        #[arg(long, default_value = "90.0")]
        identity_threshold: f64,

        /// Minimap2 alignment preset for target-vs-target comparison
        #[arg(long, default_value = "sr")]
        minimap_preset: String,

        /// Output directory
        #[arg(short, long, default_value = "./panel_qc_results")]
        outdir: PathBuf,

        /// String to prepend to every output filename
        #[arg(long, default_value = "")]
        output_prefix: String,

        /// Report output mode: full (HTML report), none (skip), rmd (editable RMarkdown file)
        #[arg(long, default_value = "full")]
        report: ReportMode,

        /// Delete intermediate files after completion
        #[arg(long)]
        cleanup: bool,
    },

    /// Call species presence/absence from multi-target detection patterns
    Identify {
        /// Detected detail TSV from metrics step (detected_detail.tsv)
        #[arg(long)]
        detected_detail: PathBuf,

        /// Sample-to-target mapping TSV (genome_id<tab>target_id)
        #[arg(long)]
        sample_target_map: PathBuf,

        /// Pre-computed target similarity TSV (from panel-qc). If not provided, --targets is required.
        #[arg(long)]
        target_similarity: Option<PathBuf>,

        /// Target sequences FASTA (used to compute similarity on-the-fly if --target-similarity not provided)
        #[arg(short, long)]
        targets: Option<PathBuf>,

        /// Minimum sequence identity percentage for target similarity (when computing on-the-fly)
        #[arg(long, default_value = "90.0")]
        identity_threshold: f64,

        /// Minimap2 alignment preset (when computing similarity on-the-fly)
        #[arg(long, default_value = "sr")]
        minimap_preset: String,

        /// Minimum unique targets detected to call a species PRESENT
        #[arg(long, default_value = "1")]
        min_unique_targets: usize,

        /// Output directory
        #[arg(short, long, default_value = "./identify_results")]
        outdir: PathBuf,

        /// String to prepend to every output filename
        #[arg(long, default_value = "")]
        output_prefix: String,
    },

    /// Build a probe set from target sequences (collapse, tile, filter, deduplicate)
    /// and optionally assess probe quality (coverage + cross-reactivity)
    BuildProbes {
        /// Target sequences FASTA
        #[arg(short, long)]
        targets: PathBuf,

        /// Probe construction method
        #[arg(long, default_value = "tile")]
        method: ProbeMethod,

        /// Probe length in bp
        #[arg(long, default_value = "120")]
        probe_length: usize,

        /// Step from the end of each probe before starting the next one.
        /// Negative = overlap (default -60 gives 60bp overlap), 0 = perfectly tiled,
        /// positive = gap between probes. Only used with --method tile.
        #[arg(long, default_value = "-60", allow_hyphen_values = true)]
        step: i64,

        /// Additional arguments passed to CATCH's design.py (only used with --method catch).
        /// Default provides: probe stride 60, 5 mismatches, no cover extension,
        /// and LSH minhash deduplication at 0.6.
        #[arg(long, default_value = "-ps 60 -m 5 -e 0 --filter-with-lsh-minhash 0.6", allow_hyphen_values = true)]
        catch_args: String,

        /// Minimum GC fraction to keep a probe (0-1)
        #[arg(long, default_value = "0.20")]
        min_gc: f64,

        /// Maximum GC fraction to keep a probe (0-1)
        #[arg(long, default_value = "0.80")]
        max_gc: f64,

        /// Maximum fraction of ambiguous (N) bases allowed in a target sequence (0-1).
        /// Sequences exceeding this threshold are removed before probe construction.
        #[arg(long, default_value = "0.05")]
        max_n_frac: f64,

        /// DUST score threshold for low-complexity filtering (Morgulis et al. 2006)
        #[arg(long, default_value = "2.0")]
        dust_threshold: f64,

        /// DUST window size in bases for low-complexity filtering
        #[arg(long, default_value = "64")]
        dust_window: usize,

        /// Maximum fraction of bases masked by DUST to keep a probe (0-1).
        /// Probes exceeding this fraction of low-complexity bases are removed.
        /// Set to 1.0 to disable complexity filtering.
        #[arg(long, default_value = "0.25")]
        max_masked_frac: f64,

        /// cd-hit-est identity threshold for initial target collapse
        #[arg(long, default_value = "0.95")]
        collapse_threshold: f64,

        /// cd-hit-est identity threshold for final probe deduplication
        #[arg(long, default_value = "0.95")]
        dedup_threshold: f64,

        /// Number of threads for cd-hit-est
        #[arg(long, default_value = "5")]
        threads: usize,

        /// Genome FASTA(s) for cross-reactivity checking during assessment (can be specified multiple times)
        #[arg(short, long, num_args = 1..)]
        genomes: Option<Vec<PathBuf>>,

        /// Minimum homology percentage for cross-reactivity detection during assessment
        #[arg(long, default_value = "80.0")]
        threshold: f64,

        /// Minimap2 alignment preset for assessment steps
        #[arg(long, default_value = "sr")]
        minimap_preset: String,

        /// Proximity distance (bp) for pull-down zone metric in coverage assessment
        #[arg(long, default_value = "50")]
        proximity: usize,

        /// Skip the probe assessment step (coverage + cross-reactivity analysis)
        #[arg(long)]
        skip_assess: bool,

        /// Output directory
        #[arg(short, long, default_value = "./build_probes_results")]
        outdir: PathBuf,

        /// String to prepend to every output filename
        #[arg(long, default_value = "")]
        output_prefix: String,

        /// Report output mode: full (HTML report), none (skip), rmd (editable RMarkdown file)
        #[arg(long, default_value = "full")]
        report: ReportMode,

        /// Delete intermediate files after completion
        #[arg(long)]
        cleanup: bool,
    },

    /// Assess probe quality: coverage analysis + cross-reactivity (self + optional genome)
    AssessProbes {
        /// Target sequences FASTA
        #[arg(short, long)]
        targets: PathBuf,

        /// Probe sequences FASTA
        #[arg(short, long)]
        probes: PathBuf,

        /// Genome FASTA(s) for cross-reactivity checking (can be specified multiple times)
        #[arg(short, long, num_args = 1..)]
        genomes: Option<Vec<PathBuf>>,

        /// Minimum homology percentage for cross-reactivity detection
        #[arg(long, default_value = "80.0")]
        threshold: f64,

        /// Minimap2 alignment preset
        #[arg(long, default_value = "sr")]
        minimap_preset: String,

        /// Proximity distance (bp) for pull-down zone metric in coverage analysis
        #[arg(long, default_value = "50")]
        proximity: usize,

        /// Output directory
        #[arg(short, long, default_value = "./assess_probes_results")]
        outdir: PathBuf,

        /// String to prepend to every output filename
        #[arg(long, default_value = "")]
        output_prefix: String,

        /// Report output mode: full (HTML report), none (skip), rmd (editable RMarkdown file)
        #[arg(long, default_value = "full")]
        report: ReportMode,

        /// Delete intermediate files after completion
        #[arg(long)]
        cleanup: bool,

        /// 1X coverage threshold (%) below which targets are re-analyzed in refinement iterations
        #[arg(long, default_value = "80.0")]
        refine_threshold: f64,

        /// Number of refinement iterations: re-run probe coverage on low-coverage targets N times
        #[arg(long, conflicts_with = "refine_until_stable")]
        refine_iterations: Option<usize>,

        /// Repeat refinement until no targets remain below --refine-threshold (or set stabilizes)
        #[arg(long, conflicts_with = "refine_iterations")]
        refine_until_stable: bool,
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

        /// String to prepend to every output filename
        #[arg(long, default_value = "")]
        output_prefix: String,

        /// Report output mode: full (HTML report), none (skip), rmd (editable RMarkdown file)
        #[arg(long, default_value = "full")]
        report: ReportMode,

        /// Delete intermediate files after completion, keeping only report inputs and final outputs
        #[arg(long)]
        cleanup: bool,
    },
}

#[derive(Clone, Copy, ValueEnum, Debug, PartialEq, Eq)]
pub enum ProbeMethod {
    /// Tile: sliding window probes with configurable overlap
    Tile,
    /// CATCH: optimization-based probe design from the Broad Institute
    Catch,
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
