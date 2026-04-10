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

        /// Simulation mode: thermodynamic (TNN Boltzmann weighting) or simple (uniform probe-site weights)
        #[arg(long, default_value = "thermodynamic")]
        simulate_mode: SimulateMode,

        /// Hybridization temperature in °C (only used with --simulate-mode thermodynamic)
        #[arg(long, default_value = "70.0")]
        hybridization_temperature: f64,

        /// Fraction of fragments derived from probe binding sites (0.0–1.0).
        /// The remaining fraction are background fragments sampled uniformly by sequence weight.
        #[arg(long, default_value = "0.5")]
        capture_fraction: f64,

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

    /// Simulate probe-biased fragments using thermodynamic or simple weighting.
    /// Aligns probes to the reference, scores each binding site, and generates
    /// fragments via multinomial sampling. Background fragments fill the remainder.
    Simulate {
        /// Combined reference FASTA
        #[arg(short, long)]
        reference: PathBuf,

        /// Weights file
        #[arg(short, long)]
        weights: PathBuf,

        /// Probe sequences FASTA
        #[arg(short, long)]
        probes: PathBuf,

        /// Number of fragments to generate
        #[arg(short, long)]
        num_fragments: usize,

        /// Fraction of fragments derived from probe binding sites (0.0–1.0)
        #[arg(long, default_value = "0.5")]
        capture_fraction: f64,

        /// Simulation mode: thermodynamic (TNN Boltzmann weighting) or simple (uniform probe-site weights)
        #[arg(long, default_value = "thermodynamic")]
        simulate_mode: SimulateMode,

        /// Hybridization temperature in °C (only used with --simulate-mode thermodynamic)
        #[arg(long, default_value = "70.0")]
        hybridization_temperature: f64,

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

        /// Number of threads for minimap2 probe alignment
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

        /// Step (stride) between candidate probes during tiling. Used with --method catch-lite
        /// and --method catch. Smaller values generate more candidates and typically fewer final probes.
        #[arg(long, default_value = "60")]
        catch_stride: usize,

        /// Maximum Hamming-distance mismatches allowed for a probe to cover a target window.
        /// Used with --method catch-lite and --method catch.
        #[arg(long, default_value = "5")]
        catch_mismatches: usize,

        /// Extend the covered interval by this many bp on each side of a probe match
        /// (only used with --method catch-lite).
        #[arg(long, default_value = "0")]
        catch_extension: usize,

        /// Minimum fraction of each target sequence that must be covered
        /// (only used with --method catch-lite).
        #[arg(long, default_value = "1.0")]
        catch_coverage: f64,

        /// Jaccard similarity threshold for MinHash near-duplicate removal of candidate probes;
        /// 0.0 disables deduplication (only used with --method catch-lite).
        #[arg(long, default_value = "0.6")]
        catch_minhash_threshold: f64,

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

        /// 1X coverage threshold (%) below which targets are re-analyzed in refinement iterations (assessment step)
        #[arg(long, default_value = "80.0")]
        refine_threshold: f64,

        /// Number of refinement iterations: re-run probe coverage on low-coverage targets N times (assessment step)
        #[arg(long, conflicts_with = "refine_until_stable")]
        refine_iterations: Option<usize>,

        /// Repeat refinement until no targets remain below --refine-threshold or set stabilizes (assessment step)
        #[arg(long, conflicts_with = "refine_iterations")]
        refine_until_stable: bool,

        /// Maximum Hamming distance for a bait to cover a reference window (only used with --method syotti-lite)
        #[arg(long, default_value = "40")]
        syotti_mismatches: usize,

        /// Seed length (k-mer size) for the syotti approximate search (only used with --method syotti-lite)
        #[arg(long, default_value = "20")]
        syotti_seed_len: usize,
    },

    /// Standalone utility tools (probe design algorithms, sequence analysis).
    /// Run `baitbench tool --help` to see available tools.
    Tool {
        #[command(subcommand)]
        command: ToolCommands,
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

    /// Generate coverage depth curves, optionally sweeping CT, capture-fraction, and/or num-sequences
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

        /// Capture-fraction values to sweep (space-separated). Conflicts with --capture-fraction.
        #[arg(long, num_args = 1.., conflicts_with = "capture_fraction")]
        capture_fraction_values: Option<Vec<f64>>,

        /// Fixed capture fraction (when not sweeping). Conflicts with --capture-fraction-values.
        #[arg(long, default_value = "0.5", conflicts_with = "capture_fraction_values")]
        capture_fraction: f64,

        /// Simulation mode: thermodynamic or simple
        #[arg(long, default_value = "thermodynamic")]
        simulate_mode: SimulateMode,

        /// Hybridization temperature in °C (thermodynamic mode only)
        #[arg(long, default_value = "70.0")]
        hybridization_temperature: f64,

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

#[derive(Subcommand)]
pub enum ToolCommands {
    /// Design probes using the Syotti greedy set-cover algorithm
    Syotti {
        /// Input target sequences FASTA
        #[arg(short, long)]
        targets: PathBuf,

        /// Output probe sequences FASTA
        #[arg(short, long)]
        output: PathBuf,

        /// Probe (bait) length in bp
        #[arg(long, default_value = "120")]
        probe_length: usize,

        /// Maximum Hamming distance for a bait to cover a reference window
        #[arg(long, default_value = "40")]
        mismatches: usize,

        /// Seed length (k-mer size) for approximate matching
        #[arg(long, default_value = "20")]
        seed_len: usize,
    },

    /// Design probes using the CATCH optimization algorithm
    Catch {
        /// Input target sequences FASTA
        #[arg(short, long)]
        targets: PathBuf,

        /// Output probe sequences FASTA
        #[arg(short, long)]
        output: PathBuf,

        /// Probe (bait) length in bp
        #[arg(long, default_value = "120")]
        probe_length: usize,

        /// Stride (tiling step) in bp
        #[arg(long, default_value = "60")]
        stride: usize,

        /// Maximum mismatches allowed for a probe to cover a window
        #[arg(long, default_value = "5")]
        mismatches: usize,

        /// Extension length on each side of candidate probes
        #[arg(long, default_value = "0")]
        extension: usize,

        /// Fraction of each target that must be covered
        #[arg(long, default_value = "1.0")]
        coverage: f64,

        /// MinHash Jaccard similarity threshold for deduplication
        #[arg(long, default_value = "0.6")]
        minhash_threshold: f64,
    },

    /// Visualize sDUST low-complexity masking on FASTA sequences
    Dustview {
        /// Input FASTA file (reads from stdin if omitted)
        input: Option<PathBuf>,

        /// DUST score threshold — positions above this are masked
        #[arg(long, default_value = "2.0")]
        dust_threshold: f64,

        /// DUST sliding window size in bases
        #[arg(long, default_value = "64")]
        dust_window: usize,
    },

    /// Collapse near-duplicate sequences using cd-hit-est
    Collapse {
        /// Input FASTA file
        #[arg(short, long)]
        input: PathBuf,

        /// Output FASTA file (clustered representatives)
        #[arg(short, long)]
        output: PathBuf,

        /// Sequence identity threshold for clustering
        #[arg(long, default_value = "0.95")]
        threshold: f64,

        /// Number of threads for cd-hit-est
        #[arg(long, default_value = "1")]
        threads: usize,

        /// Path to write cd-hit-est log output
        #[arg(long, default_value = "cdhit.log")]
        log_file: PathBuf,
    },
}

#[derive(Clone, Copy, ValueEnum, Debug, PartialEq, Eq)]
pub enum ProbeMethod {
    /// Tile: sliding window probes with configurable overlap
    Tile,
    /// catch-lite: native Rust reimplementation of CATCH optimization-based probe design
    /// (Metsky et al. 2019, doi:10.1038/s41587-018-0006-x)
    #[value(name = "catch-lite")]
    CatchLite,
    /// syotti-lite: native Rust reimplementation of Syotti greedy set-cover bait design
    /// (Alanko et al. 2022, doi:10.1093/bioinformatics/btac226)
    #[value(name = "syotti-lite")]
    SyottiLite,
    /// catch: external CATCH tool (requires catch conda package)
    /// Uses --catch-stride and --catch-mismatches; other --catch-* flags are ignored
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

#[derive(Clone, Copy, ValueEnum, Debug, PartialEq, Eq)]
pub enum SimulateMode {
    /// Thermodynamic TNN Boltzmann weighting (default)
    Thermodynamic,
    /// Uniform probe-site weights — no temperature dependency
    Simple,
}
