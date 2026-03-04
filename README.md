# BaitBench

A tool for testing probe capture efficiency via in-silico simulation.

## Overview

BaitBench evaluates how well a probe set captures target sequences while avoiding off-target (distractor) sequences. It simulates the capture process and reports metrics including sensitivity, specificity, precision, and F1 score.

A key feature is the **sample manifest**, which lets you specify a subset of targets as the "sample" — simulating which viruses are actually present in a specimen. This tests not only whether probes capture targets over distractors, but whether they can discriminate between viruses within the target set.

## Installation

### 1. Install runtime dependencies

```bash
conda env create -f environment.yml
conda activate baitbench
```

This installs minimap2, BLAST, R, and required R packages.

### 2. Build the binary

Requires the [Rust toolchain](https://rustup.rs/).

```bash
cargo build --release
```

The binary is at `target/release/baitbench`. Add it to your PATH or copy it somewhere convenient.

## Usage

### Full pipeline

```bash
baitbench run \
  --targets targets.fa \
  --distractors distractors.fa \
  --probes probes.fa \
  --num-fragments 10000 \
  --outdir results
```

### With sample targets

Specify which targets are "present" in the sample and at what abundance. `--sample` accepts either a TSV file or inline IDs directly on the command line:

```bash
# Using a TSV manifest file
baitbench run \
  --targets all_viruses.fa \
  --distractors bacteria.fa --distractors fungi.fa \
  --probes probes.fa \
  --sample sample.tsv \
  --num-fragments 10000 \
  --outdir results

# Inline IDs (all default to weight 1.0)
baitbench run \
  --targets all_viruses.fa \
  --distractors bacteria.fa \
  --probes probes.fa \
  --sample dengue_1 zika_virus chikungunya \
  --num-fragments 10000 \
  --outdir results

# Inline IDs with weights (number after an ID sets its weight)
baitbench run \
  --targets all_viruses.fa \
  --distractors bacteria.fa \
  --probes probes.fa \
  --sample dengue_1 5 zika_virus chikungunya 0.5 \
  --num-fragments 10000 \
  --outdir results
```

### Sample format

`--sample` accepts two formats:

**Inline IDs**: List target IDs directly. A number following an ID sets its weight (default: 1.0).

```
--sample dengue_1 5 zika_virus chikungunya 0.5
```
Result: `dengue_1` (weight 5.0), `zika_virus` (weight 1.0), `chikungunya` (weight 0.5).

**TSV file**: A file with one ID per line, optional tab-separated weight. If a single argument is given and it's an existing file, it's parsed as a TSV manifest.

```
# id	weight
dengue_1	5.0
zika_virus	1.0
chikungunya	0.5
```

All IDs must match a FASTA header in the targets file. IDs are taken from the first whitespace-delimited word of each FASTA header (everything after `>` up to the first space), so **sequence names must not contain spaces**. Use underscores or other delimiters instead (e.g. `>Zika_virus`, not `>Zika virus`).

Without `--sample`, all targets are treated as present with equal weight. When `--sample` is provided, only the listed targets generate fragments; remaining targets become "non-sample targets" and are treated as negatives alongside distractors (see [3-way classification](#3-way-classification)).

### Individual steps

Each pipeline step is available as a subcommand:

```bash
baitbench prepare   # Combine FASTAs, generate weights
baitbench simulate  # Generate weighted random fragments
baitbench capture   # Probe capture (minimap2 or BLAST)
baitbench enrich    # Fold enrichment adjustment (optional)
baitbench sequence  # Simulate sequencing (trim fragments to read length)
baitbench filter    # Optional host read filtering
baitbench map       # Map reads back to reference
baitbench list      # Count reads per reference from SAM
baitbench metrics   # Calculate TP/FP/FN/TN
baitbench report    # Generate HTML report (requires R)
```

Run `baitbench <command> --help` for full options.

### Probe coverage analysis

Separate from the simulation pipeline, `probe-coverage` evaluates how well your probes tile across target sequences — a probe design QC tool.

```bash
baitbench probe-coverage \
  --targets targets.fa \
  --probes probes.fa \
  --outdir probe_coverage_results
```

This maps probes to targets using minimap2 and reports per-target statistics:

| Metric | Description |
|--------|-------------|
| `pct_covered_1x` | % of target bases covered by at least one probe |
| `mean_depth` / `median_depth` | Average and median probe depth |
| `pct_covered_2x/5x/10x` | Tiered coverage thresholds |
| `max_gap_length` / `num_gaps` | Largest uncovered stretch and number of gaps |
| `pct_near_probe` | % of bases within proximity distance of a probe (default 50bp) |

Output files: `probe_coverage_summary.tsv` (per-target stats), `probe_depth.tsv` (per-position depth), and optionally `probe_coverage_report.html` (requires R).

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--targets` | required | Target sequences FASTA |
| `--probes` | required | Probe sequences FASTA |
| `--outdir` | ./probe_coverage | Output directory |
| `--minimap-preset` | sr | Minimap2 alignment preset |
| `--proximity` | 50 | Distance (bp) for pull-down zone metric |
| `--no-report` | false | Skip HTML report generation |

### CT sweep analysis

`ct-sweep` runs the full pipeline at multiple CT values and generates coverage depth curves, showing how sequencing depth translates to genome coverage at each target concentration.

```bash
# Using inline sample IDs
baitbench ct-sweep \
  --targets all_viruses.fa \
  --distractors bacteria.fa \
  --probes probes.fa \
  --sample dengue_1 zika_virus \
  --ct-values 20 25 30 35 \
  --num-fragments 10000 \
  --seed 42 \
  --outdir ct_sweep_results

# Or using a TSV manifest file
baitbench ct-sweep \
  --targets all_viruses.fa \
  --distractors bacteria.fa \
  --probes probes.fa \
  --sample sample.tsv \
  --ct-values 20 25 30 35 \
  --num-fragments 10000 \
  --seed 42 \
  --outdir ct_sweep_results
```

This runs the pipeline (prepare through mapping) for each CT value, then generates a plot of % genome covered (Y-axis) vs depth of coverage on a log10 scale (X-axis), with one line per CT value. If the sample contains multiple targets, each gets its own panel.

Output files:

```
ct_sweep_results/
├── ct_20.0/              # Pipeline intermediates for CT 20
│   ├── fragments.fa, captured.fa, reads.fa, mapped.sam, ...
├── ct_25.0/              # Pipeline intermediates for CT 25
│   └── ...
├── ct_sweep_depth_curves.tsv   # Aggregated depth curve data
└── ct_sweep_report.html        # HTML report with depth curve plots
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--targets` | required | Target sequences FASTA |
| `--distractors` | required | Distractor sequences FASTA (can be specified multiple times) |
| `--probes` | required | Probe sequences FASTA |
| `--sample` | required | Sample targets to track: TSV file or inline IDs (e.g. `--sample t1 t2`) |
| `--ct-values` | required | CT values to sweep (space-separated) |
| `--ct-baseline` | 20.0 | CT baseline value |
| `--ct-baseline-fraction` | 0.01 | Target fraction at baseline CT |
| `--num-fragments` | 10000 | Number of fragments per CT run |
| `--num-sequences` | all | Number of sequences to sample per CT run |
| `--outdir` | ./ct_sweep_results | Output directory |

All other pipeline parameters (`--read-length`, `--fold-enrichment`, `--capture-method`, etc.) are also supported.

## Input Files

| File | Description |
|------|-------------|
| `targets.fa` | FASTA of genomes your probes are designed to capture |
| `distractors.fa` | FASTA of background genomes that should NOT be captured (can specify multiple files) |
| `probes.fa` | FASTA of probe sequences to test |
| `sample.tsv` (optional) | TSV listing which targets are present in the sample, with optional weights |
| `host.fa` (optional) | Host genome for filtering captured reads |

**Note:** Sequence IDs are derived from the first word of each FASTA header (text after `>` up to the first space). Sequence names must not contain spaces. These IDs are used throughout the pipeline and must match between input files (e.g. sample manifest IDs must match target FASTA header IDs).

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--targets` | required | Path to target genomes FASTA |
| `--distractors` | required | Path to distractor genomes FASTA (can be specified multiple times) |
| `--probes` | required | Path to probe sequences FASTA |
| `--sample` | none | Sample targets: TSV manifest file or inline IDs with optional weights (e.g. `--sample t1 t2 t3 5 t4`) |
| `--num-fragments` | 10000 | Number of fragments to simulate |
| `--fragment-length-mean` | 175 | Mean fragment length (bp) |
| `--fragment-length-min` | 150 | Minimum fragment length (bp) |
| `--fragment-length-max` | 200 | Maximum fragment length (bp) |
| `--read-length` | 120 | Sequencing read length (trim captured fragments to this) |
| `--num-sequences` | all | Number of sequences to sample in sequencing step (with replacement). If not specified, all captured fragments become reads. |
| `--distractor-fraction` | 0.9 | Fraction of fragments from distractors (0-1). Mutually exclusive with `--ct`. |
| `--ct` | none | CT (cycle threshold) score — converts to distractor fraction (see below). Mutually exclusive with `--distractor-fraction`. |
| `--ct-baseline` | 20.0 | CT value at which the target fraction equals `--ct-baseline-fraction` |
| `--ct-baseline-fraction` | 0.01 | Target fraction at the baseline CT value |
| `--capture-method` | minimap2 | Capture method: `minimap2` or `blast` |
| `--min-match-bases` | 60 | Minimum matching bases for capture |
| `--fold-enrichment` | none | Fold enrichment for capture (e.g. 100 = 100x; omit for binary capture) |
| `--max-mismatches` | 10 | Maximum mismatches allowed (minimap2 only) |
| `--host-fasta` | none | Optional host genome for filtering |
| `--seed` | random | Random seed for reproducibility |
| `--no-report` | false | Skip HTML report generation |
| `--outdir` | ./results | Output directory |

### Using CT scores instead of distractor fraction

Instead of specifying `--distractor-fraction` directly, you can use a qPCR **CT (cycle threshold) score** via `--ct`. This is more intuitive for users in a diagnostic setting — a lower CT means higher viral load (more target DNA relative to background).

The conversion formula is:

```
target_fraction = ct_baseline_fraction × 2^(ct_baseline − ct)
distractor_fraction = 1 − target_fraction
```

With the defaults (`--ct-baseline 20.0`, `--ct-baseline-fraction 0.01`):

| CT | Target fraction | Distractor fraction |
|----|----------------|---------------------|
| 15 | 32% | 0.68 |
| 20 | 1% | 0.99 |
| 25 | 0.03% | 0.9997 |
| 30 | 0.001% | 0.99999 |
| 35 | 0.00003% | ~1.0 |

**Advanced calibration:** The `--ct-baseline` and `--ct-baseline-fraction` flags let you tune the mapping between CT and target fraction to match your own experimental observations. For example, if in your lab a CT of 25 corresponds to roughly 0.1% target reads:

```bash
baitbench run \
  --targets targets.fa \
  --distractors distractors.fa \
  --probes probes.fa \
  --ct 30 \
  --ct-baseline 25 \
  --ct-baseline-fraction 0.001 \
  --num-fragments 10000 \
  --outdir results
```

## Output Files

```
results/<run_name>/
├── combined_reference.fa   # Merged targets + distractors
├── weights.txt             # Sampling weights
├── targets.txt             # Target IDs
├── distractors.txt         # Distractor IDs
├── sample.txt              # Sample IDs (subset of targets)
├── fragments.fa             # Simulated fragments
├── captured.fa              # Fragments passing capture filter
├── enriched.fa              # Post-enrichment fragments (if --fold-enrichment)
├── reads.fa                 # Sequenced reads (fragments trimmed to read length)
├── mapped.sam               # Alignments to references
├── detected.list           # Reference IDs and read counts
├── run_params.tsv          # Run parameters (key-value)
├── results.tsv             # Summary metrics
├── detected_detail.tsv     # Per-reference breakdown (with coverage stats)
├── coverage.tsv            # Per-position coverage depth for detected references
├── results.json            # Machine-readable metrics
└── report.html             # HTML report with figures (if R available)
```

## Results Columns (`results.tsv`)

| Column | Description |
|--------|-------------|
| `run_name` | Name of the run (auto-generated or user-specified) |
| `timestamp` | When the run completed |
| `num_fragments` | Number of fragments requested |
| `seed` | Random seed used (or "NA" if random) |
| `fragments_generated` | Actual number of fragments generated |
| `fragments_captured` | Total fragments passing capture filter |
| `capture_rate` | `fragments_captured / fragments_generated` |
| `sample_captured` | Captured fragments originating from sample target sequences |
| `nonsample_target_captured` | Captured fragments originating from non-sample target sequences |
| `distractor_captured` | Captured fragments originating from distractor sequences |
| `reads_correctly_mapped` | Reads that map back to their source reference |
| `reads_incorrectly_mapped` | Reads that map to a different reference than their source |
| `sample_total` | Number of distinct sample genomes |
| `nonsample_target_total` | Number of distinct non-sample target genomes |
| `distractors_total` | Number of distinct distractor genomes |
| `tp_count` | True Positives: sample genomes detected |
| `fn_count` | False Negatives: sample genomes NOT detected |
| `fp_target_count` | False Positives (target): non-sample target genomes detected |
| `fp_distractor_count` | False Positives (distractor): distractor genomes detected |
| `fp_total` | Total false positives (`fp_target_count + fp_distractor_count`) |
| `tn_target_count` | True Negatives (target): non-sample target genomes NOT detected |
| `tn_distractor_count` | True Negatives (distractor): distractor genomes NOT detected |
| `tn_total` | Total true negatives (`tn_target_count + tn_distractor_count`) |
| `sensitivity` | `TP / (TP + FN)` — fraction of sample genomes detected |
| `specificity` | `TN_total / (TN_total + FP_total)` — fraction of non-sample genomes correctly rejected |
| `precision` | `TP / (TP + FP_total)` — of detected genomes, fraction that are sample targets |
| `f1_score` | `2 * (precision * sensitivity) / (precision + sensitivity)` |

### Per-Reference Detail (`detected_detail.tsv`)

Each row corresponds to one reference sequence, with detection status and coverage statistics:

| Column | Description |
|--------|-------------|
| `reference_id` | Reference sequence ID |
| `category` | `sample`, `nonsample_target`, or `distractor` |
| `expected` / `detected` | Whether detection was expected and whether it occurred |
| `fragments_generated` / `fragments_captured` | Fragment counts for this reference |
| `reads_assigned` | Number of reads mapped to this reference |
| `classification` | `TP`, `FN`, `FP_target`, `FP_distractor`, `TN_target`, or `TN_distractor` |
| `ref_length` | Reference sequence length (bp) |
| `avg_coverage` | Average coverage depth across the reference |
| `pct_covered_5x` | Percentage of positions with >= 5X coverage |
| `pct_covered_20x` | Percentage of positions with >= 20X coverage |

### 3-way classification

BaitBench uses a 3-way genome classification:

| Category | Detected | Classification |
|----------|----------|----------------|
| Sample target | Yes | TP |
| Sample target | No | FN |
| Non-sample target | Yes | FP_target |
| Non-sample target | No | TN_target |
| Distractor | Yes | FP_distractor |
| Distractor | No | TN_distractor |

This distinguishes between two types of false positives: detecting a non-sample target (cross-reactivity within the target panel) vs detecting a distractor (off-target capture). Without `--sample`, all targets are in the sample, so there are no non-sample targets and the classification reduces to the traditional 2-way TP/FP/FN/TN.

### Genome-level vs read-level metrics

BaitBench reports two levels of metrics:

- **Genome-level** (`tp_count`, `fp_*`, `fn_count`, `tn_*`, and derived rates): Was each genome detected at all? A genome is detected if at least one read maps to it after capture and mapping.

- **Fragment/Read-level** (`sample_captured`, `nonsample_target_captured`, `distractor_captured`, `reads_correctly_mapped`, `reads_incorrectly_mapped`): Since each simulated fragment is labeled with its source genome, these columns track how fragments and reads flow through the pipeline. Capture counts are at the fragment level; mapping counts are at the read level (post-sequencing). A read from virus A that maps to virus B is counted as incorrectly mapped — this catches cross-reactivity even when genome-level metrics look perfect.

## Example

```bash
# Basic run (all targets in sample)
baitbench run \
  --targets examples/minimal/targets.fa \
  --distractors examples/minimal/distractors.fa \
  --probes examples/minimal/probes.fa \
  --num-fragments 1000 \
  --seed 42 \
  --outdir example_results

# With inline sample (subset of targets, custom weights)
baitbench run \
  --targets examples/minimal/targets.fa \
  --distractors examples/minimal/distractors.fa \
  --probes examples/minimal/probes.fa \
  --sample target_virus_1 5 \
  --num-fragments 1000 \
  --seed 42 \
  --outdir example_results

# With sample manifest file
echo -e "target_virus_1\t5.0" > sample.tsv
baitbench run \
  --targets examples/minimal/targets.fa \
  --distractors examples/minimal/distractors.fa \
  --probes examples/minimal/probes.fa \
  --sample sample.tsv \
  --num-fragments 1000 \
  --seed 42 \
  --outdir example_results

# Using a CT score (CT 25 ≈ 0.03% target)
baitbench run \
  --targets examples/minimal/targets.fa \
  --distractors examples/minimal/distractors.fa \
  --probes examples/minimal/probes.fa \
  --ct 25 \
  --num-fragments 1000 \
  --seed 42 \
  --outdir example_results

# With fold enrichment (simulate 100x enrichment)
baitbench run \
  --targets examples/minimal/targets.fa \
  --distractors examples/minimal/distractors.fa \
  --probes examples/minimal/probes.fa \
  --fold-enrichment 100 \
  --num-fragments 1000 \
  --seed 42 \
  --outdir example_results

# Multiple distractor files
baitbench run \
  --targets targets.fa \
  --distractors bacteria.fa --distractors fungi.fa --distractors protozoa.fa \
  --probes probes.fa \
  --num-fragments 10000 \
  --outdir results
```

## How It Works

1. **Prepare**: Combines target and distractor genomes into a single reference. If a sample manifest is provided, only sample targets get non-zero weights; non-sample targets get weight 0. Calculates distractor weights so that a configurable fraction of fragments come from distractors.

2. **Simulate**: Generates random fragments from the reference using weighted sampling. Fragment lengths follow a normal distribution (default: mean 175bp, range 150-200bp). Sequences with weight 0 produce no fragments.

3. **Capture**: Aligns simulated fragments against probe sequences using minimap2 or BLAST. Filters by matching bases, mismatches, and indels to simulate hybridization stringency.

4. **Enrich** (optional): If `--fold-enrichment` is specified, adjusts the post-capture fragment pool to match the requested fold enrichment. Fold enrichment is defined as the ratio of target:distractor proportions post-capture vs pre-capture. A fold enrichment of 100 means the target:distractor ratio is 100x higher after capture than before. This is achieved by subsampling captured distractors or adding back uncaptured distractors to hit the target ratio.

5. **Sequence**: Simulates sequencing by trimming captured fragments to a fixed read length (default: 120bp). Fragments shorter than the read length are kept as-is. With `--num-sequences`, a specified number of reads are sampled with replacement from the fragment library (simulating PCR amplification before sequencing) and renumbered with unique IDs.

6. **Filter** (optional): Removes reads that map to a host genome.

7. **Map**: Aligns reads back to the combined reference to identify which genomes they originated from.

8. **Metrics**: Computes genome-level metrics (3-way TP/FP/FN/TN classification across sample targets, non-sample targets, and distractors), fragment/read-level metrics (how many captured fragments came from each category, and whether mapped reads return to their correct source genome), and per-reference coverage statistics (average depth, breadth of coverage at >=5X and >=20X thresholds).

9. **Report**: Generates an HTML report with ggplot2 figures (capture summary, metrics bar chart, confusion matrix, per-reference lollipop chart, and coverage depth profiles for each detected reference).

## Dependencies

- [minimap2](https://github.com/lh3/minimap2) (alignment)
- [BLAST+](https://blast.ncbi.nlm.nih.gov/) (alternative capture method)
- [R](https://www.r-project.org/) with ggplot2, rmarkdown, dplyr, tidyr (report generation, optional)

All installable via `conda env create -f environment.yml`.

## License

MIT License - see [LICENSE](LICENSE) for details.
