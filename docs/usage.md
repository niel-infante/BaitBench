# BaitBench Usage Guide

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

The binary is at `target/release/baitbench`. Add it to your PATH or symlink it somewhere convenient.

## Quick Start

```bash
baitbench run \
  --targets targets.fa \
  --distractors distractors.fa \
  --probes probes.fa \
  --num-fragments 10000 \
  --outdir results
```

## Input Files

### targets.fa

FASTA file containing sequences your probes are designed to capture. These are your positive controls — you expect these to be detected.

```fasta
>virus_A Complete genome
ATGCGTACGT...
>virus_B Partial sequence
GCTAGCTAG...
```

### distractors.fa

FASTA file(s) containing background sequences that should NOT be captured. These serve as negative controls to test specificity. You can specify multiple distractor files.

Common choices for distractors:
- Host genome sequences
- Related but non-target organisms
- Environmental/metagenomic background

### probes.fa

FASTA file containing your probe sequences to test.

```fasta
>probe_001 Target: virus_A position 100-180
ATGCGTACGTACGTACGTACGTACGTACGTACGTACGTACG
>probe_002 Target: virus_A position 500-580
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG
```

### sample.tsv (optional)

A TSV file listing which target IDs are present in the simulated sample. Each ID must match a FASTA header in the targets file.

**Important:** Sequence IDs are derived from the first whitespace-delimited word of each FASTA header (everything after `>` up to the first space). **Sequence names must not contain spaces.** Use underscores or other delimiters instead (e.g. `>Zika_virus`, not `>Zika virus`). This applies to all FASTA input files.

The weight column is optional (default: 1.0) and controls relative read abundance:

| id | weight |
|----|--------|
| `dengue_1` | `5.0` |
| `zika_virus` | `1.0` |
| `chikungunya` | `0.5` |

Example `sample.tsv`:

```
# id	weight
dengue_1	5.0
zika_virus	1.0
chikungunya	0.5
```

Without `--sample`, all targets are treated as present with equal weight. When `--sample` is provided, only the listed targets generate fragments; remaining targets become "non-sample targets" and are treated as negatives alongside distractors (see [3-way classification](#3-way-classification)).

## Parameters Reference

### Required Parameters

| Parameter | Description |
|-----------|-------------|
| `--targets` | Path to target genomes FASTA |
| `--distractors` | Path to distractor genomes FASTA (can be specified multiple times) |
| `--probes` | Path to probe sequences FASTA |

### Sample Selection

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--sample` | none | Sample manifest TSV (id and optional weight) |

### Simulation Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--num-fragments` | 10000 | Number of fragments to simulate |
| `--fragment-length-mean` | 175 | Mean fragment length (bp) |
| `--fragment-length-min` | 150 | Minimum fragment length (bp) |
| `--fragment-length-max` | 200 | Maximum fragment length (bp) |
| `--read-length` | 120 | Sequencing read length (trim captured fragments to this) |
| `--distractor-fraction` | 0.9 | Fraction of fragments from distractors (0-1) |
| `--seed` | random | Random seed for reproducibility |

### Host Filtering (Optional)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--host-fasta` | none | Host genome for filtering captured reads |

### Output

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--outdir` | ./results | Output directory |
| `--run-name` | auto | Run name (auto-generated timestamp if not specified) |
| `--report` | full | Report output: `full` (HTML), `none` (skip), `rmd` (editable RMarkdown) |

## Output Files

```
results/<run_name>/
├── combined_reference.fa   # Merged targets + distractors (+ genomes in genome mode)
├── weights.txt             # Sampling weights
├── targets.txt             # Target IDs
├── distractors.txt         # Distractor IDs
├── sample.txt              # Sample IDs (subset of targets)
├── distractor_groups.tsv   # Distractor group assignments (always written)
├── target_groups.tsv       # Target group assignments (if --groups)
├── fragments.fa            # Simulated fragments (probe-biased + background)
├── reads.fa                # Sequencing reads
├── mapped.sam              # Alignments to references
├── detected.list           # Reference IDs and read counts
├── results.tsv             # Summary metrics
├── detected_detail.tsv     # Per-reference breakdown (with group and coverage stats)
├── group_detail.tsv        # Per-group summary (if groups are present)
├── coverage.tsv            # Run-length encoded read depth intervals
├── results.json            # Machine-readable metrics
└── report.html             # HTML report with figures (if R available)
```

## results.tsv Column Reference

The `results.tsv` file contains one row with all summary metrics.

### Run Information

| Column | Description |
|--------|-------------|
| `run_name` | Name of this pipeline run (auto-generated or user-specified) |
| `timestamp` | ISO 8601 timestamp when metrics were calculated |
| `num_fragments` | Number of fragments requested to be generated |
| `seed` | Random seed used for reproducibility (or "NA" if random) |

### Fragment and Read Counts

| Column | Description |
|--------|-------------|
| `fragments_generated` | Total number of fragments actually generated from reference sequences |
| `fragments_captured` | Number of fragments that passed the capture filter (matched probes) |
| `capture_rate` | Fraction of generated fragments that were captured (`fragments_captured / fragments_generated`) |
| `sample_captured` | Captured fragments originating from **sample target** sequences |
| `nonsample_target_captured` | Captured fragments originating from **non-sample target** sequences (cross-reactivity within the target panel) |
| `distractor_captured` | Captured fragments originating from **distractor** sequences (off-target capture) |
| `reads_correctly_mapped` | Reads that map back to their source reference |
| `reads_incorrectly_mapped` | Reads that map to a different reference than their source |

### Reference Counts

| Column | Description |
|--------|-------------|
| `sample_total` | Number of distinct sample genomes |
| `nonsample_target_total` | Number of distinct non-sample target genomes |
| `distractors_total` | Number of distinct distractor genomes |

### Detection Classification (Genome-Level)

These metrics count **distinct genomes/references**, not reads. A genome is "detected" if at least one read maps back to it after capture.

| Column | Description |
|--------|-------------|
| `tp_count` | **True Positives** — sample target genomes detected |
| `fn_count` | **False Negatives** — sample target genomes NOT detected |
| `fp_target_count` | **False Positives (target)** — non-sample target genomes detected |
| `fp_distractor_count` | **False Positives (distractor)** — distractor genomes detected |
| `fp_total` | Total false positives (`fp_target_count + fp_distractor_count`) |
| `tn_target_count` | **True Negatives (target)** — non-sample target genomes NOT detected |
| `tn_distractor_count` | **True Negatives (distractor)** — distractor genomes NOT detected |
| `tn_total` | Total true negatives (`tn_target_count + tn_distractor_count`) |

### Performance Metrics

| Column | Description |
|--------|-------------|
| `sensitivity` | `TP / (TP + FN)` — fraction of sample genomes detected. Range: 0–1. |
| `specificity` | `TN_total / (TN_total + FP_total)` — fraction of non-sample genomes correctly rejected. Range: 0–1. |
| `precision` | `TP / (TP + FP_total)` — of detected genomes, fraction that are sample targets. Range: 0–1. |
| `f1_score` | `2 * (precision * sensitivity) / (precision + sensitivity)` — harmonic mean. Range: 0–1. |

## detected_detail.tsv Column Reference

The `detected_detail.tsv` file contains one row per reference (target or distractor) with detection status and coverage statistics.

| Column | Description |
|--------|-------------|
| `reference_id` | Reference sequence ID |
| `group` | Group name this sequence belongs to (sequence's own ID if no groups file provided) |
| `category` | `sample`, `nonsample_target`, `distractor`, or `untargeted` |
| `expected` | 1 if expected to be detected (sample target), 0 otherwise |
| `detected` | 1 if at least one read maps to this reference, 0 otherwise |
| `fragments_generated` | Number of fragments generated from this reference |
| `fragments_captured` | Number of fragments captured by probes |
| `reads_assigned` | Number of reads mapped to this reference |
| `classification` | `TP`, `FN`, `FP_target`, `FP_distractor`, `TN_target`, `TN_distractor`, or `untargeted` |
| `ref_length` | Reference sequence length (bp) |
| `avg_coverage` | Average read depth across reference |
| `pct_covered_5x` | % positions with >= 5x depth |
| `pct_covered_20x` | % positions with >= 20x depth |

## coverage.tsv

The `coverage.tsv` file contains run-length encoded read depth intervals for each reference. Consecutive positions with the same depth are collapsed into a single interval (1-based inclusive coordinates). This file is used by the HTML report for coverage profile plots.

| Column | Description |
|--------|-------------|
| `reference_id` | Reference sequence ID |
| `start` | 1-based start position of the interval (inclusive) |
| `end` | 1-based end position of the interval (inclusive) |
| `depth` | Coverage depth for this interval |

## Understanding the Metrics

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

This distinguishes between two types of false positives:
- **FP_target**: detecting a non-sample target (cross-reactivity within the target panel)
- **FP_distractor**: detecting a distractor (off-target capture)

Without `--sample`, all targets are in the sample, so there are no non-sample targets and the classification reduces to the traditional 2-way TP/FP/FN/TN.

### Genome-level vs read-level metrics

BaitBench reports two levels of metrics:

- **Genome-level** (`tp_count`, `fp_*`, `fn_count`, `tn_*`, and derived rates): Was each genome detected at all? A genome is detected if at least one read maps to it after capture and mapping.

- **Fragment/Read-level** (`sample_captured`, `nonsample_target_captured`, `distractor_captured`, `reads_correctly_mapped`, `reads_incorrectly_mapped`): Since each simulated fragment is labeled with its source genome, these columns track how fragments and reads flow through the pipeline. Capture counts are at the fragment level; mapping counts are at the read level (post-sequencing). A read from virus A that maps to virus B is counted as incorrectly mapped — this catches cross-reactivity even when genome-level metrics look perfect.

**Example:** If you have 2 sample genomes and 1000 fragments are captured from them:
- `sample_captured` = 1000 (fragments)
- `tp_count` = 2 (genomes), assuming both were detected

### Sensitivity (True Positive Rate)

```
Sensitivity = TP / (TP + FN)
```

Measures what fraction of sample genomes were successfully detected. High sensitivity means the probes are capturing most of their intended targets.

### Specificity (True Negative Rate)

```
Specificity = TN_total / (TN_total + FP_total)
```

Measures what fraction of non-sample genomes (both non-sample targets and distractors) were correctly NOT captured. High specificity means the probes are avoiding off-target capture.

### Precision

```
Precision = TP / (TP + FP_total)
```

Of all detected genomes, what fraction were actual sample targets? High precision means most detections are real targets, not false positives.

### F1 Score

```
F1 = 2 * (Precision * Sensitivity) / (Precision + Sensitivity)
```

Harmonic mean of precision and sensitivity. A balanced measure of overall performance.

## Individual Steps

Each pipeline step is available as a subcommand:

```bash
baitbench prepare   # Combine FASTAs, generate weights
baitbench simulate  # Probe alignment + thermodynamic scoring + fragment generation
baitbench sequence  # Simulate sequencing (trim/error-model fragments to reads)
baitbench filter    # Optional host read filtering
baitbench map       # Map reads back to reference
baitbench list      # Count reads per reference from SAM
baitbench metrics   # Calculate TP/FP/FN/TN
baitbench report    # Generate HTML report (requires R)
```

Run `baitbench <command> --help` for full options.

## Example Workflow

```bash
# 1. Basic run (all targets in sample)
baitbench run \
  --targets targets.fa \
  --distractors distractors.fa \
  --probes probes.fa \
  --num-fragments 10000 \
  --seed 42 \
  --outdir results

# 2. With sample manifest (subset of targets, custom weights)
baitbench run \
  --targets all_viruses.fa \
  --distractors bacteria.fa --distractors fungi.fa \
  --probes probes.fa \
  --sample sample.tsv \
  --num-fragments 10000 \
  --seed 42 \
  --outdir results

# 3. With CT score simulation (low viral load at CT 30)
baitbench run \
  --targets targets.fa \
  --distractors distractors.fa \
  --probes probes.fa \
  --num-fragments 10000 \
  --ct 30 \
  --seed 42 \
  --outdir results/ct30

# 4. Compare results
cat results/*/results.tsv
```

## Troubleshooting

### "Sample ID 'X' not found in targets FASTA"

The IDs in your sample manifest don't match the FASTA headers. Sequence IDs are the first word of each FASTA header (after `>`), so names with spaces are truncated. For example, `>Zika virus` becomes ID `Zika`. Use underscores: `>Zika_virus`.

### "Cannot find R scripts directory"

The HTML report requires R and the `R/` directory. Options:
- Run from the project root directory
- Set `BAITBENCH_R_DIR` to the path of the `R/` directory
- Use `--report none` to skip report generation

### Low capture rate

If very few reads are being captured:
- Check that probes actually match your target sequences
- Try increasing `--capture-fraction` to generate more probe-biased fragments
- Verify your probes FASTA is in the correct format with unique IDs

### High false positive rate

If many distractors or non-sample targets are being captured:
- Try lowering `--capture-fraction` to reduce probe-biased fragments and rely more on background
- Consider adding host filtering with `--host-fasta`
- Review which genomes appear in the `detected_detail.tsv` to distinguish FP_target from FP_distractor
