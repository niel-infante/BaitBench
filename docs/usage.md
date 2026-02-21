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
  --num-reads 10000 \
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

Without `--sample`, all targets are treated as present with equal weight. When `--sample` is provided, only the listed targets generate reads; remaining targets become "non-sample targets" and are treated as negatives alongside distractors (see [3-way classification](#3-way-classification)).

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
| `--num-reads` | 10000 | Number of reads to simulate |
| `--distractor-fraction` | 0.9 | Fraction of reads from distractors (0-1) |
| `--seed` | random | Random seed for reproducibility |

### Capture Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--capture-method` | minimap2 | Capture method: `minimap2` or `blast` |
| `--max-mismatches` | 10 | Maximum mismatches (minimap2 only) |
| `--min-match-bases` | 60 | Minimum matching bases required |

### Host Filtering (Optional)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--host-fasta` | none | Host genome for filtering captured reads |

### Output

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--outdir` | ./results | Output directory |
| `--run-name` | auto | Run name (auto-generated timestamp if not specified) |
| `--no-report` | false | Skip HTML report generation |

## Output Files

```
results/<run_name>/
├── combined_reference.fa   # Merged targets + distractors
├── weights.txt             # Sampling weights
├── targets.txt             # Target IDs
├── distractors.txt         # Distractor IDs
├── sample.txt              # Sample IDs (subset of targets)
├── reads.fa                # Simulated reads
├── captured.fa             # Reads passing capture filter
├── mapped.sam              # Alignments to references
├── detected.list           # Reference IDs and read counts
├── results.tsv             # Summary metrics
├── detected_detail.tsv     # Per-reference breakdown
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
| `num_reads` | Number of reads requested to be generated |
| `seed` | Random seed used for reproducibility (or "NA" if random) |

### Read Counts

| Column | Description |
|--------|-------------|
| `reads_generated` | Total number of reads actually generated from reference sequences |
| `reads_captured` | Number of reads that passed the capture filter (matched probes) |
| `capture_rate` | Fraction of generated reads that were captured (`reads_captured / reads_generated`) |
| `sample_captured` | Captured reads originating from **sample target** sequences |
| `nonsample_target_captured` | Captured reads originating from **non-sample target** sequences (cross-reactivity within the target panel) |
| `distractor_captured` | Captured reads originating from **distractor** sequences (off-target capture) |
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

- **Read-level** (`sample_captured`, `nonsample_target_captured`, `distractor_captured`, `reads_correctly_mapped`, `reads_incorrectly_mapped`): Since each simulated read is labeled with its source genome, these columns track how reads flow through the pipeline. A read from virus A that maps to virus B is counted as incorrectly mapped — this catches cross-reactivity even when genome-level metrics look perfect.

**Example:** If you have 2 sample genomes and 1000 reads are captured from them:
- `sample_captured` = 1000 (reads)
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
baitbench simulate  # Generate weighted random fragments
baitbench capture   # Probe capture (minimap2 or BLAST)
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
  --num-reads 10000 \
  --seed 42 \
  --outdir results

# 2. With sample manifest (subset of targets, custom weights)
baitbench run \
  --targets all_viruses.fa \
  --distractors bacteria.fa --distractors fungi.fa \
  --probes probes.fa \
  --sample sample.tsv \
  --num-reads 10000 \
  --seed 42 \
  --outdir results

# 3. More stringent capture
baitbench run \
  --targets targets.fa \
  --distractors distractors.fa \
  --probes probes.fa \
  --num-reads 10000 \
  --min-match-bases 70 \
  --max-mismatches 5 \
  --outdir results/stringent

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
- Use `--no-report` to skip report generation

### Low capture rate

If very few reads are being captured:
- Try reducing `--min-match-bases`
- Try increasing `--max-mismatches`
- Verify your probes match your targets

### High false positive rate

If many distractors or non-sample targets are being captured:
- Try increasing `--min-match-bases`
- Try decreasing `--max-mismatches`
- Consider adding host filtering with `--host-fasta`
- Review which genomes appear in the `detected_detail.tsv` to distinguish FP_target from FP_distractor
