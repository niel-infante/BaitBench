# BaitBench

A tool for testing probe capture efficiency via in-silico simulation.

## Overview

BaitBench evaluates how well a probe set captures target sequences while avoiding off-target (distractor) sequences. It simulates the capture process and reports metrics including sensitivity, specificity, precision, and F1 score.

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
  --num-reads 10000 \
  --outdir results
```

### Individual steps

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

## Input Files

| File | Description |
|------|-------------|
| `targets.fa` | FASTA of genomes your probes are designed to capture |
| `distractors.fa` | FASTA of background genomes that should NOT be captured |
| `probes.fa` | FASTA of probe sequences to test |
| `host.fa` (optional) | Host genome for filtering captured reads |

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--targets` | required | Path to target genomes FASTA |
| `--distractors` | required | Path to distractor genomes FASTA |
| `--probes` | required | Path to probe sequences FASTA |
| `--num-reads` | 10000 | Number of reads to simulate |
| `--distractor-fraction` | 0.9 | Fraction of reads from distractors (0-1) |
| `--capture-method` | minimap2 | Capture method: `minimap2` or `blast` |
| `--min-match-bases` | 60 | Minimum matching bases for capture |
| `--max-mismatches` | 10 | Maximum mismatches allowed (minimap2 only) |
| `--host-fasta` | none | Optional host genome for filtering |
| `--seed` | random | Random seed for reproducibility |
| `--no-report` | false | Skip HTML report generation |
| `--outdir` | ./results | Output directory |

## Output Files

```
results/<run_name>/
├── combined_reference.fa   # Merged targets + distractors
├── weights.txt             # Sampling weights
├── targets.txt             # Target IDs
├── distractors.txt         # Distractor IDs
├── reads.fa                # Simulated reads
├── captured.fa             # Reads passing capture filter
├── mapped.sam              # Alignments to references
├── detected.list           # Reference IDs and read counts
├── results.tsv             # Summary metrics
├── detected_detail.tsv     # Per-reference breakdown
├── results.json            # Machine-readable metrics
└── report.html             # HTML report with figures (if R available)
```

## Results Columns (`results.tsv`)

| Column | Description |
|--------|-------------|
| `run_name` | Name of the run (auto-generated or user-specified) |
| `timestamp` | When the run completed |
| `num_reads` | Number of reads requested |
| `seed` | Random seed used (or "NA" if random) |
| `reads_generated` | Actual number of reads generated |
| `reads_captured` | Total reads passing capture filter |
| `capture_rate` | `reads_captured / reads_generated` |
| `target_captured` | Captured reads originating from target sequences |
| `distractor_captured` | Captured reads originating from distractor sequences |
| `reads_correctly_mapped` | Reads that map back to their source reference |
| `reads_incorrectly_mapped` | Reads that map to a different reference than their source (read-level false positives) |
| `targets_total` | Number of distinct target genomes |
| `distractors_total` | Number of distinct distractor genomes |
| `tp_count` | True Positives: target genomes detected |
| `fp_count` | False Positives: distractor genomes detected |
| `fn_count` | False Negatives: target genomes NOT detected |
| `tn_count` | True Negatives: distractor genomes NOT detected |
| `sensitivity` | `TP / (TP + FN)` — fraction of targets detected |
| `specificity` | `TN / (TN + FP)` — fraction of distractors correctly rejected |
| `precision` | `TP / (TP + FP)` — of detected genomes, fraction that are targets |
| `f1_score` | `2 * (precision * sensitivity) / (precision + sensitivity)` |

### Genome-level vs read-level metrics

BaitBench reports two levels of metrics:

- **Genome-level** (`tp_count`, `fp_count`, `fn_count`, `tn_count`, and derived rates): Was each genome detected at all? A genome is detected if at least one read maps to it after capture and mapping.

- **Read-level** (`target_captured`, `distractor_captured`, `reads_correctly_mapped`, `reads_incorrectly_mapped`): Since each simulated read is labeled with its source genome, these columns track how reads flow through the pipeline. A read from virus A that maps to virus B is counted as incorrectly mapped — this catches cross-reactivity even when genome-level metrics look perfect.

## Example

```bash
baitbench run \
  --targets examples/minimal/targets.fa \
  --distractors examples/minimal/distractors.fa \
  --probes examples/minimal/probes.fa \
  --num-reads 1000 \
  --seed 42 \
  --outdir example_results
```

## How It Works

1. **Prepare**: Combines target and distractor genomes into a single reference. Calculates sampling weights so that a configurable fraction of reads come from distractors.

2. **Simulate**: Generates random fragments from the reference using weighted sampling. Fragment lengths follow a normal distribution (default: mean 175bp, range 150-200bp).

3. **Capture**: Aligns simulated reads against probe sequences using minimap2 or BLAST. Filters by matching bases, mismatches, and indels to simulate hybridization stringency.

4. **Filter** (optional): Removes reads that map to a host genome.

5. **Map**: Aligns captured reads back to the combined reference to identify which genomes they originated from.

6. **Metrics**: Computes genome-level metrics (TP/FP/FN/TN via set operations) and read-level metrics (how many captured reads came from targets vs distractors, and whether mapped reads return to their correct source genome).

7. **Report**: Generates an HTML report with ggplot2 figures (capture summary, metrics bar chart, confusion matrix, per-reference lollipop chart).

## Dependencies

- [minimap2](https://github.com/lh3/minimap2) (alignment)
- [BLAST+](https://blast.ncbi.nlm.nih.gov/) (alternative capture method)
- [R](https://www.r-project.org/) with ggplot2, rmarkdown, dplyr, tidyr (report generation, optional)

All installable via `conda env create -f environment.yml`.

## License

MIT License - see [LICENSE](LICENSE) for details.
