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

## Metrics

| Metric | Formula | Interpretation |
|--------|---------|----------------|
| **Sensitivity** | TP / (TP + FN) | What fraction of targets were detected? |
| **Specificity** | TN / (TN + FP) | What fraction of distractors were correctly rejected? |
| **Precision** | TP / (TP + FP) | Of detected sequences, what fraction were targets? |
| **F1 Score** | 2 * (P * S) / (P + S) | Harmonic mean of precision and sensitivity |

Where:
- **TP (True Positive)**: Target genome detected
- **FP (False Positive)**: Distractor genome detected
- **FN (False Negative)**: Target genome NOT detected
- **TN (True Negative)**: Distractor genome NOT detected

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

6. **Metrics**: Compares detected genomes against known targets and distractors using set operations to compute TP/FP/FN/TN and derived rates.

7. **Report**: Generates an HTML report with ggplot2 figures (capture summary, metrics bar chart, confusion matrix, per-reference lollipop chart).

## Dependencies

- [minimap2](https://github.com/lh3/minimap2) (alignment)
- [BLAST+](https://blast.ncbi.nlm.nih.gov/) (alternative capture method)
- [R](https://www.r-project.org/) with ggplot2, rmarkdown, dplyr, tidyr (report generation, optional)

All installable via `conda env create -f environment.yml`.

## License

MIT License - see [LICENSE](LICENSE) for details.
