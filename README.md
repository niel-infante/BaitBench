# BaitBench

A generic tool for testing probe capture efficiency via in-silico simulation.

## Overview

BaitBench allows you to evaluate how well a probe set captures target sequences while avoiding off-target (distractor) sequences. It simulates the capture process and reports metrics including sensitivity, specificity, precision, and F1 score.

## Quick Start

### Installation

```bash
# Clone the repository
git clone https://github.com/yourusername/BaitBench.git
cd BaitBench

# Create conda environment
conda env create -f environment.yml
conda activate baitbench
```

### Basic Usage

```bash
nextflow run main.nf \
  --targets data/targets.fa \
  --distractors data/background.fa \
  --probes data/probes.fa \
  --num_reads 100000 \
  --outdir results
```

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
| `--num_reads` | 10000 | Number of reads to simulate |
| `--distractor_fraction` | 0.9 | Fraction of reads from distractors (0-1) |
| `--capture_method` | minimap2 | Capture method: 'minimap2' or 'blast' |
| `--min_match_bases` | 60 | Minimum matching bases for capture |
| `--max_mismatches` | 10 | Maximum mismatches allowed (minimap2 only) |
| `--host_fasta` | null | Optional host genome for filtering |
| `--seed` | null | Random seed for reproducibility |
| `--outdir` | ./results | Output directory |

## Output Files

```
results/
├── reads.fa              # Simulated reads
├── captured.fa           # Reads passing capture filter
├── filtered.fa           # After host filtering (if enabled)
├── mapped.sam            # Alignments to references
├── detected.list         # Reference IDs and read counts
├── results.tsv           # Summary metrics
├── detected_detail.tsv   # Per-reference breakdown
└── report.html           # HTML summary report
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

## Execution Profiles

```bash
# Local execution (default)
nextflow run main.nf -profile standard ...

# SLURM cluster
nextflow run main.nf -profile slurm ...

# Docker container
nextflow run main.nf -profile docker ...

# Singularity container
nextflow run main.nf -profile singularity ...
```

## Example

Run with the included minimal example:

```bash
nextflow run main.nf \
  --targets examples/minimal/targets.fa \
  --distractors examples/minimal/distractors.fa \
  --probes examples/minimal/probes.fa \
  --num_reads 1000 \
  --outdir example_results
```

## License

MIT License - see [LICENSE](LICENSE) for details.
