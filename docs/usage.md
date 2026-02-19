# BaitBench Usage Guide

## Installation

### Using Conda (Recommended)

```bash
# Clone the repository
git clone https://github.com/yourusername/BaitBench.git
cd BaitBench

# Create and activate the conda environment
conda env create -f environment.yml
conda activate baitbench
```

### Manual Installation

Ensure you have the following tools installed:
- Nextflow (>=21.04)
- Python (>=3.9)
- minimap2
- seqtk
- BLAST+ (optional, for BLAST capture method)
- Jinja2 Python package

## Quick Start

```bash
nextflow run main.nf \
  --targets your_targets.fa \
  --distractors your_distractors.fa \
  --probes your_probes.fa \
  --num_reads 100000 \
  --outdir results
```

## Input Files

### targets.fa

FASTA file containing sequences your probes are designed to capture. These are your positive controls - you expect these to be detected.

```fasta
>virus_A Complete genome
ATGCGTACGT...
>virus_B Partial sequence
GCTAGCTAG...
```

### distractors.fa

FASTA file containing background sequences that should NOT be captured. These serve as negative controls to test specificity.

Common choices for distractors:
- Host genome sequences
- Related but non-target organisms
- Environmental/metagenomic background

### probes.fa

FASTA file containing your probe sequences to test.

```fasta
>probe_001 Target: virus_A position 100-180
ATGCGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACG
>probe_002 Target: virus_A position 500-580
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA
```

## Parameters Reference

### Required Parameters

| Parameter | Description |
|-----------|-------------|
| `--targets` | Path to target genomes FASTA |
| `--distractors` | Path to distractor genomes FASTA |
| `--probes` | Path to probe sequences FASTA |

### Simulation Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--num_reads` | 10000 | Number of reads to simulate |
| `--distractor_fraction` | 0.9 | Fraction of reads from distractors |
| `--seed` | null | Random seed for reproducibility |

### Capture Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--capture_method` | minimap2 | Capture method: 'minimap2' or 'blast' |
| `--max_mismatches` | 10 | Maximum mismatches (minimap2 only) |
| `--min_match_bases` | 60 | Minimum matching bases required |

### Host Filtering (Optional)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--host_fasta` | null | Host genome for filtering captured reads |
| `--host_minimap_preset` | sr | Minimap2 preset for host mapping |

### Output

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--outdir` | ./results | Output directory |
| `--run_name` | auto | Run name (auto-generated if not specified) |

## Using Configuration Files

Instead of command-line parameters, you can use a YAML config file:

```yaml
# config.yaml
targets: "data/targets.fa"
distractors: "data/background.fa"
probes: "data/probes.fa"

num_reads: 100000
distractor_fraction: 0.9
seed: 42

capture_method: "minimap2"
max_mismatches: 10
min_match_bases: 60

outdir: "./results"
```

Run with:
```bash
nextflow run main.nf -params-file config.yaml
```

## Output Files

```
results/{run_name}/
├── reads.fa              # Simulated reads
├── captured.fa           # Reads passing capture filter
├── filtered.fa           # After host filtering (if enabled)
├── mapped.sam            # Alignments to references
├── detected.list         # Reference IDs and read counts
├── results.tsv           # Summary metrics
├── results.json          # Metrics in JSON format
├── detected_detail.tsv   # Per-reference breakdown
├── report.html           # HTML summary report
├── weights.txt           # Generated weights file
├── targets.txt           # List of target IDs
└── distractors.txt       # List of distractor IDs
```

## results.tsv Column Reference

The `results.tsv` file contains one row with all summary metrics. Here is a description of each column:

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
| `target_captured` | Number of captured reads that originated from **target** sequences. These are reads we *want* to capture. |
| `distractor_captured` | Number of captured reads that originated from **distractor** sequences. These are reads we do *not* want to capture (off-target capture). |

### Reference Counts

| Column | Description |
|--------|-------------|
| `targets_total` | Total number of distinct target reference sequences in the input |
| `distractors_total` | Total number of distinct distractor reference sequences in the input |

### Detection Classification (Genome-Level)

These metrics count **distinct genomes/references**, not reads. A genome is "detected" if at least one read maps back to it after capture.

| Column | Description |
|--------|-------------|
| `tp_count` | **True Positives** - Number of target genomes that were detected (correctly captured) |
| `fp_count` | **False Positives** - Number of distractor genomes that were detected (incorrectly captured) |
| `fn_count` | **False Negatives** - Number of target genomes that were NOT detected (missed) |
| `tn_count` | **True Negatives** - Number of distractor genomes that were NOT detected (correctly rejected) |

### Performance Metrics

| Column | Description |
|--------|-------------|
| `sensitivity` | True Positive Rate = `TP / (TP + FN)`. What fraction of targets were detected? Range: 0-1. |
| `specificity` | True Negative Rate = `TN / (TN + FP)`. What fraction of distractors were correctly rejected? Range: 0-1. |
| `precision` | Positive Predictive Value = `TP / (TP + FP)`. Of detected genomes, what fraction were targets? Range: 0-1. |
| `f1_score` | Harmonic mean of precision and sensitivity = `2 * (precision * sensitivity) / (precision + sensitivity)`. Range: 0-1. |

### Key Distinctions

**Read-level vs Genome-level metrics:**
- `target_captured` and `distractor_captured` count **individual reads** by their source
- `tp_count`, `fp_count`, etc. count **distinct genomes/references**

**Example:** If you have 2 target genomes and 1000 reads are captured from them:
- `target_captured` = 1000 (reads)
- `tp_count` = 2 (genomes), assuming both were detected

## Understanding the Metrics

### Sensitivity (True Positive Rate)

```
Sensitivity = TP / (TP + FN)
```

Measures what fraction of target genomes were successfully detected. High sensitivity means the probes are capturing most of their intended targets.

### Specificity (True Negative Rate)

```
Specificity = TN / (TN + FP)
```

Measures what fraction of distractor genomes were correctly NOT captured. High specificity means the probes are avoiding off-target capture.

### Precision

```
Precision = TP / (TP + FP)
```

Of all detected genomes, what fraction were actual targets? High precision means most detections are real targets, not false positives.

### F1 Score

```
F1 = 2 * (Precision * Sensitivity) / (Precision + Sensitivity)
```

Harmonic mean of precision and sensitivity. A balanced measure of overall performance.

## Using BLAST for Capture

To use BLAST instead of minimap2 for capture simulation:

1. Build BLAST database for your probes:
```bash
mkdir blast_db
makeblastdb -in probes.fa -dbtype nucl -out blast_db/probes
```

2. Run with BLAST capture:
```bash
nextflow run main.nf \
  --targets targets.fa \
  --distractors distractors.fa \
  --probes probes.fa \
  --capture_method blast \
  --blast_db_dir ./blast_db \
  --outdir results
```

## Execution Profiles

### Local (Default)
```bash
nextflow run main.nf -profile standard ...
```

### SLURM Cluster
```bash
nextflow run main.nf -profile slurm ...
```

### Docker
```bash
nextflow run main.nf -profile docker ...
```

### Singularity
```bash
nextflow run main.nf -profile singularity ...
```

## Troubleshooting

### "No sequences with positive weights found"

This error means the sequence IDs in your FASTA files don't match the weights file. Check that:
- FASTA headers are formatted correctly (ID is first word after `>`)
- There are no extra spaces or special characters

### Low capture rate

If very few reads are being captured:
- Try reducing `--min_match_bases`
- Try increasing `--max_mismatches`
- Verify your probes match your targets

### High false positive rate

If many distractors are being captured:
- Try increasing `--min_match_bases`
- Try decreasing `--max_mismatches`
- Consider adding host filtering with `--host_fasta`

## Example Workflow

```bash
# 1. Run with default parameters
nextflow run main.nf \
  --targets targets.fa \
  --distractors background.fa \
  --probes probes.fa \
  --num_reads 100000 \
  --outdir results/default

# 2. Try more stringent capture
nextflow run main.nf \
  --targets targets.fa \
  --distractors background.fa \
  --probes probes.fa \
  --num_reads 100000 \
  --min_match_bases 70 \
  --max_mismatches 5 \
  --outdir results/stringent

# 3. Compare results
cat results/default/results.tsv
cat results/stringent/results.tsv
```
