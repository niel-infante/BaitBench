# BaitBench

In-silico simulation of probe capture for evaluating probe set performance.

## What It Does

BaitBench tests how well a probe set captures intended target sequences while rejecting off-target (distractor) sequences. Given a set of probes, targets, and distractors, it simulates the full capture-sequencing workflow and reports detection metrics including sensitivity, specificity, precision, and F1 score.

Key capabilities:

- **Sample discrimination** -- test whether probes can distinguish between organisms within the target panel by specifying a subset of targets as "present" in the sample
- **Genome mode** -- model bacteria and other large pathogens where the probe target region is a small part of the full genome
- **CT-based simulation** -- set target abundance using qPCR CT scores instead of abstract fractions
- **Coverage curves** -- sweep parameters (CT, fold enrichment, sequencing depth) and visualize how coverage changes across conditions
- **Probe QC** -- evaluate probe tiling coverage across targets independently of the simulation
- **Species identification** -- call species present/absent from multi-target detection patterns, accounting for cross-reactive targets
- **Panel QC** -- assess whether a target panel can discriminate between species before running simulations

## Installation

### 1. Install dependencies

```bash
conda env create -f environment.yml
conda activate baitbench
```

### 2. Build

Requires the [Rust toolchain](https://rustup.rs/).

```bash
cargo build --release
```

The binary is at `target/release/baitbench`.

## Quick Start

### Basic run

```bash
baitbench run \
  --targets targets.fa \
  --distractors distractors.fa \
  --probes probes.fa \
  --num-fragments 10000 \
  --seed 42 \
  --outdir results
```

### With a sample subset

Test detection of specific targets within a larger panel:

```bash
baitbench run \
  --targets all_viruses.fa \
  --distractors bacteria.fa \
  --probes probes.fa \
  --sample dengue_1 zika_virus chikungunya \
  --num-fragments 10000 \
  --outdir results
```

### Using CT scores

Simulate a clinical specimen at CT 25 (~0.03% target DNA):

```bash
baitbench run \
  --targets targets.fa \
  --distractors distractors.fa \
  --probes probes.fa \
  --ct 25 \
  --num-fragments 10000 \
  --outdir results
```

### Genome mode (bacteria + viruses)

When genomes are much larger than probe target regions:

```bash
baitbench run \
  --targets gene_targets.fa \
  --genomes full_genomes.fa \
  --distractors human.fa \
  --probes probes.fa \
  --sample-target-map mapping.tsv \
  --sample e_coli influenza_a \
  --num-fragments 50000 \
  --outdir results
```

### Coverage curves

Sweep CT values and visualize coverage depth:

```bash
baitbench coverage-curve \
  --targets targets.fa \
  --distractors distractors.fa \
  --probes probes.fa \
  --sample dengue_1 zika_virus \
  --ct-values 20 25 30 35 \
  --num-fragments 10000 \
  --outdir coverage_results
```

### Probe coverage QC

Check how well probes tile across targets:

```bash
baitbench probe-coverage \
  --targets targets.fa \
  --probes probes.fa \
  --outdir probe_qc
```

### Cross-reactivity analysis

Check which probes have high homology to off-target genomes:

```bash
baitbench xreact \
  --probes probes.fa \
  --against human_genome.fa other_genomes.fa \
  --threshold 80 \
  --outdir xreact_results
```

Check for probes that are too similar to each other:

```bash
baitbench xreact \
  --probes probes.fa \
  --self \
  --threshold 80 \
  --outdir xreact_self
```

### Target panel discriminability QC

Check if targets can distinguish species before running simulations:

```bash
baitbench panel-qc \
  --targets gene_targets.fa \
  --sample-target-map mapping.tsv \
  --outdir panel_qc_results
```

### Species identification

Run the pipeline with species-level calling (genome mode only):

```bash
baitbench run \
  --targets gene_targets.fa \
  --genomes full_genomes.fa \
  --distractors human.fa \
  --probes probes.fa \
  --sample-target-map mapping.tsv \
  --sample e_coli influenza_a \
  --identify \
  --num-fragments 50000 \
  --outdir results
```

Or call species from existing pipeline results:

```bash
baitbench identify \
  --detected-detail results/run/detected_detail.tsv \
  --sample-target-map mapping.tsv \
  --targets gene_targets.fa \
  --outdir identify_results
```

## Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--targets` | required | Target sequences FASTA |
| `--distractors` | required | Distractor sequences FASTA (repeatable) |
| `--probes` | required | Probe sequences FASTA |
| `--sample` | all targets | Subset of targets present in specimen (TSV file or inline IDs) |
| `--num-fragments` | 10000 | Number of fragments to simulate |
| `--ct` | -- | CT score (alternative to `--distractor-fraction`) |
| `--distractor-fraction` | 0.9 | Fraction of fragments from distractors |
| `--fold-enrichment` | -- | Post-capture enrichment factor |
| `--seed` | random | Random seed for reproducibility |
| `--outdir` | ./results | Output directory |
| `--output-prefix` | (empty) | String prepended to every auto-generated output filename |
| `--cleanup` | false | Delete intermediate files after completion, keeping only report inputs |

Run `baitbench run --help` for the full list.

## Output

Results are written to `<outdir>/<run_name>/` and include:

- `results.tsv` -- summary metrics (sensitivity, specificity, precision, F1)
- `detected_detail.tsv` -- per-reference detection status and coverage
- `results.json` -- machine-readable metrics
- `report.html` -- HTML report with figures (requires R)

## Documentation

See [MANUAL.md](MANUAL.md) for complete documentation including:

- Detailed explanations of every parameter
- CT score calculation and calibration
- Pipeline flowcharts
- Genome mode walkthrough
- All output file formats
- Coverage curve analysis guide

## Dependencies

- [minimap2](https://github.com/lh3/minimap2) -- alignment
- [BLAST+](https://blast.ncbi.nlm.nih.gov/) -- alternative capture method
- [R](https://www.r-project.org/) with ggplot2, rmarkdown -- report generation (optional)

All installable via `conda env create -f environment.yml`.

## License

MIT License -- see [LICENSE](LICENSE) for details.
