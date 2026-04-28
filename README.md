# BaitBench

In-silico simulation of probe capture for evaluating probe set performance.

## Download

[**Download the latest GUI installer**](https://github.com/niel-infante/BaitBench/releases/latest) — available for macOS (Apple Silicon), macOS (Intel), and Windows.

## What It Does

BaitBench tests how well a probe set captures intended target sequences while rejecting off-target (distractor) sequences. Given a set of probes, targets, and distractors, it simulates the full capture-sequencing workflow and reports detection metrics including sensitivity, specificity, precision, and F1 score.

Key capabilities:

- **Sample discrimination** -- test whether probes can distinguish between organisms within the target panel by specifying a subset of targets as "present" in the sample
- **Group-level metrics** -- collapse multiple sequence variants of the same organism into a single entity; all contigs from each distractor FASTA are automatically grouped by file name
- **Genome mode** -- model bacteria and other large pathogens where the probe target region is a small part of the full genome
- **CT-based simulation** -- set target abundance using qPCR CT scores instead of abstract fractions
- **Thermodynamic simulation** -- probe binding sites scored by nearest-neighbor free energy (SantaLucia 1998); fragments biased toward high-affinity sites via multinomial sampling
- **Coverage curves** -- sweep parameters (CT, capture fraction, sequencing depth) and visualize how coverage changes across conditions
- **Probe QC** -- evaluate probe tiling coverage across targets independently of the simulation
- **Species identification** -- call species present/absent from multi-target detection patterns, accounting for cross-reactive targets
- **Panel QC** -- assess whether a target panel can discriminate between species before running simulations
- **Probe building** -- construct a probe set from target sequences (collapse, tile, GC filter, complexity filter, deduplicate) with automatic quality assessment
- **Probe assessment** -- combined probe coverage + cross-reactivity analysis in a single report

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

### 3. Desktop GUI (optional)

A Tauri v2 desktop GUI is available in `gui/`. It wraps the `baitbench` binary as a sidecar and provides a point-and-click interface for all major tools, real-time log streaming, and in-app report viewing.

**Additional prerequisites:**
- [Node.js](https://nodejs.org/) v18 or later

**First-time setup:**

```bash
cd gui
npm install       # install frontend dependencies
make copy-sidecar # build the CLI and copy it into the GUI package
```

**Launch in development mode (hot-reload):**

```bash
make dev
# or equivalently:
npm run tauri:dev
```

**Build a distributable app bundle:**

```bash
make build
# output: gui/src-tauri/target/release/bundle/
```

**On first launch** the app shows a setup screen to select and validate your conda environment. Once saved it goes straight to the tool picker on subsequent launches.

---

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

### Build probes from targets

Construct a filtered, deduplicated probe set from target sequences (automatically runs probe assessment):

```bash
# Tiling method (default)
baitbench build-probes \
  --targets targets.fa \
  --probe-length 120 \
  --min-gc 0.20 \
  --max-gc 0.80 \
  --outdir probes_output

# catch-lite method (native Rust reimplementation of CATCH)
baitbench build-probes \
  --targets targets.fa \
  --method catch-lite \
  --probe-length 120 \
  --outdir probes_output

# catch method (external CATCH tool; requires catch conda package)
baitbench build-probes \
  --targets targets.fa \
  --method catch \
  --probe-length 120 \
  --outdir probes_output
```

### Assess existing probes

Run combined probe coverage + cross-reactivity analysis:

```bash
baitbench assess-probes \
  --targets targets.fa \
  --probes probes.fa \
  --genomes human_genome.fa \
  --outdir assess_results
```

## Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--targets` | required | Target sequences FASTA |
| `--distractors` | required | Distractor sequences FASTA (repeatable) |
| `--probes` | required | Probe sequences FASTA |
| `--sample` | all targets | Subset of targets present in specimen (TSV file or inline IDs) |
| `--groups` | -- | TSV mapping target sequence IDs to group names for group-level metrics |
| `--distractor-groups` | -- | TSV overriding automatic file-stem grouping of distractor sequences |
| `--num-fragments` | 10000 | Number of fragments to simulate |
| `--ct` | -- | CT score (alternative to `--distractor-fraction`) |
| `--distractor-fraction` | 0.9 | Fraction of fragments from distractors |
| `--simulate-mode` | thermodynamic | `thermodynamic` (TNN Boltzmann weighting) or `simple` (uniform probe-site weights) |
| `--hybridization-temperature` | 70.0 | Hybridization temperature in °C (thermodynamic mode only) |
| `--capture-fraction` | 0.5 | Fraction of fragments from probe binding sites (vs. background) |
| `--seed` | random | Random seed for reproducibility |
| `--outdir` | ./results | Output directory |
| `--output-prefix` | (empty) | String prepended to every auto-generated output filename |
| `--cleanup` | false | Delete intermediate files after completion, keeping only report inputs |

Run `baitbench run --help` for the full list.

## Output

Results are written to `<outdir>/<run_name>/` and include:

- `results.tsv` -- summary metrics (sensitivity, specificity, precision, F1)
- `detected_detail.tsv` -- per-reference detection status and coverage (with `group` column)
- `group_detail.tsv` -- per-group summary when group files are present
- `distractor_groups.tsv` -- distractor group assignments (auto-generated from FASTA file stems)
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
- [cd-hit](https://github.com/weizhongli/cdhit) -- sequence clustering (used by build-probes)
- [CATCH](https://github.com/broadinstitute/catch) -- optimization-based probe design (optional; required for `build-probes --method catch`)
- [R](https://www.r-project.org/) with ggplot2, rmarkdown -- report generation (optional)

All installable via `conda env create -f environment.yml`.

## License

MIT License -- see [LICENSE](LICENSE) for details.
