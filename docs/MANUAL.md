# BaitBench Manual

Complete reference for BaitBench, an in-silico probe capture simulation tool.

## Table of Contents

- [Overview](#overview)
- [Installation](#installation)
- [Concepts](#concepts)
  - [Pipeline Overview](#pipeline-overview)
  - [Standard Mode vs Genome Mode](#standard-mode-vs-genome-mode)
  - [Sample Manifest](#sample-manifest)
  - [3-Way Classification](#3-way-classification)
  - [CT Scores](#ct-scores)
  - [Capture Fraction and Thermodynamic Simulation](#capture-fraction-and-thermodynamic-simulation)
  - [Weight Calculation](#weight-calculation)
  - [Sequence ID Conventions](#sequence-id-conventions)
- [Pipeline Flowcharts](#pipeline-flowcharts)
  - [Standard Mode Flowchart](#standard-mode-flowchart)
  - [Genome Mode Flowchart](#genome-mode-flowchart)
- [Commands](#commands)
  - [run](#run)
  - [prepare](#prepare)
  - [simulate](#simulate)
  - [sequence](#sequence)
  - [filter](#filter)
  - [map](#map)
  - [list](#list)
  - [metrics](#metrics)
  - [report](#report)
  - [probe-coverage](#probe-coverage)
  - [xreact](#xreact)
  - [coverage-curve](#coverage-curve)
  - [panel-qc](#panel-qc)
  - [identify](#identify)
  - [build-probes](#build-probes)
  - [tool](#tool)
  - [assess-probes](#assess-probes)
- [Parameter Reference](#parameter-reference)
  - [Input Files](#input-files)
  - [Fragment Generation](#fragment-generation)
  - [Target Abundance](#target-abundance)
  - [CT Score Parameters](#ct-score-parameters)
  - [Simulation Parameters](#simulation-parameters)
  - [Sequencing Parameters](#sequencing-parameters)
  - [Execution Parameters](#execution-parameters)
- [CT Score Calculation](#ct-score-calculation)
  - [The Formula](#the-formula)
  - [Default Calibration](#default-calibration)
  - [Custom Calibration](#custom-calibration)
  - [CT Reference Table](#ct-reference-table)
- [Output Files](#output-files)
  - [Run Output Directory](#run-output-directory)
  - [results.tsv Columns](#resultstsv-columns)
  - [detected_detail.tsv Columns](#detected_detailtsv-columns)
  - [group_detail.tsv Columns](#group_detailtsv-columns)
  - [results.json Structure](#resultsjson-structure)
  - [coverage.tsv Format](#coveragetsv-format)
- [Usage Examples](#usage-examples)
  - [Basic Probe Evaluation](#basic-probe-evaluation)
  - [Sample Discrimination Testing](#sample-discrimination-testing)
  - [Clinical Specimen Simulation with CT](#clinical-specimen-simulation-with-ct)
  - [Genome Mode for Bacteria](#genome-mode-for-bacteria)
  - [Mixed Panels (Virus + Bacteria)](#mixed-panels-virus--bacteria)
  - [Multiple Distractor Sources](#multiple-distractor-sources)
  - [Group-Level Grouping](#group-level-grouping)
  - [Capture Fraction Sweep](#capture-fraction-sweep)
  - [Sequencing Depth Control](#sequencing-depth-control)
  - [Host Filtering](#host-filtering)
  - [Coverage Curve Analysis](#coverage-curve-analysis)
  - [Probe Design QC](#probe-design-qc)
  - [Cross-Reactivity Analysis](#cross-reactivity-analysis)
  - [Target Panel QC](#target-panel-qc)
  - [Species Identification](#species-identification)
  - [Probe Assessment](#probe-assessment)
  - [Running Individual Steps](#running-individual-steps)
  - [Reproducible Runs](#reproducible-runs)
  - [Batch Comparisons](#batch-comparisons)
- [Report Guide](#report-guide)
  - [HTML Report Sections](#html-report-sections)
  - [Coverage Curve Report](#coverage-curve-report)
  - [Probe Coverage Report](#probe-coverage-report)
  - [Panel QC Report](#panel-qc-report)
  - [Species Identification in Main Report](#species-identification-in-main-report)
  - [Probe Assessment Report](#probe-assessment-report)
- [Metrics Definitions](#metrics-definitions)
  - [Genome-Level Metrics](#genome-level-metrics)
  - [Read-Level Metrics](#read-level-metrics)
- [Input File Formats](#input-file-formats)
  - [FASTA Files](#fasta-files)
  - [Sample Manifest](#sample-manifest-format)
  - [Sample-Target Map](#sample-target-map-format)
  - [Groups File Format](#groups-file-format)
- [Dependencies](#dependencies)

---

## Overview

BaitBench simulates a probe capture and sequencing workflow to evaluate how well a probe set performs. It answers questions like:

- Does the probe set capture all target sequences?
- Does it reject background (distractor) sequences?
- Can it discriminate between organisms within the target panel?
- How does performance change at different target abundances (CT values)?
- What sequencing depth is needed for adequate genome coverage?

The tool aligns probes to reference sequences, scores each binding site using thermodynamic nearest-neighbor free energy (SantaLucia 1998), and generates fragments biased toward high-affinity binding sites. Background fragments fill the remainder. Reads are then mapped back to references and detection metrics are computed.

## Installation

### Prerequisites

- [Rust toolchain](https://rustup.rs/) (for building)
- [Conda](https://docs.conda.io/) or [Mamba](https://mamba.readthedocs.io/) (for runtime dependencies)

### Steps

```bash
# 1. Install runtime dependencies
conda env create -f environment.yml
conda activate baitbench

# 2. Build
cargo build --release

# 3. Verify
./target/release/baitbench --help
```

The binary is at `target/release/baitbench`. Copy it to a location on your PATH or use the full path.

### Runtime Dependencies

Installed via `environment.yml`:

| Tool | Version | Purpose |
|------|---------|---------|
| minimap2 | >= 2.24 | Sequence alignment (simulate, mapping, filtering) |
| BLAST+ | >= 2.12 | Cross-reactivity analysis (xreact) |
| R | >= 4.2 | Report generation (optional) |
| r-ggplot2 | >= 3.4 | Figures |
| r-rmarkdown | >= 2.20 | HTML report rendering |
| r-dplyr | >= 1.1 | Data manipulation |
| r-tidyr | >= 1.3 | Data reshaping |
| r-scales | >= 1.2 | Axis formatting |
| r-knitr | >= 1.40 | Report rendering |
| r-optparse | >= 1.7 | R script CLI parsing |
| r-DT | >= 0.27 | Interactive tables |
| pandoc | >= 2.19 | Document conversion |

R and its packages are only required for full HTML report generation (`--report full`). Use `--report none` to skip report generation entirely, or `--report rmd` to produce an editable RMarkdown file that can be rendered later without requiring R at pipeline run time.

---

## Concepts

### Pipeline Overview

BaitBench runs a multi-step simulation pipeline:

1. **Prepare** -- Combine target and distractor sequences; generate sampling weights
2. **Simulate** -- Align probes to reference; score binding sites by TNN thermodynamics; generate probe-biased fragments + background (controlled by `--capture-fraction`)
3. **Sequence** -- Trim fragments to read length; optionally sample to model sequencing depth
4. **Filter** (optional) -- Remove reads mapping to a host genome
5. **Map** -- Align reads back to reference sequences
6. **List** -- Count reads per reference
7. **Metrics** -- Classify each reference as TP/FP/FN/TN; compute summary statistics
8. **Report** (optional) -- Generate HTML report with figures

The `baitbench run` command chains all steps automatically. Each step is also available as a standalone subcommand for custom workflows.

### Standard Mode vs Genome Mode

**Standard mode** (default): Fragments are generated from target sequences and distractors. Reads are mapped back to the same sequences. Use this for viruses and other small genomes where the probe target IS the genome.

**Genome mode** (`--genomes`): Fragments are generated from full genome sequences, but reads are mapped back to probe target subsequences. Use this for bacteria and other large pathogens where probes target specific gene regions (e.g., 16S rRNA) rather than the whole genome. A `--sample-target-map` links genome IDs to their target IDs.

### Sample Manifest

The `--sample` flag specifies which targets (or genomes) are "present" in the simulated specimen. Only sample entries generate fragments; non-sample targets become negatives that should NOT be detected.

Without `--sample`, all targets (or genomes) are treated as present with equal weight. This tests basic capture efficiency. With `--sample`, the tool tests discrimination -- can the probes detect sample targets while rejecting non-sample targets within the same panel?

See [Sample Manifest Format](#sample-manifest-format) for syntax details.

### 3-Way Classification

BaitBench classifies each reference sequence into one of three categories, then evaluates detection:

| Category | Detected | Classification | Meaning |
|----------|----------|----------------|---------|
| Sample target | Yes | **TP** | Correctly detected |
| Sample target | No | **FN** | Missed detection |
| Non-sample target | Yes | **FP_target** | Cross-reactive within panel |
| Non-sample target | No | **TN_target** | Correctly rejected |
| Distractor | Yes | **FP_distractor** | Off-target capture |
| Distractor | No | **TN_distractor** | Correctly rejected |
| Untargeted genome | -- | **untargeted** | No expected target (genome mode only) |

This distinguishes two types of false positives:
- **FP_target**: Cross-reactivity within the target panel (e.g., probe for virus A captures virus B)
- **FP_distractor**: True off-target capture (e.g., probe captures bacterial DNA)

Without `--sample`, all targets are in the sample, reducing to a 2-way classification (TP/FP/FN/TN with no FP_target).

### Group-Level Metrics

When a target panel contains multiple sequence variants of the same organism (e.g., `West_Nile_virus_0001`, `West_Nile_virus_0002`, `West_Nile_virus_0003`), it may be desirable to treat all variants as a single entity for metrics purposes. Similarly, a distractor FASTA with thousands of contigs (e.g., host genome) should count as a single FP entity rather than thousands.

BaitBench supports **group-level metrics** via two optional flags:

- `--groups <groups.tsv>` -- maps target sequence IDs to group names. Sequences not mentioned form singleton groups using their own ID. When absent, each target sequence is its own group (backward-compatible behavior).
- `--distractor-groups <distractor_groups.tsv>` -- explicit distractor grouping. When absent (default), all contigs from each `--distractors` FASTA file are automatically grouped under the file stem name (e.g., all contigs in `Aaegypti.fa` → group `"Aaegypti"`).

A group is **detected** if any member sequence has at least one read mapped to it. Classification (TP/FN/FP/TN) operates on groups rather than individual sequences. A read mapping from a sequence in group A to any other sequence in group A is counted as **correctly mapped** (not as a cross-mapping error).

Results include a `group_detail.tsv` file with one row per group, showing group name, category, detection status, member count, detected member count, and total reads. The `detected_detail.tsv` gains a `group` column showing each sequence's group assignment.

See [Groups File Format](#groups-file-format) for syntax. See [Group-Level Grouping Example](#group-level-grouping) for a complete usage example.

### CT Scores

CT (cycle threshold) scores from qPCR provide a natural way to express target abundance. BaitBench converts CT values to distractor fractions using a calibrated exponential formula. Lower CT = more target DNA = easier to detect.

See [CT Score Calculation](#ct-score-calculation) for the formula, default calibration, and how to customize it.

### Capture Fraction and Thermodynamic Simulation

`--capture-fraction` (default 0.5) controls what fraction of simulated fragments come from probe binding sites. The remaining fraction are background fragments drawn uniformly by sequence weight × length.

Probe binding sites are scored using the SantaLucia (1998) nearest-neighbor thermodynamic model: ΔG is computed from consecutive Watson-Crick stacking interactions along each probe-reference alignment, and the Boltzmann factor `exp(-ΔG / RT)` weights sampling toward high-affinity sites. Use `--simulate-mode simple` to skip TNN scoring and use uniform weights instead (no temperature required).

Target enrichment is emergent from TNN affinity × sequence weights rather than being imposed post-hoc — sequences with weight 0.0 (non-sample targets) never generate probe-biased fragments. Fold enrichment is no longer a parameter.

### Weight Calculation

Sampling weights determine how many fragments each sequence generates. The number of fragments from a sequence is proportional to `weight * sequence_length`.

**Standard mode:**
- Sample targets: weight from sample manifest (default 1.0)
- Non-sample targets: weight = 0 (no fragments)
- Distractors: calculated to achieve the requested distractor fraction

**Genome mode:**
- Sample genomes: weight from sample manifest (default 1.0)
- Non-sample genomes: weight = 0
- Distractors: same formula as standard mode

The distractor weight formula ensures the requested fraction of total fragments come from distractors:

```
distractor_weight = (distractor_fraction * total_sample_weight) / (n_distractors * (1 - distractor_fraction))
```

### Sequence ID Conventions

Sequence IDs are taken from the first whitespace-delimited word of each FASTA header (everything after `>` up to the first space). These IDs must be unique within each file and consistent across input files.

**Sequence names must not contain spaces.** Use underscores or other delimiters: `>Zika_virus` not `>Zika virus`.

Fragment names follow the pattern `{seq_id}_fragment_{n}`, using the last occurrence of `_fragment_` as the delimiter. This allows sequence IDs to contain the substring `_fragment_` without ambiguity.

---

## Pipeline Flowcharts

### Standard Mode Flowchart

```
INPUT FILES                    STEP                          OUTPUT FILES
=============                  ====                          ============

targets.fa ──────────┐
                     │
distractors.fa ──────┤
                     ├──── 1. PREPARE ──────────────────── combined_reference.fa
sample (optional) ───┤         │                            weights.txt
                     │         │                            targets.txt
--distractor-fraction│         │                            distractors.txt
  or --ct ───────────┘         │                            sample.txt
                               │
                               ▼
combined_reference.fa ──┐
probes.fa ──────────────┤
weights.txt ────────────┤ 2. SIMULATE ────────────────── fragments.fa
--num-fragments ────────┤      │      (probe-biased + background)
--capture-fraction ─────┤      │
--simulate-mode ────────┤      │
--hybridization-temp ───┤      │
--fragment-length-* ────┤      │
--seed ─────────────────┘      │
                               │
                               ▼
fragments.fa ───────────┐
                        ├─ 3. SEQUENCE ──────────────────── reads.fa
--read-length ──────────┤      │
--num-sequences ────────┤      │
--seed ─────────────────┘      │
                             │
                    ┌────────┴────────┐
                    │  --host-fasta   │
                    │   specified?    │
                    └──┬──────────┬───┘
                   yes │          │ no
                       ▼          │
reads.fa ───────┐                 │
host.fa ────────┤ 4. FILTER       │
--host-minimap- ┤    │            │
  preset ───────┘    │            │
                     ▼            │
              filtered.fa         │
                     │            │
                     ▼            ▼
              (filtered or reads).fa
                     │
combined_            ├─ 5. MAP ────────────────────────── mapped.sam
  reference.fa ──────┤      │
--minimap-preset ────┘      │
                            │
                            ▼
mapped.sam ──────────── 6. LIST ───────────────────────── detected.list
                            │
                            ▼
targets.txt ─────────┐
distractors.txt ─────┤
sample.txt ──────────┤
detected.list ───────┤ 7. METRICS ────────────────────── results.tsv
fragments.fa ────────┤                                    detected_detail.tsv
fragments.fa ────────┤                                    results.json
mapped.sam ──────────┘                                    coverage.tsv
                            │
                            ▼
results.tsv ─────────┐
detected_detail.tsv ─┤ 8. REPORT (optional) ──────────── report.html
run_params.tsv ──────┤
coverage.tsv ────────┘
```

### Genome Mode Flowchart

Genome mode adds a separate mapping reference and genome-aware classification:

```
INPUT FILES                    STEP                          OUTPUT FILES
=============                  ====                          ============

targets.fa ──────────┐
                     │
genomes.fa ──────────┤
                     │
distractors.fa ──────┤
                     ├──── 1. PREPARE ──────────────────── combined_reference.fa
sample ──────────────┤                                        (genomes + distractors)
                     │                                      mapping_reference.fa
sample-target-map ───┤                                        (targets + distractors)
                     │                                      weights.txt
--distractor-fraction│                                      targets.txt
  or --ct ───────────┘                                      distractors.txt
                                                            genomes.txt
                                                            sample.txt
                                                            sample_target_map.txt

    Steps 2-4 are identical to standard mode:
      - Simulate uses combined_reference.fa (genomes + distractors); probes align to genomes

                     ... (steps 2-4) ...

              (filtered or reads).fa
                     │
mapping_             ├─ 5. MAP ────────────────────────── mapped.sam
  reference.fa ──────┤       (targets + distractors)
                     │
                            │
                            ▼
                     ... (step 6 same) ...
                            │
                            ▼
targets.txt ─────────┐
distractors.txt ─────┤
sample.txt ──────────┤
sample_target_map ───┤ 7. METRICS ────────────────────── results.tsv
detected.list ───────┤   (genome-aware classification)    detected_detail.tsv
fragments.fa ────────┤                                    results.json
fragments.fa ────────┤                                    coverage.tsv
mapped.sam ──────────┘

                     ... (step 8 same) ...
```

Key differences in genome mode:
- **combined_reference.fa** = genomes + distractors (fragments generated from full genomes)
- **mapping_reference.fa** = targets + distractors (reads mapped to target regions)
- A read from genome G mapping to target T is correct if T is linked to G in the sample-target-map
- Untargeted genomes (no target mapping) are tracked separately and do not affect TP/FP/FN/TN

---

## Commands

### run

Runs the complete pipeline from input files to metrics and report.

```bash
baitbench run [OPTIONS]
```

This is the primary command for most use cases. It chains all pipeline steps (prepare through report) automatically. Use `--cleanup` to delete intermediate files (FASTA, SAM, logs) after completion, keeping only report inputs and final outputs. See [Parameter Reference](#parameter-reference) for all options.

In genome mode with `--sample-target-map`, use `--identify` to add species-level calling after metrics. This computes target similarity, calls species PRESENT/ABSENT/AMBIGUOUS, and includes the results in the HTML report with ground-truth comparison against the `--sample` manifest.

### prepare

Combines target and distractor sequences into a single reference, generates sampling weights, and writes ID lists.

```bash
baitbench prepare \
  --targets targets.fa \
  --distractors distractors.fa \
  [--genomes genomes.fa] \
  [--sample manifest.tsv] \
  [--sample-target-map mapping.tsv] \
  [--groups target_groups.tsv] \
  [--distractor-groups distractor_groups.tsv] \
  [--distractor-fraction 0.9 | --ct 25] \
  [--ct-baseline 20.0] \
  [--ct-baseline-fraction 0.01] \
  --outdir prep_output
```

**Output files:**
- `combined_reference.fa` -- merged sequences for fragment generation
- `weights.txt` -- per-sequence sampling weights (TSV: `id\tweight`)
- `targets.txt` -- target sequence IDs (one per line)
- `distractors.txt` -- distractor sequence IDs (one per line)
- `sample.txt` -- sample sequence IDs (one per line)
- `target_groups.tsv` -- target group assignments (written only if `--groups` provided)
- `distractor_groups.tsv` -- distractor group assignments (always written; from `--distractor-groups` or auto-generated from FASTA file stems)
- `mapping_reference.fa` -- targets + distractors for read mapping (genome mode only)
- `genomes.txt` -- genome IDs (genome mode only)
- `sample_target_map.txt` -- genome-to-target mapping (genome mode only)

### simulate

Generates weighted random fragments from a reference.

```bash
baitbench simulate \
  --reference combined_reference.fa \
  --weights weights.txt \
  --num-fragments 10000 \
  --output fragments.fa \
  [--fragment-length-mean 175] \
  [--fragment-length-min 150] \
  [--fragment-length-max 200] \
  [--seed 42]
```

Fragment lengths follow a truncated normal distribution clamped to [min, max]. Fragments are named `{seq_id}_fragment_{n} start={pos} length={len}`.

**Output files:**
- `fragments.fa` -- simulated DNA fragments

### sequence

Simulates sequencing of captured fragments. Three modes are available:

#### Perfect (default)

Trims each fragment to `--read-length` bp. No errors. One read per fragment.

```bash
baitbench sequence \
  --input fragments.fa \
  --output reads.fa \
  [--read-length 120] \
  [--num-sequences 5000] \
  [--seed 42]
```

Fragments shorter than `--read-length` are kept as-is. With `--num-sequences`, reads are sampled with replacement from the fragment pool.

#### Illumina (ART-modern)

Generates error-realistic Illumina reads using [ART-modern](https://github.com/YU-Zhejian/art_modern). Requires `art_modern` on PATH (`conda install -c bioconda art_modern`).

```bash
baitbench sequence \
  --input fragments.fa \
  --output reads.fa \
  --read-simulator art \
  --sequencer-profile HiSeq2500_150bp \
  --read-length 150 \
  [--coverage-depth 1.0] \
  [--paired-end] \
  [--pe-frag-len-mean 200 --pe-frag-len-sd 50]
```

Common profiles: `HiSeq2500_150bp` (default), `HiSeq2500_100bp`, `MiSeq_250bp`. Run `art_modern --list-profiles` for the full list. `--coverage-depth` controls how many reads are generated per fragment (total reads ≈ num_fragments × mean_fragment_len / read_length × coverage_depth). `--paired-end` produces both `reads.fa` (R1) and `reads_R2.fa` (R2).

#### Long reads (badread)

Generates ONT or PacBio CLR long reads using [badread](https://github.com/rrwick/Badread). Requires `badread` on PATH (`conda install -c conda-forge badread`).

```bash
baitbench sequence \
  --input fragments.fa \
  --output reads.fa \
  --read-simulator badread \
  --sequencer-profile ont \
  [--coverage-depth 1.0]
```

`--sequencer-profile` selects the chemistry:

| Profile | Platform / Chemistry | Notes |
|---------|---------------------|-------|
| `ont` (default) | ONT R10.4.1 / Kit14 | nanopore2023 error model |
| `ont-2020` | ONT R9.4.1 | nanopore2020 error model |
| `pacbio` | PacBio CLR | pacbio2016 error model |

`--read-length` is not used — read length is bounded by fragment length. `--coverage-depth 1` produces approximately one read per captured fragment. Paired-end is not supported for long reads.

#### Choosing a simulator

| Scenario | Recommended simulator |
|----------|-----------------------|
| Fast development / debugging | `perfect` |
| Comparing against Illumina data | `art` with matching profile |
| Comparing against ONT data | `badread --sequencer-profile ont` or `ont-2020` |
| Comparing against PacBio CLR data | `badread --sequencer-profile pacbio` |
| Paired-end Illumina panel evaluation | `art --paired-end` |

When `--num-sequences` is set, the final read count is capped by sampling, regardless of simulator. This lets you model a fixed sequencing depth even when `art`/`badread` generate more reads than needed.

When using `baitbench run`, `--minimap-preset` and `--host-minimap-preset` are auto-selected to match the simulator: `sr` for `perfect`/`art`, `map-ont` for badread `ont`/`ont-2020`, `map-pb` for badread `pacbio`. Pass either flag explicitly to override.

**Output files:**
- `reads.fa` -- sequencing reads (R1 for paired-end)
- `reads_R2.fa` -- R2 reads (paired-end only)

### filter

Removes reads that map to a host genome.

```bash
baitbench filter \
  --host host_genome.fa \
  --reads reads.fa \
  --output filtered.fa \
  [--minimap-preset sr]
```

Uses minimap2 to align reads against the host genome. Reads that map are removed; unmapped reads are kept.

**Output files:**
- `filtered.fa` -- reads after host depletion

### map

Maps reads back to a reference using minimap2.

```bash
baitbench map \
  --reference combined_reference.fa \
  --reads reads.fa \
  --output mapped.sam \
  [--minimap-preset sr]
```

In standard mode, reads are mapped to `combined_reference.fa`. In genome mode, reads are mapped to `mapping_reference.fa` (targets + distractors).

**Output files:**
- `mapped.sam` -- SAM alignment file

### list

Counts reads per reference from a SAM file.

```bash
baitbench list \
  --sam mapped.sam \
  --output detected.list
```

**Output files:**
- `detected.list` -- TSV: `reference_id\tcount` (sorted ascending by count)

### metrics

Computes classification metrics and coverage statistics.

```bash
baitbench metrics \
  --targets targets.txt \
  --distractors distractors.txt \
  --sample sample.txt \
  --detected detected.list \
  --fragments fragments.fa \
  --captured fragments.fa \
  --sam mapped.sam \
  --run-name "my_run" \
  --num-fragments 10000 \
  --output-summary results.tsv \
  --output-detail detected_detail.tsv \
  [--output-json results.json] \
  [--output-coverage coverage.tsv] \
  [--sample-target-map sample_target_map.txt] \
  [--target-groups target_groups.tsv] \
  [--distractor-groups distractor_groups.tsv] \
  [--seed 42]
```

**Output files:**
- `results.tsv` -- genome-level summary metrics (group-level if groups provided)
- `detected_detail.tsv` -- per-reference detection and coverage detail (includes `group` column)
- `group_detail.tsv` -- per-group summary (written when group files are provided)
- `results.json` -- structured JSON output (optional)
- `coverage.tsv` -- run-length encoded read depth intervals (optional)

### report

Generates an HTML report with figures, or outputs an editable RMarkdown file.

```bash
baitbench report \
  --summary results.tsv \
  --detail detected_detail.tsv \
  --params run_params.tsv \
  --output report.html \
  [--coverage coverage.tsv] \
  [--run-name "BaitBench Run"] \
  [--report full|rmd]
```

**Output files:**
- `report.html` -- HTML report with ggplot2 visualizations (`--report full`)
- `report.Rmd` -- editable RMarkdown file with parameters pre-filled (`--report rmd`)

### probe-coverage

Standalone probe design QC tool. Not part of the main simulation pipeline.

```bash
baitbench probe-coverage \
  --targets targets.fa \
  --probes probes.fa \
  [--outdir probe_coverage] \
  [--minimap-preset sr] \
  [--proximity 50] \
  [--report full|none|rmd]
```

Maps probes to targets and computes per-target tiling statistics.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--targets` | required | Target sequences FASTA |
| `--probes` | required | Probe sequences FASTA |
| `--outdir` | ./probe_coverage | Output directory |
| `--minimap-preset` | sr | Minimap2 alignment preset |
| `--proximity` | 50 | Pull-down zone distance in bp |
| `--report` | full | Report mode: `full` (HTML), `none` (skip), `rmd` (editable RMarkdown) |
| `--cleanup` | false | Delete intermediate files (SAM, logs) after completion |

**Output files:**
- `probe_depth.tsv` -- run-length encoded probe depth intervals (TSV: `reference_id\tstart\tend\tdepth`)
- `probe_coverage_summary.tsv` -- per-target coverage statistics
- `multi_mapping_probes.tsv` -- probes mapping to multiple targets
- `probe_coverage_report.html` -- HTML report (`--report full`, requires R)
- `probe_coverage_report.Rmd` -- editable RMarkdown file (`--report rmd`)

**Coverage summary columns:**

| Column | Description |
|--------|-------------|
| `reference_id` | Target sequence ID |
| `pct_covered_1x` | % bases with >= 1 probe |
| `pct_covered_2x` | % bases with >= 2 probes |
| `pct_covered_5x` | % bases with >= 5 probes |
| `pct_covered_10x` | % bases with >= 10 probes |
| `mean_depth` | Average probe depth across target |
| `median_depth` | Median probe depth |
| `max_gap_length` | Longest uncovered stretch (bp) |
| `num_gaps` | Number of gaps with no probe coverage |
| `pct_near_probe` | % bases within `--proximity` distance of a probe alignment |

### xreact

Standalone cross-reactivity analysis tool. Checks whether probes have high homology to off-target genomes or to each other. Not part of the main simulation pipeline.

```bash
baitbench xreact \
  --probes probes.fa \
  [--against genome1.fa genome2.fa ...] \
  [--self] \
  [--threshold 80.0] \
  [--minimap-preset sr] \
  [--outdir xreact_results]
```

At least one of `--against` or `--self` must be specified; both can be used together.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--probes` | required | Probe sequences FASTA |
| `--against` | none | Reference genome FASTA(s) to check cross-reactivity against (repeatable) |
| `--self` | false | Check probe-vs-probe cross-reactivity (self-hits excluded) |
| `--threshold` | 80.0 | Minimum homology % to report: `matching_bases / probe_length * 100` |
| `--minimap-preset` | sr | Minimap2 alignment preset |
| `--outdir` | ./xreact_results | Output directory |
| `--cleanup` | false | Delete intermediate files (logs) after completion |

**Homology metric:** `matching_bases / probe_length * 100`. This single number captures both alignment identity and query coverage -- a probe with 90% identity over 90% of its length scores ~81%.

**Self-mode filtering:** In `--self` mode, self-hits (probeA mapping to probeA) are excluded from all output. Only cross-probe hits (probeA mapping to probeB where A != B) are reported.

**Output files:**

- `hits.tsv` -- All alignments above the threshold
- `summary.tsv` -- Per-probe summary (every probe gets a row, even with zero hits)

**hits.tsv columns:**

| Column | Description |
|--------|-------------|
| `probe_id` | Query probe ID |
| `target_id` | Reference sequence the probe maps to (genome ID or other probe ID) |
| `homology_pct` | `matching_bases / probe_length * 100` |
| `identity_pct` | `matching_bases / alignment_block_length * 100` |
| `query_coverage_pct` | `aligned_query_span / probe_length * 100` |
| `matching_bases` | Number of matching bases in the alignment |
| `alignment_length` | Alignment block length |
| `probe_length` | Total probe sequence length |
| `mode` | `against` (probe-to-genome) or `self` (probe-to-probe) |

**summary.tsv columns:**

| Column | Description |
|--------|-------------|
| `probe_id` | Probe ID |
| `mode` | `against` or `self` |
| `max_homology_pct` | Highest homology % across all hits (0.0 if no hits) |
| `best_hit` | Target ID with highest homology (NA if no hits) |
| `num_hits_above_threshold` | Number of distinct alignments above threshold |

### coverage-curve

Runs the pipeline at multiple parameter combinations and generates coverage depth curves.

```bash
baitbench coverage-curve \
  --targets targets.fa \
  --distractors distractors.fa \
  --probes probes.fa \
  --sample dengue_1 zika_virus \
  [--ct-values 20 25 30 | --ct 25] \
  [--hybridization-temperature-values 55 65 70 75 | --hybridization-temperature 70] \
  [--capture-fraction-values 0.3 0.5 0.8 | --capture-fraction 0.5] \
  [--num-sequences-values 100 500 | --num-sequences 500] \
  [--outdir coverage_curve_results] \
  [--cleanup] \
  [... other pipeline parameters ...]
```

Four parameters can be swept (each has a singular fixed form and a plural sweep form):

| Sweep flag | Fixed flag | Default | Description |
|-----------|------------|---------|-------------|
| `--ct-values 20 25 30` | `--ct 25` | — | CT values (converted to distractor fractions) |
| `--hybridization-temperature-values 55 65 70 75` | `--hybridization-temperature 70` | 70 °C | Hybridization temperature; thermodynamic mode only |
| `--capture-fraction-values 0.3 0.5 0.8` | `--capture-fraction 0.5` | 0.5 | Capture fraction (probe-biased fragment proportion) |
| `--num-sequences-values 100 500` | `--num-sequences 500` | all | Number of sequences to sample |

Sweep and fixed forms of the same parameter are mutually exclusive. `--ct-values` and `--distractor-fraction` are also mutually exclusive.

`--sample` is **required** for coverage-curve (must specify which targets to track).

The pipeline shares intermediate files across combinations for efficiency: prepare is shared per CT value; simulate is shared per CT × temperature × capture-fraction combination.

**Output files:**
- Combo subdirectories named by swept params (e.g., `ct_20/`, `ct_20_temp_65_cf_0.50/`, `ct_20_temp_65_cf_0.50_ns_500/`)
- `coverage_curve_depth_curves.tsv` -- aggregated depth data (columns: ct, hybridization_temperature, capture_fraction, num_sequences, ...)
- `coverage_curve_report.html` -- HTML report with depth curves (`--report full`)
- `coverage_curve_report.Rmd` -- editable RMarkdown file (`--report rmd`)

### panel-qc

Standalone pre-experiment QC tool that assesses whether a target panel can discriminate between species. This evaluates target uniqueness before running simulations.

```bash
baitbench panel-qc \
  --targets gene_targets.fa \
  --sample-target-map mapping.tsv \
  [--identity-threshold 90.0] \
  [--minimap-preset sr] \
  [--outdir panel_qc_results] \
  [--report full] \
  [--cleanup]
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--targets` | required | Target sequences FASTA |
| `--sample-target-map` | required | Mapping TSV linking species/genome IDs to target IDs |
| `--identity-threshold` | 90.0 | Minimum sequence identity % to consider two targets "similar" |
| `--minimap-preset` | sr | Minimap2 alignment preset for target-vs-target comparison |
| `--outdir` | ./panel_qc_results | Output directory |
| `--report` | full | Report mode: `full` (HTML), `none` (skip), `rmd` (editable RMarkdown) |
| `--cleanup` | false | Delete intermediate files after completion |

**Algorithm:**

1. All targets are aligned against all targets using minimap2 (`--minimap-preset`)
2. Pairwise similarity is computed as `matching_bases / min(len_a, len_b) * 100`
3. Pairs above `--identity-threshold` are reported as similar
4. Using the sample-target-map, targets are classified as "unique" (no cross-species similarity) or "shared" (has similar targets in other species)
5. Per-species discriminability score is `unique_targets / total_targets`
6. A species confusion matrix shows which species pairs share similar targets

**Interpreting results:**

- A species with discriminability score 0.0 has **zero** unique targets -- it cannot be reliably distinguished from other species. Consider adding more targets.
- The confusion matrix highlights species pairs that share similar targets, indicating potential cross-reactivity in identification.
- High discriminability (close to 1.0) means most targets are unique to that species -- identification should be reliable.

**Output files:**

- `target_similarity.tsv` -- pairwise target similarities above threshold
- `species_discriminability.tsv` -- per-species discriminability scores
- `species_confusion_matrix.tsv` -- species-by-species shared target counts
- `panel_qc_report.html` -- HTML report with heatmap and discriminability charts (`--report full`)
- `panel_qc_report.Rmd` -- editable RMarkdown file (`--report rmd`)

**target_similarity.tsv columns:**

| Column | Description |
|--------|-------------|
| `target_a` | First target ID |
| `target_b` | Second target ID |
| `identity_pct` | `matching_bases / min(len_a, len_b) * 100` |
| `matching_bases` | Number of matching bases |
| `len_a` | Length of target A |
| `len_b` | Length of target B |

**species_discriminability.tsv columns:**

| Column | Description |
|--------|-------------|
| `species_id` | Species/genome ID |
| `total_targets` | Total targets assigned to this species |
| `unique_targets` | Targets with no cross-species similarity |
| `shared_targets` | Targets similar to targets in other species |
| `discriminability_score` | `unique_targets / total_targets` (0.0–1.0) |
| `confusable_species` | Comma-separated species IDs with shared targets |

### identify

Call species presence/absence from multi-target detection patterns. Can be run standalone on existing pipeline results or integrated into `baitbench run` with `--identify`.

```bash
# Using pre-computed similarity from panel-qc
baitbench identify \
  --detected-detail results/run/detected_detail.tsv \
  --sample-target-map mapping.tsv \
  --target-similarity panel_qc/target_similarity.tsv \
  [--min-unique-targets 1] \
  [--outdir identify_results]

# Computing similarity on-the-fly
baitbench identify \
  --detected-detail results/run/detected_detail.tsv \
  --sample-target-map mapping.tsv \
  --targets gene_targets.fa \
  [--identity-threshold 90.0] \
  [--minimap-preset sr] \
  [--min-unique-targets 1] \
  [--outdir identify_results]
```

Either `--target-similarity` or `--targets` must be provided (not both).

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--detected-detail` | required | `detected_detail.tsv` from metrics step |
| `--sample-target-map` | required | Mapping TSV linking species/genome IDs to target IDs |
| `--target-similarity` | none | Pre-computed similarity TSV from `panel-qc` |
| `--targets` | none | Target sequences FASTA (computes similarity on-the-fly) |
| `--identity-threshold` | 90.0 | Similarity threshold (only with `--targets`) |
| `--minimap-preset` | sr | Minimap2 preset (only with `--targets`) |
| `--min-unique-targets` | 1 | Minimum unique target detections to call a species PRESENT |
| `--outdir` | ./identify_results | Output directory |

**Calling algorithm (ordered-explanation approach):**

1. **Classify targets**: Each target is labeled "unique" (no similar targets in other species) or "shared" (has cross-species similarity), using the similarity data.

2. **Collect evidence**: For each species, count how many unique vs shared targets were detected, and how many total reads were observed.

3. **Sort by evidence strength**: Species are ranked by `(unique_detected DESC, total_reads DESC)`.

4. **Process in order**: Each species is assigned one of three calls:
   - **PRESENT**: `≥ min_unique_targets` unique targets detected
   - **ABSENT** (no detections): zero targets detected for this species
   - **ABSENT** (cross-reactivity explained): all detected targets are "shared" AND every one can be explained by a species already called PRESENT earlier in the ordered processing
   - **AMBIGUOUS** (no unique markers): species has zero unique targets in the panel -- cannot confirm or deny
   - **AMBIGUOUS** (insufficient evidence): species has some unique targets but not enough detected (< `min_unique_targets`), and not all shared detections are explained

This creates a natural parsimony effect: the species with the strongest unique evidence is called first, then its presence "explains away" shared target hits in subsequent species.

**Output files:**

- `species_calls.tsv` -- per-species call with evidence breakdown
- `species_calls.json` -- structured JSON format

**species_calls.tsv columns:**

| Column | Description |
|--------|-------------|
| `species_id` | Species/genome ID |
| `call` | PRESENT, ABSENT, or AMBIGUOUS |
| `total_targets` | Total targets for this species in the panel |
| `unique_targets` | Targets unique to this species |
| `shared_targets` | Targets shared with other species |
| `unique_detected` | Unique targets that were detected |
| `shared_detected` | Shared targets that were detected |
| `total_detected` | Total targets detected |
| `total_reads` | Total reads across all detected targets |
| `explained_by` | Comma-separated species IDs that explain shared hits |
| `reason` | Call reason: `unique_markers_detected`, `no_detections`, `cross_reactivity_explained`, `no_unique_markers`, `insufficient_unique_evidence` |

**Integration with `baitbench run`:**

When `--identify` is passed to `baitbench run` (genome mode with `--sample-target-map` required), species identification runs automatically after the metrics step. The species calls are included in the HTML report and compared against ground truth (the `--sample` manifest) to compute species-level sensitivity and specificity.

### build-probes

Build a probe set from target sequences. Runs a multi-step pipeline: collapse redundant targets, construct probes, filter by GC content and sequence complexity, and deduplicate. After building, automatically chains into probe assessment (probe coverage + cross-reactivity analysis) unless `--skip-assess` is specified.

Five probe construction methods are available: `tile` (sliding window, default), `catch-lite` (native Rust reimplementation of CATCH optimization-based design), `catch` (external CATCH tool from the Broad Institute; requires the `catch` conda package), `syotti-lite` (native Rust reimplementation of Syotti greedy set-cover design), and `probetools-lite` (native Rust reimplementation of ProbeTools iterative k-mer clustering design; requires cd-hit-est).

```bash
baitbench build-probes \
  --targets targets.fa \
  [--method tile|catch-lite|syotti-lite|catch|probetools-lite] \
  [--probe-length 120] \
  [--step -60] \
  [--catch-probe-stride 60] \
  [--catch-mismatches 5] \
  [--catch-extension 0] \
  [--catch-coverage 1.0] \
  [--catch-minhash-threshold 0.6] \
  [--syotti-mismatches 40] \
  [--syotti-seed-len 20] \
  [--pt-step 1] \
  [--pt-identity 0.9] \
  [--pt-coverage 0.9] \
  [--pt-batch-size 100] \
  [--pt-max-panel-size N] \
  [--pt-min-depth 1] \
  [--pt-max-iterations 20] \
  [--pt-min-coverage-gain 0.001] \
  [--min-gc 0.20] \
  [--max-gc 0.80] \
  [--max-n-frac 0.05] \
  [--dust-threshold 2.0] \
  [--dust-window 64] \
  [--max-masked-frac 0.25] \
  [--collapse-threshold 0.95] \
  [--dedup-threshold 0.95] \
  [--threads 5] \
  [--genomes genome1.fa genome2.fa ...] \
  [--threshold 80.0] \
  [--minimap-preset sr] \
  [--proximity 50] \
  [--skip-assess] \
  [--outdir build_probes_results] \
  [--report full|none|rmd] \
  [--refine-iterations N | --refine-until-stable] \
  [--refine-threshold 80.0]
```

**Pipeline steps:**

1. **N filter**: Remove target sequences with more than `--max-n-frac` fraction of ambiguous (non-ACGT) bases. Sequences with excessive N content are poor probe sources and would generate uninformative probes.
2. **Collapse**: cd-hit-est clusters targets at `--collapse-threshold` identity to remove near-duplicates
3. **Build**: Construct probes from collapsed sequences. Method `tile` generates sliding-window probes of `--probe-length` bp across each sequence with `--step` controlling overlap/gap. A final probe is anchored to the end of each sequence to ensure full coverage. Method `catch-lite` uses BaitBench's native Rust reimplementation of the CATCH algorithm. Method `catch` calls the external CATCH tool (`design_probes.py`; requires `catch` conda package). Method `syotti-lite` uses BaitBench's native Rust reimplementation of the Syotti greedy set-cover algorithm. Method `probetools-lite` uses BaitBench's native Rust reimplementation of the ProbeTools iterative k-mer clustering algorithm (requires cd-hit-est).
4. **GC filter**: Remove probes with GC content outside `--min-gc` to `--max-gc` range
5. **Complexity filter**: Remove low-complexity probes using the sDUST algorithm (Morgulis et al. 2006). Probes where more than `--max-masked-frac` of bases are identified as low-complexity (e.g., homopolymers, dinucleotide repeats) are removed. Set `--max-masked-frac 1.0` to disable.
6. **Deduplicate**: cd-hit-est clusters probes at `--dedup-threshold` identity to remove redundant probes

**Tiling geometry (`--step`):**

The stride between consecutive probes is `probe_length + step`. The step is measured from the end of the previous probe:

- `--step -60` (default): stride = 60, probes overlap by 60bp (50% overlap with 120bp probes)
- `--step 0`: stride = 120, probes are perfectly tiled (no overlap, no gap)
- `--step 10`: stride = 130, 10bp gap between probes

Probes are named `probe_{target_id}|tile_{n}`. A final probe is always placed at the sequence end regardless of overlap.

**catch-lite method (`--method catch-lite`):**

BaitBench includes a native Rust reimplementation of the CATCH algorithm (Metsky et al. 2019, Nature Biotechnology). Unlike tiling, CATCH minimizes the number of probes needed while guaranteeing a configurable fraction of each target sequence is covered. The algorithm tiles candidate probes at a configurable stride, removes near-duplicates via MinHash LSH, then runs a greedy set-cover to select the minimum probe set that covers all targets to the required depth.

Parameters:

| Flag | Default | Description |
|------|---------|-------------|
| `--catch-probe-stride` | 60 | Step between candidate probes (bp) |
| `--catch-mismatches` | 5 | Mismatches tolerated for a probe to cover a target window |
| `--catch-extension` | 0 | Flanking bp beyond probe boundaries counted as covered |
| `--catch-coverage` | 1.0 | Fraction of each target that must be covered (0.0–1.0) |
| `--catch-minhash-threshold` | 0.6 | Jaccard similarity threshold for near-deduplication; set to 0.0 to disable |

Probes are named `probe_{source_id}|catch_{n}`.

Example with custom parameters:

```bash
baitbench build-probes \
  --targets targets.fa \
  --method catch-lite \
  --probe-length 120 \
  --catch-probe-stride 30 \
  --catch-mismatches 3 \
  --catch-extension 10 \
  --catch-coverage 0.95 \
  --outdir probes_output
```

#### Differences from external CATCH

| Aspect | External CATCH | catch-lite (native) |
|--------|---------------|---------------------|
| Coverage model | LCF-k (full algorithm) | Hamming distance (= LCF-k at default threshold); LCF-k planned |
| Greedy tie-breaking | Python dict order | Rust iteration order (different, equally valid) |
| MinHash hash functions | Python random | Fixed-seed RNG (deterministic) |
| Multi-taxa optimization | Supported | Not implemented (not needed by BaitBench) |
| Blacklisting | Supported | Planned extension |
| Probe output | Equivalent coverage | Equivalent coverage; specific probes may differ |

**catch method (`--method catch`):**

Calls the external CATCH tool (`design.py`) from the Broad Institute. Requires the `catch` conda package (`conda install -c bioconda catch`). All `--catch-*` flags apply.

Example:

```bash
baitbench build-probes \
  --targets targets.fa \
  --method catch \
  --probe-length 120 \
  --catch-probe-stride 30 \
  --catch-mismatches 3 \
  --outdir probes_output
```

**syotti-lite method (`--method syotti-lite`):**

[Syotti](https://github.com/jnalanko/syotti) (Alanko et al. 2022) is a greedy set-cover bait designer. It scans the input sequences; at every uncovered position, it extracts a bait of `--probe-length` bp and marks all reference windows within `--syotti-mismatches` Hamming distance as covered (checking both strands). This is more targeted than tiling — probes are only generated where coverage is not already achieved by an earlier probe, yielding a smaller set while guaranteeing full coverage. Probes are named `probe_{target_id}|syotti_{n}`.

The BaitBench implementation replaces the original FM-index (SDSL/Divsufsort C++ libraries) with a k-mer hash index, which is well-suited to the MB-scale inputs typical in BaitBench and requires no additional dependencies.

Key design decisions:

| Decision | Choice | Rationale |
|----------|--------|-----------|
| Index type | k-mer HashMap | No external dependencies; correct and fast for MB-scale inputs |
| N handling | N always mismatches (not even N≡N) | Matches Syotti paper semantics |
| RC matching | Query RC of bait against forward index | Equivalent to bidirectional index; simpler |
| Seed guarantee | All overlapping seeds checked | Correct for mismatches ≤ probe_length − seed_len; robust beyond that |
| Output format | `probe_{id}\|syotti_{n}` | Consistent with tile format; compatible with downstream pipeline |

Example:

```bash
baitbench build-probes \
  --targets targets.fa \
  --method syotti-lite \
  --probe-length 120 \
  --syotti-mismatches 40 \
  --syotti-seed-len 20 \
  --outdir probes_output
```

Memory note: the k-mer index stores one entry per seed position per input base. For a 10 MB input with seed_len=20, expect ~500 MB peak memory. For very large inputs, reduce `--syotti-seed-len` (shorter seeds use more memory) or pre-filter with collapse.

> Alanko JN, Slizovskiy IB, Lokshtanov D, Gagie T, Noyes NR, Boucher C. "Syotti: scalable bait design for DNA enrichment." *Bioinformatics.* 2022;38(Supplement_1):i177–i184. doi:10.1093/bioinformatics/btac226

**probetools-lite method (`--method probetools-lite`):**

BaitBench includes a native Rust reimplementation of the ProbeTools algorithm (Kuchinski et al. 2022, BMC Genomics). ProbeTools takes a coverage-first, diversity-aware approach: rather than starting from every position and pruning (tiling, Syotti) or optimizing a coverage model (CATCH), it identifies which k-mers in the target space are the most *representative* of sequence diversity — then greedily builds a panel that iterates on remaining coverage gaps.

This makes it well-suited for highly variable targets (e.g., diverse virus families) where each probe needs to cover many slightly different sequence variants.

**Algorithm:**

1. **K-mer enumeration** — Extract all k-mers of length `--probe-length` from the input sequences using a sliding window with step `--pt-step`. K-mers with more than 50% N bases are discarded. With the default `--pt-step 1`, every overlapping window is enumerated; larger values reduce candidates and run time.

2. **Clustering** — Cluster all k-mers by sequence identity using cd-hit-est at `--pt-identity` threshold. Each cluster groups similar k-mers together; the centroid represents the most "central" sequence in that cluster.

3. **Ranking** — Sort cluster centroids by cluster size (number of members) in descending order. A large cluster means many similar k-mers exist across the input, so its centroid is a highly representative probe candidate that will hybridize broadly across variants.

4. **Batch selection** — Add the top `--pt-batch-size` ranked centroids to the probe panel.

5. **Coverage assessment** — Align the full accumulated panel against all original targets using minimap2 (`probe_align`), then compute per-position depth from the SAM output. For each target, calculate the fraction of positions with depth ≥ `--pt-min-depth`. The stopping metric is the **10th-percentile** of these fractions across all targets — meaning 90% of targets must reach the coverage goal before the algorithm considers itself done.

6. **Low-coverage extraction** — Find runs of positions where depth < `--pt-min-depth`. Expand any run shorter than `--probe-length` bp bidirectionally to exactly `--probe-length` bp (centred on the run, clipped to sequence boundaries). Merge overlapping expanded regions. Write these sub-sequences to a new FASTA.

7. **Iterate** — Repeat steps 1–6 on the under-covered sub-sequences. In each subsequent iteration the k-mer pool is drawn only from remaining problem regions, so probes become progressively more targeted to gaps.

**Algorithm substitutions vs. original ProbeTools:**

The original ProbeTools tool uses VSEARCH for k-mer clustering and BLAST+ for coverage assessment. BaitBench replaces both with tools already integrated into the binary:

| Step | Original ProbeTools | probetools-lite |
|------|--------------------|-----------------------|
| K-mer clustering | VSEARCH (`--cluster_fast`) | cd-hit-est (equivalent identity-based clustering; already used by build-probes) |
| Coverage assessment | blastn (`-task blastn-short`) | minimap2 (`probe_align`, via the embedded rammap library; no external process needed) |

Both substitutions produce equivalent results for probe design. cd-hit-est and VSEARCH use the same sequence identity model. minimap2 and BLAST both produce ungapped short-read alignments against targets of this size, and minimap2 is already embedded in the BaitBench binary as a compiled library (no external `minimap2` process is spawned).

**Termination conditions:**

The iterative loop exits as soon as *any* of the following six conditions is met:

| # | Condition | Log level | Notes |
|---|-----------|-----------|-------|
| 1 | 10th-pct coverage ≥ `--pt-coverage` | INFO | Normal success path |
| 2 | Panel size ≥ `--pt-max-panel-size` | INFO | Only applies when `--pt-max-panel-size` is set |
| 3 | Iteration count ≥ `--pt-max-iterations` | WARN | Hard safety cap; default 20 iterations |
| 4 | Coverage gain < `--pt-min-coverage-gain` | WARN | Stagnation: each new batch adds negligible coverage |
| 5 | Low-coverage extraction yields no sequences | INFO | All targets fully covered before goal was reached |
| 6 | Under-covered regions too short or all-N to yield k-mers | INFO | Rare; happens with very short or degenerate gap regions |

Conditions 3 and 4 are safety guards against non-termination. If the loop exits with a stagnation warning on real data, it typically means the remaining under-covered regions are too divergent from each other to form clusters at the current `--pt-identity` threshold. Solutions: lower `--pt-identity`, lower `--pt-coverage`, or increase `--pt-batch-size`.

**Parameters:**

| Flag | Default | Description |
|------|---------|-------------|
| `--pt-step` | 1 | Sliding window step between consecutive k-mers during enumeration. `1` = every position (densest, matching original ProbeTools default). Increasing to `30` or `60` speeds up clustering significantly and is reasonable for lower-diversity panels. |
| `--pt-identity` | 0.9 | cd-hit-est sequence identity threshold for k-mer clustering (0.0–1.0). Lower values merge more divergent k-mers into fewer, broader clusters, producing probes with wider cross-variant coverage. Higher values keep clusters tighter, preserving probe specificity. |
| `--pt-coverage` | 0.9 | Coverage goal (0.0–1.0). The loop continues until the 10th-percentile of per-target coverage fractions reaches this value. Setting 0.9 means the loop stops when 90% of targets have ≥ 90% of their positions covered at depth ≥ `--pt-min-depth`. |
| `--pt-batch-size` | 100 | Probes added per iteration. Larger batches converge faster but may overshoot the goal; smaller batches are more precise but require more iterations and minimap2 alignments. For large, diverse panels a value of 200–500 is practical. |
| `--pt-max-panel-size` | (none) | Hard cap on total probes. Stops once this many probes are selected even if the coverage goal is not met. Use when targeting a fixed array format with a known probe budget. |
| `--pt-min-depth` | 1 | Minimum per-position depth to count a position as covered. `1` (default) counts any alignment as sufficient. Set to `2` or higher to require redundant probe coverage at every position. |
| `--pt-max-iterations` | 20 | Hard iteration cap. Always terminates after this many iterations regardless of coverage progress. Prevents infinite loops on pathological inputs such as highly repetitive regions or very divergent outlier sequences. A WARN is logged when this limit fires. |
| `--pt-min-coverage-gain` | 0.001 | Stagnation threshold (0.0–1.0). The loop stops if the 10th-percentile coverage improves by less than this fraction between two consecutive iterations. At the default of `0.001`, the loop stops if a full batch of probes moves coverage by less than 0.1 percentage points. A WARN is logged when stagnation is detected. |

Example for a diverse viral panel (dense enumeration, lower identity threshold):

```bash
baitbench build-probes \
  --targets targets.fa \
  --method probetools-lite \
  --probe-length 120 \
  --pt-step 1 \
  --pt-identity 0.85 \
  --pt-coverage 0.9 \
  --pt-batch-size 200 \
  --pt-max-panel-size 5000 \
  --pt-max-iterations 15 \
  --outdir probes_output
```

Example for a lower-diversity panel (sparser enumeration, faster):

```bash
baitbench build-probes \
  --targets targets.fa \
  --method probetools-lite \
  --probe-length 120 \
  --pt-step 30 \
  --pt-identity 0.9 \
  --pt-coverage 0.95 \
  --pt-batch-size 50 \
  --outdir probes_output
```

Temporary working files are written to `<outdir>/probetools_work/` during the run and are removed automatically if `--cleanup` is specified.

> Kuchinski KS, Christropher-Hennings J, Bhide K, Bhide M. "ProbeTools: designing hybridization probes for targeted genomic sequencing of diverse and hypervariable viral taxa." *BMC Genomics.* 2022;23(1):579. doi:10.1186/s12864-022-08790-4

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--targets` | required | Input target sequences FASTA |
| `--method` | tile | Probe construction method: `tile`, `catch-lite`, `syotti-lite`, `catch`, or `probetools-lite` |
| `--probe-length` | 120 | Probe length in bp |
| `--step` | -60 | Step from end of previous probe. Negative = overlap, 0 = tiled, positive = gap. Only used with `--method tile`. |
| `--catch-probe-stride` | 60 | Step between candidate probes (bp). Used with `--method catch-lite` and `--method catch`. |
| `--catch-mismatches` | 5 | Mismatches tolerated for a probe to cover a target window. Used with `--method catch-lite` and `--method catch`. |
| `--catch-extension` | 0 | Flanking bp beyond probe boundaries counted as covered. Used with `--method catch-lite` and `--method catch`. |
| `--catch-coverage` | 1.0 | Fraction of each target that must be covered (0.0–1.0). Used with `--method catch-lite` and `--method catch`. |
| `--catch-minhash-threshold` | 0.6 | Jaccard similarity threshold for near-deduplication; 0.0 disables. Used with `--method catch-lite` and `--method catch`. |
| `--syotti-mismatches` | 40 | Maximum Hamming distance for a bait to cover a reference window. Only used with `--method syotti-lite`. |
| `--syotti-seed-len` | 20 | K-mer seed length for Syotti approximate matching. Only used with `--method syotti-lite`. |
| `--pt-step` | 1 | Sliding window step between k-mer enumeration windows. Only used with `--method probetools-lite`. |
| `--pt-identity` | 0.9 | cd-hit-est identity threshold for k-mer clustering (0.0–1.0). Only used with `--method probetools-lite`. |
| `--pt-coverage` | 0.9 | 10th-percentile coverage goal (0.0–1.0). Only used with `--method probetools-lite`. |
| `--pt-batch-size` | 100 | Probes added per iteration. Only used with `--method probetools-lite`. |
| `--pt-max-panel-size` | (none) | Hard cap on total probes; no limit if omitted. Only used with `--method probetools-lite`. |
| `--pt-min-depth` | 1 | Minimum per-position depth to count as covered. Only used with `--method probetools-lite`. |
| `--pt-max-iterations` | 20 | Hard iteration cap (termination guard). Only used with `--method probetools-lite`. |
| `--pt-min-coverage-gain` | 0.001 | Stagnation threshold; stops if coverage improvement per iteration falls below this (termination guard). Only used with `--method probetools-lite`. |
| `--min-gc` | 0.20 | Minimum GC fraction (0–1) |
| `--max-gc` | 0.80 | Maximum GC fraction (0–1) |
| `--max-n-frac` | 0.05 | Maximum fraction of ambiguous (non-ACGT) bases in a target sequence (0–1). Targets exceeding this are removed before collapse. |
| `--dust-threshold` | 2.0 | sDUST score threshold *T* for low-complexity detection |
| `--dust-window` | 64 | sDUST window size *W* in bases |
| `--max-masked-frac` | 0.25 | Maximum fraction of bases masked by sDUST to keep a probe (0–1). Set to 1.0 to disable. |
| `--collapse-threshold` | 0.95 | cd-hit-est identity threshold for initial collapse |
| `--dedup-threshold` | 0.95 | cd-hit-est identity threshold for final dedup |
| `--threads` | 5 | Threads for cd-hit-est |
| `--genomes` | none | Genome FASTA(s) to check cross-reactivity against (assessment step) |
| `--threshold` | 80.0 | Homology threshold for cross-reactivity (assessment step) |
| `--minimap-preset` | sr | Minimap2 alignment preset (assessment step) |
| `--proximity` | 50 | Pull-down zone distance in bp (assessment step) |
| `--skip-assess` | false | Skip automatic probe assessment after building |
| `--outdir` | ./build_probes_results | Output directory |
| `--report` | full | Report mode (full, none, rmd) |
| `--cleanup` | false | Delete intermediate files |
| `--refine-iterations` | none | Number of refinement iterations on low-coverage targets (assessment step; mutually exclusive with `--refine-until-stable`) |
| `--refine-until-stable` | false | Repeat refinement until no targets remain below the threshold or set stabilizes (assessment step; mutually exclusive with `--refine-iterations`) |
| `--refine-threshold` | 80.0 | 1X coverage threshold (%) for refinement iterations (assessment step) |

**Auto-assessment:**

After building probes, `build-probes` automatically chains into `assess-probes` which runs probe coverage analysis and self-homology cross-reactivity. If `--genomes` is specified, cross-reactivity is also checked against those genomes. The combined report includes both build pipeline statistics and assessment results. Use `--skip-assess` to produce only the build pipeline output.

**Complexity filtering (sDUST):**

Low-complexity sequences (homopolymer runs, dinucleotide repeats, etc.) make poor probes because they hybridize non-specifically. The sDUST algorithm identifies low-complexity regions by computing a score based on triplet (3-mer) frequencies within sliding windows. Regions where a single triplet dominates receive high scores. The `--dust-threshold` parameter controls the sensitivity (lower = more aggressive masking). The default threshold of 2.0 and window size of 64 match NCBI's dustmasker defaults.

> Morgulis A, Gertz EM, Schäffer AA, Agarwala R. "A Fast and Symmetric DUST Implementation to Mask Low-Complexity DNA Sequences." *J Comput Biol.* 2006;13(5):1028-1040. doi:10.1089/cmb.2006.13.1028

**Output files:**

- `probes_final.fa` -- final deduplicated probe set
- `build_probes_stats.tsv` -- sequence/base counts at each pipeline step
- `assess_probes_report.html` -- combined HTML report with build stats + assessment (unless `--skip-assess`)
- `cov_probe_coverage_summary.tsv` -- per-target coverage statistics (from assessment)
- `cov_probe_depth.tsv` -- probe depth intervals (from assessment)
- `xreact_hits.tsv` -- cross-reactivity hits (from assessment)
- `xreact_summary.tsv` -- cross-reactivity summary (from assessment)

With `--skip-assess`, only produces `probes_final.fa`, `build_probes_stats.tsv`, and optionally `build_probes_report.html`.

### tool

Standalone utility tools grouped under a single subcommand. Run `baitbench tool --help` to list available tools.

```bash
baitbench tool <TOOL> [OPTIONS]
```

#### tool syotti

Run the Syotti greedy bait design algorithm directly, without the `build-probes` pipeline (no collapse, GC filter, complexity filter, or deduplication). Useful for testing Syotti parameters in isolation or when you want to bypass the pipeline.

```bash
baitbench tool syotti \
  --targets targets.fa \
  --output probes.fa \
  [--probe-length 120] \
  [--mismatches 40] \
  [--seed-len 20]
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--targets` | required | Input target sequences FASTA |
| `--output` | required | Output probe sequences FASTA |
| `--probe-length` | 120 | Probe (bait) length in bp |
| `--mismatches` | 40 | Maximum Hamming distance for a bait to cover a reference window. N never matches. |
| `--seed-len` | 20 | K-mer seed length for approximate matching. Matching is guaranteed correct when mismatches ≤ probe_length − seed_len. |

#### tool catch

Run the CATCH optimization probe design algorithm directly, without the `build-probes` pipeline.

```bash
baitbench tool catch \
  --targets targets.fa \
  --output probes.fa \
  [--probe-length 120] \
  [--stride 60] \
  [--mismatches 5] \
  [--extension 0] \
  [--coverage 1.0] \
  [--minhash-threshold 0.6]
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--targets` | required | Input target sequences FASTA |
| `--output` | required | Output probe sequences FASTA |
| `--probe-length` | 120 | Probe (bait) length in bp |
| `--stride` | 60 | Tiling step in bp |
| `--mismatches` | 5 | Maximum mismatches for a probe to cover a window |
| `--extension` | 0 | Extension length on each side of candidate probes |
| `--coverage` | 1.0 | Fraction of each target that must be covered |
| `--minhash-threshold` | 0.6 | MinHash Jaccard similarity threshold for deduplication |

#### tool dustview

Visualize sDUST low-complexity masking on FASTA sequences. Outputs to stdout: original sequence, masked sequence (X marks low-complexity regions), and per-sequence statistics.

```bash
baitbench tool dustview [input.fa] [--dust-threshold 2.0] [--dust-window 64]
# or from stdin:
cat sequences.fa | baitbench tool dustview
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `input` | stdin | Input FASTA file (positional, optional) |
| `--dust-threshold` | 2.0 | DUST score threshold — positions above this are masked |
| `--dust-window` | 64 | DUST sliding window size in bases |

#### tool collapse

Cluster near-duplicate sequences using cd-hit-est and write cluster representatives to a FASTA file.

```bash
baitbench tool collapse \
  --input sequences.fa \
  --output collapsed.fa \
  [--threshold 0.95] \
  [--threads 1] \
  [--log-file cdhit.log]
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input` | required | Input FASTA file |
| `--output` | required | Output FASTA file (cluster representatives) |
| `--threshold` | 0.95 | Sequence identity threshold for clustering |
| `--threads` | 1 | Number of threads for cd-hit-est |
| `--log-file` | cdhit.log | Path to write cd-hit-est log output |


### assess-probes

Standalone combined probe assessment. Runs probe coverage analysis and cross-reactivity analysis (self-homology always; against genomes if `--genomes` provided), producing a single combined HTML report.

```bash
baitbench assess-probes \
  --targets targets.fa \
  --probes probes.fa \
  [--genomes genome1.fa genome2.fa ...] \
  [--threshold 80.0] \
  [--minimap-preset sr] \
  [--proximity 50] \
  [--outdir assess_probes_results] \
  [--output-prefix ""] \
  [--report full|none|rmd] \
  [--cleanup] \
  [--all-individual-targets] \
  [--refine-iterations N | --refine-until-stable] \
  [--refine-threshold 80.0]
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--targets` | required | Target sequences FASTA |
| `--probes` | required | Probe sequences FASTA |
| `--genomes` | none | Genome FASTA(s) to check cross-reactivity against (repeatable) |
| `--threshold` | 80.0 | Minimum homology % to report cross-reactive hits |
| `--minimap-preset` | sr | Minimap2 alignment preset |
| `--proximity` | 50 | Pull-down zone distance in bp |
| `--outdir` | ./assess_probes_results | Output directory |
| `--output-prefix` | (empty) | String prepended to every output filename |
| `--report` | full | Report mode: `full` (HTML), `none` (skip), `rmd` (editable RMarkdown) |
| `--cleanup` | false | Delete intermediate files (SAM, logs) after completion |
| `--all-individual-targets` | false | Also compute probe coverage for each target individually. Runs minimap2 once per target against that target alone, eliminating all probe competition from similar targets. Produces `individual_target_coverage_summary.tsv` and adds an **Individual Target Coverage** section to the report. |
| `--refine-iterations` | none | Number of refinement iterations (mutually exclusive with `--refine-until-stable`) |
| `--refine-until-stable` | false | Repeat refinement until no targets remain below the threshold or the set stops changing (mutually exclusive with `--refine-iterations`) |
| `--refine-threshold` | 80.0 | 1X coverage threshold (%) used to identify low-coverage targets for refinement |

**Output files:**

- `cov_probe_coverage_summary.tsv` -- per-target coverage statistics (pangenome alignment)
- `cov_probe_depth.tsv` -- run-length encoded probe depth intervals
- `cov_multi_mapping_probes.tsv` -- probes mapping to multiple targets
- `xreact_hits.tsv` -- cross-reactivity hits above threshold
- `xreact_summary.tsv` -- per-probe cross-reactivity summary
- `assess_run_params.tsv` -- run parameters
- `individual_target_coverage_summary.tsv` -- per-target coverage without probe competition (only with `--all-individual-targets`)
- `assess_probes_report.html` -- combined HTML report (`--report full`)
- `assess_probes_report.Rmd` -- editable RMarkdown file (`--report rmd`)
- `refine_N_targets.fa` -- filtered targets for refinement iteration N (when `--refine-iterations` or `--refine-until-stable`)
- `refine_N_cov_probe_coverage_summary.tsv` -- coverage statistics for refinement iteration N
- `refine_N_probe_coverage_report.html` -- probe coverage report for refinement iteration N

**Report sections:**

1. **Probe Coverage** -- summary table, coverage breadth bar charts, individual target coverage (if `--all-individual-targets`), tiered coverage, gap analysis, pangenome depth (subtitle shows % pangenome ≥1X), depth profiles, proximity coverage, multi-mapping probes
2. **Self-Homology** -- heatmap (≤1000 probes), density plots, hits table
3. **Cross-Reactivity vs Genomes** (if `--genomes` provided) -- heatmap, per-genome bar chart, density plots, hits table
4. **Parameters** -- run configuration under a collapsible fold

**Refinement iterations:**

Many target panels contain highly similar sequences (e.g. closely related viruses or gene variants) that are unlikely to occur together in the same sample. When all targets are assessed together, probe coverage for any one target may appear low because the probes covering it also tile many similar targets — but those similar targets would not be present in real samples. Refinement iterations address this by re-running probe coverage on only the targets that showed poor coverage (below `--refine-threshold`), so you can assess how well the probes cover each subset in isolation.

- **`--refine-iterations N`** runs exactly N additional probe-coverage-only analyses after the initial full assessment. Each iteration filters to targets with `pct_covered_1x < --refine-threshold` from the previous iteration's summary, runs probe coverage on that subset, and produces a separate `refine_N_probe_coverage_report.html`. Stops early if no targets remain below the threshold.
- **`--refine-until-stable`** repeats automatically until no targets fall below the threshold, or until the set of low-coverage targets stops changing between iterations (indicating no further improvement is possible).
- **`--refine-threshold`** (default 80.0) sets the 1X coverage percentage below which a target is considered poorly covered and included in the next refinement iteration. Applies to both modes.

The refinement reports are probe-coverage-only (no cross-reactivity re-analysis, since that does not depend on the target subset).

---

## Parameter Reference

### Input Files

| Parameter | Flag | Default | Applies to | Description |
|-----------|------|---------|------------|-------------|
| Targets | `--targets`, `-t` | required | run, prepare, probe-coverage, coverage-curve | FASTA of target sequences the probes are designed to capture |
| Genomes | `--genomes`, `-g` | none | run, prepare, coverage-curve | FASTA of full genomes for fragment generation (genome mode) |
| Distractors | `--distractors`, `-d` | required | run, prepare, coverage-curve | FASTA of background sequences that should not be captured. Can be specified multiple times to provide multiple distractor files |
| Probes | `--probes`, `-p` | required | run, capture, probe-coverage, coverage-curve | FASTA of probe sequences |
| Sample | `--sample` | all targets | run, coverage-curve | Sample targets or genomes: TSV file path OR inline IDs with optional weights. See [Sample Manifest Format](#sample-manifest-format) |
| Sample-target map | `--sample-target-map` | none | run, prepare, coverage-curve | TSV mapping genome IDs to target IDs (genome mode). See [Sample-Target Map Format](#sample-target-map-format) |
| Groups | `--groups` | none | run, prepare | TSV mapping target sequence IDs to group names for group-level metrics. See [Groups File Format](#groups-file-format) |
| Distractor groups | `--distractor-groups` | none | run, prepare | TSV mapping distractor sequence IDs to group names (overrides default file-stem grouping). See [Groups File Format](#groups-file-format) |
| Host FASTA | `--host-fasta` | none | run, coverage-curve | Host genome for read filtering |

### Fragment Generation

| Parameter | Flag | Default | Description |
|-----------|------|---------|-------------|
| Num fragments | `--num-fragments` | 10000 | Number of fragments to generate. More fragments = better statistical power but slower |
| Fragment length mean | `--fragment-length-mean` | 175 | Mean fragment length in bp. Center of the normal distribution |
| Fragment length min | `--fragment-length-min` | 150 | Minimum fragment length in bp. Fragments shorter than this are discarded |
| Fragment length max | `--fragment-length-max` | 200 | Maximum fragment length in bp. Fragments longer than this are truncated |

### Target Abundance

| Parameter | Flag | Default | Description |
|-----------|------|---------|-------------|
| Distractor fraction | `--distractor-fraction`, `-f` | 0.9 | Fraction of fragments from distractor sequences (0-1). Higher = lower target abundance. **Mutually exclusive with `--ct`** |
| CT score | `--ct` | none | qPCR CT (cycle threshold) score. Converted to distractor fraction via calibration formula. Lower CT = more target. **Mutually exclusive with `--distractor-fraction`** |

If neither `--distractor-fraction` nor `--ct` is specified, defaults to a distractor fraction of 0.9 (10% target).

### CT Score Parameters

These parameters calibrate the CT-to-fraction conversion. Only relevant when using `--ct`.

| Parameter | Flag | Default | Description |
|-----------|------|---------|-------------|
| CT baseline | `--ct-baseline` | 20.0 | The CT value at which the target fraction equals the baseline fraction. Blocked by `--ct-calibration` |
| CT baseline fraction | `--ct-baseline-fraction` | 0.01 | The target fraction at the baseline CT value. Blocked by `--ct-calibration` |
| CT efficiency | `--ct-efficiency` | 1.0 | PCR amplification efficiency (0–1). 1.0 = perfect doubling per cycle; typical assays run at 0.90–0.98. Blocked by `--ct-calibration` |
| CT calibration | `--ct-calibration` | — | Two `"CT,fraction"` reference points (e.g. `"20.0,0.01" "30.0,0.00001"`). Derives efficiency automatically from the slope; replaces `--ct-baseline`, `--ct-baseline-fraction`, and `--ct-efficiency` |

See [CT Score Calculation](#ct-score-calculation) for details.

### Simulation Parameters

| Parameter | Flag | Default | Description |
|-----------|------|---------|-------------|
| Simulate mode | `--simulate-mode` | thermodynamic | `thermodynamic` (TNN Boltzmann weighting) or `simple` (uniform probe-site weights) |
| Hybridization temperature | `--hybridization-temperature` | 70.0 | Fixed hybridization temperature in °C; only used in thermodynamic mode. Use `--hybridization-temperature-values` to sweep |
| Hybridization temperature sweep | `--hybridization-temperature-values` | — | Space-separated temperatures to sweep in `coverage-curve` (e.g. `55 65 70 75`). Conflicts with `--hybridization-temperature`. Thermodynamic mode only |
| Capture fraction | `--capture-fraction` | 0.5 | Fraction of fragments from probe binding sites (0.0–1.0); remainder are background |

### Sequencing Parameters

| Parameter | Flag | Default | Description |
|-----------|------|---------|-------------|
| Read length | `--read-length` | 120 | Trim fragments to this length (bp). Used by `perfect` and `art`. Not applicable for `badread` |
| Num sequences | `--num-sequences` | all | Number of reads to sample with replacement. If not set, all fragments become reads. Models sequencing depth control |
| Read simulator | `--read-simulator` | `perfect` | `perfect` (trim, no errors), `art` (Illumina via ART-modern), `badread` (long reads — ONT or PacBio CLR) |
| Sequencer profile | `--sequencer-profile` | `HiSeq2500_150bp` / `ont` | Chemistry / error model. Required for `art` and `badread`. See [sequence command docs](#sequence) for details |
| Coverage depth | `--coverage-depth` | 1.0 | Reads generated per fragment for `art`/`badread`. With `badread`, depth=1 ≈ 1 read per captured fragment |
| Paired-end | `--paired-end` | false | Paired-end output (art only). Produces reads.fa + reads_R2.fa |
| PE fragment mean | `--pe-frag-len-mean` | 200 | Mean insert size for paired-end (`art` + `--paired-end`) |
| PE fragment SD | `--pe-frag-len-sd` | 50 | Insert size std-dev for paired-end (`art` + `--paired-end`) |

#### Minimap2 Preset Auto-Selection

When `--minimap-preset` (and `--host-minimap-preset`) are not specified, `baitbench run` automatically picks the appropriate preset based on `--read-simulator` and `--sequencer-profile`:

| Simulator | Profile | Auto preset |
|-----------|---------|-------------|
| `perfect` | — | `sr` |
| `art` | any | `sr` |
| `badread` | `ont`, `ont-2020` | `map-ont` |
| `badread` | `pacbio` | `map-pb` |

You can always override by passing `--minimap-preset` explicitly.

### Execution Parameters

| Parameter | Flag | Default | Description |
|-----------|------|---------|-------------|
| Threads | `--threads` | 1 | Number of threads for external tools (minimap2, BLAST) |
| Output dir | `--outdir`, `-o` | ./results | Output directory. A timestamped subdirectory is created for each run |
| Output prefix | `--output-prefix` | (empty) | String prepended to every auto-generated output filename. Available on `run`, `prepare`, `probe-coverage`, `coverage-curve`, `xreact`, `panel-qc`, `identify`. E.g., `--output-prefix myrun_` produces `myrun_results.tsv` instead of `results.tsv` |
| Run name | `--run-name` | auto | Custom name for the run. Default: `run_YYYYMMDD_HHMMSS` |
| Report mode | `--report` | full | Report output: `full` (render HTML), `none` (skip), `rmd` (editable RMarkdown file) |
| Seed | `--seed`, `-s` | random | Random seed for reproducibility. If not set, results vary between runs |
| Verbose | `--verbose` | false | Enable debug logging (global flag) |
| Minimap preset | `--minimap-preset` | auto | Minimap2 preset for read mapping. Auto-selected based on `--read-simulator` and `--sequencer-profile` (see below) |
| Host minimap preset | `--host-minimap-preset` | auto | Minimap2 preset for host read filtering. Same auto-selection as `--minimap-preset` |
| Cleanup | `--cleanup` | false | Delete intermediate files after completion, keeping only report inputs and final outputs. Available on `run`, `coverage-curve`, `probe-coverage`, and `xreact` |
| Identify | `--identify` | false | Enable species-level identification after metrics (genome mode only, requires `--sample-target-map`). Available on `run` |
| Identity threshold | `--identity-threshold` | 90.0 | Minimum sequence identity % to consider targets "similar" for species identification. Available on `run`, `panel-qc`, `identify` |
| Min unique targets | `--min-unique-targets` | 1 | Minimum unique target detections required to call a species PRESENT. Available on `run`, `identify` |

---

## CT Score Calculation

CT (cycle threshold) scores from qPCR provide an intuitive way to express target abundance. In qPCR, each cycle amplifies the DNA by a factor of `(1 + E)`, where E is the amplification efficiency. At 100% efficiency, each cycle doubles the DNA. Lower CT = more target DNA.

### The Formula

```
target_fraction = ct_baseline_fraction * (1 + efficiency)^(ct_baseline - ct)
distractor_fraction = 1 - target_fraction
```

Where:
- `ct_baseline` is a known CT value (default: 20.0)
- `ct_baseline_fraction` is the target fraction at that CT (default: 0.01)
- `efficiency` is the PCR amplification efficiency (default: 1.0 = 100%)
- `ct` is the CT value you want to simulate

### Default Calibration

With defaults (`--ct-baseline 20.0`, `--ct-baseline-fraction 0.01`, `--ct-efficiency 1.0`), the interpretation is: "at CT 20, 1% of DNA is from targets, assuming 100% PCR efficiency."

### CT Reference Table

Values below use default parameters (efficiency = 1.0, baseline CT 20 = 1% target):

| CT | Target fraction | Distractor fraction | Interpretation |
|----|-----------------|---------------------|----------------|
| 10 | 100%* | 0% | Pure target (capped at 100%) |
| 15 | 32% | 68% | Very high abundance |
| 18 | 4% | 96% | High abundance |
| 20 | 1% | 99% | Moderate (baseline) |
| 22 | 0.25% | 99.75% | Low-moderate |
| 25 | 0.031% | 99.97% | Low abundance |
| 28 | 0.004% | 99.996% | Very low |
| 30 | 0.001% | 99.999% | Near limit of detection |
| 35 | 0.00003% | ~100% | Extremely low |
| 40 | 0.000001% | ~100% | At qPCR detection limit |

*Target fractions above 100% are capped at 100% (distractor fraction = 0).

At lower efficiency (e.g. 95%), the same CT maps to a slightly higher target fraction because each cycle amplifies less — reaching CT 25 requires more starting material.

### PCR Efficiency

Real qPCR assays typically run at 90–98% efficiency. The default assumption of 100% efficiency (`--ct-efficiency 1.0`) is an idealisation that can overestimate how much the target is diluted. Specify a measured efficiency with `--ct-efficiency`:

```bash
baitbench run \
  --targets targets.fa \
  --distractors distractors.fa \
  --probes probes.fa \
  --ct 25 \
  --ct-efficiency 0.95 \
  --num-fragments 10000 \
  --outdir results
```

Assay efficiency is usually reported in kit documentation or can be measured from a standard curve: `E = 10^(-1/slope) - 1`, where slope is the slope of CT vs. log10(concentration).

### One-Point Calibration

The default calibration assumes CT 20 = 1% target. Shift the entire curve with `--ct-baseline` and `--ct-baseline-fraction`:

**Example: Your lab data shows CT 25 = 0.1% target reads:**

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

This shifts the entire curve:
- CT 25 = 0.1% target (your calibration point)
- CT 30 = 0.003% target (5 CT higher = 32× less)
- CT 20 = 3.2% target (5 CT lower = 32× more)

### Two-Point Calibration

If you have two reference samples with known target fractions and their CT values (e.g. from ddPCR-quantified standards), use `--ct-calibration` to derive the efficiency automatically. This eliminates all modelling assumptions about efficiency and the baseline:

```bash
baitbench run \
  --targets targets.fa \
  --distractors distractors.fa \
  --probes probes.fa \
  --ct 25 \
  --ct-calibration "20.0,0.01" "30.0,0.00001" \
  --num-fragments 10000 \
  --outdir results
```

The formula used to derive efficiency from the two points is:

```
E = (f1 / f2)^(1 / (ct2 - ct1)) - 1
```

Where `(ct1, f1)` and `(ct2, f2)` are the two calibration points. The first point also serves as the baseline. The derived efficiency is logged at run time so you can inspect it.

`--ct-calibration` conflicts with `--ct-baseline`, `--ct-baseline-fraction`, and `--ct-efficiency` (it replaces all three).

**Which calibration method to use:**

| Situation | Recommended approach |
|-----------|---------------------|
| Quick simulation with reasonable defaults | `--ct` only |
| Known assay efficiency from kit docs | `--ct --ct-efficiency 0.95` |
| One reference sample with known fraction | `--ct --ct-baseline` + `--ct-baseline-fraction` |
| Two ddPCR-quantified reference standards | `--ct --ct-calibration "CT1,frac1" "CT2,frac2"` |

### Tips for Using CT Scores

- **Match your experimental system.** If you have empirical data linking CT to target fraction, use calibration flags to match your curve rather than relying on defaults.
- **Use coverage-curve to sweep.** The `coverage-curve` command with `--ct-values` lets you visualize performance across a range of CT values in a single analysis. Calibration flags apply to every value in the sweep.
- **Remember the log scale.** Each CT unit represents a `(1 + E)`-fold change. At 100% efficiency a 10-CT range spans ~1000-fold differences in abundance; at 95% efficiency it spans ~614-fold.

---

## Output Files

### Run Output Directory

Each `baitbench run` creates a timestamped subdirectory:

```
results/run_20250101_120000/
├── combined_reference.fa       # All sequences merged for fragment generation
├── mapping_reference.fa        # Targets + distractors for mapping (genome mode only)
├── weights.txt                 # Per-sequence sampling weights
├── targets.txt                 # Target sequence IDs
├── distractors.txt             # Distractor sequence IDs
├── sample.txt                  # Sample sequence IDs
├── genomes.txt                 # Genome IDs (genome mode only)
├── sample_target_map.txt       # Genome-to-target mapping (genome mode only)
├── fragments.fa                # Simulated DNA fragments (probe-biased + background)
├── reads.fa                    # Sequencing reads (trimmed to read length)
├── filtered.fa                 # Host-filtered reads (if --host-fasta)
├── mapped.sam                  # Read alignments to reference
├── detected.list               # Read counts per reference
├── run_params.tsv              # Run configuration (used by report)
├── target_groups.tsv           # Target group assignments (if --groups)
├── distractor_groups.tsv       # Distractor group assignments (always; auto or from --distractor-groups)
├── results.tsv                 # Summary metrics
├── detected_detail.tsv         # Per-reference detection and coverage detail
├── group_detail.tsv            # Per-group summary (if groups are present)
├── results.json                # Machine-readable JSON metrics
├── coverage.tsv                # Run-length encoded read depth intervals
├── report.html                 # HTML report (--report full, requires R)
├── report.Rmd                  # Editable RMarkdown (--report rmd)
├── species_calls.tsv           # Species-level calls (if --identify)
├── species_calls.json          # Species calls JSON (if --identify)
├── target_similarity.tsv       # Target pairwise similarity (if --identify)
├── capture.log                 # Capture alignment log
├── mapping.log                 # Read mapping log
└── host_filter.log             # Host filtering log (if --host-fasta)
```

### results.tsv Columns

| Column | Description |
|--------|-------------|
| `run_name` | Run identifier |
| `timestamp` | Completion time |
| `num_fragments` | Fragments requested |
| `seed` | Random seed (or "NA") |
| `fragments_generated` | Fragments actually generated |
| `fragments_captured` | Fragments passing capture |
| `capture_rate` | fragments_captured / fragments_generated |
| `sample_captured` | Captured fragments from sample targets |
| `nonsample_target_captured` | Captured fragments from non-sample targets |
| `distractor_captured` | Captured fragments from distractors |
| `untargeted_captured` | Captured fragments from untargeted genomes (genome mode) |
| `reads_correctly_mapped` | Reads mapping to their source reference |
| `reads_incorrectly_mapped` | Reads mapping to a different reference |
| `sample_total` | Number of distinct sample targets |
| `nonsample_target_total` | Number of non-sample targets |
| `distractors_total` | Number of distractor sequences |
| `tp_count` | True Positives: sample targets detected |
| `fn_count` | False Negatives: sample targets not detected |
| `fp_target_count` | False Positives: non-sample targets detected |
| `fp_distractor_count` | False Positives: distractors detected |
| `fp_total` | fp_target_count + fp_distractor_count |
| `tn_target_count` | True Negatives: non-sample targets not detected |
| `tn_distractor_count` | True Negatives: distractors not detected |
| `tn_total` | tn_target_count + tn_distractor_count |
| `sensitivity` | TP / (TP + FN) |
| `specificity` | TN_total / (TN_total + FP_total) |
| `precision` | TP / (TP + FP_total) |
| `f1_score` | 2 * (precision * sensitivity) / (precision + sensitivity) |
| `reads_sequenced` | Number of reads after the sequencing step (0 if not tracked) |
| `reads_after_filter` | Number of reads after host filtering (0 if filter not applied) |
| `reads_mapped` | reads_correctly_mapped + reads_incorrectly_mapped |
| `reads_unmapped` | Reads that entered mapping but did not map to any reference |

### detected_detail.tsv Columns

One row per reference sequence:

| Column | Description |
|--------|-------------|
| `reference_id` | Sequence ID |
| `group` | Group name this sequence belongs to (sequence's own ID if no groups file provided) |
| `category` | `sample`, `nonsample_target`, `distractor`, or `untargeted` |
| `expected` | 1 if expected to be detected (sample target), 0 otherwise |
| `detected` | 1 if at least one read maps to this reference, 0 otherwise |
| `fragments_generated` | Number of fragments generated from this sequence |
| `fragments_captured` | Number of fragments captured by probes |
| `reads_assigned` | Number of reads mapped to this reference |
| `classification` | `TP`, `FN`, `FP_target`, `FP_distractor`, `TN_target`, `TN_distractor`, or `untargeted` |
| `ref_length` | Reference sequence length (bp) |
| `avg_coverage` | Average read depth across reference |
| `pct_covered_5x` | % positions with >= 5x depth |
| `pct_covered_20x` | % positions with >= 20x depth |

### group_detail.tsv Columns

Written when group files are present (`target_groups.tsv` or `distractor_groups.tsv`). One row per group:

| Column | Description |
|--------|-------------|
| `group_name` | Group identifier |
| `category` | `sample`, `nonsample_target`, or `distractor` |
| `expected` | `true` if the group is expected to be detected (sample group) |
| `detected` | `true` if at least one member sequence has reads mapped to it |
| `classification` | `TP`, `FN`, `FP_target`, `FP_distractor`, `TN_target`, or `TN_distractor` |
| `member_count` | Number of sequences in this group |
| `detected_member_count` | Number of member sequences individually detected |
| `total_reads` | Sum of reads assigned to all members of this group |

### results.json Structure

Structured JSON output with nested sections:

```json
{
  "run_info": {
    "run_name": "...",
    "timestamp": "...",
    "num_fragments": 10000,
    "seed": "42"
  },
  "capture_stats": {
    "fragments_generated": 10000,
    "fragments_captured": 3500,
    "capture_rate": 0.35
  },
  "read_level": {
    "reads_correctly_mapped": 3400,
    "reads_incorrectly_mapped": 100,
    "reads_mapped": 3500,
    "reads_unmapped": 0,
    "reads_sequenced": 3500,
    "reads_after_filter": 0
  },
  "metrics": {
    "sensitivity": 1.0,
    "specificity": 0.95,
    "precision": 0.8,
    "f1_score": 0.89,
    "tp": 5, "fn": 0,
    "fp_target": 2, "fp_distractor": 1,
    "tn_target": 10, "tn_distractor": 50
  },
  "details": [ ... ]
}
```

### coverage.tsv Format

Run-length encoded read depth intervals. Consecutive positions with the same depth are collapsed into a single interval (1-based inclusive coordinates):

```
reference_id	start	end	depth
dengue_1	1	50	0
dengue_1	51	100	3
dengue_1	101	200	5
...
```

This format is typically 100-1000x smaller than per-position output, making it feasible for large target panels.

---

## Usage Examples

### Basic Probe Evaluation

Test whether probes capture all targets and reject distractors:

```bash
baitbench run \
  --targets targets.fa \
  --distractors distractors.fa \
  --probes probes.fa \
  --num-fragments 10000 \
  --seed 42 \
  --outdir results
```

All targets are treated as "present" (no `--sample`). The default distractor fraction is 0.9 (90% background, 10% target).

### Sample Discrimination Testing

Test whether probes can detect specific targets while rejecting others in the panel:

```bash
# Inline sample IDs
baitbench run \
  --targets all_viruses.fa \
  --distractors bacteria.fa \
  --probes probes.fa \
  --sample dengue_1 zika_virus \
  --num-fragments 10000 \
  --outdir results

# With custom weights (dengue at 5x abundance)
baitbench run \
  --targets all_viruses.fa \
  --distractors bacteria.fa \
  --probes probes.fa \
  --sample dengue_1 5 zika_virus \
  --num-fragments 10000 \
  --outdir results

# Using a TSV manifest file
baitbench run \
  --targets all_viruses.fa \
  --distractors bacteria.fa \
  --probes probes.fa \
  --sample sample.tsv \
  --num-fragments 10000 \
  --outdir results
```

Non-sample targets will have FP_target classification if detected, testing cross-reactivity within the panel.

### Clinical Specimen Simulation with CT

Simulate specimens at different viral loads:

```bash
# High viral load (CT 20)
baitbench run \
  --targets targets.fa \
  --distractors distractors.fa \
  --probes probes.fa \
  --ct 20 \
  --num-fragments 10000 \
  --outdir results_ct20

# Low viral load (CT 30)
baitbench run \
  --targets targets.fa \
  --distractors distractors.fa \
  --probes probes.fa \
  --ct 30 \
  --num-fragments 10000 \
  --outdir results_ct30
```

### Genome Mode for Bacteria

When probe targets are sub-regions of large genomes:

```bash
# targets.fa: 16S gene sequences
# genomes.fa: full bacterial genomes
# mapping.tsv links genome IDs to target gene IDs

baitbench run \
  --targets 16S_targets.fa \
  --genomes bacteria_genomes.fa \
  --distractors human.fa \
  --probes probes.fa \
  --sample-target-map mapping.tsv \
  --sample e_coli s_aureus \
  --num-fragments 50000 \
  --outdir results
```

Use higher `--num-fragments` for bacteria since genomes are much larger than target regions, requiring more fragments to achieve adequate target coverage.

### Mixed Panels (Virus + Bacteria)

Genome mode handles mixed panels naturally. Virus genomes that match their target IDs auto-link:

```bash
# genomes.fa: influenza_a (13kb), e_coli (5Mb)
# targets.fa: influenza_a (same seq), e_coli_16S (1.5kb subsequence)
# mapping.tsv only needs the e_coli entry (influenza_a auto-links)

baitbench run \
  --targets targets.fa \
  --genomes genomes.fa \
  --distractors human.fa \
  --probes probes.fa \
  --sample-target-map mapping.tsv \
  --sample influenza_a e_coli \
  --num-fragments 50000 \
  --outdir results
```

### Multiple Distractor Sources

Provide multiple distractor FASTA files:

```bash
baitbench run \
  --targets targets.fa \
  --distractors bacteria.fa \
  --distractors fungi.fa \
  --distractors protozoa.fa \
  --probes probes.fa \
  --num-fragments 10000 \
  --outdir results
```

All distractor sequences are concatenated and share the same per-sequence weight. Each distinct `--distractors` file automatically forms its own distractor group (by file stem), so multiple distractor files are still counted separately in the group-level metrics.

### Group-Level Grouping

When your target panel contains multiple sequence variants of the same organism, use `--groups` to collapse them into a single entity for metrics:

```bash
# Create groups file
cat > groups.tsv <<'EOF'
West_Nile_virus_0001	West_Nile_virus
West_Nile_virus_0002	West_Nile_virus
West_Nile_virus_0003	West_Nile_virus
Dengue_virus_1_6275	Dengue_virus_1
Dengue_virus_1_2274	Dengue_virus_1
Dengue_virus_2_8773	Dengue_virus_2
Dengue_virus_2_1822	Dengue_virus_2
EOF

baitbench run \
  --targets all_variants.fa \
  --distractors Aaegypti.fa \
  --probes probes.fa \
  --sample West_Nile_virus_0001 \
  --groups groups.tsv \
  --num-fragments 10000 \
  --outdir results
```

With this configuration:
- Detection of **any** WNV variant counts as a TP for the `West_Nile_virus` group
- `Dengue_virus_1_6275` and `Dengue_virus_1_2274` are grouped as `Dengue_virus_1` (one entity for FP/TN counting)
- `Dengue_virus_2_*` is a separate group
- All ~2300 Aaegypti contigs are automatically grouped as `"Aaegypti"` (one FP_distractor)
- Cross-mapping between variants of the same group (e.g., WNV_0001 reads mapping to WNV_0002) is counted as **correctly mapped**

For distractors, the default file-stem grouping is automatic. To override it (e.g., if multiple organisms share one FASTA file), provide an explicit distractor groups file:

```bash
# Override automatic distractor grouping
cat > distractor_groups.tsv <<'EOF'
# seq_id	group_name
contig_001	Aedes_aegypti
contig_002	Aedes_aegypti
bacterial_16S_1	E_coli
bacterial_16S_2	E_coli
EOF

baitbench run \
  --targets targets.fa \
  --distractors mixed_distractors.fa \
  --probes probes.fa \
  --distractor-groups distractor_groups.tsv \
  --num-fragments 10000 \
  --outdir results
```

### Capture Fraction Sweep

Control what fraction of simulated fragments come from probe binding sites. With thermodynamic mode (default), high-affinity probe-reference alignments receive higher weight.

```bash
# Default: 50% probe-biased, 50% background (thermodynamic mode)
baitbench run \
  --targets targets.fa \
  --distractors distractors.fa \
  --probes probes.fa \
  --num-fragments 10000 \
  --capture-fraction 0.5 \
  --outdir results_thermo

# High capture fraction with simple (uniform) weighting
baitbench run \
  --targets targets.fa \
  --distractors distractors.fa \
  --probes probes.fa \
  --simulate-mode simple \
  --capture-fraction 0.8 \
  --num-fragments 10000 \
  --outdir results_simple

# Sweep capture fractions with coverage-curve
baitbench coverage-curve \
  --targets targets.fa \
  --distractors distractors.fa \
  --probes probes.fa \
  --sample target_1 target_2 \
  --capture-fraction-values 0.2 0.4 0.6 0.8 \
  --ct 25 \
  --outdir cf_sweep
```

### Sequencing Depth Control

Control the number of reads output by the sequencing step:

```bash
# Sample 5000 reads with replacement (models limited sequencing)
baitbench run \
  --targets targets.fa \
  --distractors distractors.fa \
  --probes probes.fa \
  --num-fragments 50000 \
  --num-sequences 5000 \
  --seed 42 \
  --outdir results
```

### Host Filtering

Remove host reads before mapping:

```bash
baitbench run \
  --targets targets.fa \
  --distractors distractors.fa \
  --probes probes.fa \
  --host-fasta human_genome.fa \
  --num-fragments 10000 \
  --outdir results
```

### Coverage Curve Analysis

Sweep parameters to understand how conditions affect coverage:

```bash
# Sweep CT values only
baitbench coverage-curve \
  --targets targets.fa \
  --distractors distractors.fa \
  --probes probes.fa \
  --sample dengue_1 zika_virus \
  --ct-values 20 25 30 35 \
  --num-fragments 10000 \
  --seed 42 \
  --outdir coverage_ct

# Sweep CT and capture fraction (combinatorial)
baitbench coverage-curve \
  --targets targets.fa \
  --distractors distractors.fa \
  --probes probes.fa \
  --sample dengue_1 \
  --ct-values 20 25 30 \
  --capture-fraction-values 0.3 0.5 0.7 \
  --num-fragments 10000 \
  --outdir coverage_ct_cf

# Sweep hybridization temperature (thermodynamic mode)
baitbench coverage-curve \
  --targets targets.fa \
  --distractors distractors.fa \
  --probes probes.fa \
  --sample dengue_1 \
  --ct 25 \
  --hybridization-temperature-values 55 60 65 70 75 \
  --num-fragments 10000 \
  --outdir coverage_temp

# Sweep all four parameters
baitbench coverage-curve \
  --targets targets.fa \
  --distractors distractors.fa \
  --probes probes.fa \
  --sample dengue_1 \
  --ct-values 20 25 30 \
  --hybridization-temperature-values 65 70 75 \
  --capture-fraction-values 0.3 0.5 0.7 \
  --num-sequences-values 500 1000 5000 \
  --num-fragments 10000 \
  --outdir coverage_full

# Fixed CT with capture fraction sweep
baitbench coverage-curve \
  --targets targets.fa \
  --distractors distractors.fa \
  --probes probes.fa \
  --sample dengue_1 \
  --ct 25 \
  --capture-fraction-values 0.1 0.3 0.5 0.7 0.9 \
  --num-fragments 10000 \
  --outdir coverage_cf
```

### Probe Design QC

Evaluate probe tiling independently of the simulation:

```bash
baitbench probe-coverage \
  --targets targets.fa \
  --probes probes.fa \
  --proximity 100 \
  --outdir probe_qc
```

### Cross-Reactivity Analysis

Check whether probes have off-target homology to specific genomes:

```bash
# Probe-to-genome: which probes match the human genome?
baitbench xreact \
  --probes probes.fa \
  --against human_genome.fa \
  --threshold 80 \
  --outdir xreact_human

# Probe-to-genome: check against multiple references
baitbench xreact \
  --probes probes.fa \
  --against human_genome.fa mouse_genome.fa \
  --threshold 80 \
  --outdir xreact_hosts

# Probe-to-probe: find probes that are too similar to each other
baitbench xreact \
  --probes probes.fa \
  --self \
  --threshold 80 \
  --outdir xreact_self

# Both modes together
baitbench xreact \
  --probes probes.fa \
  --against human_genome.fa \
  --self \
  --threshold 80 \
  --outdir xreact_full
```

### Target Panel QC

Assess whether a target panel can distinguish between species before running simulations:

```bash
# Basic panel QC
baitbench panel-qc \
  --targets gene_targets.fa \
  --sample-target-map mapping.tsv \
  --outdir panel_qc_results

# Stricter similarity threshold (95% instead of default 90%)
baitbench panel-qc \
  --targets gene_targets.fa \
  --sample-target-map mapping.tsv \
  --identity-threshold 95 \
  --outdir panel_qc_strict

# Skip HTML report (just produce TSV files)
baitbench panel-qc \
  --targets gene_targets.fa \
  --sample-target-map mapping.tsv \
  --report none \
  --outdir panel_qc_tsv
```

The HTML report includes a species discriminability chart, confusion matrix heatmap, and target composition breakdown.

### Species Identification

Call species from existing pipeline results or as part of `baitbench run`:

```bash
# Standalone: using pre-computed similarity from panel-qc
baitbench identify \
  --detected-detail results/run/detected_detail.tsv \
  --sample-target-map mapping.tsv \
  --target-similarity panel_qc_results/target_similarity.tsv \
  --outdir identify_results

# Standalone: compute similarity on-the-fly from target FASTA
baitbench identify \
  --detected-detail results/run/detected_detail.tsv \
  --sample-target-map mapping.tsv \
  --targets gene_targets.fa \
  --outdir identify_results

# Integrated into pipeline (genome mode)
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

# With stricter calling threshold (require 2 unique markers)
baitbench run \
  --targets gene_targets.fa \
  --genomes full_genomes.fa \
  --distractors human.fa \
  --probes probes.fa \
  --sample-target-map mapping.tsv \
  --sample e_coli influenza_a \
  --identify \
  --min-unique-targets 2 \
  --num-fragments 50000 \
  --outdir results
```

When `--identify` is used with `baitbench run`, the species calls are compared against the ground-truth `--sample` manifest and included in the HTML report.

### Probe Assessment

Run combined probe coverage + cross-reactivity analysis on an existing probe set:

```bash
# Basic assessment (probe coverage + self-homology)
baitbench assess-probes \
  --targets targets.fa \
  --probes probes.fa \
  --outdir assess_results

# With cross-reactivity against genomes
baitbench assess-probes \
  --targets targets.fa \
  --probes probes.fa \
  --genomes human_genome.fa other_genomes.fa \
  --threshold 80 \
  --outdir assess_results

# Skip HTML report (produce only TSV outputs)
baitbench assess-probes \
  --targets targets.fa \
  --probes probes.fa \
  --report none \
  --outdir assess_results

# Refinement: re-run coverage on low-coverage targets 3 times
baitbench assess-probes \
  --targets targets.fa \
  --probes probes.fa \
  --refine-iterations 3 \
  --refine-threshold 80 \
  --outdir assess_results

# Refinement: repeat until no targets remain below 80% 1X coverage
baitbench assess-probes \
  --targets targets.fa \
  --probes probes.fa \
  --refine-until-stable \
  --refine-threshold 80 \
  --outdir assess_results
```

Build probes and automatically assess them:

```bash
# Build + assess (default behavior)
baitbench build-probes \
  --targets targets.fa \
  --outdir probes_output

# Build + assess with cross-reactivity against genomes
baitbench build-probes \
  --targets targets.fa \
  --genomes human_genome.fa \
  --outdir probes_output

# Build only, skip assessment
baitbench build-probes \
  --targets targets.fa \
  --skip-assess \
  --outdir probes_output

# Build using external CATCH tool (requires catch conda package)
baitbench build-probes \
  --targets targets.fa \
  --method catch \
  --probe-length 120 \
  --outdir probes_output

# Build using catch-lite (native) with custom parameters
baitbench build-probes \
  --targets targets.fa \
  --method catch-lite \
  --catch-probe-stride 30 \
  --catch-mismatches 3 \
  --catch-extension 10 \
  --outdir probes_output
```

### Running Individual Steps

Run pipeline steps independently for custom workflows:

```bash
# 1. Prepare
baitbench prepare \
  --targets targets.fa \
  --distractors distractors.fa \
  --distractor-fraction 0.95 \
  --outdir prep

# 2. Simulate (probe alignment + TNN scoring + multinomial sampling)
baitbench simulate \
  --reference prep/combined_reference.fa \
  --weights prep/weights.txt \
  --probes probes.fa \
  --num-fragments 50000 \
  --capture-fraction 0.5 \
  --seed 42 \
  --output prep/fragments.fa

# 3. Sequence
baitbench sequence \
  --input prep/fragments.fa \
  --read-length 150 \
  --output prep/reads.fa

# 4. Map
baitbench map \
  --reference prep/combined_reference.fa \
  --reads prep/reads.fa \
  --output prep/mapped.sam

# 5. List
baitbench list \
  --sam prep/mapped.sam \
  --output prep/detected.list

# 6. Metrics
baitbench metrics \
  --targets prep/targets.txt \
  --distractors prep/distractors.txt \
  --sample prep/sample.txt \
  --detected prep/detected.list \
  --fragments prep/fragments.fa \
  --captured prep/fragments.fa \
  --sam prep/mapped.sam \
  --run-name custom_run \
  --num-fragments 50000 \
  --output-summary prep/results.tsv \
  --output-detail prep/detected_detail.tsv
```

### Reproducible Runs

Use `--seed` for reproducibility. The same seed with the same inputs produces identical results:

```bash
baitbench run \
  --targets targets.fa \
  --distractors distractors.fa \
  --probes probes.fa \
  --num-fragments 10000 \
  --seed 42 \
  --outdir results
```

### Batch Comparisons

Run multiple configurations and compare results:

```bash
# Compare different CT values
for ct in 20 25 30 35; do
  baitbench run \
    --targets targets.fa \
    --distractors distractors.fa \
    --probes probes.fa \
    --ct $ct \
    --num-fragments 10000 \
    --seed 42 \
    --report none \
    --outdir "results_ct${ct}"
done

# Aggregate results
head -1 results_ct20/*/results.tsv > comparison.tsv
for ct in 20 25 30 35; do
  tail -1 "results_ct${ct}"/*/results.tsv >> comparison.tsv
done
```

---

## Report Guide

### Report Modes

The `--report` flag controls report output for `run`, `probe-coverage`, `coverage-curve`, and `report` commands:

| Mode | Description |
|------|-------------|
| `full` (default) | Render the full HTML report using R/rmarkdown. Requires R and pandoc. |
| `none` | Skip report generation entirely. All other outputs (TSV, JSON) are still produced. |
| `rmd` | Write a parameterized RMarkdown (`.Rmd`) file with all file paths and parameters pre-filled. Does not require R at run time. |

**Using `--report rmd`:**

The `rmd` mode produces an `.Rmd` file in the output directory with all parameters baked in. You can then:

1. Open the `.Rmd` file in RStudio or any text editor
2. Customize the report -- add sections, change figures, adjust formatting
3. Render it when ready:

```bash
Rscript -e 'rmarkdown::render("results/run_20250101_120000/report.Rmd")'
```

This is useful when you want to:
- Customize the report before rendering
- Render on a different machine that has R installed
- Add project-specific analysis sections
- Iterate on the report without re-running the pipeline

### HTML Report Sections

The main pipeline report (`report.html`) includes:

1. **Run Parameters** -- Table of all configuration values. Also shows the reconstructed command line for reproducibility.

2. **Capture Summary** -- Bar chart comparing fragments generated vs captured, broken down by source (sample, non-sample target, distractor, untargeted).

3. **Detection Performance** -- Bar chart of sensitivity, specificity, precision, and F1 score.

4. **Read Mapping Accuracy** -- Correctly vs incorrectly mapped reads. Incorrect mapping indicates cross-reactivity (e.g., virus A reads mapping to virus B).

5. **Confusion Matrix** -- Heatmap showing TP, FN, FP, and TN counts.

6. **Detection Detail** -- Table of every reference sequence with detection status, fragment counts, read counts, and coverage statistics.

7. **Detection Lollipop** -- Reads per detected reference, colored by classification (TP, FP_target, FP_distractor).

8. **Coverage Plots** (if coverage data available) -- Per-position read depth plots for each detected reference, with faceted overview and expandable per-reference detail views.

### Coverage Curve Report

The coverage curve report (`coverage_curve_report.html`) shows:

- **Depth curves** -- % genome covered (Y-axis) vs depth of coverage threshold on log10 scale (X-axis), with one line per parameter combination
- **Faceting** -- With < 10 combinations, all lines on one plot. With >= 10 combinations, faceted by the parameter with the fewest levels
- **Per-target panels** -- If the sample contains multiple targets, each gets its own panel
- **Summary table** -- Key depth thresholds (1x, 5x, 10x, 20x, 50x, 100x) for each combination

### Probe Coverage Report

The probe coverage report (`probe_coverage_report.html`) shows:

- **Summary table** -- Per-target coverage statistics (adapts to dataset size: simple table for <= 20 targets, interactive DT table for > 20)
- **Coverage bar charts** -- Per-target 1x/2x/5x/10x coverage (switches to histograms/boxplots for > 100 targets)
- **Depth profiles** -- Per-position probe depth plots for each target (omitted for > 100 targets)
- **Gap analysis** -- Uncovered regions and gap statistics
- **Multi-mapping probes** -- Probes that align to multiple targets (specificity concerns)

### Panel QC Report

The panel QC report (`panel_qc_report.html`) shows:

- **Panel Summary** -- Total species, targets, similar pairs, and species with zero unique targets
- **Species Discriminability** -- Bar chart (≤50 species) or histogram (>50) of discriminability scores
- **Target Composition** -- Stacked bar chart of unique vs shared targets per species
- **Species Confusion Matrix** -- Heatmap (≤30 species) or distribution statistics (>30) of shared target counts
- **Discriminability Table** -- Full per-species discriminability data (simple table ≤20, interactive DT table >20)
- **Target Similarity Pairs** -- All pairwise target similarities above the threshold

### Species Identification in Main Report

When species calls are available (from `--identify` or standalone `baitbench identify`), the main HTML report includes a "Species Identification" section with:

- **Summary table** -- Species-level sensitivity and specificity (when ground truth is available via `--sample`)
- **Species call chart** -- Bar chart of PRESENT/ABSENT/AMBIGUOUS calls per species
- **Evidence detail table** -- Full breakdown with unique/shared detected counts, reads, and explanation

### Probe Assessment Report

The probe assessment report (`assess_probes_report.html`) combines coverage and cross-reactivity analysis into a single document:

- **Build Pipeline** (conditional, when chained from build-probes) -- Pipeline stats table, sequence/base count bar charts
- **Probe Coverage** -- Summary table, coverage breadth bar charts, individual target coverage (conditional, when `--all-individual-targets` was used), tiered coverage, gap analysis, pangenome depth (subtitle includes % pangenome ≥1X), per-target depth profiles, proximity coverage, multi-mapping probes
- **Self-Homology** -- Plotly heatmap (≤1000 probes), density plots, hits table
- **Cross-Reactivity vs Genomes** (conditional, when `--genomes` provided) -- Plotly heatmap, per-genome bar chart, density plots, hits table
- **Parameters** -- Run configuration under a collapsible fold

---

## Metrics Definitions

### Genome-Level Metrics

These answer: "Was each genome detected?"

| Metric | Formula | Meaning |
|--------|---------|---------|
| **Sensitivity** | TP / (TP + FN) | Fraction of sample targets that were detected |
| **Specificity** | TN / (TN + FP) | Fraction of non-sample references correctly not detected |
| **Precision** | TP / (TP + FP) | Of detected references, fraction that are sample targets |
| **F1 Score** | 2 * (Precision * Sensitivity) / (Precision + Sensitivity) | Harmonic mean of precision and sensitivity |

Where:
- **TP** = sample targets detected (at least one read maps)
- **FN** = sample targets not detected
- **FP** = FP_target + FP_distractor (non-sample references incorrectly detected)
- **TN** = TN_target + TN_distractor (non-sample references correctly not detected)

### Read-Level Metrics

These track how fragments and reads flow through the pipeline:

| Metric | Description |
|--------|-------------|
| `sample_captured` | Fragments from sample targets that were captured by probes |
| `nonsample_target_captured` | Fragments from non-sample targets that were captured |
| `distractor_captured` | Fragments from distractors that were captured |
| `untargeted_captured` | Fragments from untargeted genomes that were captured (genome mode) |
| `reads_correctly_mapped` | Reads that map back to their source reference |
| `reads_incorrectly_mapped` | Reads that map to a different reference than their source |
| `reads_sequenced` | Number of reads after the sequencing step (may differ from captured if `--num-sequences` is used) |
| `reads_after_filter` | Number of reads after host filtering (0 if `--host-fasta` not provided) |
| `reads_mapped` | Total reads that mapped to any reference (= correctly + incorrectly mapped) |
| `reads_unmapped` | Reads that entered the mapping step but did not map to any reference |

In genome mode, a read from genome G mapping to target T is considered correctly mapped if T is linked to G in the sample-target-map.

Read source is determined from the fragment naming pattern `{seq_id}_fragment_{n}`, using the last `_fragment_` occurrence as the delimiter.

---

## Input File Formats

### FASTA Files

Standard FASTA format. Sequence IDs are the first whitespace-delimited word of the header:

```
>dengue_1 Dengue virus type 1
ATGCTAGCTAGCTAGC...
>zika_virus
GCTAGCTAGCTAGCTA...
```

**Requirements:**
- Sequence IDs must be unique within each file
- Sequence IDs must not contain spaces (use underscores)
- IDs must be consistent across input files (sample manifest IDs must match FASTA headers)

### Sample Manifest Format

The `--sample` flag accepts two formats:

**Inline IDs** (on the command line):

```bash
--sample id1 id2 id3
```

All IDs default to weight 1.0. A number following an ID sets that ID's weight:

```bash
--sample dengue_1 5 zika_virus chikungunya 0.5
# Result: dengue_1=5.0, zika_virus=1.0, chikungunya=0.5
```

**TSV file** (if a single argument that is an existing file):

```
# Optional comment lines starting with #
dengue_1	5.0
zika_virus
chikungunya	0.5
```

- First column: sequence ID (required)
- Second column: weight (optional, defaults to 1.0)
- Empty lines and lines starting with `#` are ignored

In standard mode, IDs must match target FASTA headers. In genome mode, IDs must match genome FASTA headers.

### Sample-Target Map Format

TSV file mapping genome IDs to target IDs (used with `--genomes` via `--sample-target-map`):

```
# genome_id	target_id
e_coli	e_coli_16S
e_coli	e_coli_gyrB
influenza_a	influenza_a
```

- One mapping per line: `genome_id<TAB>target_id`
- Multiple targets per genome supported (one line per mapping)
- Lines starting with `#` are ignored
- Empty lines are ignored

**Auto-linking:** When `--sample-target-map` is omitted (or for genomes not listed in the map), BaitBench auto-links genomes to targets by:

1. **Exact match**: genome ID equals a target ID (e.g., genome `influenza_a` → target `influenza_a`)
2. **Prefix match**: target ID starts with `{genome_id}|` (e.g., genome `Bartonella_grahamii` → targets `Bartonella_grahamii|ompB`, `Bartonella_grahamii|16S`)

This means you can name targets using the `organism|gene` convention and genomes using just `organism`, and they will auto-link without needing an explicit map file:

```
# genomes.fa
>Bartonella_grahamii
ATGC...
>Rickettsia_montanensis
ATGC...

# targets.fa (organism|gene naming)
>Bartonella_grahamii|ompB
ATGC...
>Rickettsia_montanensis|ompA
ATGC...
>Rickettsia_montanensis|gltA
ATGC...
```

With this naming, `Bartonella_grahamii` auto-links to `Bartonella_grahamii|ompB`, and `Rickettsia_montanensis` auto-links to both `Rickettsia_montanensis|ompA` and `Rickettsia_montanensis|gltA`.

**Using `--sample-target-map` for non-standard naming:** If your genome and target IDs don't follow either naming convention (exact match or `organism|gene`), provide an explicit mapping file:

```
# mapping.tsv — needed when genome IDs don't match target IDs
NC_012846.1	bartonella_ompB
NC_012846.1	bartonella_16S
GCF_000022725.1	rickettsia_gltA
```

Explicit mappings take precedence over auto-linking for the same genome ID.

**Untargeted genomes:** Sample genomes with no target mapping (explicit or auto-linked) become "untargeted" -- they generate fragments but have no expected target to detect. This models unknown organisms.

**Validation:** BaitBench errors if the map references genome or target IDs not found in their respective FASTA files.

### Groups File Format

TSV file mapping sequence IDs to group names. Used by `--groups` (target grouping) and `--distractor-groups` (distractor grouping):

```
# Optional comment lines starting with #
# seq_id	group_name
West_Nile_virus_0001	West_Nile_virus
West_Nile_virus_0002	West_Nile_virus
West_Nile_virus_0003	West_Nile_virus
Dengue_virus_1_6275	Dengue_virus_1
Dengue_virus_1_2274	Dengue_virus_1
Dengue_virus_2_8773	Dengue_virus_2
```

- One mapping per line: `seq_id<TAB>group_name`
- Lines starting with `#` are ignored
- Empty lines are ignored
- Leading `>` characters on sequence IDs are stripped automatically

**Target groups (`--groups`):**
- Sequence IDs must exist in the targets FASTA (BaitBench errors on unknown IDs)
- Sequences not listed in the file form singleton groups (their own ID = their group name)
- Without `--groups`, each target sequence is its own singleton group (no behavioral change)

**Distractor groups (`--distractor-groups`):**
- Overrides the default automatic grouping by FASTA file stem
- Sequence IDs must exist in the distractor sequences (BaitBench errors on unknown IDs)
- Without `--distractor-groups`, all contigs from each `--distractors` FASTA file are automatically grouped together using the file stem as the group name (e.g., all contigs in `Aaegypti.fa` → group `"Aaegypti"`)

---

## Dependencies

### Rust (build-time, managed by Cargo)

| Crate | Purpose |
|-------|---------|
| clap | CLI argument parsing (derive macros) |
| anyhow | Error handling |
| serde, serde_json | Serialization (JSON output) |
| rand, rand_distr | Random sampling and normal distribution |
| chrono | Timestamps |
| log, env_logger | Logging |

### External (runtime, installed via conda)

| Tool | Purpose | Required? |
|------|---------|-----------|
| minimap2 | Alignment (simulate, mapping, filtering) | Yes |
| BLAST+ | Cross-reactivity analysis (xreact) | Only if `baitbench xreact` is used |
| R + packages | HTML report generation | Only if reports are enabled |

Install all via:

```bash
conda env create -f environment.yml
conda activate baitbench
```
