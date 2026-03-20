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
  - [Fold Enrichment](#fold-enrichment)
  - [Weight Calculation](#weight-calculation)
  - [Sequence ID Conventions](#sequence-id-conventions)
- [Pipeline Flowcharts](#pipeline-flowcharts)
  - [Standard Mode Flowchart](#standard-mode-flowchart)
  - [Genome Mode Flowchart](#genome-mode-flowchart)
- [Commands](#commands)
  - [run](#run)
  - [prepare](#prepare)
  - [simulate](#simulate)
  - [capture](#capture)
  - [enrich](#enrich)
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
- [Parameter Reference](#parameter-reference)
  - [Input Files](#input-files)
  - [Fragment Generation](#fragment-generation)
  - [Target Abundance](#target-abundance)
  - [CT Score Parameters](#ct-score-parameters)
  - [Capture Parameters](#capture-parameters)
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
  - [results.json Structure](#resultsjson-structure)
  - [coverage.tsv Format](#coveragetsv-format)
- [Usage Examples](#usage-examples)
  - [Basic Probe Evaluation](#basic-probe-evaluation)
  - [Sample Discrimination Testing](#sample-discrimination-testing)
  - [Clinical Specimen Simulation with CT](#clinical-specimen-simulation-with-ct)
  - [Genome Mode for Bacteria](#genome-mode-for-bacteria)
  - [Mixed Panels (Virus + Bacteria)](#mixed-panels-virus--bacteria)
  - [Multiple Distractor Sources](#multiple-distractor-sources)
  - [Fold Enrichment Modeling](#fold-enrichment-modeling)
  - [Sequencing Depth Control](#sequencing-depth-control)
  - [Host Filtering](#host-filtering)
  - [Coverage Curve Analysis](#coverage-curve-analysis)
  - [Probe Design QC](#probe-design-qc)
  - [Cross-Reactivity Analysis](#cross-reactivity-analysis)
  - [Target Panel QC](#target-panel-qc)
  - [Species Identification](#species-identification)
  - [Running Individual Steps](#running-individual-steps)
  - [Reproducible Runs](#reproducible-runs)
  - [Batch Comparisons](#batch-comparisons)
- [Report Guide](#report-guide)
  - [HTML Report Sections](#html-report-sections)
  - [Coverage Curve Report](#coverage-curve-report)
  - [Probe Coverage Report](#probe-coverage-report)
  - [Panel QC Report](#panel-qc-report)
  - [Species Identification in Main Report](#species-identification-in-main-report)
- [Metrics Definitions](#metrics-definitions)
  - [Genome-Level Metrics](#genome-level-metrics)
  - [Read-Level Metrics](#read-level-metrics)
- [Input File Formats](#input-file-formats)
  - [FASTA Files](#fasta-files)
  - [Sample Manifest](#sample-manifest-format)
  - [Sample-Target Map](#sample-target-map-format)
- [Dependencies](#dependencies)

---

## Overview

BaitBench simulates a probe capture and sequencing workflow to evaluate how well a probe set performs. It answers questions like:

- Does the probe set capture all target sequences?
- Does it reject background (distractor) sequences?
- Can it discriminate between organisms within the target panel?
- How does performance change at different target abundances (CT values)?
- What sequencing depth is needed for adequate genome coverage?

The tool generates weighted random fragments from target and distractor sequences, simulates probe capture using minimap2 or BLAST, maps the resulting reads back to references, and computes detection and coverage metrics.

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
| minimap2 | >= 2.24 | Sequence alignment (capture, mapping, filtering) |
| BLAST+ | >= 2.12 | Alternative capture method |
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
2. **Simulate** -- Generate random DNA fragments proportional to weights
3. **Capture** -- Align fragments against probes; keep those that hybridize
4. **Enrich** (optional) -- Adjust target:distractor ratio to model fold enrichment
5. **Sequence** -- Trim fragments to read length; optionally sample to model sequencing depth
6. **Filter** (optional) -- Remove reads mapping to a host genome
7. **Map** -- Align reads back to reference sequences
8. **List** -- Count reads per reference
9. **Metrics** -- Classify each reference as TP/FP/FN/TN; compute summary statistics
10. **Report** (optional) -- Generate HTML report with figures

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

### CT Scores

CT (cycle threshold) scores from qPCR provide a natural way to express target abundance. BaitBench converts CT values to distractor fractions using a calibrated exponential formula. Lower CT = more target DNA = easier to detect.

See [CT Score Calculation](#ct-score-calculation) for the formula, default calibration, and how to customize it.

### Fold Enrichment

Post-capture enrichment adjusts the ratio of target to distractor fragments after the capture step. A fold enrichment of 100 means the target:distractor ratio is 100x higher after enrichment than before capture.

When `--fold-enrichment` is specified, BaitBench adjusts the captured fragment pool by either subsampling captured distractors or adding back uncaptured distractors to achieve the requested enrichment ratio.

When omitted, the capture step operates as binary (captured or not) with no ratio adjustment.

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
                        ├─ 2. SIMULATE ────────────────── fragments.fa
weights.txt ────────────┤      │
--num-fragments ────────┤      │
--fragment-length-* ────┤      │
--seed ─────────────────┘      │
                               │
                               ▼
fragments.fa ───────────┐
                        ├─ 3. CAPTURE ─────────────────── captured.fa
probes.fa ──────────────┤      │
--capture-method ───────┤      │
--min-match-bases ──────┤      │
--max-mismatches ───────┘      │
                               │
                      ┌────────┴────────┐
                      │ --fold-enrichment│
                      │   specified?     │
                      └──┬───────────┬───┘
                     yes │           │ no
                         ▼           │
captured.fa ─────┐                   │
fragments.fa ────┤                   │
targets.txt ─────┤ 3b. ENRICH       │
distractors.txt ─┤    │             │
--fold-enrichment┤    │             │
--seed ──────────┘    │             │
                      ▼             │
               enriched.fa          │
                      │             │
                      ▼             ▼
                 (enriched or captured).fa
                      │
                      ├─ 4. SEQUENCE ──────────────────── reads.fa
--read-length ────────┤      │
--num-sequences ──────┤      │
--seed ───────────────┘      │
                             │
                    ┌────────┴────────┐
                    │  --host-fasta   │
                    │   specified?    │
                    └──┬──────────┬───┘
                   yes │          │ no
                       ▼          │
reads.fa ───────┐                 │
host.fa ────────┤ 5. FILTER       │
--host-minimap- ┤    │            │
  preset ───────┘    │            │
                     ▼            │
              filtered.fa         │
                     │            │
                     ▼            ▼
              (filtered or reads).fa
                     │
combined_            ├─ 6. MAP ────────────────────────── mapped.sam
  reference.fa ──────┤      │
--minimap-preset ────┘      │
                            │
                            ▼
mapped.sam ──────────── 7. LIST ───────────────────────── detected.list
                            │
                            ▼
targets.txt ─────────┐
distractors.txt ─────┤
sample.txt ──────────┤
detected.list ───────┤ 8. METRICS ────────────────────── results.tsv
fragments.fa ────────┤                                    detected_detail.tsv
captured.fa ─────────┤                                    results.json
mapped.sam ──────────┘                                    coverage.tsv
                            │
                            ▼
results.tsv ─────────┐
detected_detail.tsv ─┤ 9. REPORT (optional) ──────────── report.html
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

    Steps 2-5 are identical to standard mode, except:
      - Simulate uses combined_reference.fa (genomes + distractors)
      - Enrich uses genomes.txt to classify fragment sources

                     ... (steps 2-5) ...

              (filtered or reads).fa
                     │
mapping_             ├─ 6. MAP ────────────────────────── mapped.sam
  reference.fa ──────┤       (targets + distractors)
                     │
                            │
                            ▼
                     ... (step 7 same) ...
                            │
                            ▼
targets.txt ─────────┐
distractors.txt ─────┤
sample.txt ──────────┤
sample_target_map ───┤ 8. METRICS ────────────────────── results.tsv
detected.list ───────┤   (genome-aware classification)    detected_detail.tsv
fragments.fa ────────┤                                    results.json
captured.fa ─────────┤                                    coverage.tsv
mapped.sam ──────────┘

                     ... (step 9 same) ...
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

### capture

Aligns fragments against probes to simulate hybridization capture.

```bash
baitbench capture \
  --probes probes.fa \
  --fragments fragments.fa \
  --output captured.fa \
  [--method minimap2|blast] \
  [--min-match-bases 60] \
  [--max-mismatches 10] \
  [--blast-db path] \
  [--threads 1]
```

minimap2 mode: Aligns fragments to probes using minimap2, then filters the PAF output by matching bases, mismatches, and indels.

BLAST mode: Runs blastn, then filters by identity and gaps. Requires `--blast-db` pointing to a pre-built BLAST database.

**Output files:**
- `captured.fa` -- fragments that pass the capture filter

### enrich

Adjusts the post-capture fragment pool to achieve a target fold enrichment.

```bash
baitbench enrich \
  --captured captured.fa \
  --fragments fragments.fa \
  --targets targets.txt \
  --distractors distractors.txt \
  --fold-enrichment 100 \
  --output enriched.fa \
  [--seed 42]
```

Fold enrichment is the ratio of target:distractor proportions post-enrichment vs pre-capture. A value of 100 means the target:distractor ratio improves 100-fold after enrichment.

**Output files:**
- `enriched.fa` -- enriched fragment pool

### sequence

Simulates sequencing by trimming fragments to read length.

```bash
baitbench sequence \
  --input captured.fa \
  --output reads.fa \
  [--read-length 120] \
  [--num-sequences 5000] \
  [--seed 42]
```

Fragments shorter than `--read-length` are kept as-is. With `--num-sequences`, reads are sampled with replacement from the fragment pool (modeling PCR amplification before sequencing) and given unique IDs.

**Output files:**
- `reads.fa` -- sequencing reads

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
  --captured captured.fa \
  --sam mapped.sam \
  --run-name "my_run" \
  --num-fragments 10000 \
  --output-summary results.tsv \
  --output-detail detected_detail.tsv \
  [--output-json results.json] \
  [--output-coverage coverage.tsv] \
  [--sample-target-map sample_target_map.txt] \
  [--seed 42]
```

**Output files:**
- `results.tsv` -- genome-level summary metrics
- `detected_detail.tsv` -- per-reference detection and coverage detail
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
  [--fold-enrichment-values 10 100 | --fold-enrichment 100] \
  [--num-sequences-values 100 500 | --num-sequences 500] \
  [--outdir coverage_curve_results] \
  [--cleanup] \
  [... other pipeline parameters ...]
```

Three parameters can be swept (each has a singular fixed form and a plural sweep form):

| Sweep flag | Fixed flag | Description |
|-----------|------------|-------------|
| `--ct-values 20 25 30` | `--ct 25` | CT values |
| `--fold-enrichment-values 10 100` | `--fold-enrichment 100` | Fold enrichment values |
| `--num-sequences-values 100 500` | `--num-sequences 500` | Number of sequences to sample |

Sweep and fixed forms of the same parameter are mutually exclusive. `--ct-values` and `--distractor-fraction` are also mutually exclusive.

`--sample` is **required** for coverage-curve (must specify which targets to track).

The pipeline shares intermediate files across combinations for efficiency: prepare/simulate/capture are shared per CT value; enrich is shared per CT x fold-enrichment combination.

**Output files:**
- Combo subdirectories named by swept params (e.g., `ct_20/`, `ct_20_fe_100/`, `ct_20_fe_100_ns_500/`)
- `coverage_curve_depth_curves.tsv` -- aggregated depth data
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
| CT baseline | `--ct-baseline` | 20.0 | The CT value at which the target fraction equals the baseline fraction |
| CT baseline fraction | `--ct-baseline-fraction` | 0.01 | The target fraction at the baseline CT value |

See [CT Score Calculation](#ct-score-calculation) for details.

### Capture Parameters

| Parameter | Flag | Default | Description |
|-----------|------|---------|-------------|
| Capture method | `--capture-method` | minimap2 | Alignment tool for capture: `minimap2` or `blast` |
| Min match bases | `--min-match-bases` | 60 | Minimum number of matching bases for a fragment to be considered captured |
| Max mismatches | `--max-mismatches` | 10 | Maximum mismatches allowed (minimap2 only) |
| BLAST database | `--blast-db` | none | Path to pre-built BLAST database (required if `--capture-method blast`) |
| Fold enrichment | `--fold-enrichment` | none | Post-capture enrichment factor (>= 1.0). Omit for binary capture |

### Sequencing Parameters

| Parameter | Flag | Default | Description |
|-----------|------|---------|-------------|
| Read length | `--read-length` | 120 | Trim captured fragments to this length (bp). Fragments shorter than this are kept as-is |
| Num sequences | `--num-sequences` | all | Number of reads to sample with replacement. If not set, all captured fragments become reads. Models sequencing depth control |

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
| Minimap preset | `--minimap-preset` | sr | Minimap2 preset for read mapping |
| Host minimap preset | `--host-minimap-preset` | sr | Minimap2 preset for host read filtering |
| Cleanup | `--cleanup` | false | Delete intermediate files after completion, keeping only report inputs and final outputs. Available on `run`, `coverage-curve`, `probe-coverage`, and `xreact` |
| Identify | `--identify` | false | Enable species-level identification after metrics (genome mode only, requires `--sample-target-map`). Available on `run` |
| Identity threshold | `--identity-threshold` | 90.0 | Minimum sequence identity % to consider targets "similar" for species identification. Available on `run`, `panel-qc`, `identify` |
| Min unique targets | `--min-unique-targets` | 1 | Minimum unique target detections required to call a species PRESENT. Available on `run`, `identify` |

---

## CT Score Calculation

CT (cycle threshold) scores from qPCR provide an intuitive way to express target abundance. In qPCR, each cycle doubles the DNA, so a CT difference of 1 corresponds to a 2-fold change in DNA quantity. Lower CT = more target DNA.

### The Formula

```
target_fraction = ct_baseline_fraction * 2^(ct_baseline - ct)
distractor_fraction = 1 - target_fraction
```

Where:
- `ct_baseline` is a known CT value (default: 20.0)
- `ct_baseline_fraction` is the target fraction at that CT (default: 0.01)
- `ct` is the CT value you want to simulate

### Default Calibration

With defaults (`--ct-baseline 20.0`, `--ct-baseline-fraction 0.01`), the interpretation is: "at CT 20, 1% of DNA is from targets."

### CT Reference Table

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

### Custom Calibration

The default calibration assumes CT 20 = 1% target. If your experimental system has different characteristics, use `--ct-baseline` and `--ct-baseline-fraction` to calibrate:

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
- CT 30 = 0.003% target (5 CT higher = 32x less)
- CT 20 = 3.2% target (5 CT lower = 32x more)

**Example: Calibrate with a strong-positive sample (CT 15 = 50% target):**

```bash
baitbench run \
  --targets targets.fa \
  --distractors distractors.fa \
  --probes probes.fa \
  --ct 25 \
  --ct-baseline 15 \
  --ct-baseline-fraction 0.5 \
  --num-fragments 10000 \
  --outdir results
```

### Tips for Using CT Scores

- **Match your experimental system.** If you have empirical data linking CT to target fraction, use `--ct-baseline` and `--ct-baseline-fraction` to match your curve.
- **Use coverage-curve to sweep.** The `coverage-curve` command with `--ct-values` lets you visualize performance across a range of CT values in a single analysis.
- **Remember the log scale.** Each CT unit represents a 2-fold change. A 10-CT range spans ~1000-fold differences in abundance.

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
├── fragments.fa                # Simulated DNA fragments
├── captured.fa                 # Fragments passing capture filter
├── enriched.fa                 # Post-enrichment fragments (if --fold-enrichment)
├── reads.fa                    # Sequencing reads (trimmed to read length)
├── filtered.fa                 # Host-filtered reads (if --host-fasta)
├── mapped.sam                  # Read alignments to reference
├── detected.list               # Read counts per reference
├── run_params.tsv              # Run configuration (used by report)
├── results.tsv                 # Summary metrics
├── detected_detail.tsv         # Per-reference detection and coverage detail
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

All distractor sequences are concatenated and share the same per-sequence weight.

### Fold Enrichment Modeling

Model post-capture enrichment:

```bash
# Binary capture (no enrichment adjustment)
baitbench run \
  --targets targets.fa \
  --distractors distractors.fa \
  --probes probes.fa \
  --num-fragments 10000 \
  --outdir results_binary

# 100x fold enrichment
baitbench run \
  --targets targets.fa \
  --distractors distractors.fa \
  --probes probes.fa \
  --fold-enrichment 100 \
  --num-fragments 10000 \
  --outdir results_enriched
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

# Sweep CT and fold enrichment (combinatorial)
baitbench coverage-curve \
  --targets targets.fa \
  --distractors distractors.fa \
  --probes probes.fa \
  --sample dengue_1 \
  --ct-values 20 25 30 \
  --fold-enrichment-values 10 100 \
  --num-fragments 10000 \
  --outdir coverage_ct_fe

# Sweep all three parameters
baitbench coverage-curve \
  --targets targets.fa \
  --distractors distractors.fa \
  --probes probes.fa \
  --sample dengue_1 \
  --ct-values 20 25 30 \
  --fold-enrichment-values 10 100 \
  --num-sequences-values 500 1000 5000 \
  --num-fragments 10000 \
  --outdir coverage_full

# Fixed CT with fold enrichment sweep
baitbench coverage-curve \
  --targets targets.fa \
  --distractors distractors.fa \
  --probes probes.fa \
  --sample dengue_1 \
  --ct 25 \
  --fold-enrichment-values 1 10 100 1000 \
  --num-fragments 10000 \
  --outdir coverage_fe
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

### Running Individual Steps

Run pipeline steps independently for custom workflows:

```bash
# 1. Prepare
baitbench prepare \
  --targets targets.fa \
  --distractors distractors.fa \
  --distractor-fraction 0.95 \
  --outdir prep

# 2. Simulate
baitbench simulate \
  --reference prep/combined_reference.fa \
  --weights prep/weights.txt \
  --num-fragments 50000 \
  --seed 42 \
  --output prep/fragments.fa

# 3. Capture
baitbench capture \
  --probes probes.fa \
  --fragments prep/fragments.fa \
  --min-match-bases 80 \
  --output prep/captured.fa

# 4. Sequence
baitbench sequence \
  --input prep/captured.fa \
  --read-length 150 \
  --output prep/reads.fa

# 5. Map
baitbench map \
  --reference prep/combined_reference.fa \
  --reads prep/reads.fa \
  --output prep/mapped.sam

# 6. List
baitbench list \
  --sam prep/mapped.sam \
  --output prep/detected.list

# 7. Metrics
baitbench metrics \
  --targets prep/targets.txt \
  --distractors prep/distractors.txt \
  --sample prep/sample.txt \
  --detected prep/detected.list \
  --fragments prep/fragments.fa \
  --captured prep/captured.fa \
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
| minimap2 | Alignment (capture, mapping, filtering) | Yes (for default capture method) |
| BLAST+ | Alternative capture method | Only if `--capture-method blast` |
| R + packages | HTML report generation | Only if reports are enabled |

Install all via:

```bash
conda env create -f environment.yml
conda activate baitbench
```
