# BaitBench - Claude Code Guide

## Project Overview

BaitBench is a generic tool for testing probe capture efficiency via in-silico simulation. Users provide probesets, target genomes, and distractor genomes to evaluate how well probes capture intended targets while avoiding off-target sequences.

A key feature is the **sample manifest** (`--sample`), which specifies a subset of targets as the "sample" with optional weights, enabling testing of discrimination between viruses within the target panel.

**Genome mode** (`--genomes`) adds support for bacteria and other large pathogens where the sample genome differs from probe targets (e.g., full bacterial genome vs 16S gene target). Fragments are generated from full genomes, but reads are mapped back to targets for evaluation. An optional `--sample-target-map` links genome IDs to their corresponding target IDs.

## Architecture

BaitBench is a Rust CLI binary with R/ggplot2 for visualization.

### Pipeline Flow (Standard Mode)
```
targets.fa + distractors.fa [+ sample.tsv]
         ↓
   baitbench prepare   (combine, generate weights, write sample.txt)
         ↓
   baitbench simulate  (weighted random fragments → fragments.fa)
         ↓
   baitbench capture   (minimap2 or BLAST → captured.fa)
         ↓
   baitbench enrich    (optional fold enrichment → enriched.fa)
         ↓
   baitbench sequence  (optional sampling + trim to read length → reads.fa)
         ↓
   baitbench filter    (optional host filtering)
         ↓
   baitbench map       (back to combined_reference.fa)
         ↓
   baitbench list      (count reads per reference)
         ↓
   baitbench metrics   (3-way TP/FP/FN/TN)
         ↓
   baitbench report    (HTML with ggplot2 figures)
```

### Pipeline Flow (Genome Mode — `--genomes`)
```
genomes.fa + targets.fa + distractors.fa [+ sample.tsv] [+ mapping.tsv]
         ↓
   baitbench prepare   (build combined_reference.fa [genomes+distractors],
                         mapping_reference.fa [targets+distractors],
                         genomes.txt, sample_target_map.txt, weights)
         ↓
   baitbench simulate  (fragments from combined_reference.fa)
         ↓
   baitbench capture   (minimap2 or BLAST → captured.fa)
         ↓
   baitbench enrich    (uses genomes.txt to classify fragment sources)
         ↓
   baitbench sequence  (optional sampling + trim to read length → reads.fa)
         ↓
   baitbench filter    (optional host filtering)
         ↓
   baitbench map       (back to mapping_reference.fa — targets+distractors)
         ↓
   baitbench list      (count reads per reference)
         ↓
   baitbench metrics   (genome-aware classification with sample-target-map)
         ↓
   baitbench report    (HTML with ggplot2 figures)
```

`baitbench run` chains all steps automatically.

`baitbench ct-sweep` runs the pipeline at multiple CT values and produces coverage depth curve plots.

`baitbench xreact` checks probe cross-reactivity against genomes and/or other probes (standalone, not part of the pipeline).

### Key Files

| File | Purpose |
|------|---------|
| `src/main.rs` | CLI entry point, clap dispatch |
| `src/cli.rs` | Subcommand and argument definitions |
| `src/commands/run.rs` | Full pipeline orchestrator |
| `src/commands/prepare.rs` | Combines references, generates weights, writes ID lists; genome mode: two references + sample-target-map |
| `src/commands/simulate.rs` | Weighted random fragment generation |
| `src/commands/capture.rs` | minimap2 or BLAST probe capture |
| `src/commands/enrich.rs` | Fold enrichment adjustment (post-capture target:distractor ratio tuning) |
| `src/commands/sequence.rs` | Simulate sequencing (trim fragments to read length) |
| `src/commands/filter.rs` | Optional host read filtering |
| `src/commands/map_reads.rs` | Map reads back to reference |
| `src/commands/generate_list.rs` | SAM parsing → per-reference counts |
| `src/commands/metrics.rs` | 3-way classification (genome-aware with --sample-target-map), TSV/JSON output |
| `src/commands/report.rs` | Invokes Rscript for HTML report |
| `src/commands/xreact.rs` | Cross-reactivity analysis (probes vs genomes, probes vs probes) |
| `src/commands/ct_sweep.rs` | CT sweep: pipeline at multiple CT values → depth curves |
| `src/fasta/` | FASTA parsing, writing, extract-by-ID (replaces seqtk) |
| `src/alignment/paf.rs` | PAF format parser for minimap2 output |
| `src/alignment/sam.rs` | SAM format parser |
| `src/sampling/` | Weights calculation and fragment sampling |
| `src/io_utils.rs` | ID set parsing, sample manifest parsing, source ID extraction, sample-target-map I/O |
| `src/external/` | minimap2, blastn, Rscript process wrappers |
| `R/report.Rmd` | RMarkdown template with ggplot2 figures |
| `R/report.R` | R script entry point for report generation |
| `R/ct_sweep.R` | R script entry point for CT sweep report |
| `R/ct_sweep.Rmd` | RMarkdown template for coverage depth curves |
| `environment.yml` | Conda environment (minimap2, blast, R packages) |

### Metrics Definitions

**3-way genome-level classification** (was each genome detected at all?):

| Category | Detected | Classification |
|----------|----------|----------------|
| Sample target | Yes | TP |
| Sample target | No | FN |
| Non-sample target | Yes | FP_target |
| Non-sample target | No | TN_target |
| Distractor | Yes | FP_distractor |
| Distractor | No | TN_distractor |
| Untargeted genome | — | untargeted |

Without `--sample`, all targets are in the sample, reducing to the traditional 2-way classification.

In **genome mode**, the "sample targets" are derived from the sample-target-map: the union of all target IDs linked to sample genome IDs. Untargeted genomes (sample genomes with no target mapping) are tracked separately and do not participate in TP/FP/FN/TN classification.

**Read-level** (how reads flow through the pipeline):
- **sample_captured**: Captured fragments originating from sample target sequences (or sample genome sequences in genome mode)
- **nonsample_target_captured**: Captured fragments originating from non-sample target sequences
- **distractor_captured**: Captured fragments originating from distractor sequences
- **untargeted_captured**: Captured fragments originating from untargeted genome sequences (genome mode only)
- **reads_correctly_mapped**: Reads that map back to their source reference. In genome mode, a read from genome G mapping to target T is correct if T is in genome_to_targets[G]
- **reads_incorrectly_mapped**: Reads that map to a different reference (e.g., virus A read maps to virus B)

Read source is extracted from the fragment name pattern `{seq_id}_fragment_{n}` using the last occurrence of `_fragment_` as the delimiter.

### Weight Generation

**Standard mode:**
- Sample targets: use weight from sample manifest (default 1.0)
- Non-sample targets: weight = 0.0 (no reads generated)
- Distractors: `distractor_weight = (distractor_fraction * total_sample_weight) / (n_distractors * (1 - distractor_fraction))`
- Multiple distractor FASTA files are concatenated; all distractor sequences share the same per-sequence weight

**Genome mode:**
- Sample genomes: use weight from sample manifest (default 1.0)
- Non-sample genomes: weight = 0.0 (untargeted genomes not in sample also get 0)
- Distractors: same formula as standard mode
- Weights are assigned to genome IDs + distractor IDs (not target IDs)

### CT Score Support

Instead of `--distractor-fraction`, users can specify a qPCR CT (cycle threshold) score via `--ct`. The conversion is:

```
target_fraction = ct_baseline_fraction * 2^(ct_baseline - ct)
distractor_fraction = 1 - target_fraction
```

- `--ct` and `--distractor-fraction` are mutually exclusive (enforced by clap)
- `--ct-baseline` (default 20.0) and `--ct-baseline-fraction` (default 0.01) set the calibration point
- Default: CT 20 → 1% target (0.99 distractor), CT 25 → 0.03% target, CT 30 → 0.001% target
- If neither `--ct` nor `--distractor-fraction` is specified, defaults to distractor fraction 0.9

### Sample Manifest

TSV format: `id<tab>weight` (weight optional, defaults to 1.0). In standard mode, all IDs must exist in the targets FASTA. In genome mode, IDs refer to genome sequences. Without `--sample`, all targets (or all genomes in genome mode) are treated as sample with weight 1.0.

### Sample-Target Map

Optional TSV format: `genome_id<tab>target_id` (one mapping per line, `#` comments). Links genome IDs to their corresponding target IDs. Supports one-to-one, one-to-many, and many-to-one relationships. If omitted, auto-linking matches genome IDs to target IDs by (1) exact name match or (2) prefix match where a target ID starts with `{genome_id}|` (e.g., genome `Bartonella_grahamii` auto-links to `Bartonella_grahamii|ompB`). Genome IDs with no matching target become "untargeted" — they generate fragments but aren't expected to produce reads mapping to any target. Errors if map references IDs not found in genomes or targets FASTA.

## Important Rules

- **Read `ARCHITECTURE.md` first.** Before exploring the codebase, read `ARCHITECTURE.md` for a complete map of every source file, its public types/functions, data flow between modules, and intermediate file formats. This avoids redundant exploration.
- **Keep `ARCHITECTURE.md` up to date.** After adding, removing, or modifying source files, structs, public functions, intermediate files, or report sections, update `ARCHITECTURE.md` to reflect the changes.
- **Always update documentation when making changes.** When adding, removing, or modifying CLI flags, subcommands, metrics, output files, or any user-facing behavior, update the relevant documentation files: `README.md`, `ARCHITECTURE.md`, `MANUAL.md`, this file (`CLAUDE.md`), and CLI help text in `src/cli.rs`.

## Development Guidelines

### Environment Setup
```bash
conda activate baitbench
```
This ensures pandoc and other dependencies (R, minimap2, blast) are available for report generation.

Cargo is located at `/Users/niel/.cargo/bin/cargo` — ensure it's on `PATH` before building.

### Building
```bash
cargo build --release
```

### Rust Conventions
- Modules in `src/` follow a commands/library split
- Each command module exposes an `execute()` function taking an args struct
- External tools (minimap2, blastn) are called via `std::process::Command`
- FASTA operations are done natively in Rust (no seqtk dependency)
- Use `anyhow` for error handling, `log`/`env_logger` for logging
- Use `clap` derive macros for CLI argument definitions
- `--distractors` accepts multiple values via `num_args = 1..`

### R Scripts
- Located in `R/` directory
- `report.Rmd` is the parameterized RMarkdown template
- Uses ggplot2, dplyr, tidyr for visualization
- Called via `Rscript R/report.R --summary ... --detail ... --output ...`

### Testing Changes
```bash
# Build
cargo build --release

# Run with minimal example (all targets in sample)
./target/release/baitbench run \
  --targets examples/minimal/targets.fa \
  --distractors examples/minimal/distractors.fa \
  --probes examples/minimal/probes.fa \
  --num-fragments 1000 \
  --seed 42 \
  --report none \
  --outdir test_results

# Run with sample manifest (subset of targets)
echo "target_virus_1" > /tmp/sample.tsv
./target/release/baitbench run \
  --targets examples/minimal/targets.fa \
  --distractors examples/minimal/distractors.fa \
  --probes examples/minimal/probes.fa \
  --sample /tmp/sample.tsv \
  --num-fragments 1000 \
  --seed 42 \
  --report none \
  --outdir test_results_sample

# Run with genomes mode (bacteria + virus mix)
# genomes.fa contains full genomes; targets.fa contains probe target subsequences
# sample-target-map links genome IDs to their target IDs
./target/release/baitbench run \
  --targets targets.fa \
  --distractors distractors.fa \
  --genomes genomes.fa \
  --sample-target-map mapping.tsv \
  --sample genome_id_1 genome_id_2 \
  --probes probes.fa \
  --num-fragments 1000 \
  --seed 42 \
  --report none \
  --outdir test_results_genomes

# Check outputs
cat test_results/*/results.tsv
cat test_results_sample/*/results.tsv
cat test_results_genomes/*/detected_detail.tsv
```

### Common Modifications

**Adding a new metric**: Edit `src/commands/metrics.rs`, update `calculate_metrics()` and the TSV/JSON output.

**Adding a new figure**: Add to `R/report.Rmd` or create a new R script in `R/`.

**Changing capture parameters**: Pass via CLI flags (e.g., `--max-mismatches`, `--min-match-bases`).

**Adding a new subcommand**: Add to `src/cli.rs` (clap definition), create `src/commands/new_cmd.rs`, wire into `main.rs`.

**Modifying fragment generation**: Edit `src/sampling/fragment.rs`.

**Modifying sequencing**: Edit `src/commands/sequence.rs` (currently trims to read length; future: paired-end, error models, nanopore).

## Dependencies

### External (installed via conda)
- minimap2 (alignment)
- blastn (alternative capture)
- R + ggplot2 + rmarkdown (report generation, optional)

### Rust (managed by Cargo)
- clap (CLI), anyhow (errors), serde/serde_json (serialization)
- rand/rand_distr (sampling), chrono (timestamps)
- log/env_logger (logging)
