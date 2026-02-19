# BaitBench - Claude Code Guide

## Project Overview

BaitBench is a generic tool for testing probe capture efficiency via in-silico simulation. Users provide probesets, target genomes, and distractor genomes to evaluate how well probes capture intended targets while avoiding off-target sequences.

## Architecture

BaitBench is a Rust CLI binary with R/ggplot2 for visualization.

### Pipeline Flow
```
targets.fa + distractors.fa
         ↓
   baitbench prepare  (combine, generate weights)
         ↓
   baitbench simulate (weighted random fragments)
         ↓
   baitbench capture  (minimap2 or BLAST)
         ↓
   baitbench filter   (optional host filtering)
         ↓
   baitbench map      (back to references)
         ↓
   baitbench list     (count reads per reference)
         ↓
   baitbench metrics  (TP/FP/FN/TN)
         ↓
   baitbench report   (HTML with ggplot2 figures)
```

`baitbench run` chains all steps automatically.

### Key Files

| File | Purpose |
|------|---------|
| `src/main.rs` | CLI entry point, clap dispatch |
| `src/cli.rs` | Subcommand and argument definitions |
| `src/commands/run.rs` | Full pipeline orchestrator |
| `src/commands/prepare.rs` | Combines targets/distractors, generates weights |
| `src/commands/simulate.rs` | Weighted random fragment generation |
| `src/commands/capture.rs` | minimap2 or BLAST probe capture |
| `src/commands/filter.rs` | Optional host read filtering |
| `src/commands/map_reads.rs` | Map reads back to reference |
| `src/commands/generate_list.rs` | SAM parsing → per-reference counts |
| `src/commands/metrics.rs` | TP/FP/FN/TN calculation, TSV/JSON output |
| `src/commands/report.rs` | Invokes Rscript for HTML report |
| `src/fasta/` | FASTA parsing, writing, extract-by-ID (replaces seqtk) |
| `src/alignment/paf.rs` | PAF format parser for minimap2 output |
| `src/alignment/sam.rs` | SAM format parser |
| `src/sampling/` | Weights calculation and fragment sampling |
| `src/external/` | minimap2, blastn, Rscript process wrappers |
| `R/report.Rmd` | RMarkdown template with ggplot2 figures |
| `R/report.R` | R script entry point for report generation |
| `environment.yml` | Conda environment (minimap2, blast, R packages) |

### Metrics Definitions

**Genome-level** (was each genome detected at all?):
- **TP (True Positive)**: Target genome detected
- **FP (False Positive)**: Distractor genome detected
- **FN (False Negative)**: Target genome NOT detected
- **TN (True Negative)**: Distractor genome NOT detected

**Read-level** (how reads flow through the pipeline):
- **target_captured**: Captured reads originating from target sequences
- **distractor_captured**: Captured reads originating from distractor sequences
- **reads_correctly_mapped**: Reads that map back to their source reference
- **reads_incorrectly_mapped**: Reads that map to a different reference (e.g., virus A read maps to virus B)

Read source is extracted from the fragment name pattern `{seq_id}_fragment_{n}` using the last occurrence of `_fragment_` as the delimiter.

## Important Rules

- **Always update documentation when making changes.** When adding, removing, or modifying CLI flags, subcommands, metrics, output files, or any user-facing behavior, update the relevant documentation files: `README.md`, this file (`CLAUDE.md`), and CLI help text in `src/cli.rs`.

## Development Guidelines

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

### R Scripts
- Located in `R/` directory
- `report.Rmd` is the parameterized RMarkdown template
- Uses ggplot2, dplyr, tidyr for visualization
- Called via `Rscript R/report.R --summary ... --detail ... --output ...`

### Testing Changes
```bash
# Build
cargo build --release

# Run with minimal example
./target/release/baitbench run \
  --targets examples/minimal/targets.fa \
  --distractors examples/minimal/distractors.fa \
  --probes examples/minimal/probes.fa \
  --num-reads 1000 \
  --seed 42 \
  --outdir test_results

# Check outputs
ls test_results/
cat test_results/results.tsv
```

### Common Modifications

**Adding a new metric**: Edit `src/commands/metrics.rs`, update `calculate_metrics()` and the TSV/JSON output.

**Adding a new figure**: Add to `R/report.Rmd` or create a new R script in `R/`.

**Changing capture parameters**: Pass via CLI flags (e.g., `--max-mismatches`, `--min-match-bases`).

**Adding a new subcommand**: Add to `src/cli.rs` (clap definition), create `src/commands/new_cmd.rs`, wire into `main.rs`.

**Modifying read generation**: Edit `src/sampling/fragment.rs`.

## Dependencies

### External (installed via conda)
- minimap2 (alignment)
- blastn (alternative capture)
- R + ggplot2 + rmarkdown (report generation, optional)

### Rust (managed by Cargo)
- clap (CLI), anyhow (errors), serde/serde_json (serialization)
- rand/rand_distr (sampling), chrono (timestamps)
- log/env_logger (logging)
