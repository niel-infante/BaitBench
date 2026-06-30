# BaitBench - Claude Code Guide

## Project Overview

BaitBench is a generic tool for testing probe capture efficiency via in-silico simulation. Users provide probesets, target genomes, and distractor genomes to evaluate how well probes capture intended targets while avoiding off-target sequences.

A key feature is the **sample manifest** (`--sample`), which specifies a subset of targets as the "sample" with optional weights, enabling testing of discrimination between viruses within the target panel.

**Genome mode** (`--genomes`) adds support for bacteria and other large pathogens where the sample genome differs from probe targets (e.g., full bacterial genome vs 16S gene target). Fragments are generated from full genomes, but reads are mapped back to targets for evaluation. An optional `--sample-target-map` links genome IDs to their corresponding target IDs.

## Architecture

BaitBench is a Rust CLI binary with R/ggplot2 for visualization.

### Pipeline Flow (Standard Mode)
```
targets.fa + distractors.fa + probes.fa [+ sample.tsv]
         ↓
   baitbench prepare   (combine, generate weights, write sample.txt)
         ↓
   baitbench simulate  (probes→ref align + TNN scoring + multinomial sampling → fragments.fa)
         ↓
   baitbench sequence  (optional sampling + trim to read length → reads.fa / reads.fq)
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
   baitbench simulate  (probe-biased fragments from combined_reference.fa, TNN-scored)
         ↓
   baitbench sequence  (optional sampling + trim to read length → reads.fa / reads.fq)
         ↓
   baitbench filter    (optional host filtering)
         ↓
   baitbench map       (back to mapping_reference.fa — targets+distractors)
         ↓
   baitbench list      (count reads per reference)
         ↓
   baitbench metrics   (genome-aware classification with sample-target-map)
         ↓
   baitbench identify  (optional: species-level calling from multi-target patterns)
         ↓
   baitbench report    (HTML with ggplot2 figures + optional species ID section)
```

`baitbench run` chains all steps automatically. With `--identify`, adds species-level calling.

`baitbench coverage-curve` runs the pipeline at multiple parameter combinations (CT × hybridization temperature × capture fraction × num-sequences) and produces coverage depth curve plots. Use `--ct-values`, `--hybridization-temperature-values`, `--capture-fraction-values`, and `--num-sequences-values` to sweep each dimension independently or in combination.

`baitbench xreact` checks probe cross-reactivity against genomes and/or other probes (standalone, not part of the pipeline).

`baitbench panel-qc` assesses whether a target panel can discriminate between species by computing target-vs-target similarity and per-species discriminability scores (standalone, pre-experiment QC).

`baitbench identify` calls species PRESENT/ABSENT/AMBIGUOUS from multi-target detection patterns, using cross-reactivity knowledge to explain away false positives (standalone or as pipeline step via `--identify`).

`baitbench build-probes` builds a probeset from target sequences: filter high-N targets, collapse redundant targets (cd-hit-est), filter sequences shorter than probe length, construct probes (`--method tile`: sliding window with configurable overlap via `--step`; `--method catch-lite`: native Rust reimplementation of CATCH optimization-based probe design (Metsky et al. 2019) with `--catch-probe-stride`, `--catch-mismatches`, `--catch-extension`, `--catch-coverage`, `--catch-minhash-threshold`; `--method catch`: external CATCH tool (requires `catch` conda package) with `--catch-probe-stride` and `--catch-mismatches`; `--method syotti-lite`: native Rust reimplementation of Syotti greedy set-cover bait design (Alanko et al. 2022) with `--syotti-mismatches` and `--syotti-seed-len`; `--method probetools-lite`: native Rust reimplementation of ProbeTools iterative k-mer clustering probe design (Kuchinski et al. 2022) with `--pt-step`, `--pt-identity`, `--pt-coverage`, `--pt-batch-size`, `--pt-max-panel-size`, `--pt-min-depth`, `--pt-max-iterations`, `--pt-min-coverage-gain`; requires cd-hit-est), fix N bases in probes (`--no-n-in-probes`: replaces each N with T unless adjacent to T, then A/C/G), filter by GC content, filter by sequence complexity (sDUST; Morgulis et al. 2006), deduplicate (cd-hit-est). Auto-chains into `assess-probes` unless `--skip-assess` is specified. Standalone, not part of the simulation pipeline.

`baitbench tool <TOOL>` groups standalone utility tools under a single subcommand to keep the main help clean. Run `baitbench tool --help` to list available tools. Current tools: `syotti` (greedy set-cover probe design), `catch` (CATCH optimization probe design), `dustview` (sDUST masking visualization), `collapse` (cd-hit-est sequence clustering). More tools may be added here without cluttering top-level help.


`baitbench assess-probes` runs combined probe assessment: probe coverage analysis + cross-reactivity (self-homology always, against genomes if `--genomes` provided), producing a single combined HTML report. Can include build pipeline stats when chained from `build-probes`. Standalone, not part of the simulation pipeline. `--all-individual-targets` reruns coverage once per target in isolation (minimap2 called separately for each target) and adds an **Individual Target Coverage** section to the report, useful when the panel contains many similar strain variants.

### Key Files

| File | Purpose |
|------|---------|
| `src/main.rs` | CLI entry point, clap dispatch |
| `src/cli.rs` | Subcommand and argument definitions |
| `src/commands/run.rs` | Full pipeline orchestrator |
| `src/commands/prepare.rs` | Combines references, generates weights, writes ID lists; genome mode: two references + sample-target-map |
| `src/thermodynamics.rs` | SantaLucia (1998) nearest-neighbor TNN model: `ThermoModel` struct (temp_c, na_conc_m), `delta_g()` (stacking + initiation terms + Owczarzy salt correction), `boltzmann_score()` |
| `src/commands/simulate.rs` | Thermodynamic/simple probe-biased fragment simulation (replaces simulate+capture+enrich) |
| `src/sampling/thermo_sim.rs` | ProbeHit, SimulateMode, load_probe_hits, sample_capture_fragments, sample_background_fragments, write_fragments |
| `src/commands/sequence.rs` | Simulate sequencing: dispatches `ReadSimulator` enum (perfect/art/badread); `art_modern.rs` handles read renaming via SAM RNAME field; `badread.rs` via the `{ref_name},{start}-{end}` field in read descriptions |
| `src/commands/filter.rs` | Optional host read filtering |
| `src/commands/map_reads.rs` | Map reads back to reference |
| `src/commands/generate_list.rs` | SAM parsing → per-reference counts |
| `src/commands/metrics.rs` | 3-way classification (genome-aware with --sample-target-map), TSV/JSON output |
| `src/commands/report.rs` | Invokes Rscript for HTML report |
| `src/commands/xreact.rs` | Cross-reactivity analysis (probes vs genomes, probes vs probes) |
| `src/commands/panel_qc.rs` | Target panel discriminability QC (target-vs-target similarity, species discrimination) |
| `src/commands/identify.rs` | Species-level calling from multi-target detection patterns |
| `src/target_similarity.rs` | Shared library: target similarity computation, discriminability scoring, confusion matrices |
| `src/commands/coverage_curve.rs` | Coverage curve: pipeline at multiple param combos (CT × temp × CF × NS) → depth curves |
| `src/commands/build_probes.rs` | Build probes: N filter → collapse → tile/CATCH/Syotti → N-fix (`--no-n-in-probes`) → GC filter → complexity filter (sDUST) → deduplicate; auto-chains to assess-probes |
| `src/commands/assess_probes.rs` | Combined probe assessment: orchestrates probe_coverage + xreact + optional individual-target coverage, generates combined report |
| `src/sdust.rs` | sDUST low-complexity sequence detection (Morgulis et al. 2006) |
| `src/syotti.rs` | Syotti greedy bait design: design_probes() — k-mer hash index, seed-and-extend, greedy set-cover (Alanko et al. 2022) |
| `src/catch.rs` | Native CATCH probe design: design_probes() — tiling → MinHash dedup → greedy set cover (reimplementation of Metsky et al. 2019) |
| `src/external/cdhit.rs` | cd-hit-est wrapper: check_available, cluster |
| `src/fasta/` | FASTA parsing, writing, extract-by-ID (replaces seqtk) |
| `src/alignment/paf.rs` | PAF format parser for minimap2 output |
| `src/alignment/sam.rs` | SAM format parser |
| `src/sampling/` | Weights calculation and fragment sampling |
| `src/cleanup.rs` | Post-pipeline cleanup: delete intermediate files/dirs, keep report inputs |
| `src/io_utils.rs` | `prefixed_join` helper, ID set parsing, sample manifest parsing, source ID extraction, sample-target-map I/O, groups file I/O |
| `src/external/` | minimap2, blastn, catch, Rscript process wrappers |
| `R/report.Rmd` | RMarkdown template with ggplot2 figures |
| `R/report.R` | R script entry point for report generation |
| `R/ct_sweep.R` | R script entry point for CT sweep report |
| `R/ct_sweep.Rmd` | RMarkdown template for coverage depth curves |
| `R/panel_qc.R` | R script entry point for panel QC report |
| `R/panel_qc.Rmd` | RMarkdown template for panel discriminability report |
| `R/build_probes.R` | R script entry point for build probes report |
| `R/build_probes.Rmd` | RMarkdown template for probe building pipeline stats |
| `R/assess_probes.R` | R script entry point for combined probe assessment report |
| `R/assess_probes.Rmd` | RMarkdown template: build stats (optional) + probe coverage (incl. individual target coverage section) + cross-reactivity |
| `environment.yml` | Conda environment (minimap2, blast, cd-hit, R packages) |

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

When `--groups` is provided, all classification is at **group level**: a group is detected if any of its member sequences has mapped reads. Multiple variant sequences of the same organism can be collapsed into one logical entity. Distractor contigs are always grouped by source FASTA file (e.g., all contigs in `Aaegypti.fa` → one `"Aaegypti"` entity); `--distractor-groups` overrides this with an explicit mapping.

**Read-level** (how reads flow through the pipeline):
- **sample_captured**: Captured fragments originating from sample target sequences (or sample genome sequences in genome mode)
- **nonsample_target_captured**: Captured fragments originating from non-sample target sequences
- **distractor_captured**: Captured fragments originating from distractor sequences
- **untargeted_captured**: Captured fragments originating from untargeted genome sequences (genome mode only)
- **reads_correctly_mapped**: Reads that map back to their source reference. In genome mode, a read from genome G mapping to target T is correct if T is in genome_to_targets[G]. With `--groups`, a read mapping to any sequence in the same group as the source counts as correctly mapped.
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
target_fraction = ct_baseline_fraction * (1 + efficiency)^(ct_baseline - ct)
distractor_fraction = 1 - target_fraction
```

- `--ct` and `--distractor-fraction` are mutually exclusive (enforced by clap)
- `--ct-baseline` (default 20.0) and `--ct-baseline-fraction` (default 0.01) set the one-point calibration
- `--ct-efficiency` (default 1.0 = 100%) sets the PCR amplification efficiency; typical real assays run at 0.90–0.98
- Default: CT 20 → 1% target (0.99 distractor), CT 25 → ~0.03% target, CT 30 → ~0.001% target (at 100% efficiency)
- If neither `--ct` nor `--distractor-fraction` is specified, defaults to distractor fraction 0.9

**Two-point calibration** (`--ct-calibration "CT1,FRAC1" "CT2,FRAC2"`): provide two (CT, target-fraction) reference points; the tool derives PCR efficiency automatically via `E = (f1/f2)^(1/(ct2-ct1)) - 1` and uses the first point as the baseline. Replaces `--ct-baseline`, `--ct-baseline-fraction`, and `--ct-efficiency` (all three are blocked when `--ct-calibration` is active). The derived efficiency is logged.

CT conversion logic lives in `src/main.rs`: `ct_to_distractor_fraction()`, `parse_calibration_point()`, `derive_calibration_params()`, `resolve_ct_params()`, `resolve_distractor_fraction()`.

### Sample Manifest

TSV format: `id<tab>weight` (weight optional, defaults to 1.0). In standard mode, all IDs must exist in the targets FASTA. In genome mode, IDs refer to genome sequences. Without `--sample`, all targets (or all genomes in genome mode) are treated as sample with weight 1.0.

### Sample-Target Map

Optional TSV format: `genome_id<tab>target_id` (one mapping per line, `#` comments). Links genome IDs to their corresponding target IDs. Supports one-to-one, one-to-many, and many-to-one relationships. If omitted, auto-linking matches genome IDs to target IDs by (1) exact name match or (2) prefix match where a target ID starts with `{genome_id}|` (e.g., genome `Bartonella_grahamii` auto-links to `Bartonella_grahamii|ompB`). Genome IDs with no matching target become "untargeted" — they generate fragments but aren't expected to produce reads mapping to any target. Errors if map references IDs not found in genomes or targets FASTA.

### Target Groups (`--groups`)

Optional TSV format: `seq_id<tab>group_name` (one mapping per line, `#` comments). Maps individual target sequence IDs to logical group names. When provided, all classification metrics (TP/FP/FN/TN) are computed at the group level: a group is detected if **any** of its member sequences has mapped reads.

**Use case**: Multiple variant sequences of the same organism (e.g., `West_Nile_virus_0001`, `_0002`, `_0003`) should be treated as one entity. Reads from any variant count as detecting the group.

Sequences not present in the groups file are treated as singleton groups (their own ID is the group name), so ungrouped sequences behave exactly as before.

### Distractor Groups (`--distractor-groups`)

Optional TSV format: `contig_id<tab>group_name`. By default (when omitted), all contigs from each `--distractors` FASTA file are automatically grouped by file stem (e.g., all contigs in `Aaegypti.fa` → group `"Aaegypti"`). This means FP_distractor / TN_distractor counts are at the file level, not per-contig.

Use `--distractor-groups` when multiple organisms are mixed into one FASTA file and you want per-organism distractor counting.

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

### Report Conventions
When adding or modifying any RMarkdown report (`R/*.Rmd`):

- **Scalability guards**: Adapt visualizations and tables based on data size. Use tiered strategies:
  - **Small** (≤20 items): Full detail — named axis labels, `knitr::kable()` tables, per-item plots
  - **Medium** (21–100 items): Compressed — `DT::datatable()` with pagination/filtering, smaller axis text, limit faceted plots
  - **Large** (>100 items): Distribution-based — histograms/boxplots instead of per-item bars, omit individual detail plots, downsample data for rendering speed
  - Reference `R/probe_coverage.Rmd` as the canonical example of this pattern
- **Parameters under fold**: Every report must include a `<details><summary>Parameters</summary>` section showing the parameters used to generate the report and a reconstructed CLI command. Pass a `run_params.tsv` (3-column: `parameter`, `flag`, `value`) from the Rust side. Use the data-driven command reconstruction pattern from `R/report.Rmd`.
- **Interactive tables**: Use `DT::datatable()` for any table that could exceed ~20 rows. Always set `scrollX = TRUE` and a reasonable `pageLength`.
- **Self-contained HTML**: Always use `self_contained: true` in the YAML front matter so reports are portable single-file HTML.

### Testing Changes
```bash
# Build
cargo build --release

# Run with minimal example — thermodynamic mode (default)
./target/release/baitbench run \
  --targets examples/minimal/targets.fa \
  --distractors examples/minimal/distractors.fa \
  --probes examples/minimal/probes.fa \
  --num-fragments 1000 \
  --capture-fraction 0.5 \
  --seed 42 \
  --report none \
  --outdir test_results

# Run with simple mode (no TNN scoring)
./target/release/baitbench run \
  --targets examples/minimal/targets.fa \
  --distractors examples/minimal/distractors.fa \
  --probes examples/minimal/probes.fa \
  --num-fragments 1000 \
  --simulate-mode simple \
  --seed 42 \
  --report none \
  --outdir test_results_simple

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
  --capture-fraction 0.6 \
  --hybridization-temperature 70 \
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

**Adding a new CLI flag (GUI)**: Add the corresponding input field in `gui/src/lib/views/RunView.svelte` — either in the main form or inside the `<AdvancedOptions>` block. The Rust GUI backend passes all args through as plain strings via a `HashMap<String, String>`, so no Rust changes to the GUI are needed. See the GUI section below for the full architecture.

**Modifying fragment generation**: Edit `src/sampling/fragment.rs`.

**Modifying sequencing**: Edit `src/commands/sequence.rs`. To add a new simulator: add a variant to `ReadSimulator` in `sequence.rs`, create a wrapper in `src/external/`, add profiles to `badread.rs` if needed, update `from_str()` and `execute()` dispatch, then add the new flag default in `main.rs`'s `Commands::Run` profile resolver.

**Changing paired-end behavior**: paired-end flows through `SequenceArgs.output_r2`, `FilterArgs.reads_r2`/`output_r2`, `MapArgs.reads_r2`, and `minimap2::map_reads`/`host_align`'s `reads_r2` parameter. All accept `Option<&Path>` — `None` is the default for non-PE mode.

## GUI (Tauri v2 + Svelte)

The `gui/` directory is a separate Tauri v2 + Svelte desktop app that wraps the CLI as a sidecar binary.

### GUI Development Workflow

```bash
cd gui

# First-time setup
make install          # npm install
make copy-sidecar     # cargo build --release + copy binary to src-tauri/binaries/

# Development (hot-reload frontend + Tauri window)
make dev

# Production build + package into dist-release/
make release
```

The sidecar binary must be rebuilt and copied any time the Rust CLI changes: `make copy-sidecar`.

### GUI Architecture

**View routing** (`src/App.svelte`): on mount, reads `conda_env_path` from the Tauri store (`settings.json`). If set → `RunView`; if missing → `SetupView`. After launching a pipeline → `LogView`.

**Three views:**
- `SetupView.svelte` — conda env picker; calls `detect_conda_envs` (auto-scans common paths) and `validate_conda_env` (checks for required binaries). Saves path to Tauri store on confirm.
- `RunView.svelte` — tool selector + form for all BaitBench subcommands. Builds a `PipelineConfig` object and calls `run_pipeline`. Add new CLI flags here only (no Rust changes needed).
- `LogView.svelte` — real-time log display via `listen('pipeline-log', ...)` events. Detects the report path in log output and calls `open_report` to open it in a new webview window.

**Shared state** (`src/lib/stores.ts`): `currentView`, `pipelineStatus`, `logLines`, `reportPath`, `condaEnvPath` — all Svelte writable stores.

**Tauri command layer** (`src-tauri/src/commands/`):
- `setup.rs` — `detect_conda_envs`, `validate_conda_env`
- `pipeline.rs` — `run_pipeline` (spawns sidecar with conda PATH prepended), `get_pipeline_status`, `cancel_pipeline`
- `report.rs` — `open_report` (opens HTML in a new webview window)

**Sidecar invocation**: `run_pipeline` prepends the conda env's bin directories to PATH (platform-aware: `bin/` on Unix, `Scripts/` + `Library/bin/` etc. on Windows), then spawns the `baitbench` sidecar with the subcommand and all flags as argv. stdout/stderr stream back as `pipeline-log` events.

### Adding a New Tauri Command

1. Add the `#[tauri::command]` function to the appropriate file in `src-tauri/src/commands/`
2. Register it in `src-tauri/src/lib.rs` inside `generate_handler![...]`
3. Call it from the frontend with `invoke('command_name', { args })` from `@tauri-apps/api/core`

### GUI Distribution

`.github/workflows/gui-release.yml` builds the GUI for macOS ARM64 (`macos-14`), macOS x64 (`macos-13`), and Windows (`windows-latest`) on push of a `v*` tag or manual dispatch. Each job compiles the CLI sidecar natively then runs `npm run tauri:build`. Artifacts (`.dmg`, `.msi`, `.exe`) are uploaded and attached to a GitHub Release.

## Dependencies

### External (installed via conda)
- minimap2 (alignment)
- blastn (alternative capture)
- cd-hit (sequence clustering, used by build-probes)
- R + ggplot2 + rmarkdown (report generation, optional)

### Rust (managed by Cargo)
- clap (CLI), anyhow (errors), serde/serde_json (serialization)
- rand/rand_distr (sampling), chrono (timestamps)
- log/env_logger (logging)
