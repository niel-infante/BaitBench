# BaitBench Architecture

Quick reference for the codebase structure, module responsibilities, and key data types.

## Pipeline Flow

```
Step 1: prepare    → combined_reference.fa, weights.txt, targets.txt, distractors.txt, sample.txt
Step 2: simulate   → fragments.fa
Step 3: capture    → captured.fa
Step 3b: enrich    → enriched.fa (optional, if --fold-enrichment)
Step 4: sequence   → reads.fa
Step 5: filter     → filtered.fa (optional, if --host-fasta)
Step 6: map_reads  → mapped.sam
Step 7: list       → detected.list
Step 8: metrics    → results.tsv, detected_detail.tsv, results.json, coverage.tsv
Report: report     → report.html (optional, requires R)
```

`run.rs` orchestrates all steps. Each step is also available as a standalone subcommand.

## Source Layout

```
src/
├── main.rs              # Entry point: CLI parsing → command dispatch
├── cli.rs               # clap argument definitions (Commands enum, all flags)
├── io_utils.rs           # Shared helpers: parse_id_set, extract_source_id, parse_sample_manifest
├── alignment/
│   ├── coverage.rs      # CIGAR-based per-position coverage from SAM
│   ├── paf.rs           # PAF record filtering (mismatch/indel criteria)
│   └── sam.rs           # SAM parsing: read counts, mappings, mapped IDs
├── commands/
│   ├── run.rs           # Full pipeline orchestrator (steps 1-8 + report)
│   ├── prepare.rs       # Combine FASTAs, generate weights, write ID lists
│   ├── simulate.rs      # Generate weighted random fragments
│   ├── capture.rs       # Probe capture via minimap2 or BLAST
│   ├── enrich.rs        # Fold enrichment adjustment (post-capture ratio tuning)
│   ├── sequence.rs      # Trim fragments to read length
│   ├── filter.rs        # Remove host-mapping reads
│   ├── map_reads.rs     # Map reads to reference (minimap2)
│   ├── generate_list.rs # Count reads per reference from SAM → detected.list
│   ├── metrics.rs       # TP/FP/FN/TN classification, coverage stats, TSV/JSON output
│   ├── report.rs        # Invoke Rscript to render HTML report
│   └── probe_coverage.rs # Standalone probe tiling QC (maps probes to targets)
├── external/
│   ├── minimap2.rs      # minimap2 wrapper: capture_align (PAF), map_reads (SAM), host_align, probe_align
│   ├── blastn.rs        # BLAST+ wrapper: capture_align, filter_blast_results
│   └── rscript.rs       # Rscript discovery (BAITBENCH_R_DIR, binary walk, ./R/) and execution
├── fasta/
│   ├── reader.rs        # parse_fasta (id→seq), parse_fasta_ids, count_sequences
│   └── writer.rs        # write_fasta_record, extract_by_ids (streaming), concatenate_fastas
└── sampling/
    ├── fragment.rs      # generate_fragments: weighted sampling, normal-dist length
    └── weights.rs       # parse_weights, generate_weights (sample/distractor fraction)

R/
├── report.R             # CLI wrapper: parse args, invoke rmarkdown::render
├── report.Rmd           # RMarkdown template: tables, ggplot2 figures, coverage plots
├── probe_coverage.R     # CLI wrapper for probe coverage report
└── probe_coverage.Rmd   # RMarkdown template: probe tiling depth, gap analysis, proximity
```

## Key Data Types

### CLI (`cli.rs`)

- **`Commands`** enum — one variant per subcommand (Run, Prepare, Simulate, Capture, Enrich, Sequence, Filter, Map, List, Metrics, ProbeCoverage, Report), each with its own fields
- **`CaptureMethodArg`** — ValueEnum: Minimap2 | Blast

### Command Args Pattern

Every command module exports an `Args` struct and an `execute(&Args) -> Result<()>` function:

| Command | Args struct | Key inputs | Key outputs |
|---------|-------------|------------|-------------|
| `prepare` | `PrepareArgs` | targets, distractors, sample, distractor_fraction | combined_reference.fa, weights.txt, ID lists |
| `simulate` | `SimulateArgs` | reference, weights, num_fragments, seed, fragment_length_* | fragments.fa |
| `capture` | `CaptureArgs` | method, probes, fragments, max_mismatches, min_match_bases | captured.fa |
| `enrich` | `EnrichArgs` | captured, fragments, targets, distractors, fold_enrichment, seed | enriched.fa |
| `sequence` | `SequenceArgs` | input, read_length | reads.fa (trimmed) |
| `filter` | `FilterArgs` | host, reads, minimap_preset | filtered.fa |
| `map_reads` | `MapArgs` | reference, reads, minimap_preset | mapped.sam |
| `generate_list` | `ListArgs` | sam | detected.list |
| `metrics` | `MetricsArgs` | targets, distractors, sample, detected, fragments, captured, sam | results.tsv, detected_detail.tsv, results.json, coverage.tsv |
| `report` | `ReportArgs` | summary, detail, params, coverage, run_name | report.html |
| `probe_coverage` | `ProbeCoverageArgs` | targets, probes, minimap_preset, proximity | probe_depth.tsv, probe_coverage_summary.tsv, probe_coverage_report.html |
| `run` | `RunArgs` | all pipeline inputs | all of the above |

### Metrics (`metrics.rs`)

- **`MetricsResult`** — genome-level classification: tp/fn/fp_target/fp_distractor/tn_target/tn_distractor counts, sensitivity/specificity/precision/f1, plus ID lists (false_negatives, etc.)
- **`ReadLevelMetrics`** — fragment/read counts: sample_captured, nonsample_target_captured, distractor_captured, reads_correctly_mapped, reads_incorrectly_mapped
- **`DetailRow`** (Serialize) — per-reference row: reference_id, category, expected, detected, fragments_generated, fragments_captured, reads_assigned, classification, ref_length, avg_coverage, pct_covered_5x, pct_covered_20x
- **`JsonOutput`** (Serialize) — wraps RunInfo, CaptureStats, ReadLevelStats, JsonMetrics, JsonDetails

### Coverage (`alignment/coverage.rs`)

- **`CoverageResult`** — ref_lengths: HashMap<String, usize>, coverage: HashMap<String, Vec<u32>>
- **`CoverageStats`** — ref_length, avg_coverage, pct_covered_5x, pct_covered_20x (used by pipeline read coverage)
- **`ProbeCoverageStats`** — ref_length, covered_bases, pct_covered_1x, mean_depth, median_depth, pct_covered_2x/5x/10x, max_gap_length, num_gaps, pct_near_probe (used by probe-coverage QC)
- `compute_coverage(sam)` — pipeline read coverage (skips secondary alignments)
- `compute_probe_coverage(sam)` — probe tiling depth (includes secondary alignments)

### FASTA (`fasta/`)

- `parse_fasta(path) → HashMap<String, String>` (id → sequence, full load)
- `parse_fasta_ids(path) → Vec<String>` (IDs only, streaming)
- `extract_by_ids(fasta, ids, output) → usize` (streaming extraction)
- `concatenate_fastas(inputs, output)` (streaming concatenation)

### Sampling (`sampling/`)

- `generate_fragments(sequences, weights, num, output, seed, length_params) → usize`
- Fragment naming: `{seq_id}_fragment_{n} start={pos} length={len}`
- `generate_weights(target_ids, distractor_ids, sample_weights, distractor_fraction, output)`

### IO Utilities (`io_utils.rs`)

- `parse_id_set(path) → HashSet<String>` — one ID per line, # comments
- `extract_source_id(read_name) → Option<&str>` — extract `seq_id` from `{seq_id}_fragment_{n}...`
- `parse_sample_manifest(path) → HashMap<String, f64>` — TSV id\tweight

## External Tool Wrappers (`external/`)

All wrappers follow the pattern: `check_available() → bool/Result`, then specific invocation functions.

| Tool | Functions | Output format |
|------|-----------|---------------|
| minimap2 | `capture_align` (PAF), `map_reads` (SAM), `host_align` (SAM), `probe_align` (SAM, with secondary) | PAF or SAM |
| blastn | `capture_align` (TSV outfmt 6), `filter_blast_results` | TSV |
| rscript | `check_available`, `find_r_dir`, `run_rscript` | HTML |

## R Report (`R/`)

`report.R` accepts CLI args and calls `rmarkdown::render()` on `report.Rmd`.

### RMarkdown Parameters

| Param | Source |
|-------|--------|
| `summary_file` | results.tsv |
| `detail_file` | detected_detail.tsv |
| `params_file` | run_params.tsv |
| `coverage_file` | coverage.tsv (optional) |
| `run_name` | string |

### Report Sections

1. **Run parameters** — curated table from run_params.tsv
2. **Command** — reconstructed from run_params.tsv (auto-adapts to new params)
3. **Capture Summary** — generated vs captured bar chart, captured by source
4. **Detection Performance** — sensitivity/specificity/precision/F1 bar chart
5. **Read Mapping Accuracy** — correct vs incorrect mapped reads
6. **Confusion Matrix** — TP/FN/FP/TN heatmap
7. **Detection Detail** — per-reference table with coverage stats
8. **Detection Lollipop** — reads per detected reference, colored by classification
9. **Coverage** (conditional) — faceted overview + expandable per-reference detail plots

## Intermediate Files

| File | Format | Written by | Read by |
|------|--------|------------|---------|
| `combined_reference.fa` | FASTA | prepare | simulate, map_reads |
| `weights.txt` | TSV (id weight) | prepare | simulate |
| `targets.txt` | ID list | prepare | metrics |
| `distractors.txt` | ID list | prepare | metrics |
| `sample.txt` | ID list | prepare | metrics |
| `fragments.fa` | FASTA | simulate | capture, metrics |
| `captured.fa` | FASTA | capture | enrich (if --fold-enrichment), sequence, metrics |
| `enriched.fa` | FASTA | enrich | sequence, metrics (only if --fold-enrichment) |
| `reads.fa` | FASTA | sequence | filter/map_reads |
| `filtered.fa` | FASTA | filter | map_reads |
| `mapped.sam` | SAM | map_reads | generate_list, metrics |
| `detected.list` | TSV (id count) | generate_list | metrics |
| `run_params.tsv` | TSV (parameter flag value) | run | report |
| `results.tsv` | TSV | metrics | report |
| `detected_detail.tsv` | TSV | metrics | report |
| `results.json` | JSON | metrics | — |
| `coverage.tsv` | TSV (reference_id position depth) | metrics | report |
| `report.html` | HTML | report | — |

### Probe Coverage (standalone, not part of pipeline)

| File | Format | Written by | Read by |
|------|--------|------------|---------|
| `probe_alignment.sam` | SAM | probe_coverage | probe_coverage (then deleted) |
| `probe_depth.tsv` | TSV (reference_id position depth) | probe_coverage | probe_coverage report |
| `probe_coverage_summary.tsv` | TSV (per-target stats) | probe_coverage | probe_coverage report |
| `multi_mapping_probes.tsv` | TSV (probe_id num_targets targets) | probe_coverage | probe_coverage report |
| `probe_coverage_report.html` | HTML | probe_coverage (via R) | — |

## Key Conventions

- **Error handling**: `anyhow::Result` throughout; commands return `Result<()>`
- **Streaming I/O**: Large files (FASTA, SAM) processed via BufReader/BufWriter, not loaded fully
- **Fragment naming**: `{source_id}_fragment_{n} start={pos} length={len}` — source ID extracted via `io_utils::extract_source_id`
- **Sequence IDs**: First whitespace-delimited word of FASTA header (no spaces allowed in names)
- **Weights**: `explicit_weight * sequence_length` for sampling probability; weight 0 = no fragments
- **Capture filtering**: minimap2 → PAF → filter by mismatches/indels/match-bases; BLAST → outfmt 6 → filter by gaps/nident
- **Coverage**: Single-pass SAM parsing, CIGAR ops M/=/X increment depth, D/N advance position only
- **Report extensibility**: run_params.tsv drives command reconstruction — add new params there and the report picks them up automatically
