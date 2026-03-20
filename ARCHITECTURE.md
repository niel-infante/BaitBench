# BaitBench Architecture

Quick reference for the codebase structure, module responsibilities, and key data types.

## Pipeline Flow

### Standard Mode
```
Step 1: prepare    → combined_reference.fa, weights.txt, targets.txt, distractors.txt, sample.txt
Step 2: simulate   → fragments.fa
Step 3: capture    → captured.fa
Step 3b: enrich    → enriched.fa (optional, if --fold-enrichment)
Step 4: sequence   → reads.fa
Step 5: filter     → filtered.fa (optional, if --host-fasta)
Step 6: map_reads  → mapped.sam  (against combined_reference.fa)
Step 7: list       → detected.list
Step 8: metrics    → results.tsv, detected_detail.tsv, results.json, coverage.tsv
Report: report     → report.html (optional, requires R)
```

### Genome Mode (`--genomes`)
```
Step 1: prepare    → combined_reference.fa (genomes+distractors), mapping_reference.fa (targets+distractors),
                     weights.txt, targets.txt, distractors.txt, genomes.txt, sample.txt, sample_target_map.txt
Step 2: simulate   → fragments.fa  (from combined_reference.fa — genomes+distractors)
Step 3: capture    → captured.fa
Step 3b: enrich    → enriched.fa  (uses genomes.txt to classify fragment sources)
Step 4: sequence   → reads.fa
Step 5: filter     → filtered.fa
Step 6: map_reads  → mapped.sam  (against mapping_reference.fa — targets+distractors)
Step 7: list       → detected.list
Step 8: metrics    → results.tsv, detected_detail.tsv, results.json, coverage.tsv  (genome-aware classification)
Step 9: identify   → species_calls.tsv, species_calls.json  (optional, if --identify)
Report: report     → report.html  (includes Species Identification section if species_calls.tsv exists)
```

`run.rs` orchestrates all steps. Each step is also available as a standalone subcommand.

## Source Layout

```
src/
├── main.rs              # Entry point: CLI parsing → command dispatch
├── cli.rs               # clap argument definitions (Commands enum, all flags)
├── cleanup.rs           # Post-pipeline cleanup: delete intermediate files/dirs, keep report inputs
├── io_utils.rs           # Shared helpers: prefixed_join, parse_id_set, extract_source_id, parse_sample_manifest, parse_sample_target_map
├── target_similarity.rs  # Shared library: target-vs-target similarity computation, discriminability scoring, confusion matrices
├── alignment/
│   ├── coverage.rs      # CIGAR-based per-position coverage from SAM
│   ├── paf.rs           # PAF record filtering (mismatch/indel criteria) + structured PafRecord parsing
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
│   ├── report.rs        # Report generation: HTML (Rscript), RMarkdown (template substitution), or skip; shared substitute_rmd_params utility
│   ├── probe_coverage.rs # Standalone probe tiling QC (maps probes to targets)
│   ├── xreact.rs        # Standalone cross-reactivity analysis (probes vs genomes, probes vs probes)
│   ├── panel_qc.rs      # Standalone target panel discriminability QC (target-vs-target similarity, species discrimination)
│   ├── identify.rs      # Species-level calling from multi-target detection patterns (standalone or pipeline step)
│   └── coverage_curve.rs # Coverage curve: pipeline at multiple param combos → depth curves
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
├── probe_coverage.Rmd   # RMarkdown template: probe tiling depth, gap analysis, proximity
├── coverage_curve.R     # CLI wrapper for coverage curve report
├── coverage_curve.Rmd   # RMarkdown template: coverage depth curves (multi-param sweep)
├── xreact.R             # CLI wrapper for cross-reactivity report
├── xreact.Rmd           # RMarkdown template: plotly heatmaps, density plots, DT hit tables
├── panel_qc.R           # CLI wrapper for panel QC report
└── panel_qc.Rmd         # RMarkdown template: discriminability charts, confusion heatmaps, target tables
```

## Key Data Types

### CLI (`cli.rs`)

- **`Commands`** enum — one variant per subcommand (Run, Prepare, Simulate, Capture, Enrich, Sequence, Filter, Map, List, Metrics, ProbeCoverage, Xreact, PanelQc, Identify, Report, CoverageCurve), each with its own fields
- **`CaptureMethodArg`** — ValueEnum: Minimap2 | Blast
- **`ReportMode`** — ValueEnum: Full | None | Rmd — controls report output (HTML, skip, or editable RMarkdown)
- **CT score flags** — `--ct`, `--ct-baseline`, `--ct-baseline-fraction` on Run and Prepare; `--ct` conflicts with `--distractor-fraction`
- **Genome mode flags** — `--genomes` (optional genome FASTA for fragment generation), `--sample-target-map` (optional genome-to-target mapping TSV) on Run, Prepare, and CoverageCurve
- **`--output-prefix`** — string prepended to every auto-generated output filename; available on Run, Prepare, ProbeCoverage, CoverageCurve, Xreact, PanelQc, Identify (default: empty string)

### Command Args Pattern

Every command module exports an `Args` struct and an `execute(&Args) -> Result<()>` function:

| Command | Args struct | Key inputs | Key outputs |
|---------|-------------|------------|-------------|
| `prepare` | `PrepareArgs` | targets, distractors, sample, distractor_fraction, genomes, sample_target_map | combined_reference.fa, weights.txt, ID lists; genome mode also: mapping_reference.fa, genomes.txt, sample_target_map.txt |
| `simulate` | `SimulateArgs` | reference, weights, num_fragments, seed, fragment_length_* | fragments.fa |
| `capture` | `CaptureArgs` | method, probes, fragments, max_mismatches, min_match_bases | captured.fa |
| `enrich` | `EnrichArgs` | captured, fragments, targets, distractors, fold_enrichment, seed | enriched.fa |
| `sequence` | `SequenceArgs` | input, read_length, num_sequences, seed | reads.fa (trimmed, optionally sampled) |
| `filter` | `FilterArgs` | host, reads, minimap_preset | filtered.fa |
| `map_reads` | `MapArgs` | reference, reads, minimap_preset | mapped.sam |
| `generate_list` | `ListArgs` | sam | detected.list |
| `metrics` | `MetricsArgs` | targets, distractors, sample, detected, fragments, captured, sam, sample_target_map, reads_sequenced, reads_after_filter | results.tsv, detected_detail.tsv, results.json, coverage.tsv |
| `report` | `ReportArgs` | summary, detail, params, coverage, run_name, report (ReportMode) | report.html or report.Rmd |
| `probe_coverage` | `ProbeCoverageArgs` | targets, probes, minimap_preset, proximity | probe_depth.tsv, probe_coverage_summary.tsv, probe_coverage_report.html |
| `xreact` | `XreactArgs` | probes, against (genome FASTAs), self_mode, threshold, minimap_preset, report (ReportMode) | hits.tsv, summary.tsv, xreact_report.html |
| `panel_qc` | `PanelQcArgs` | targets, sample_target_map, identity_threshold, minimap_preset, report (ReportMode) | target_similarity.tsv, species_discriminability.tsv, species_confusion_matrix.tsv, panel_qc_report.html |
| `identify` | `IdentifyArgs` | detected_detail, sample_target_map, target_similarity (or targets for on-the-fly), identity_threshold, min_unique_targets | species_calls.tsv, species_calls.json |
| `run` | `RunArgs` | all pipeline inputs + ct, ct_baseline, ct_baseline_fraction, num_sequences, genomes, sample_target_map, identify, identity_threshold, min_unique_targets | all of the above |
| `coverage_curve` | `CoverageCurveArgs` | targets, distractors, probes, sample (required), ct/fe/ns values (sweep or fixed), all pipeline params, genomes, sample_target_map | coverage_curve_depth_curves.tsv, coverage_curve_report.html, combo subdirs |

### Metrics (`metrics.rs`)

- **`MetricsResult`** — genome-level classification: tp/fn/fp_target/fp_distractor/tn_target/tn_distractor counts, sensitivity/specificity/precision/f1, plus ID lists (false_negatives, etc.), untargeted_genomes list
- **`ReadLevelMetrics`** — fragment/read counts: sample_captured, nonsample_target_captured, distractor_captured, untargeted_captured, reads_correctly_mapped, reads_incorrectly_mapped. Summary TSV also includes derived pipeline flow counts: reads_sequenced, reads_after_filter, reads_mapped, reads_unmapped
- **`GenomeContext`** — derived from sample_target_map: sample_targets, genome_ids, sample_genome_ids, genome_to_targets, target_to_genomes, untargeted_genomes (genomes with no target mapping)
- **`DetailRow`** (Serialize) — per-reference row: reference_id, category, expected, detected, fragments_generated, fragments_captured, reads_assigned, classification, ref_length, avg_coverage, pct_covered_5x, pct_covered_20x
- **`JsonOutput`** (Serialize) — wraps RunInfo, CaptureStats, ReadLevelStats, JsonMetrics, JsonDetails

### Coverage (`alignment/coverage.rs`)

- **`CoverageResult`** — ref_lengths: HashMap<String, usize>, coverage: HashMap<String, Vec<u32>>
- **`CoverageStats`** — ref_length, avg_coverage, pct_covered_5x, pct_covered_20x (used by pipeline read coverage)
- **`ProbeCoverageStats`** — ref_length, covered_bases, pct_covered_1x, mean_depth, median_depth, pct_covered_2x/5x/10x, max_gap_length, num_gaps, pct_near_probe (used by probe-coverage QC)
- `compute_coverage(sam)` — pipeline read coverage (skips secondary alignments)
- `compute_probe_coverage(sam)` — probe tiling depth (includes secondary alignments)
- `write_coverage_intervals(path, coverage)` — run-length encode Vec<u32> depth vectors into interval TSV (reference_id, start, end, depth; 1-based inclusive)

### PAF (`alignment/paf.rs`)

- **`PafRecord`** — structured PAF alignment record: query_name, query_length, query_start, query_end, target_name, target_length, target_start, target_end, matching_bases, block_length, mapq
- `filter_paf(paf, max_mismatches, min_match_bases) → HashSet<String>` — filter for capture
- `parse_paf_records(paf) → Vec<PafRecord>` — parse all records (no filtering)

### FASTA (`fasta/`)

- `parse_fasta(path) → HashMap<String, String>` (id → sequence, full load)
- `parse_fasta_ids(path) → Vec<String>` (IDs only, streaming)
- `extract_by_ids(fasta, ids, output) → usize` (streaming extraction)
- `concatenate_fastas(inputs, output)` (streaming concatenation)

### Sampling (`sampling/`)

- `generate_fragments(sequences, weights, num, output, seed, length_params) → usize`
- Fragment naming: `{seq_id}_fragment_{n} start={pos} length={len}`
- `generate_weights(target_ids, distractor_ids, sample_weights, distractor_fraction, output)`

### Target Similarity (`target_similarity.rs`)

- **`TargetSimilarity`** — pairwise record: target_a, target_b, identity_pct, matching_bases, len_a, len_b
- **`SimilarityContext`** — cross_species_similar (target→set of similar targets from other species), target_is_unique, species_targets, target_to_species
- **`SpeciesDiscriminability`** — per-species: total/unique/shared targets, discriminability_score, confusable_with
- `compute_target_similarity(fasta, preset, threshold, work_dir)` — minimap2 all-vs-all, filter by identity
- `build_similarity_context(similarities, genome_to_targets)` — classify targets as unique/shared
- `compute_discriminability(ctx)` — per-species discriminability scores
- `build_confusion_matrix(ctx)` — species×species shared target count matrix

### Species Identification (`commands/identify.rs`)

- **`SpeciesCall`** enum — Present | Absent | Ambiguous
- **`SpeciesCallResult`** — species_id, call, target counts (total/unique/shared/detected), reads, explained_by, reason
- `call_species(ctx, detail_rows, min_unique_targets)` — ordered-explanation algorithm: sort by evidence strength, call PRESENT for unique marker hits, ABSENT when all hits explained by cross-reactivity, AMBIGUOUS when indeterminate

### IO Utilities (`io_utils.rs`)

- `prefixed_join(dir, prefix, filename) → PathBuf` — join directory with optionally prefixed filename; used by all directory-based commands to support `--output-prefix`
- `parse_id_set(path) → HashSet<String>` — one ID per line, # comments
- `extract_source_id(read_name) → Option<&str>` — extract `seq_id` from `{seq_id}_fragment_{n}...`
- `parse_sample_manifest(path) → HashMap<String, f64>` — TSV id\tweight
- `resolve_sample_arg(tokens) → HashMap<String, f64>` — resolves `--sample` CLI arg: single token that is an existing file → parse as TSV; otherwise → parse as inline ID list
- `parse_sample_inline(tokens) → HashMap<String, f64>` — inline format: IDs default to weight 1.0; a number after an ID sets that ID's weight (e.g. `t1 t2 t3 5 t4`)
- `format_sample_display(samples) → String` — human-readable display of sample HashMap (e.g. `t1 t2 t3(5.0) t4`)
- `parse_sample_target_map(path) → HashMap<String, Vec<String>>` — TSV genome_id\ttarget_id, # comments; returns genome→[targets] mapping
- `write_sample_target_map(map, path)` — writes genome→target mapping as commented TSV

## External Tool Wrappers (`external/`)

All wrappers follow the pattern: `check_available() → bool/Result`, then specific invocation functions.

| Tool | Functions | Output format |
|------|-----------|---------------|
| minimap2 | `capture_align` (PAF), `map_reads` (SAM), `host_align` (SAM), `probe_align` (SAM, with secondary), `xreact_align` (PAF, with secondary) | PAF or SAM |
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
| `species_calls_file` | species_calls.tsv (optional, from --identify) |
| `run_name` | string |

### Report Sections

1. **Run parameters** — curated table from run_params.tsv
2. **Command** — reconstructed from run_params.tsv (auto-adapts to new params)
3. **Capture Summary** — captured vs not-captured pie chart, captured by source pie chart
4. **Pipeline Flow** — custom ggplot2 flow diagram (sigmoid-curved bands) showing: [source bars] → generated → captured/not-captured → [sequenced] → [filtered] → correctly/incorrectly mapped. Features a zoom expansion effect where the small captured flow widens to ~85% of generated height for downstream detail. Defined inline via `pipeline_flow_diagram()` function.
5. **Detection Performance** — sensitivity/specificity/precision/F1 bar chart
6. **Read Mapping Accuracy** — correct vs incorrect mapped reads
7. **Confusion Matrix** — TP/FN/FP/TN heatmap (distractor row hidden when no distractors present)
8. **Detection Detail** — per-reference table with coverage stats
9. **Detection Lollipop** — reads per detected reference, colored by classification
10. **Coverage** (conditional) — faceted overview + expandable per-reference detail plots
11. **Species Identification** (conditional, if species_calls.tsv exists) — species calls summary, bar chart, evidence table

## Intermediate Files

| File | Format | Written by | Read by |
|------|--------|------------|---------|
| `combined_reference.fa` | FASTA | prepare | simulate; map_reads (standard mode) |
| `mapping_reference.fa` | FASTA | prepare (genome mode) | map_reads (genome mode) |
| `weights.txt` | TSV (id weight) | prepare | simulate |
| `targets.txt` | ID list | prepare | metrics |
| `distractors.txt` | ID list | prepare | metrics |
| `genomes.txt` | ID list | prepare (genome mode) | enrich (genome mode), metrics |
| `sample.txt` | ID list | prepare | metrics |
| `sample_target_map.txt` | TSV (genome_id target_id) | prepare (genome mode) | metrics, coverage_curve |
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
| `coverage.tsv` | TSV intervals (reference_id start end depth) | metrics | report |
| `species_calls.tsv` | TSV | identify (optional) | report |
| `species_calls.json` | JSON | identify (optional) | — |
| `target_similarity.tsv` | TSV | identify / panel_qc | identify |
| `report.html` | HTML | report | — |

### Panel QC (standalone, not part of pipeline)

| File | Format | Written by | Read by |
|------|--------|------------|---------|
| `target_similarity.tsv` | TSV (target_a, target_b, identity_pct, matching_bases, len_a, len_b) | panel_qc | panel_qc report, identify |
| `species_discriminability.tsv` | TSV (species_id, total/unique/shared targets, score, confusable) | panel_qc | panel_qc report |
| `species_confusion_matrix.tsv` | TSV (species × species shared target counts) | panel_qc | panel_qc report |
| `panel_qc_report.html` | HTML | panel_qc (via R) | — |

### Panel QC Report (`R/panel_qc.Rmd`)

`panel_qc.R` accepts CLI args and calls `rmarkdown::render()` on `panel_qc.Rmd`.

| Param | Source |
|-------|--------|
| `discriminability_file` | species_discriminability.tsv |
| `matrix_file` | species_confusion_matrix.tsv |
| `similarity_file` | target_similarity.tsv |
| `params_file` | run_params.tsv |

Report sections:
1. **Panel Summary** — species count, target count, unique/partial/zero discriminability counts
2. **Species Discriminability** — bar chart (≤50 species) or histogram (>50), colored by tier
3. **Target Composition** — stacked bar of unique vs shared targets per species
4. **Species Confusion Matrix** — heatmap (≤30 species) or summary stats (>30)
5. **Discriminability Table** — full per-species detail (DT::datatable for >20)
6. **Target Similarity Pairs** — pairwise similarity hits above threshold

### Species Identification (standalone or optional pipeline step)

| File | Format | Written by | Read by |
|------|--------|------------|---------|
| `species_calls.tsv` | TSV (species_id, call, targets/detected counts, reads, explained_by, reason) | identify | report |
| `species_calls.json` | JSON (structured species calls array) | identify | — |

### Probe Coverage Report (`R/probe_coverage.Rmd`)

`probe_coverage.R` accepts CLI args and calls `rmarkdown::render()` on `probe_coverage.Rmd`.

| Param | Source |
|-------|--------|
| `summary_file` | probe_coverage_summary.tsv |
| `depth_file` | probe_depth.tsv |
| `multi_mapping_file` | multi_mapping_probes.tsv (optional) |
| `proximity` | `--proximity` CLI value (integer, default 50) |

Report sections adapt dynamically based on target count: tables switch from kable to DT::datatable for >20 targets, bar charts switch to histograms/boxplots for >100 targets, and individual depth plots are omitted for >100 targets.

### Probe Coverage (standalone, not part of pipeline)

| File | Format | Written by | Read by |
|------|--------|------------|---------|
| `probe_alignment.sam` | SAM | probe_coverage | probe_coverage |
| `probe_depth.tsv` | TSV intervals (reference_id start end depth) | probe_coverage | probe_coverage report |
| `probe_coverage_summary.tsv` | TSV (per-target stats) | probe_coverage | probe_coverage report |
| `multi_mapping_probes.tsv` | TSV (probe_id num_targets targets) | probe_coverage | probe_coverage report |
| `probe_coverage_report.html` | HTML | probe_coverage (via R) | — |

### Cross-Reactivity (standalone, not part of main pipeline)

| File | Format | Written by | Read by |
|------|--------|------------|---------|
| `hits.tsv` | TSV (probe_id, target_id, homology_pct, identity_pct, query_coverage_pct, matching_bases, alignment_length, probe_length, mode) | xreact | xreact report |
| `summary.tsv` | TSV (probe_id, mode, max_homology_pct, best_hit, num_hits_above_threshold) | xreact | xreact report |
| `xreact_report.html` | HTML | xreact (via R) | — |

### Cross-Reactivity Report (`R/xreact.Rmd`)

`xreact.R` accepts CLI args and calls `rmarkdown::render()` on `xreact.Rmd`.

| Param | Source |
|-------|--------|
| `hits_file` | hits.tsv |
| `summary_file` | summary.tsv |
| `threshold` | `--threshold` CLI value (float, default 80.0) |

Report sections (conditional on mode):

1. **Summary** — threshold, probe count, per-mode hit counts
2. **Self-Homology** (if `--self`) — plotly heatmap (<1000 probes), density plot, DT hits table
3. **Cross-Reactivity** (if `--against`) — plotly heatmap (<1000 probes), per-genome bar chart, density plot, DT hits table

Heatmaps show axis labels when ≤20 items on that axis.

### Coverage Curve (standalone, not part of main pipeline)

Runs the pipeline for each parameter combination (CT × fold-enrichment × num-sequences) and computes depth curves. Nested loop optimization: prepare/simulate/capture shared per CT, enrich shared per CT×FE, sequence/filter/map per combo.

| File | Format | Written by | Read by |
|------|--------|------------|---------|
| `_prep_ct_N/combined_reference.fa` | FASTA | prepare (per CT) | simulate, map_reads |
| `_prep_ct_N/fragments.fa` | FASTA | simulate (per CT) | capture, enrich |
| `_prep_ct_N/captured.fa` | FASTA | capture (per CT) | enrich, sequence |
| `_prep_ct_N/_enrich_fe_X/enriched.fa` | FASTA | enrich (per CT×FE) | sequence |
| `{combo}/reads.fa` | FASTA | sequence (per combo) | filter/map_reads |
| `{combo}/mapped.sam` | SAM | map_reads (per combo) | coverage_curve (coverage) |
| `coverage_curve_depth_curves.tsv` | TSV (ct, fold_enrichment, num_sequences, reference_id, depth_threshold, pct_covered) | coverage_curve | coverage_curve report |
| `coverage_curve_report.html` | HTML | coverage_curve (via R) | — |

Combo directory names use only swept params: `ct_20`, `ct_20_fe_100`, `ct_20_fe_100_ns_500`. Single combo uses `run/`.

### Coverage Curve Report (`R/coverage_curve.Rmd`)

`coverage_curve.R` accepts CLI args and calls `rmarkdown::render()` on `coverage_curve.Rmd`.

| Param | Source |
|-------|--------|
| `sweep_file` | coverage_curve_depth_curves.tsv |
| `sample_ids` | comma-separated sample target IDs |
| `swept_params` | comma-separated swept parameter names |

Report logic: detects swept params from data, builds combo labels. <10 combos: single plot colored by combo. ≥10 combos: faceted by param with fewest levels, colored by remaining. Sample name in title. Summary table at key depth thresholds.

## Key Conventions

- **Error handling**: `anyhow::Result` throughout; commands return `Result<()>`
- **Streaming I/O**: Large files (FASTA, SAM) processed via BufReader/BufWriter, not loaded fully
- **Fragment naming**: `{source_id}_fragment_{n} start={pos} length={len}` — source ID extracted via `io_utils::extract_source_id`
- **Sequence IDs**: First whitespace-delimited word of FASTA header (no spaces allowed in names)
- **Weights**: `explicit_weight * sequence_length` for sampling probability; weight 0 = no fragments
- **Sample resolution**: `--sample` accepts either a TSV file path or inline IDs with optional weights; resolved to `HashMap<String, f64>` in `main.rs` via `io_utils::resolve_sample_arg` before pipeline runs. In genome mode, sample IDs refer to genome IDs (not target IDs)
- **Genome mode**: When `--genomes` is provided, fragments are generated from full genomes (combined_reference.fa = genomes+distractors), but reads are mapped back to targets (mapping_reference.fa = targets+distractors). The sample-target-map links genome IDs to target IDs. Auto-linking: genome IDs are matched to targets by (1) exact name match or (2) prefix match where target ID starts with `{genome_id}|` (e.g., genome `Bartonella_grahamii` auto-links to targets `Bartonella_grahamii|ompB` and `Bartonella_grahamii|16S`). Explicit mappings take precedence. Untargeted genomes (no target mapping) are tracked separately
- **Correct mapping (genome mode)**: A read from genome G mapping to target T is correct if T is in the genome_to_targets mapping for G
- **CT conversion**: `target_fraction = ct_baseline_fraction * 2^(ct_baseline - ct)`; resolved to `distractor_fraction` in `main.rs` before pipeline runs
- **Capture filtering**: minimap2 → PAF → filter by mismatches/indels/match-bases; BLAST → outfmt 6 → filter by gaps/nident
- **Coverage**: Single-pass SAM parsing, CIGAR ops M/=/X increment depth, D/N advance position only
- **Report extensibility**: run_params.tsv drives command reconstruction — add new params there and the report picks them up automatically
