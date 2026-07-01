# BaitBench Architecture

Quick reference for the codebase structure, module responsibilities, and key data types.

## Pipeline Flow

### Standard Mode
```
Step 1: prepare    → combined_reference.fa, weights.txt, targets.txt, distractors.txt, sample.txt
Step 2: simulate   → fragments.fa  (probe-biased via TNN thermodynamics + background, controlled by --capture-fraction)
Step 3: sequence   → reads.fa / reads.fq  (extension matches --output-format; default fasta)
Step 4: filter     → filtered.fa / filtered.fq (optional, if --host-fasta)
Step 5: map_reads  → mapped.sam  (against combined_reference.fa)
Step 6: list       → detected.list
Step 7: metrics    → results.tsv, detected_detail.tsv, results.json, coverage.tsv
Report: report     → report.html (optional, requires R)
```

### Genome Mode (`--genomes`)
```
Step 1: prepare    → combined_reference.fa (genomes+distractors), mapping_reference.fa (targets+distractors),
                     weights.txt, targets.txt, distractors.txt, genomes.txt, sample.txt, sample_target_map.txt
Step 2: simulate   → fragments.fa  (probe-biased fragments from combined_reference.fa — genomes+distractors)
Step 3: sequence   → reads.fa / reads.fq  (extension matches --output-format; default fasta)
Step 4: filter     → filtered.fa / filtered.fq (optional, if --host-fasta)
Step 5: map_reads  → mapped.sam  (against mapping_reference.fa — targets+distractors)
Step 6: list       → detected.list
Step 7: metrics    → results.tsv, detected_detail.tsv, results.json, coverage.tsv  (genome-aware classification)
Step 8: identify   → species_calls.tsv, species_calls.json  (optional, if --identify)
Report: report     → report.html  (includes Species Identification section if species_calls.tsv exists)
```

The `simulate` step replaces the old separate simulate → capture → enrich pipeline.
It aligns probes to `combined_reference.fa` internally, scores each binding site using
the SantaLucia (1998) thermodynamic nearest-neighbor model, and generates two pools of
fragments: probe-site-biased (fraction = `--capture-fraction`) and background (remainder).
Fold enrichment is no longer a pipeline parameter — target enrichment is emergent from
probe TNN affinity × sequence weights.

`run.rs` orchestrates all steps. Each step is also available as a standalone subcommand.

## Source Layout

```
src/
├── main.rs              # Entry point: CLI parsing → command dispatch
├── cli.rs               # clap argument definitions (Commands enum, all flags)
├── cleanup.rs           # Post-pipeline cleanup: delete intermediate files/dirs, keep report inputs
├── io_utils.rs           # Shared helpers: prefixed_join, parse_id_set, extract_source_id, parse_sample_manifest, parse_sample_target_map, parse_groups_file, write_groups_file
├── target_similarity.rs  # Shared library: target-vs-target similarity computation, discriminability scoring, confusion matrices
├── alignment/
│   ├── coverage.rs      # CIGAR-based per-position coverage from SAM
│   ├── paf.rs           # PAF record filtering (mismatch/indel criteria) + structured PafRecord parsing
│   └── sam.rs           # SAM parsing: read counts, mappings, mapped IDs
├── thermodynamics.rs    # SantaLucia (1998) nearest-neighbor TNN model: delta_g(), boltzmann_score()
├── commands/
│   ├── run.rs           # Full pipeline orchestrator (steps 1-7 + report)
│   ├── prepare.rs       # Combine FASTAs, generate weights, write ID lists
│   ├── simulate.rs      # Thermodynamic/simple probe-biased fragment simulation (replaces simulate+capture+enrich)
│   ├── sequence.rs      # Simulate sequencing: perfect (trim), art (Illumina via ART-modern), badread (long reads — ONT/PacBio CLR); dispatches ReadSimulator enum; art renaming via SAM RNAME; badread renaming via ref_name in read description
│   ├── filter.rs        # Remove host-mapping reads
│   ├── map_reads.rs     # Map reads to reference (minimap2)
│   ├── generate_list.rs # Count reads per reference from SAM → detected.list
│   ├── metrics.rs       # TP/FP/FN/TN classification, coverage stats, TSV/JSON output
│   ├── report.rs        # Report generation: HTML (Rscript), RMarkdown (template substitution), or skip; shared substitute_rmd_params utility
│   ├── probe_coverage.rs # Standalone probe tiling QC (maps probes to targets)
│   ├── xreact.rs        # Standalone cross-reactivity analysis (probes vs genomes, probes vs probes)
│   ├── panel_qc.rs      # Standalone target panel discriminability QC (target-vs-target similarity, species discrimination)
│   ├── identify.rs      # Species-level calling from multi-target detection patterns (standalone or pipeline step)
│   ├── coverage_curve.rs # Coverage curve: pipeline at multiple param combos (CT × temp × CF × NS) → depth curves
│   ├── build_probes.rs  # Standalone probe building: N filter → collapse → length filter → tile/CATCH/Syotti/ProbeTools → GC filter → complexity filter (sDUST) → deduplicate
│   ├── assess_probes.rs # Combined probe assessment: probe coverage + cross-reactivity (self + optional genomes), orchestrates sub-commands
│   └── tool_dustview.rs # `baitbench tool dustview` handler: execute() — sDUST masking visualization on FASTA sequences
├── sdust.rs             # sDUST low-complexity detection: sdust(), masked_fraction() (Morgulis et al. 2006)
├── syotti.rs            # Syotti greedy bait design: design_probes() — k-mer hash index, seed-and-extend, greedy set-cover (Alanko et al. 2022)
├── catch.rs             # Native CATCH probe design: design_probes() — tiling → MinHash dedup → greedy set cover (reimplementation of Metsky et al. 2019)
├── probetools.rs        # ProbeTools-Lite probe design: design_probes() — iterative k-mer enumeration → cd-hit-est clustering → greedy batch selection → minimap2 coverage → low-cov extraction (reimplementation of Kuchinski et al. 2022)
├── external/
│   ├── minimap2.rs      # minimap2 wrapper: capture_align (PAF), map_reads (SAM, optional reads_r2 for PE), host_align (optional reads_r2), probe_align(max_secondary), xreact_align
│   ├── art_modern.rs    # ART-modern wrapper: check_available, run_simulation (Illumina SE/PE reads + SAM for renaming)
│   ├── badread.rs       # badread wrapper: check_available, profile_params (ont/ont-2020/pacbio → error+qscore models), run_simulation (long reads; FASTQ to stdout)
│   ├── blastn.rs        # BLAST+ wrapper: capture_align, filter_blast_results
│   ├── cdhit.rs         # cd-hit-est wrapper: check_available, cluster (sequence clustering by identity)
│   └── rscript.rs       # Rscript discovery (BAITBENCH_R_DIR, binary walk, ./R/) and execution
├── fasta/
│   ├── reader.rs        # parse_fasta (id→seq), parse_fasta_ids, count_sequences
│   └── writer.rs        # write_fasta_record, extract_by_ids (streaming), concatenate_fastas
└── sampling/
    ├── fragment.rs      # generate_fragments: weighted sampling, normal-dist length
    ├── thermo_sim.rs    # Thermodynamic simulation: load_probe_hits (SAM→ProbeHit), sample_capture_fragments, sample_background_fragments, write_fragments; SimulateMode enum
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
├── panel_qc.Rmd         # RMarkdown template: discriminability charts, confusion heatmaps, target tables
├── build_probes.R       # CLI wrapper for build probes report
├── build_probes.Rmd     # RMarkdown template: pipeline stats table, sequence/base count bar charts
├── assess_probes.R      # CLI wrapper for combined probe assessment report
└── assess_probes.Rmd    # RMarkdown template: build stats (optional) + probe coverage + cross-reactivity
```

## Key Data Types

### CLI (`cli.rs`)

- **`Commands`** enum — one variant per subcommand (Run, Prepare, Simulate, Sequence, Filter, Map, List, Metrics, ProbeCoverage, Xreact, PanelQc, Identify, Report, CoverageCurve, BuildProbes, AssessProbes, Syotti), each with its own fields
- **`SimulateMode`** — ValueEnum: Thermodynamic | Simple — controls probe-site weighting in simulate
- **`ProbeMethod`** — ValueEnum: Tile | Catch | Syotti
- **`ReportMode`** — ValueEnum: Full | None | Rmd | BothR — controls report output (HTML, skip, editable RMarkdown, or both HTML and RMarkdown)
- **CT score flags** — `--ct`, `--ct-baseline`, `--ct-baseline-fraction` on Run and Prepare; `--ct` conflicts with `--distractor-fraction`
- **Simulate flags** — `--probes` (probe FASTA), `--simulate-mode` (thermodynamic/simple), `--hybridization-temperature` (°C, default 70), `--capture-fraction` (0–1, default 0.5) on Run and Simulate; `--hybridization-temperature-values` on CoverageCurve for temperature sweep
- **Genome mode flags** — `--genomes` (optional genome FASTA for fragment generation), `--sample-target-map` (optional genome-to-target mapping TSV) on Run, Prepare, and CoverageCurve
- **`--output-prefix`** — string prepended to every auto-generated output filename; available on Run, Prepare, ProbeCoverage, CoverageCurve, Xreact, PanelQc, Identify (default: empty string)

### Command Args Pattern

Every command module exports an `Args` struct and an `execute(&Args) -> Result<()>` function:

| Command | Args struct | Key inputs | Key outputs |
|---------|-------------|------------|-------------|
| `prepare` | `PrepareArgs` | targets, distractors, sample, distractor_fraction, genomes, sample_target_map, groups, distractor_groups | combined_reference.fa, weights.txt, ID lists, distractor_groups.tsv (always); target_groups.tsv (if --groups); genome mode also: mapping_reference.fa, genomes.txt, sample_target_map.txt |
| `simulate` | `SimulateArgs` | reference, weights, probes, num_fragments, capture_fraction, simulate_mode, hybridization_temperature, seed, fragment_length_*, threads | fragments.fa (probe-biased + background) |
| `sequence` | `SequenceArgs` | input, output, output_r2, read_length, num_sequences, seed, simulator (perfect/art/badread), sequencer_profile, coverage_depth, paired_end, pe_frag_len_mean, pe_frag_len_sd, output_format | reads.fa/.fq (SE) or reads.fa/.fq + reads_R2.fa/.fq (PE); extension determined by output_format; errors injected by art_modern or badread when simulator ≠ perfect |
| `filter` | `FilterArgs` | host, reads, reads_r2, minimap_preset | filtered.fa/.fq (+ filtered_R2.fa/.fq for PE); format detected automatically from input |
| `map_reads` | `MapArgs` | reference, reads, reads_r2, minimap_preset | mapped.sam |
| `generate_list` | `ListArgs` | sam | detected.list |
| `metrics` | `MetricsArgs` | targets, distractors, sample, detected, fragments, captured, sam, sample_target_map, target_groups, distractor_groups, reads_sequenced, reads_after_filter, output_group_detail | results.tsv, detected_detail.tsv (with group col), group_detail.tsv (when groups present), results.json, coverage.tsv |
| `report` | `ReportArgs` | summary, detail, params, coverage, group_detail, run_name, report (ReportMode) | report.html or report.Rmd |
| `probe_coverage` | `ProbeCoverageArgs` | targets, probes, minimap_preset, proximity | probe_depth.tsv, probe_coverage_summary.tsv, probe_coverage_report.html |
| `xreact` | `XreactArgs` | probes, against (genome FASTAs), self_mode, threshold, minimap_preset, report (ReportMode) | hits.tsv, summary.tsv, xreact_report.html |
| `panel_qc` | `PanelQcArgs` | targets, sample_target_map, identity_threshold, minimap_preset, report (ReportMode) | target_similarity.tsv, species_discriminability.tsv, species_confusion_matrix.tsv, panel_qc_report.html |
| `identify` | `IdentifyArgs` | detected_detail, sample_target_map, target_similarity (or targets for on-the-fly), identity_threshold, min_unique_targets | species_calls.tsv, species_calls.json |
| `run` | `RunArgs` | all pipeline inputs + ct, ct_baseline, ct_baseline_fraction, simulate_mode, hybridization_temperature, capture_fraction, num_sequences, simulator, sequencer_profile, coverage_depth, paired_end, pe_frag_len_mean/sd, genomes, sample_target_map, identify, identity_threshold, min_unique_targets | all of the above |
| `coverage_curve` | `CoverageCurveArgs` | targets, distractors, probes, sample (required), ct/df/cf/ns values (sweep or fixed), simulate_mode, hybridization_temperature, all pipeline params, genomes, sample_target_map | coverage_curve_depth_curves.tsv, coverage_curve_report.html, combo subdirs |
| `build_probes` | `BuildProbesArgs` | targets, method (tile/catch-lite/syotti-lite/catch), probe_length, step, catch_probe_stride/mismatches/extension/coverage/minhash_threshold, syotti_mismatches, syotti_seed_len, max_n_frac, min/max_gc, dust_threshold/dust_window/max_masked_frac, collapse/dedup thresholds, threads, genomes, threshold, skip_assess | probes_final.fa, build_probes_stats.tsv; filters sequences shorter than probe_length after collapse; auto-chains to assess_probes unless --skip-assess |
| `tool syotti` | — (standalone) | targets, output, probe_length, mismatches, seed_len | output FASTA of probes; direct access to Syotti algorithm |
| `tool catch` | — (standalone) | targets, output, probe_length, stride, mismatches, extension, coverage, minhash_threshold | output FASTA of probes; direct access to CATCH algorithm |
| `tool dustview` | — (standalone) | input (optional, defaults stdin), dust_threshold, dust_window | stdout: per-sequence masked view + score stats |
| `tool collapse` | — (standalone) | input, output, threshold, threads, log_file | output FASTA of cd-hit-est cluster representatives |
| `assess_probes` | `AssessProbesArgs` | targets, probes, genomes (optional), threshold, minimap_preset, proximity, build_stats_file (optional), build_params_file (optional), all_individual_targets | cov_probe_coverage_summary.tsv, cov_probe_depth.tsv, xreact_hits.tsv, xreact_summary.tsv, assess_run_params.tsv, individual_target_coverage_summary.tsv (optional), assess_probes_report.html |

### Metrics (`metrics.rs`)

- **`MetricsResult`** — group-level classification: tp/fn/fp_target/fp_distractor/tn_target/tn_distractor counts, sensitivity/specificity/precision/f1, plus group name lists (true_positives, false_negatives, etc.), untargeted_genomes list
- **`ReadLevelMetrics`** — fragment/read counts: sample_captured, nonsample_target_captured, distractor_captured, untargeted_captured, reads_correctly_mapped, reads_incorrectly_mapped. Summary TSV also includes derived pipeline flow counts: reads_sequenced, reads_after_filter, reads_mapped, reads_unmapped
- **`GenomeContext`** — derived from sample_target_map: sample_targets, genome_ids, sample_genome_ids, genome_to_targets, target_to_genomes, untargeted_genomes (genomes with no target mapping)
- **`DetailRow`** (Serialize) — per-reference row: reference_id, group, category, expected, detected, fragments_generated, fragments_captured, reads_assigned, classification, ref_length, avg_coverage, pct_covered_5x, pct_covered_20x
- **`GroupDetailRow`** — per-group row for `group_detail.tsv`: group_name, category, expected, detected, classification, member_count, detected_member_count, total_reads
- **`JsonOutput`** (Serialize) — wraps RunInfo, CaptureStats, ReadLevelStats, JsonMetrics, JsonDetails
- Group maps are loaded via `load_group_map(path)` → `HashMap<String, String>` (seq_id → group_name). Missing entries default to identity mapping.
- `collapse_to_groups(detected, group_map)` — collapses per-sequence detection counts to per-group counts

### Coverage (`alignment/coverage.rs`)

- **`CoverageResult`** — ref_lengths: HashMap<String, usize>, coverage: HashMap<String, Vec<u32>>
- **`CoverageStats`** — ref_length, avg_coverage, pct_covered_5x, pct_covered_20x (used by pipeline read coverage)
- **`ProbeCoverageStats`** — ref_length, covered_bases, pct_covered_1x, mean_depth, median_depth, pct_covered_2x/5x/10x, max_gap_length, num_gaps, pct_near_probe (used by probe-coverage QC)
- `compute_coverage(sam)` — pipeline read coverage (skips secondary alignments)
- `compute_probe_coverage(sam)` — probe tiling depth (includes secondary alignments; used for pangenome coverage and for each per-target SAM in individual-target mode)
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

### Thermodynamics (`thermodynamics.rs`)

- `ThermoModel { temp_c: f64, na_conc_m: f64 }` — parameter struct for hybridization temperature (°C) and Na+ concentration (molar); constructed via `ThermoModel::new(temp_c, na_conc_m)` from `--hybridization-temperature` and `--salt-concentration` (mM→M)
- `delta_g(aligned_pairs: &[(u8,u8)], model: &ThermoModel) -> f64` — SantaLucia (1998) nearest-neighbor free energy in kcal/mol; three contributions:
  1. NN stacking energy (SkipStacking: mismatches break stacking chain)
  2. Initiation terms (SantaLucia 1998 Table 2): AT (+2.3/+4.1) or GC (+0.1/−2.8) applied once per terminal WC pair, **only when at least one stacking step exists** (see design note below)
  3. Salt correction (Owczarzy et al. 1997): `ΔS += 0.368 × (n_wc−1) × ln([Na+])`; at 1 M Na+ ln(1)=0 so no correction; at 50 mM (~typical hybridization buffer) makes ΔG more positive (weaker binding)
- `boltzmann_score(dg: f64, model: &ThermoModel) -> f64` — `exp(-ΔG / (R × T_K))` where R = 1.987e-3 kcal/(mol·K); score = 1 means neutral (no hybridization), > 1 means favorable binding, < 1 is not produced by this model

#### Thermodynamic model design notes

**Score semantics.** `boltzmann_score` is a Boltzmann factor proportional to the hybridization equilibrium constant. It is used as a relative sampling weight — only ratios between scores matter, not their absolute values. Score = 1 (ΔG = 0) is the neutral baseline: no thermodynamic tendency to bind.

**Initiation terms.** SantaLucia (1998) was designed to predict the *absolute* ΔG of duplex formation from two free strands in solution. Initiation parameters capture two physical costs: (1) the loss of translational/rotational entropy when two separate molecules associate, and (2) the instability of terminal base pairs, which are less constrained than interior ones. AT termini (+2.3/+4.1 kcal/mol) are more costly than GC termini (+0.1/−2.8) because A-T pairs have only two hydrogen bonds and are weaker duplex anchors. Initiation always adds positive ΔG (destabilising), reducing the Boltzmann score relative to stacking alone.

**Why initiation is conditional on stacking.** Initiation terms are only applied when at least one NN stacking step exists (i.e., ≥ 2 consecutive WC pairs). An isolated single WC pair flanked by mismatches on both sides cannot sustain a stable duplex — applying the initiation nucleation cost to it would produce ΔG > 0 and score < 1, implying the probe is actively repelled from that locus, which is unphysical. When no stacking exists, ΔG = 0 and score = 1 (neutral), matching the behaviour of RAmpSim (see below).

**Comparison with RAmpSim.** RAmpSim (Rooney et al. 2025, biorxiv 2025.12.05.692407) implements the same SantaLucia NN stacking table and SkipStacking mismatch strategy but omits initiation and salt correction entirely. For a relative-scoring use case, the constant part of the initiation penalty cancels across all probes when normalising sampling weights. The residual effect — AT vs GC terminal distinction — is a secondary correction (~2 kcal/mol) that BaitBench preserves because it is physically real: probes with GC-rich termini form more stable duplexes than AT-rich ones, and this difference affects relative capture probabilities. The salt correction has a larger practical impact: at 50 mM Na+ (typical hybridization buffer), it reduces scores by an additional factor that scales with duplex length, meaningfully changing the relative ranking of long vs short probe hits.

### Sampling (`sampling/`)

- `generate_fragments(sequences, weights, num, output, seed, length_params) → usize`
- Fragment naming: `{seq_id}_fragment_{n} start={pos} length={len}`
- `generate_weights(target_ids, distractor_ids, sample_weights, distractor_fraction, output)`

**`thermo_sim.rs`:**
- `SimulateMode` — enum: `Thermodynamic` | `Simple`
- `ProbeHit` — struct: probe_name, seq_id, start (0-based), end (exclusive), score
- `load_probe_hits(sam_path, weights, mode, model: &ThermoModel) → HashMap<String, Vec<ProbeHit>>`
  - Parses SAM; scores each hit as `boltzmann_score(delta_g(..., model), model) × seq_weight` (thermodynamic) or `seq_weight` (simple)
  - Skips unmapped, secondary-flag-dropped, and weight-0 sequences
- `sample_capture_fragments(hits_by_probe, sequences, n, length_params, seed, counter) → Vec<(header, seq)>`
  - Two-level multinomial: probe (uniform) → hit (weighted by score) → fragment center (probe_center ± frag_len/4 jitter) → length (Normal)
- `sample_background_fragments(sequences, weights, n, length_params, seed, counter) → Vec<(header, seq)>`
- `write_fragments(fragments, output_path)`

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
- `parse_groups_file(path) → HashMap<String, String>` — TSV seq_id\tgroup_name, # comments; strips leading `>` from IDs
- `write_groups_file(map, path)` — writes sorted seq_id→group_name TSV with header comment

## External Tool Wrappers (`external/`)

All wrappers follow the pattern: `check_available() → bool/Result`, then specific invocation functions.

| Tool | Functions | Output format |
|------|-----------|---------------|
| minimap2 | `capture_align` (PAF), `map_reads` (SAM), `host_align` (SAM), `probe_align` (SAM, with secondary), `xreact_align` (PAF, with secondary) | PAF or SAM |
| blastn | `capture_align` (TSV outfmt 6), `filter_blast_results` | TSV |
| cd-hit-est | `check_available`, `cluster` (identity-based sequence clustering) | FASTA + .clstr |
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
3. **Pipeline Flow — Capture** — interactive plotly Sankey diagram showing: [source nodes] → generated → captured/not-captured → captured breakdown by source type (sample, distractor, non-sample, untargeted)
4. **Pipeline Flow — Read Processing** — interactive plotly Sankey diagram showing: sequences → [sequenced/sampling loss] → [filtered/host removed] → correctly/incorrectly mapped/unmapped. Displayed at full size since it is separate from the capture diagram.
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
| `target_groups.tsv` | TSV (seq_id group_name) | prepare (only if `--groups`) | metrics, report |
| `distractor_groups.tsv` | TSV (seq_id group_name) | prepare (always; auto or `--distractor-groups`) | metrics, report |
| `fragments.fa` | FASTA | simulate | capture, metrics |
| `captured.fa` | FASTA | capture | enrich (if --fold-enrichment), sequence, metrics |
| `enriched.fa` | FASTA | enrich | sequence, metrics (only if --fold-enrichment) |
| `reads.fa` / `reads.fq` | FASTA or FASTQ | sequence | filter/map_reads |
| `filtered.fa` / `filtered.fq` | FASTA or FASTQ | filter | map_reads |
| `mapped.sam` | SAM | map_reads | generate_list, metrics |
| `detected.list` | TSV (id count) | generate_list | metrics |
| `run_params.tsv` | TSV (parameter flag value) | run | report |
| `results.tsv` | TSV | metrics | report |
| `detected_detail.tsv` | TSV (includes `group` column) | metrics | report |
| `group_detail.tsv` | TSV (one row per group) | metrics (when groups present) | report |
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

Runs the pipeline for each parameter combination (CT or distractor fraction × temp × CF × NS) and computes depth curves. Nested loop optimization: prepare/simulate shared per CT/DF × temp × CF combo; sequence/filter/map per full combo.

| File | Format | Written by | Read by |
|------|--------|------------|---------|
| `_prep_ct_N/combined_reference.fa` | FASTA | prepare (per CT/DF) | simulate, map_reads |
| `_prep_ct_N/fragments.fa` | FASTA | simulate (per CT/DF) | sequence |
| `{combo}/reads.fa` | FASTA | sequence (per combo) | filter/map_reads |
| `{combo}/mapped.sam` | SAM | map_reads (per combo) | coverage_curve (coverage) |
| `coverage_curve_depth_curves.tsv` | TSV (ct, hybridization_temperature, capture_fraction, num_sequences, reference_id, depth_threshold, pct_covered) | coverage_curve | coverage_curve report |
| `coverage_curve_report.html` | HTML | coverage_curve (via R) | — |

Combo directory names use only swept params: `ct_20`, `df_0.50`, `ct_20_temp_65`, `ct_20_temp_65_cf_0.50`, `ct_20_temp_65_cf_0.50_ns_500`. Single combo uses `run/`. `df_{:.2}` prefix is used when `--distractor-fraction-values` is swept directly.

### Coverage Curve Report (`R/coverage_curve.Rmd`)

`coverage_curve.R` accepts CLI args and calls `rmarkdown::render()` on `coverage_curve.Rmd`.

| Param | Source |
|-------|--------|
| `sweep_file` | coverage_curve_depth_curves.tsv |
| `sample_ids` | comma-separated sample target IDs |
| `swept_params` | comma-separated swept parameter names |

Report logic: detects swept params from data, builds combo labels. <10 combos: single plot colored by combo. ≥10 combos: faceted by param with fewest levels, colored by remaining. Sample name in title. Summary table at key depth thresholds.

### Build Probes (standalone, not part of pipeline)

| File | Format | Written by | Read by |
|------|--------|------------|---------|
| `targets_clean.fa` | FASTA | build_probes (N filter) | build_probes (cd-hit-est collapse) |
| `collapsed.fa` | FASTA | build_probes (cd-hit-est) | build_probes (length filter) |
| `length_filtered.fa` | FASTA | build_probes (length filter) | build_probes (tile or CATCH) |
| `probes_raw.fa` | FASTA | build_probes (tile, CATCH, or Syotti) | build_probes (GC filter) |
| `probes_gc.fa` | FASTA | build_probes (GC filter) | build_probes (complexity filter) |
| `probes_complexity.fa` | FASTA | build_probes (sDUST filter) | build_probes (cd-hit-est dedup) |
| `probes_final.fa` | FASTA | build_probes (cd-hit-est dedup) | user (final output) |
| `build_probes_stats.tsv` | TSV (step, sequences, bases) | build_probes | build_probes report |
| `build_probes_report.html` | HTML | build_probes (via R) | — |

### Assess Probes (standalone, or chained from build-probes)

| File | Format | Written by | Read by |
|------|--------|------------|---------|
| `cov_probe_coverage_summary.tsv` | TSV | assess_probes (via probe_coverage) | assess_probes report |
| `cov_probe_depth.tsv` | TSV intervals | assess_probes (via probe_coverage) | assess_probes report |
| `cov_multi_mapping_probes.tsv` | TSV | assess_probes (via probe_coverage) | assess_probes report |
| `xreact_hits.tsv` | TSV | assess_probes (via xreact) | assess_probes report |
| `xreact_summary.tsv` | TSV | assess_probes (via xreact) | assess_probes report |
| `assess_run_params.tsv` | TSV (parameter, flag, value) | assess_probes | assess_probes report |
| `individual_target_coverage_summary.tsv` | TSV (same columns as cov_probe_coverage_summary.tsv) | assess_probes (`--all-individual-targets` only) | assess_probes report |
| `refine_summary.tsv` | TSV (step, num_targets, num_below_threshold, pct_covered) | assess_probes (when `--refine-iterations` or `--refine-until-stable` used) | assess_probes report |
| `assess_probes_report.html` | HTML | assess_probes (via R) | — |

### Assess Probes Report (`R/assess_probes.Rmd`)

`assess_probes.R` accepts CLI args and calls `rmarkdown::render()` on `assess_probes.Rmd`.

| Param | Source |
|-------|--------|
| `build_stats_file` | build_probes_stats.tsv (optional, from build-probes chain) |
| `build_params_file` | run_params.tsv (optional, from build-probes chain) |
| `xreact_hits_file` | xreact_hits.tsv |
| `xreact_summary_file` | xreact_summary.tsv |
| `xreact_threshold` | `--threshold` CLI value (float, default 80.0) |
| `coverage_summary_file` | cov_probe_coverage_summary.tsv |
| `coverage_depth_file` | cov_probe_depth.tsv |
| `coverage_multi_mapping_file` | cov_multi_mapping_probes.tsv (optional) |
| `coverage_proximity` | `--proximity` CLI value (integer, default 50) |
| `params_file` | assess_run_params.tsv |
| `individual_coverage_file` | individual_target_coverage_summary.tsv (optional, from `--all-individual-targets`) |
| `refine_summary_file` | refine_summary.tsv (optional, from `--refine-iterations` / `--refine-until-stable`) |

Report sections (conditionally rendered):
1. **Build Pipeline** (if build_stats_file provided) — stats table, sequences/bases bar charts
2. **Probe Coverage** (always) — summary table, coverage breadth, individual target coverage (if `individual_coverage_file` provided), tiered coverage, gap analysis, pangenome depth (subtitle shows % pangenome ≥1X), depth profiles, proximity, multi-mapping probes
3. **Refinement Summary** (if refine_summary_file provided) — inline stats (iterations ran, initial/final below-threshold counts), step-by-step table, stacked bar chart of covered vs below-threshold per step
4. **Cross-Reactivity Summary** (always) — summary table
5. **Self-Homology** (if self-mode hits present) — heatmap (≤1000 probes), density plots, hits table
6. **Cross-Reactivity vs Genomes** (if against-mode hits present) — heatmap, per-genome bar chart, density plots, hits table
7. **Parameters** (under `<details>` fold) — assess params + optional build params

### Build Probes Report (`R/build_probes.Rmd`)

`build_probes.R` accepts CLI args and calls `rmarkdown::render()` on `build_probes.Rmd`.

| Param | Source |
|-------|--------|
| `stats_file` | build_probes_stats.tsv |
| `params_file` | run_params.tsv |

Report sections:
1. **Pipeline Summary** — table showing sequences/bases at each step with drop counts and retention %
2. **Sequences per Step** — bar chart of sequence counts through pipeline
3. **Total Bases per Step** — bar chart of total bases through pipeline

## Desktop GUI (`gui/`)

A Tauri v2 desktop application that wraps the `baitbench` CLI as a sidecar. The GUI does **not** share library code with the CLI — it invokes the compiled binary as a subprocess and communicates via stdin/stdout/stderr.

### Structure

```
gui/
├── src-tauri/               # Tauri Rust backend
│   ├── Cargo.toml
│   ├── tauri.conf.json      # sidecar registration, window config, plugin config
│   ├── capabilities/
│   │   └── default.json     # permission grants (shell, dialog, fs, notification, store)
│   ├── icons/
│   │   └── icon.png         # app icon (RGBA PNG)
│   ├── binaries/
│   │   └── baitbench-<triple>  # compiled CLI binary (gitignored; copied by make copy-sidecar)
│   └── src/
│       ├── main.rs          # Tauri entry point
│       ├── lib.rs           # plugin registration, AppState, command registration
│       └── commands/
│           ├── setup.rs     # detect_conda_envs, validate_conda_env
│           ├── pipeline.rs  # run_pipeline, get_pipeline_status, cancel_pipeline, AppState
│           └── report.rs    # open_report (opens HTML in a new WebviewWindow)
├── src/                     # Svelte + TypeScript frontend
│   ├── App.svelte           # root component; loads persisted env on mount, routes to correct view
│   ├── app.css              # design tokens (light/dark), reset, drag-region rules
│   ├── main.ts
│   └── lib/
│       ├── stores.ts        # currentView, pipelineStatus, logLines, reportPath, condaEnvPath
│       ├── types.ts         # CondaEnv, ValidationResult, LogLine, PipelineConfig, ToolDef
│       ├── views/
│       │   ├── SetupView.svelte  # conda env detection, validation badge, persist to store
│       │   ├── RunView.svelte    # tool selector sidebar + per-tool forms + run button
│       │   └── LogView.svelte   # real-time log panel, cancel, view-report button
│       └── components/
│           ├── FilePicker.svelte      # text input + browse button wrapping dialog:open
│           └── AdvancedOptions.svelte # <details> collapsible section
├── scripts/
│   └── copy-sidecar.mjs     # copies target/release/baitbench → binaries/baitbench-<triple>
├── Makefile                 # dev / build / copy-sidecar targets
├── package.json
├── vite.config.ts           # Svelte + svelte-preprocess (TypeScript support)
└── tsconfig.json
```

### Tauri Commands (Rust → Frontend)

| Command | Signature | Description |
|---------|-----------|-------------|
| `detect_conda_envs` | `() → Vec<CondaEnv>` | Scans common conda paths + `$CONDA_PREFIX`; returns name+path pairs |
| `validate_conda_env` | `(path) → ValidationResult` | Checks for `minimap2`, `cd-hit-est` (required) and `R`, `pandoc` (optional warnings) |
| `run_pipeline` | `(PipelineConfig) → ()` | Spawns sidecar with conda bin prepended to PATH; streams stdout/stderr as `pipeline-log` events; emits `pipeline-complete` / `pipeline-error` + OS notification on exit |
| `get_pipeline_status` | `() → PipelineStatus` | Returns current state: `idle \| running \| complete \| failed` |
| `cancel_pipeline` | `() → ()` | Sends SIGTERM to the running child process |
| `open_report` | `(path) → ()` | Opens an HTML file in a new 1200×900 WebviewWindow |

### Frontend Events (Backend → Frontend)

| Event | Payload | Emitted when |
|-------|---------|--------------|
| `pipeline-log` | `{ text: string, stream: "stdout" \| "stderr" }` | Each line of output from the sidecar |
| `pipeline-complete` | — | Sidecar exits with code 0 |
| `pipeline-error` | — | Sidecar exits with non-zero code or spawn error |

### View Flow

```
App mount
  └─ load conda_env_path from tauri-plugin-store
       ├─ found → RunView
       └─ not found → SetupView
            └─ save → RunView
                  └─ Run button clicked → LogView
                        └─ pipeline-complete → "View Report" button → report WebviewWindow
```

### Sidecar Invocation

`PipelineConfig` holds `{ tool, args: HashMap<flag, value>, conda_env }`. The Rust backend converts it to argv:
- Boolean flags (value `""`) → append flag only
- Single-value flags → `[flag, value]`
- Multi-value flags (tab-separated value) → `[flag, v1, v2, ...]`

The sidecar is invoked as `baitbench <subcommand> [flags...]` with `PATH` prepended by `{conda_env}/bin`.

### Build / Dev Workflow

```bash
# One-time setup
cd gui && npm install

# Every time the CLI changes
make copy-sidecar   # cargo build --release (root) + copy binary

# Development (hot-reload frontend, recompiles Rust backend on change)
make dev            # = make copy-sidecar && npm run tauri:dev

# Production bundle (macOS .app / .dmg, Windows .msi, Linux .AppImage)
make build          # = make copy-sidecar && npm run tauri:build
```

### Required Plugins

| Plugin | Purpose |
|--------|---------|
| `tauri-plugin-shell` | Spawn sidecar, stream output |
| `tauri-plugin-dialog` | File and directory pickers |
| `tauri-plugin-fs` | File existence checks |
| `tauri-plugin-store` | Persist conda env path across launches |
| `tauri-plugin-notification` | OS notification on pipeline completion |

---

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
