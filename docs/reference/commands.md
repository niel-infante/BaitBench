# Commands

All BaitBench subcommands with their flags and output files.

## Quick Reference

| Subcommand | Purpose |
|-----------|---------|
| [`run`](#run) | Full pipeline orchestrator |
| [`prepare`](#prepare) | Combine references, generate weights |
| [`simulate`](#simulate) | Generate probe-biased fragments |
| [`sequence`](#sequence) | Trim fragments to read length |
| [`filter`](#filter) | Optional host read filtering |
| [`map`](#map) | Align reads back to reference |
| [`list`](#list) | Count reads per reference |
| [`metrics`](#metrics) | 3-way classification and coverage stats |
| [`report`](#report) | Generate HTML report |
| [`coverage-curve`](#coverage-curve) | Sweep parameters and plot depth curves |
| [`build-probes`](#build-probes) | Design a probe set from target sequences |
| [`assess-probes`](#assess-probes) | Coverage + cross-reactivity assessment |
| [`probe-coverage`](#probe-coverage) | Standalone probe coverage analysis |
| [`xreact`](#xreact) | Probe cross-reactivity analysis |
| [`panel-qc`](#panel-qc) | Target panel discriminability QC |
| [`identify`](#identify) | Species-level calling from detection patterns |
| [`tool`](#tool) | Standalone utility tools |

---

## run

Runs the complete pipeline from input files to metrics and report.

```bash
baitbench run [OPTIONS]
```

The primary command for most use cases. Chains all pipeline steps (prepare through report) automatically. Use `--cleanup` to delete intermediate files after completion, keeping only report inputs and final outputs.

In genome mode with `--sample-target-map`, use `--identify` to add species-level calling after metrics. Species calls are included in the HTML report with ground-truth comparison against `--sample`.

`--report` accepts: `full` (HTML, default), `none` (skip), `rmd` (editable RMarkdown), `both-r` (HTML + RMarkdown).

See [Parameters](parameters.md) for all options.

---

## prepare

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

| File | Description |
|------|-------------|
| `combined_reference.fa` | Merged sequences for fragment generation |
| `weights.txt` | Per-sequence sampling weights (TSV: `id<TAB>weight`) |
| `targets.txt` | Target sequence IDs (one per line) |
| `distractors.txt` | Distractor sequence IDs (one per line) |
| `sample.txt` | Sample sequence IDs (one per line) |
| `target_groups.tsv` | Target group assignments (only if `--groups` provided) |
| `distractor_groups.tsv` | Distractor group assignments (always written) |
| `mapping_reference.fa` | Targets + distractors for read mapping (genome mode only) |
| `genomes.txt` | Genome IDs (genome mode only) |
| `sample_target_map.txt` | Genome-to-target mapping (genome mode only) |

---

## simulate

Generates weighted random fragments from a reference, biased toward probe binding sites.

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

**Output:** `fragments.fa` — simulated DNA fragments.

---

## sequence

Simulates sequencing by trimming fragments to read length.

```bash
baitbench sequence \
  --input fragments.fa \
  --output reads.fa \
  [--read-length 120] \
  [--num-sequences 5000] \
  [--seed 42]
```

Fragments shorter than `--read-length` are kept as-is. With `--num-sequences`, reads are sampled with replacement from the fragment pool, modelling PCR amplification before sequencing.

**Output:** `reads.fa` — sequencing reads.

---

## filter

Removes reads that map to a host genome.

```bash
baitbench filter \
  --host host_genome.fa \
  --reads reads.fa \
  --output filtered.fa \
  [--minimap-preset sr]
```

Uses the embedded aligner to map reads against the host genome. Reads that map are removed; unmapped reads are retained.

**Output:** `filtered.fa` — reads after host depletion.

---

## map

Maps reads back to a reference.

```bash
baitbench map \
  --reference combined_reference.fa \
  --reads reads.fa \
  --output mapped.sam \
  [--minimap-preset sr]
```

In standard mode, reads are mapped to `combined_reference.fa`. In genome mode, reads are mapped to `mapping_reference.fa` (targets + distractors).

**Output:** `mapped.sam` — SAM alignment file.

---

## list

Counts reads per reference from a SAM file.

```bash
baitbench list \
  --sam mapped.sam \
  --output detected.list
```

**Output:** `detected.list` — TSV: `reference_id<TAB>count` (sorted ascending by count).

---

## metrics

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
  [--seed 42]
```

**Output files:** `results.tsv`, `detected_detail.tsv`, optionally `results.json` and `coverage.tsv`.

See [Output Formats](output-formats.md) for column definitions.

---

## report

Generates an HTML report, or outputs an editable RMarkdown file.

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

**Output:**
- `report.html` — HTML report with ggplot2 visualisations (`--report full`)
- `report.Rmd` — editable RMarkdown with parameters pre-filled (`--report rmd`)

---

## coverage-curve

Runs the pipeline at multiple parameter combinations and generates coverage depth curves.

```bash
baitbench coverage-curve \
  --targets targets.fa \
  --distractors distractors.fa \
  --probes probes.fa \
  --sample dengue_1 zika_virus \
  [--ct-values 20 25 30 | --ct 25] \
  [--distractor-fraction-values 0.9 0.99 | --distractor-fraction 0.9] \
  [--hybridization-temperature-values 55 65 70 75 | --hybridization-temperature 70] \
  [--capture-fraction-values 0.3 0.5 0.8 | --capture-fraction 0.5] \
  [--num-sequences-values 100 500 | --num-sequences 500] \
  [--outdir coverage_curve_results] \
  [--cleanup] \
  [... other pipeline parameters ...]
```

`--sample` is **required**. Five parameters can be swept independently or in combination:

| Sweep flag | Fixed flag | Default | Description |
|-----------|------------|---------|-------------|
| `--ct-values 20 25 30` | `--ct 25` | — | CT values (converted to distractor fractions) |
| `--distractor-fraction-values 0.9 0.99` | `--distractor-fraction 0.9` | 0.9 | Distractor fractions directly |
| `--hybridization-temperature-values 55 65 70 75` | `--hybridization-temperature 70` | 70 °C | Hybridisation temperature |
| `--capture-fraction-values 0.3 0.5 0.8` | `--capture-fraction 0.5` | 0.5 | Capture fraction |
| `--num-sequences-values 100 500` | `--num-sequences 500` | all | Number of sequences to sample |

`--ct-values` and `--distractor-fraction-values` are mutually exclusive. Sweep and fixed forms of the same parameter are also mutually exclusive.

The pipeline shares intermediates across combinations: prepare is shared per CT/distractor-fraction; simulate is shared per CT × temperature × capture-fraction.

**Output:**
- Combo subdirectories named by swept params (e.g., `ct_20/`, `ct_20_temp_65_cf_0.50/`)
- `coverage_curve_depth_curves.tsv` — aggregated depth data
- `coverage_curve_report.html` — HTML report with depth curves (`--report full`)
- `coverage_curve_report.Rmd` — editable RMarkdown (`--report rmd`)

---

## build-probes

Build a probe set from target sequences. Runs a multi-step pipeline: collapse redundant targets, construct probes, filter by GC content and sequence complexity, deduplicate. Automatically chains into `assess-probes` unless `--skip-assess` is specified.

```bash
baitbench build-probes \
  --targets targets.fa \
  [--method tile|catch-lite|syotti-lite|catch] \
  [--probe-length 120] \
  [--step -60] \
  [--catch-probe-stride 60] \
  [--catch-mismatches 5] \
  [--catch-extension 0] \
  [--catch-coverage 1.0] \
  [--catch-minhash-threshold 0.6] \
  [--syotti-mismatches 40] \
  [--syotti-seed-len 20] \
  [--no-n-in-probes] \
  [--min-gc 0.20] \
  [--max-gc 0.80] \
  [--max-n-frac 0.05] \
  [--dust-threshold 2.0] \
  [--dust-window 64] \
  [--max-masked-frac 0.25] \
  [--collapse-threshold 0.95] \
  [--dedup-threshold 0.95] \
  [--threads 5] \
  [--genomes genome1.fa ...] \
  [--threshold 80.0] \
  [--skip-assess] \
  [--outdir build_probes_results] \
  [--report full|none|rmd] \
  [--refine-iterations N | --refine-until-stable] \
  [--refine-threshold 80.0]
```

### Pipeline steps

1. **N filter**: Remove target sequences with more than `--max-n-frac` fraction of ambiguous bases.
2. **Collapse**: cd-hit-est clusters targets at `--collapse-threshold` identity.
3. **Build**: Construct probes using the selected method.
4. **N-fix** (if `--no-n-in-probes`): Replace each N in probes with T (or A/C/G if adjacent to T).
5. **GC filter**: Remove probes outside `--min-gc` to `--max-gc`.
6. **Complexity filter**: Remove low-complexity probes via sDUST.
7. **Deduplicate**: cd-hit-est at `--dedup-threshold`.

### Probe design methods

| Method | Flag | Description |
|--------|------|-------------|
| Tiling | `--method tile` (default) | Sliding window of `--probe-length` bp with `--step` controlling overlap/gap |
| catch-lite | `--method catch-lite` | Native Rust reimplementation of CATCH (Metsky et al. 2019) — greedy set-cover |
| syotti-lite | `--method syotti-lite` | Native Rust reimplementation of Syotti (Alanko et al. 2022) — Hamming-distance greedy |
| catch | `--method catch` | External CATCH tool (requires `catch` conda package) |
| probetools-lite | `--method probetools-lite` | Native Rust reimplementation of ProbeTools (Kuchinski et al. 2022) — iterative k-mer clustering. Use `--pt-*` flags. Requires cd-hit-est. |

### Tiling geometry (`--step`)

The stride between consecutive probes is `probe_length + step`:

| `--step` | Stride | Effect |
|---------|--------|--------|
| `-60` (default) | 60 bp | 60 bp overlap (50% overlap with 120 bp probes) |
| `0` | 120 bp | End-to-end, no overlap |
| `10` | 130 bp | 10 bp gap between probes |

A final probe is always anchored to the sequence end regardless of step.

### catch-lite parameters

| Flag | Default | Description |
|------|---------|-------------|
| `--catch-probe-stride` | 60 | Step between candidate probes (bp) |
| `--catch-mismatches` | 5 | Mismatches tolerated for a probe to cover a window |
| `--catch-extension` | 0 | Flanking bp beyond probe boundaries counted as covered |
| `--catch-coverage` | 1.0 | Fraction of each target that must be covered (0.0–1.0) |
| `--catch-minhash-threshold` | 0.6 | Jaccard similarity threshold for near-deduplication; 0.0 disables |

### syotti-lite parameters

| Flag | Default | Description |
|------|---------|-------------|
| `--syotti-mismatches` | 40 | Maximum Hamming distance for a probe to cover a window |
| `--syotti-seed-len` | 20 | K-mer seed length for approximate matching |

### Full parameter table

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--targets` | required | Input target sequences FASTA |
| `--method` | tile | Probe design method: `tile`, `catch-lite`, `syotti-lite`, or `catch` |
| `--probe-length` | 120 | Probe length in bp |
| `--step` | -60 | Tiling step (negative = overlap, 0 = tiled, positive = gap). Tile method only. |
| `--no-n-in-probes` | false | Replace N bases in probes with real nucleotides |
| `--min-gc` | 0.20 | Minimum GC fraction |
| `--max-gc` | 0.80 | Maximum GC fraction |
| `--max-n-frac` | 0.05 | Maximum N fraction to keep a target sequence |
| `--dust-threshold` | 2.0 | sDUST score threshold for low-complexity detection |
| `--dust-window` | 64 | sDUST window size in bases |
| `--max-masked-frac` | 0.25 | Max sDUST-masked fraction to keep a probe. Set to 1.0 to disable. |
| `--collapse-threshold` | 0.95 | cd-hit-est identity for initial target collapse |
| `--dedup-threshold` | 0.95 | cd-hit-est identity for final probe deduplication |
| `--threads` | 5 | Threads for cd-hit-est |
| `--genomes` | none | Genome FASTA(s) for cross-reactivity check (assessment step) |
| `--threshold` | 80.0 | Homology threshold for cross-reactivity (assessment step) |
| `--syotti-mismatches` | 40 | Max Hamming distance for coverage (`syotti-lite` only) |
| `--syotti-seed-len` | 20 | K-mer seed length for approximate matching (`syotti-lite` only) |
| `--pt-step` | 1 | K-mer enumeration step (`probetools-lite` only) |
| `--pt-identity` | 0.9 | cd-hit-est clustering threshold for k-mer clustering (`probetools-lite` only) |
| `--pt-coverage` | 0.9 | Target coverage fraction (10th-percentile) to reach (`probetools-lite` only) |
| `--pt-batch-size` | 100 | Probes added per iteration (`probetools-lite` only) |
| `--pt-max-panel-size` | none | Hard cap on total panel size (`probetools-lite` only) |
| `--pt-min-depth` | 1 | Minimum depth to count as covered (`probetools-lite` only) |
| `--pt-max-iterations` | 20 | Maximum iterations (`probetools-lite` only) |
| `--pt-min-coverage-gain` | 0.001 | Stagnation guard: stop if 10th-percentile improvement falls below this (`probetools-lite` only) |
| `--skip-assess` | false | Skip automatic probe assessment after building |
| `--outdir` | ./build_probes_results | Output directory |
| `--report` | full | Report mode: `full`, `none`, `rmd`, or `both-r` |
| `--cleanup` | false | Delete intermediate files |
| `--refine-iterations` | none | Number of refinement iterations on low-coverage targets |
| `--refine-until-stable` | false | Repeat until no targets remain below the threshold or set stabilises |
| `--refine-threshold` | 80.0 | 1× coverage threshold (%) for refinement |

**Output files:**

| File | Description |
|------|-------------|
| `probes_final.fa` | Final deduplicated probe set |
| `build_probes_stats.tsv` | Sequence/base counts at each pipeline step |
| `assess_probes_report.html` | Combined HTML report (unless `--skip-assess`) |
| `cov_probe_coverage_summary.tsv` | Per-target coverage statistics |
| `cov_probe_depth.tsv` | Probe depth intervals |
| `xreact_hits.tsv` | Cross-reactivity hits (assessment) |
| `xreact_summary.tsv` | Per-probe cross-reactivity summary (assessment) |

---

## assess-probes

Standalone combined probe assessment: probe coverage + cross-reactivity (self-homology always; against genomes if `--genomes` provided).

```bash
baitbench assess-probes \
  --targets targets.fa \
  --probes probes.fa \
  [--genomes genome1.fa ...] \
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
| `--genomes` | none | Genome FASTA(s) for cross-reactivity (repeatable) |
| `--threshold` | 80.0 | Minimum homology % to report cross-reactive hits |
| `--minimap-preset` | sr | Alignment preset |
| `--proximity` | 50 | Pull-down zone distance in bp |
| `--outdir` | ./assess_probes_results | Output directory |
| `--output-prefix` | (empty) | String prepended to every output filename |
| `--report` | full | Report mode |
| `--cleanup` | false | Delete intermediate files |
| `--all-individual-targets` | false | Rerun coverage per target in isolation; adds individual target coverage section |
| `--refine-iterations` | none | Number of refinement iterations |
| `--refine-until-stable` | false | Repeat until stable |
| `--refine-threshold` | 80.0 | 1× coverage threshold for refinement |

**Output files:** `cov_probe_coverage_summary.tsv`, `cov_probe_depth.tsv`, `cov_multi_mapping_probes.tsv`, `xreact_hits.tsv`, `xreact_summary.tsv`, `assess_probes_report.html`.

---

## probe-coverage

Standalone probe coverage analysis tool. Not part of the main pipeline.

```bash
baitbench probe-coverage \
  --targets targets.fa \
  --probes probes.fa \
  [--outdir probe_coverage] \
  [--minimap-preset sr] \
  [--proximity 50] \
  [--report full|none|rmd]
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--targets` | required | Target sequences FASTA |
| `--probes` | required | Probe sequences FASTA |
| `--outdir` | ./probe_coverage | Output directory |
| `--minimap-preset` | sr | Alignment preset |
| `--proximity` | 50 | Pull-down zone distance in bp |
| `--report` | full | Report mode |
| `--cleanup` | false | Delete intermediate files |

**Output files:** `probe_depth.tsv`, `probe_coverage_summary.tsv`, `multi_mapping_probes.tsv`, `probe_coverage_report.html`.

**probe_coverage_summary.tsv columns:**

| Column | Description |
|--------|-------------|
| `reference_id` | Target sequence ID |
| `pct_covered_1x` | % bases with ≥ 1 probe |
| `pct_covered_2x` | % bases with ≥ 2 probes |
| `pct_covered_5x` | % bases with ≥ 5 probes |
| `pct_covered_10x` | % bases with ≥ 10 probes |
| `mean_depth` | Average probe depth across target |
| `median_depth` | Median probe depth |
| `max_gap_length` | Longest uncovered stretch (bp) |
| `num_gaps` | Number of gaps with no probe coverage |
| `pct_near_probe` | % bases within `--proximity` distance of a probe alignment |

---

## xreact

Standalone cross-reactivity analysis. Checks whether probes have high homology to off-target genomes or to each other. Not part of the main pipeline.

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
| `--against` | none | Reference genome FASTA(s) (repeatable) |
| `--self` | false | Check probe-vs-probe cross-reactivity (self-hits excluded) |
| `--threshold` | 80.0 | Minimum homology %: `matching_bases / probe_length × 100` |
| `--minimap-preset` | sr | Alignment preset |
| `--outdir` | ./xreact_results | Output directory |
| `--cleanup` | false | Delete intermediate files |

**Homology metric:** `matching_bases / probe_length × 100`. Captures both alignment identity and query coverage in one number — a probe with 90% identity over 90% of its length scores ~81%.

**Output files:** `hits.tsv`, `summary.tsv`.

**hits.tsv columns:**

| Column | Description |
|--------|-------------|
| `probe_id` | Query probe ID |
| `target_id` | Reference sequence the probe hits |
| `homology_pct` | `matching_bases / probe_length × 100` |
| `identity_pct` | `matching_bases / alignment_block_length × 100` |
| `query_coverage_pct` | `aligned_query_span / probe_length × 100` |
| `matching_bases` | Number of matching bases |
| `alignment_length` | Alignment block length |
| `probe_length` | Total probe length |
| `mode` | `against` or `self` |

**summary.tsv columns:**

| Column | Description |
|--------|-------------|
| `probe_id` | Probe ID |
| `mode` | `against` or `self` |
| `max_homology_pct` | Highest homology across all hits (0.0 if no hits) |
| `best_hit` | Target with highest homology (`NA` if no hits) |
| `num_hits_above_threshold` | Distinct alignments above threshold |

---

## panel-qc

Standalone pre-experiment QC: assesses whether a target panel can discriminate between species by computing target-vs-target similarity.

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
| `--sample-target-map` | required | Mapping TSV linking species IDs to target IDs |
| `--identity-threshold` | 90.0 | Minimum sequence identity % to consider two targets "similar" |
| `--minimap-preset` | sr | Alignment preset |
| `--outdir` | ./panel_qc_results | Output directory |
| `--report` | full | Report mode |
| `--cleanup` | false | Delete intermediate files |

**Algorithm:**

1. All-vs-all target alignment
2. Pairs above `--identity-threshold` are marked "similar"
3. Targets are classified as "unique" (no cross-species similarity) or "shared"
4. Per-species discriminability score = `unique_targets / total_targets`

**Output files:** `target_similarity.tsv`, `species_discriminability.tsv`, `species_confusion_matrix.tsv`, `panel_qc_report.html`.

**target_similarity.tsv columns:**

| Column | Description |
|--------|-------------|
| `target_a` | First target ID |
| `target_b` | Second target ID |
| `identity_pct` | `matching_bases / min(len_a, len_b) × 100` |
| `matching_bases` | Number of matching bases |
| `len_a` | Length of target A |
| `len_b` | Length of target B |

**species_discriminability.tsv columns:**

| Column | Description |
|--------|-------------|
| `species_id` | Species/genome ID |
| `total_targets` | Total targets for this species |
| `unique_targets` | Targets unique to this species |
| `shared_targets` | Targets similar to those in other species |
| `discriminability_score` | `unique_targets / total_targets` (0.0–1.0) |
| `confusable_species` | Species with shared targets (comma-separated) |

---

## identify

Call species PRESENT/ABSENT/AMBIGUOUS from multi-target detection patterns. Can be run standalone or integrated into `baitbench run` with `--identify`.

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
  [--min-unique-targets 1] \
  [--outdir identify_results]
```

Either `--target-similarity` or `--targets` must be provided (not both).

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--detected-detail` | required | `detected_detail.tsv` from metrics step |
| `--sample-target-map` | required | Mapping TSV linking species IDs to target IDs |
| `--target-similarity` | none | Pre-computed similarity TSV from `panel-qc` |
| `--targets` | none | Target sequences FASTA (computes similarity on-the-fly) |
| `--identity-threshold` | 90.0 | Similarity threshold (only with `--targets`) |
| `--min-unique-targets` | 1 | Minimum unique target detections to call a species PRESENT |
| `--outdir` | ./identify_results | Output directory |

**Calling algorithm:**

Species are ranked by evidence strength (unique targets detected, then total reads). For each species in order:
- **PRESENT**: ≥ `--min-unique-targets` unique targets detected
- **ABSENT** (no detections): zero targets detected
- **ABSENT** (explained): all detected targets are shared AND every one is explained by a species already called PRESENT
- **AMBIGUOUS** (no unique markers): species has zero unique targets — cannot confirm or deny
- **AMBIGUOUS** (insufficient evidence): some unique targets but fewer than threshold; not all shared hits explained

**Output files:** `species_calls.tsv`, `species_calls.json`.

**species_calls.tsv columns:**

| Column | Description |
|--------|-------------|
| `species_id` | Species/genome ID |
| `call` | PRESENT, ABSENT, or AMBIGUOUS |
| `total_targets` | Total targets for this species |
| `unique_targets` | Targets unique to this species |
| `shared_targets` | Targets shared with other species |
| `unique_detected` | Unique targets that were detected |
| `shared_detected` | Shared targets that were detected |
| `total_detected` | Total targets detected |
| `total_reads` | Total reads across all detected targets |
| `explained_by` | Species that explain shared hits (comma-separated) |
| `reason` | `unique_markers_detected`, `no_detections`, `cross_reactivity_explained`, `no_unique_markers`, or `insufficient_unique_evidence` |

---

## tool

Standalone utility tools grouped under a single subcommand. Run `baitbench tool --help` to list available tools.

```bash
baitbench tool <TOOL> [OPTIONS]
```

### tool syotti

Run the Syotti greedy bait design algorithm directly.

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
| `--probe-length` | 120 | Probe length in bp |
| `--mismatches` | 40 | Maximum Hamming distance for coverage |
| `--seed-len` | 20 | K-mer seed length for approximate matching |

### tool catch

Run the CATCH optimization probe design algorithm directly.

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
| `--probe-length` | 120 | Probe length in bp |
| `--stride` | 60 | Tiling step in bp |
| `--mismatches` | 5 | Maximum mismatches for a probe to cover a window |
| `--extension` | 0 | Extension length on each side of candidate probes |
| `--coverage` | 1.0 | Fraction of each target that must be covered |
| `--minhash-threshold` | 0.6 | MinHash Jaccard similarity threshold for deduplication |

### tool dustview

Visualize sDUST low-complexity masking on FASTA sequences. Outputs original sequence, masked sequence (X marks low-complexity regions), and per-sequence statistics to stdout.

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

### tool collapse

Cluster near-duplicate sequences using cd-hit-est and write cluster representatives.

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
| `--input` | required | Input FASTA |
| `--output` | required | Output FASTA (cluster representatives) |
| `--threshold` | 0.95 | Identity threshold for clustering |
| `--threads` | 1 | Threads for cd-hit-est |
| `--log-file` | cdhit.log | Path for cd-hit-est log |
