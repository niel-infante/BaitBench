# Commands

## run

Runs the complete pipeline from input files to metrics and report.

```bash
baitbench run [OPTIONS]
```

This is the primary command for most use cases. It chains all pipeline steps (prepare through report) automatically. Use `--cleanup` to delete intermediate files (FASTA, SAM, logs) after completion, keeping only report inputs and final outputs. See [Parameter Reference](parameters.md) for all options.

In genome mode with `--sample-target-map`, use `--identify` to add species-level calling after metrics. This computes target similarity, calls species PRESENT/ABSENT/AMBIGUOUS, and includes the results in the HTML report with ground-truth comparison against the `--sample` manifest.

## prepare

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

## simulate

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

Fragments shorter than `--read-length` are kept as-is. With `--num-sequences`, reads are sampled with replacement from the fragment pool (modeling PCR amplification before sequencing) and given unique IDs.

**Output files:**
- `reads.fa` -- sequencing reads

## filter

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

## map

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

## list

Counts reads per reference from a SAM file.

```bash
baitbench list \
  --sam mapped.sam \
  --output detected.list
```

**Output files:**
- `detected.list` -- TSV: `reference_id\tcount` (sorted ascending by count)

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

**Output files:**
- `results.tsv` -- genome-level summary metrics
- `detected_detail.tsv` -- per-reference detection and coverage detail
- `results.json` -- structured JSON output (optional)
- `coverage.tsv` -- run-length encoded read depth intervals (optional)

## report

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

## probe-coverage

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

## xreact

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

## coverage-curve

Runs the pipeline at multiple parameter combinations and generates coverage depth curves.

```bash
baitbench coverage-curve \
  --targets targets.fa \
  --distractors distractors.fa \
  --probes probes.fa \
  --sample dengue_1 zika_virus \
  [--ct-values 20 25 30 | --ct 25] \
  [--capture-fraction-values 0.3 0.5 0.8 | --capture-fraction 0.5] \
  [--num-sequences-values 100 500 | --num-sequences 500] \
  [--outdir coverage_curve_results] \
  [--cleanup] \
  [... other pipeline parameters ...]
```

Three parameters can be swept (each has a singular fixed form and a plural sweep form):

| Sweep flag | Fixed flag | Description |
|-----------|------------|-------------|
| `--ct-values 20 25 30` | `--ct 25` | CT values |
| `--capture-fraction-values 0.3 0.5 0.8` | `--capture-fraction 0.5` | Capture fraction (probe-biased fragment proportion) |
| `--num-sequences-values 100 500` | `--num-sequences 500` | Number of sequences to sample |

Sweep and fixed forms of the same parameter are mutually exclusive. `--ct-values` and `--distractor-fraction` are also mutually exclusive.

`--sample` is **required** for coverage-curve (must specify which targets to track).

The pipeline shares intermediate files across combinations for efficiency: prepare is shared per CT value; simulate is shared per CT x capture-fraction combination.

**Output files:**
- Combo subdirectories named by swept params (e.g., `ct_20/`, `ct_20_cf_0.50/`, `ct_20_cf_0.50_ns_500/`)
- `coverage_curve_depth_curves.tsv` -- aggregated depth data (columns: ct, capture_fraction, num_sequences, ...)
- `coverage_curve_report.html` -- HTML report with depth curves (`--report full`)
- `coverage_curve_report.Rmd` -- editable RMarkdown file (`--report rmd`)

## panel-qc

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

## identify

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

## build-probes

Build a probe set from target sequences. Runs a multi-step pipeline: collapse redundant targets, construct probes, filter by GC content and sequence complexity, and deduplicate. After building, automatically chains into probe assessment (probe coverage + cross-reactivity analysis) unless `--skip-assess` is specified.

Four probe construction methods are available: `tile` (sliding window, default), `catch-lite` (native Rust reimplementation of CATCH optimization-based design), `catch` (external CATCH tool from the Broad Institute; requires the `catch` conda package), and `syotti-lite` (native Rust reimplementation of Syotti greedy set-cover design).

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

1. **N filter**: Remove target sequences with more than `--max-n-frac` fraction of ambiguous (non-ACGT) bases.
2. **Collapse**: cd-hit-est clusters targets at `--collapse-threshold` identity to remove near-duplicates
3. **Build**: Construct probes from collapsed sequences. Method `tile` generates sliding-window probes of `--probe-length` bp across each sequence with `--step` controlling overlap/gap. A final probe is anchored to the end of each sequence to ensure full coverage. Method `catch-lite` uses BaitBench's native Rust reimplementation of CATCH. Method `catch` calls the external CATCH tool (`design_probes.py`). Method `syotti-lite` uses BaitBench's native Rust reimplementation of the Syotti greedy set-cover algorithm.
4. **GC filter**: Remove probes with GC content outside `--min-gc` to `--max-gc` range
5. **Complexity filter**: Remove low-complexity probes using the sDUST algorithm (Morgulis et al. 2006). Probes where more than `--max-masked-frac` of bases are identified as low-complexity are removed. Set `--max-masked-frac 1.0` to disable.
6. **Deduplicate**: cd-hit-est clusters probes at `--dedup-threshold` identity to remove redundant probes

**Tiling geometry (`--step`):**

The stride between consecutive probes is `probe_length + step`. The step is measured from the end of the previous probe:

- `--step -60` (default): stride = 60, probes overlap by 60bp (50% overlap with 120bp probes)
- `--step 0`: stride = 120, probes are perfectly tiled (no overlap, no gap)
- `--step 10`: stride = 130, 10bp gap between probes

Probes are named `probe_{target_id}|tile_{n}`. A final probe is always placed at the sequence end regardless of overlap.

**catch-lite method (`--method catch-lite`):**

BaitBench includes a native Rust reimplementation of the CATCH algorithm (Metsky et al. 2019, Nature Biotechnology). Unlike tiling, CATCH minimizes the number of probes needed while guaranteeing a configurable fraction of each target sequence is covered. The algorithm tiles candidate probes at a configurable stride, removes near-duplicates via MinHash LSH, then runs a greedy set-cover to select the minimum probe set that covers all targets to the required depth.

| Flag | Default | Description |
|------|---------|-------------|
| `--catch-probe-stride` | 60 | Step between candidate probes (bp) |
| `--catch-mismatches` | 5 | Mismatches tolerated for a probe to cover a target window |
| `--catch-extension` | 0 | Flanking bp beyond probe boundaries counted as covered |
| `--catch-coverage` | 1.0 | Fraction of each target that must be covered (0.0–1.0) |
| `--catch-minhash-threshold` | 0.6 | Jaccard similarity threshold for near-deduplication; set to 0.0 to disable |

Probes are named `probe_{source_id}|catch_{n}`.

**catch method (`--method catch`):**

Calls the external CATCH tool (`design.py`) from the Broad Institute. Requires the `catch` conda package (`conda install -c bioconda catch`). All `--catch-*` flags apply.

**syotti-lite method (`--method syotti-lite`):**

[Syotti](https://github.com/jnalanko/syotti) (Alanko et al. 2022) is a greedy set-cover bait designer. It scans the input sequences; at every uncovered position, it extracts a bait of `--probe-length` bp and marks all reference windows within `--syotti-mismatches` Hamming distance as covered (checking both strands). This is more targeted than tiling — probes are only generated where coverage is not already achieved by an earlier probe, yielding a smaller set while guaranteeing full coverage.

The BaitBench implementation uses a k-mer hash index (no external dependencies required).

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--targets` | required | Input target sequences FASTA |
| `--method` | tile | Probe construction method: `tile`, `catch-lite`, `syotti-lite`, or `catch` |
| `--probe-length` | 120 | Probe length in bp |
| `--step` | -60 | Step from end of previous probe. Negative = overlap, 0 = tiled, positive = gap. Only used with `--method tile`. |
| `--catch-probe-stride` | 60 | Step between candidate probes (bp). Used with `--method catch-lite` and `--method catch`. |
| `--catch-mismatches` | 5 | Mismatches tolerated for a probe to cover a target window. Used with `--method catch-lite` and `--method catch`. |
| `--catch-extension` | 0 | Flanking bp beyond probe boundaries counted as covered. Used with `--method catch-lite` and `--method catch`. |
| `--catch-coverage` | 1.0 | Fraction of each target that must be covered (0.0–1.0). Used with `--method catch-lite` and `--method catch`. |
| `--catch-minhash-threshold` | 0.6 | Jaccard similarity threshold for near-deduplication; 0.0 disables. Used with `--method catch-lite` and `--method catch`. |
| `--syotti-mismatches` | 40 | Maximum Hamming distance for a bait to cover a reference window. Only used with `--method syotti-lite`. |
| `--syotti-seed-len` | 20 | K-mer seed length for Syotti approximate matching. Only used with `--method syotti-lite`. |
| `--min-gc` | 0.20 | Minimum GC fraction (0–1) |
| `--max-gc` | 0.80 | Maximum GC fraction (0–1) |
| `--max-n-frac` | 0.05 | Maximum fraction of ambiguous (non-ACGT) bases in a target sequence (0–1) |
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
| `--refine-iterations` | none | Number of refinement iterations on low-coverage targets |
| `--refine-until-stable` | false | Repeat refinement until no targets remain below the threshold or set stabilizes |
| `--refine-threshold` | 80.0 | 1X coverage threshold (%) for refinement iterations |

**Output files:**

- `probes_final.fa` -- final deduplicated probe set
- `build_probes_stats.tsv` -- sequence/base counts at each pipeline step
- `assess_probes_report.html` -- combined HTML report with build stats + assessment (unless `--skip-assess`)
- `cov_probe_coverage_summary.tsv` -- per-target coverage statistics (from assessment)
- `cov_probe_depth.tsv` -- probe depth intervals (from assessment)
- `xreact_hits.tsv` -- cross-reactivity hits (from assessment)
- `xreact_summary.tsv` -- cross-reactivity summary (from assessment)

With `--skip-assess`, only produces `probes_final.fa`, `build_probes_stats.tsv`, and optionally `build_probes_report.html`.

## tool

Standalone utility tools grouped under a single subcommand. Run `baitbench tool --help` to list available tools.

```bash
baitbench tool <TOOL> [OPTIONS]
```

### tool syotti

Run the Syotti greedy bait design algorithm directly, without the `build-probes` pipeline.

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
| `--seed-len` | 20 | K-mer seed length for approximate matching. |

### tool catch

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

### tool dustview

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

### tool collapse

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

## assess-probes

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
| `--refine-iterations` | none | Number of refinement iterations (mutually exclusive with `--refine-until-stable`) |
| `--refine-until-stable` | false | Repeat refinement until no targets remain below the threshold or the set stops changing |
| `--refine-threshold` | 80.0 | 1X coverage threshold (%) used to identify low-coverage targets for refinement |

**Output files:**

- `cov_probe_coverage_summary.tsv` -- per-target coverage statistics
- `cov_probe_depth.tsv` -- run-length encoded probe depth intervals
- `cov_multi_mapping_probes.tsv` -- probes mapping to multiple targets
- `xreact_hits.tsv` -- cross-reactivity hits above threshold
- `xreact_summary.tsv` -- per-probe cross-reactivity summary
- `assess_run_params.tsv` -- run parameters
- `assess_probes_report.html` -- combined HTML report (`--report full`)
- `assess_probes_report.Rmd` -- editable RMarkdown file (`--report rmd`)
- `refine_N_targets.fa` -- filtered targets for refinement iteration N
- `refine_N_cov_probe_coverage_summary.tsv` -- coverage statistics for refinement iteration N
- `refine_N_probe_coverage_report.html` -- probe coverage report for refinement iteration N

**Report sections:**

1. **Probe Coverage** -- summary table, coverage breadth bar charts, tiered coverage, gap analysis, depth profiles, proximity coverage, multi-mapping probes
2. **Self-Homology** -- heatmap (≤1000 probes), density plots, hits table
3. **Cross-Reactivity vs Genomes** (if `--genomes` provided) -- heatmap, per-genome bar chart, density plots, hits table
4. **Parameters** -- run configuration under a collapsible fold

**Refinement iterations:**

Many target panels contain highly similar sequences that are unlikely to occur together in the same sample. Refinement iterations address this by re-running probe coverage on only the targets that showed poor coverage (below `--refine-threshold`), so you can assess how well the probes cover each subset in isolation.

- **`--refine-iterations N`** runs exactly N additional probe-coverage-only analyses after the initial full assessment.
- **`--refine-until-stable`** repeats automatically until no targets fall below the threshold, or until the set of low-coverage targets stops changing.
- **`--refine-threshold`** (default 80.0) sets the 1X coverage percentage below which a target is considered poorly covered.
