# Analyze Cross-Reactivity

Three tools address different aspects of cross-reactivity:

| Tool | Question | When to use |
|------|----------|-------------|
| `xreact` | Do probes bind to off-target sequences? | Evaluating probe specificity against host or environmental DNA |
| `panel-qc` | Can the target panel distinguish between species? | Before running simulations; pre-experiment QC |
| `identify` | Which species are present, accounting for cross-reactivity? | Calling species from detection results |

---

## xreact: probe-to-genome and probe-to-probe

### Probe cross-reactivity against external genomes

Check how much your probes bind to sequences outside your target panel:

```bash
baitbench xreact \
  --probes probes.fa \
  --against human_genome.fa \
  --threshold 80 \
  --outdir xreact_results
```

Check multiple organisms at once:

```bash
baitbench xreact \
  --probes probes.fa \
  --against human.fa mosquito.fa bacteria.fa \
  --threshold 80 \
  --outdir xreact_results
```

**The homology metric** is `matching_bases / probe_length × 100`. A probe with 90% identity over 90% of its length scores ~81% — this single number captures both alignment identity and query coverage.

**`xreact_hits.tsv`** — every alignment above the threshold:

```
probe_id        target_id     homology_pct  identity_pct  query_coverage_pct
probe_flu_001   Human_chr1    82.5          91.7          90.0
probe_flu_015   Human_chr22   78.2          86.9          90.0
```

**`xreact_summary.tsv`** — worst-case per probe:

```
probe_id        mode     max_homology_pct  best_hit       num_hits_above_threshold
probe_flu_001   against  82.5              Human_chr1     1
probe_flu_002   against  0.0              NA              0
```

Probes with `max_homology_pct` above your threshold are likely to produce false positive reads in real experiments. Consider removing or redesigning them.

### Probe self-homology

Check whether probes are too similar to each other (which can cause cross-capture between panel members):

```bash
baitbench xreact \
  --probes probes.fa \
  --self \
  --threshold 80 \
  --outdir self_homology
```

Self-hits (probe A vs itself) are automatically excluded. High self-homology between probes from different targets indicates they will compete for the same binding site — useful for diagnosing FP_target detections in simulations.

### Combined probe assessment

You can run both `--against` and `--self` together. Or use `assess-probes`, which combines probe coverage and cross-reactivity in one report:

```bash
baitbench assess-probes \
  --targets targets.fa \
  --probes probes.fa \
  --genomes human.fa \
  --threshold 80 \
  --outdir assessment
```

---

## panel-qc: target panel discriminability

`panel-qc` asks whether your target sequences are distinct enough to identify species uniquely. It does this before any simulation — purely from sequence similarity.

```bash
baitbench panel-qc \
  --targets gene_targets.fa \
  --sample-target-map mapping.tsv \
  --identity-threshold 90 \
  --outdir panel_qc
```

`--sample-target-map` links species/genome IDs to their target sequence IDs (same format as in genome mode).

**What it computes:**

1. Aligns all targets against all targets (all-vs-all)
2. Marks pairs with ≥ `--identity-threshold` % similarity as "shared"
3. For each species, counts how many of its targets are unique (no cross-species similarity) vs shared
4. Reports a discriminability score: `unique_targets / total_targets`

**`species_discriminability.tsv`:**

```
species_id              total_targets  unique_targets  shared_targets  discriminability_score  confusable_species
Dengue_virus_1          3              2               1               0.667                   Dengue_virus_2
Dengue_virus_2          3              2               1               0.667                   Dengue_virus_1
Zika_virus              2              2               0               1.000
```

A discriminability score of 0.0 means the species has **zero unique targets** — it cannot be reliably distinguished from others using this panel. Add more targets for that species, or redesign targets in a less conserved gene region.

**`species_confusion_matrix.tsv`** shows which species pairs share targets — a direct readout of which identifications are ambiguous.

---

## identify: calling species from detection results

`identify` converts a `detected_detail.tsv` from the simulation pipeline into species-level PRESENT/ABSENT/AMBIGUOUS calls.

### Standalone on existing results

```bash
# Using target similarity pre-computed by panel-qc
baitbench identify \
  --detected-detail results/run/detected_detail.tsv \
  --sample-target-map mapping.tsv \
  --target-similarity panel_qc/target_similarity.tsv \
  --outdir identify_results

# Computing similarity on-the-fly (slower; skips panel-qc step)
baitbench identify \
  --detected-detail results/run/detected_detail.tsv \
  --sample-target-map mapping.tsv \
  --targets gene_targets.fa \
  --identity-threshold 90 \
  --outdir identify_results
```

### Integrated into the pipeline

Add `--identify` to `baitbench run` (requires genome mode with `--sample-target-map`):

```bash
baitbench run \
  --targets targets.fa \
  --genomes genomes.fa \
  --distractors distractors.fa \
  --probes probes.fa \
  --sample-target-map mapping.tsv \
  --sample species_a species_b \
  --identify \
  --outdir results
```

Species calls are included in the HTML report with ground-truth comparison.

### How the calling algorithm works

Species are processed in order of evidence strength (most unique target detections first). For each species:

- **PRESENT**: ≥ `--min-unique-targets` unique targets detected
- **ABSENT** (no detections): zero targets detected
- **ABSENT** (explained): all detected targets are "shared" AND every shared detection is explained by a species already called PRESENT
- **AMBIGUOUS** (no unique markers): species has zero unique targets in the panel — cannot confirm or deny
- **AMBIGUOUS** (insufficient evidence): has unique targets but not enough detected, and not all shared detections are explained

The ordered approach creates a parsimony effect: the species with the strongest unique evidence is called first, then its presence explains away shared target hits in weaker candidates.

**`species_calls.tsv`:**

```
species_id          call      unique_detected  shared_detected  total_reads  explained_by  reason
Dengue_virus_1      PRESENT   2               1                1450         —             unique_markers_detected
Dengue_virus_2      ABSENT    0               1                320          Dengue_virus_1  cross_reactivity_explained
Zika_virus          PRESENT   2               0                890          —             unique_markers_detected
```

### Tuning the calling threshold

`--min-unique-targets` (default 1) sets how many unique targets must be detected to call a species PRESENT. Increase it to reduce false positives at the cost of sensitivity:

```bash
--min-unique-targets 2   # require at least 2 unique targets
```

---

## Typical cross-reactivity workflow

```bash
# 1. Pre-experiment: assess panel discriminability
baitbench panel-qc \
  --targets targets.fa \
  --sample-target-map mapping.tsv \
  --outdir panel_qc
# Review: panel_qc/panel_qc_report.html

# 2. Check probe cross-reactivity against host sequences
baitbench xreact \
  --probes probes.fa \
  --against host.fa \
  --threshold 80 \
  --outdir xreact

# 3. Run simulation
baitbench run \
  --targets targets.fa --genomes genomes.fa \
  --distractors host.fa \
  --probes probes.fa \
  --sample-target-map mapping.tsv \
  --sample org_a org_b \
  --outdir sim_results

# 4. Call species, using pre-computed similarity
baitbench identify \
  --detected-detail sim_results/*/detected_detail.tsv \
  --sample-target-map mapping.tsv \
  --target-similarity panel_qc/target_similarity.tsv \
  --outdir species_calls
```
