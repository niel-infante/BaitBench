# Design and Assess Probes

This guide covers building a probe set from target sequences and evaluating it with `assess-probes`.

---

## Overview

```
Target sequences
       │
  build-probes    ← collapse → design → GC filter → complexity filter → deduplicate
       │
  probes_final.fa
       │
  assess-probes   ← probe coverage + cross-reactivity → combined report
       │
  assess_probes_report.html
```

`build-probes` automatically runs `assess-probes` at the end unless you pass `--skip-assess`.

---

## Choosing a design method

Four methods are available:

| Method | Flag | Approach | Best for |
|--------|------|----------|----------|
| Tiling | `--method tile` (default) | Sliding window across each target | General use; predictable, uniform coverage |
| catch-lite | `--method catch-lite` | Greedy set-cover, minimises probe count | Reducing panel size while guaranteeing coverage |
| syotti-lite | `--method syotti-lite` | Greedy set-cover with Hamming distance | Panels with high mismatch tolerance |
| catch | `--method catch` | External CATCH tool (Broad Institute) | Reproducing published CATCH designs |

**Start with tiling.** It is the simplest and most predictable. Move to catch-lite or syotti-lite if you need to minimise the number of probes while maintaining coverage.

---

## Basic tiling design

```bash
baitbench build-probes \
  --targets targets.fa \
  --probe-length 120 \
  --outdir probes_output
```

Default step is `--step -60` (60 bp overlap = 50% overlap for 120 bp probes). Increase overlap to improve robustness to sequence variation; reduce it to cut probe count.

```bash
--step -60    # 60 bp overlap (default, 2× coverage per position)
--step -30    # 30 bp overlap (1.25× coverage)
--step 0      # no overlap, tiled end-to-end
--step 10     # 10 bp gap between probes
```

---

## catch-lite design

```bash
baitbench build-probes \
  --targets targets.fa \
  --method catch-lite \
  --probe-length 120 \
  --catch-probe-stride 60 \
  --catch-mismatches 5 \
  --catch-coverage 1.0 \
  --outdir probes_output
```

`--catch-coverage 1.0` requires 100% of each target to be covered; lower it (e.g., 0.9) to reduce probe count at the cost of gaps.

---

## syotti-lite design

```bash
baitbench build-probes \
  --targets targets.fa \
  --method syotti-lite \
  --probe-length 120 \
  --syotti-mismatches 40 \
  --outdir probes_output
```

`--syotti-mismatches` sets the Hamming distance within which a probe is considered to "cover" a reference window. Higher values generate fewer probes but with less stringent matching.

---

## Filtering options

All methods apply the same post-design filters:

| Filter | Flags | Default | Purpose |
|--------|-------|---------|---------|
| GC content | `--min-gc` / `--max-gc` | 0.20 / 0.80 | Remove probes that won't hybridize well |
| Complexity | `--dust-threshold`, `--max-masked-frac` | 2.0, 0.25 | Remove low-complexity probes (repetitive sequences) |
| Deduplication | `--dedup-threshold` | 0.95 | Remove near-identical probes |

Tighten GC bounds for demanding hybridization conditions:

```bash
--min-gc 0.40 --max-gc 0.60   # stricter GC for high-stringency capture
```

---

## Checking for N bases

If your targets contain ambiguous (N) bases, probes derived from them will also contain Ns. Use `--no-n-in-probes` to replace each N with a real nucleotide (T preferred, or A/C/G if T is adjacent):

```bash
baitbench build-probes \
  --targets targets.fa \
  --no-n-in-probes \
  --outdir probes_output
```

Without this flag, probes with Ns pass through normally. The downstream simulation handles Ns in probe sequences but probe efficiency may be reduced.

---

## Running assess-probes on an existing probe set

To assess a probe set you already have (or one built outside BaitBench):

```bash
baitbench assess-probes \
  --targets targets.fa \
  --probes probes.fa \
  --outdir assess_results
```

To also check cross-reactivity against off-target genomes:

```bash
baitbench assess-probes \
  --targets targets.fa \
  --probes probes.fa \
  --genomes human_genome.fa \
  --threshold 80 \
  --outdir assess_results
```

`--threshold 80` reports any probe with ≥80% homology to the human genome (matching bases / probe length × 100).

---

## Reading the coverage report

The assessment report shows probe depth across each target. Key statistics in `cov_probe_coverage_summary.tsv`:

| Column | Meaning |
|--------|---------|
| `pct_covered_1x` | % of target bases with at least one probe |
| `pct_covered_2x` | % with at least two probes (robustness to one probe failure) |
| `mean_depth` | Average probe depth across the target |
| `max_gap_length` | Longest uncovered stretch (bp) — flag if > probe length |
| `num_gaps` | Number of gaps with no probe coverage |

A well-designed panel should have `pct_covered_1x` close to 100% and `max_gap_length` below the read length (120 bp by default).

---

## Refinement for highly similar targets

If your panel contains many similar strain variants, probes designed for variant A may fail to cover unique regions of variant B because they preferentially align to A. Use `--refine-until-stable` to iterate:

```bash
baitbench assess-probes \
  --targets targets.fa \
  --probes probes.fa \
  --refine-until-stable \
  --refine-threshold 80 \
  --outdir assess_results
```

Each refinement iteration:
1. Identifies targets with < `--refine-threshold`% coverage at 1×
2. Re-runs probe coverage on those targets alone (removing competition from similar sequences)
3. Repeats until no targets fall below the threshold or the set stabilises

The report includes one coverage section per iteration. Low-coverage targets in later iterations indicate genuine gaps in the probe panel that require additional probes.

---

## Interpreting cross-reactivity results

`xreact_hits.tsv` lists every probe-to-genome alignment above the homology threshold:

| Column | Meaning |
|--------|---------|
| `probe_id` | The probe that cross-reacts |
| `target_id` | The off-target genome it hits |
| `homology_pct` | matching_bases / probe_length × 100 |

`xreact_summary.tsv` gives the worst-case homology per probe:

```
probe_id           mode     max_homology_pct  best_hit       num_hits_above_threshold
probe_flu_001      against  42.5              Human_chr1     1
probe_flu_002      against  0.0               NA             0
```

Probes with `max_homology_pct` above 80% are likely to produce false positive reads in real data. Consider redesigning those probes or filtering them from the panel.

---

## Typical design workflow

```bash
# 1. Build probes (tiling, with assessment)
baitbench build-probes \
  --targets targets.fa \
  --probe-length 120 \
  --outdir probes_v1

# 2. Review the report: probes_v1/assess_probes_report.html
#    Check coverage gaps and cross-reactivity

# 3. If coverage gaps exist, try catch-lite with higher coverage requirement
baitbench build-probes \
  --targets targets.fa \
  --method catch-lite \
  --catch-coverage 1.0 \
  --outdir probes_v2

# 4. If cross-reactivity against host is too high, tighten GC or filter specific probes
baitbench build-probes \
  --targets targets.fa \
  --genomes human.fa \
  --threshold 70 \
  --min-gc 0.40 --max-gc 0.65 \
  --outdir probes_v3

# 5. Run the simulation to evaluate actual capture performance
baitbench run \
  --targets targets.fa \
  --distractors human.fa \
  --probes probes_v3/probes_final.fa \
  --num-fragments 10000 \
  --outdir sim_results
```
