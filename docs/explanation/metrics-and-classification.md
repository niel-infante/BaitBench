# Metrics and Classification

BaitBench uses a 3-way classification system that distinguishes between different types of false positives and tracks detections at both the genome (or group) level and the read level. This page explains the design rationale and how to interpret each metric.

---

## The 3-Way Classification System

A typical 2-way classification (TP/FP/FN/TN) conflates two different types of false positive: detecting a sequence that is in the probe panel but not in the specimen, and detecting a background (distractor) sequence that is outside the panel entirely. These have different causes and implications:

- **FP_target**: A non-sample target was detected. This means probes designed for a panel member accidentally captured a sequence that was not in the specimen. Caused by cross-reactivity between panel members (similar sequences sharing probes).
- **FP_distractor**: A distractor (background, off-target) sequence was detected. This means probes bind to background DNA. Caused by probe similarity to host or environmental sequences.

BaitBench tracks both separately so you can diagnose whether specificity failures are within-panel (probe redesign needed) or against background (probe filtering or host-depletion needed).

---

## Classification Table

A detection decision is made for every genome (or group) in the analysis. A genome is "detected" if at least one read maps to it after capture and mapping.

| Category | In sample? | Detected? | Classification |
|----------|-----------|-----------|----------------|
| Sample target | Yes | Yes | **TP** |
| Sample target | Yes | No | **FN** |
| Non-sample target | No | Yes | **FP_target** |
| Non-sample target | No | No | **TN_target** |
| Distractor | No | Yes | **FP_distractor** |
| Distractor | No | No | **TN_distractor** |
| Untargeted genome *(genome mode only)* | — | — | untargeted |

Without `--sample`, all targets are treated as in the sample, reducing to the traditional 2-way classification (only TP/FN; no FP_target or TN_target).

Untargeted genomes (sample genomes with no corresponding target in the sample-target-map) are tracked separately and excluded from TP/FP/FN/TN counts.

---

## Genome-Level vs Read-Level Metrics

BaitBench reports metrics at two levels:

### Genome level (group_detail.tsv, results.tsv)

Each organism (or distractor group) is called detected or not. The threshold is any reads at all: ≥ 1 read mapped = detected.

This is the primary evaluation metric because it answers the question the probe panel is designed to answer: "Is this pathogen present in the sample?"

### Read level (results.tsv)

Read-level metrics track how individual reads behave:

| Metric | What it counts |
|--------|----------------|
| `sample_captured` | Fragments from sample targets that passed the capture filter |
| `nonsample_target_captured` | Fragments from non-sample targets (should be 0 — they have weight 0) |
| `distractor_captured` | Fragments from distractor sequences that passed the capture filter |
| `reads_correctly_mapped` | Reads that mapped back to their source reference |
| `reads_incorrectly_mapped` | Reads that mapped to a different reference |

`reads_incorrectly_mapped` is a sensitive indicator of cross-mapping: even when genome-level metrics show all TP and no FP, a panel may have high cross-mapping, where reads from virus A map to virus B. This doesn't affect genome-level detection (since both A and B may still be correctly detected) but indicates the panel cannot quantify individual strain abundances reliably.

---

## Sensitivity (Recall)

```
sensitivity = TP / (TP + FN)
```

What fraction of sample targets were detected?

- **1.0**: Every sample target had at least one read mapping to it — perfect recall.
- **0.5**: Half the sample targets were missed — probes may not cover those targets, or target abundance is too low.
- **0.0**: No sample targets were detected — likely a configuration problem (wrong probes, mismatched IDs).

Low sensitivity points to probe coverage gaps, mismatches, or insufficient sequencing depth.

---

## Specificity

```
specificity = (TN_target + TN_distractor) / (TN_target + TN_distractor + FP_target + FP_distractor)
```

What fraction of non-sample entities were correctly not detected?

- **1.0**: No false positives of any kind.
- **0.0**: Every non-sample entity was detected.

Specificity is reported separately for targets and distractors in `results.tsv` so you can distinguish within-panel cross-reactivity from background noise:

```
specificity_target:     TN_target / (TN_target + FP_target)
specificity_distractor: TN_distractor / (TN_distractor + FP_distractor)
```

---

## Precision

```
precision = TP / (TP + FP_target + FP_distractor)
```

Of all detected entities, what fraction are genuine sample targets?

High sensitivity with low precision means you're detecting your targets but also picking up a lot of noise. High precision with low sensitivity means your detections are reliable but you're missing some targets.

---

## F1 Score

```
F1 = 2 × (precision × sensitivity) / (precision + sensitivity)
```

The harmonic mean of precision and sensitivity. The harmonic mean penalises extreme imbalance: a test that detects everything (sensitivity = 1.0) but has very low precision still gets a low F1.

F1 is useful when you want a single number to compare probe panels, but it flattens the FP_target vs FP_distractor distinction. For detailed diagnostics, look at `group_detail.tsv` directly.

---

## Groups and How They Affect Classification

By default, each individual sequence (target or distractor contig) is its own entity for classification. This can produce misleading results when:
- Multiple variant sequences of the same organism are in the targets FASTA (e.g., 10 dengue strains) — detecting any 1 out of 10 would count as 1 TP
- Multiple contigs of a distractor genome are in the distractors FASTA — each contig would be its own FP_distractor

`--groups` solves the first problem: a groups TSV maps sequence IDs to group names, and classification happens at the group level. A group is detected if any of its member sequences has ≥ 1 read.

For distractors, all contigs from the same FASTA file are automatically grouped by file stem (e.g., all contigs from `Aaegypti.fa` → one group `"Aaegypti"`). `--distractor-groups` overrides this with an explicit mapping.

---

## Genome Mode Classification

In genome mode, TP/FN classification is driven by the sample-target-map. A read from genome G mapping to target T is:
- **Correctly mapped** if T is listed as a target for G in the map
- **Incorrectly mapped** if T is a target for a different genome

The genome-level detection decision (detected or not) is based on whether any reads map to any of G's corresponding targets. A genome G is TP if any of its linked targets received reads; FN if none did.

This correctly handles the case where a bacterium has multiple target genes (e.g., both 16S and *rpoB*): detecting reads at either gene counts as detecting the bacterium.
