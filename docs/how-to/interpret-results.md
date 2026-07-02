# Interpret Results

After a run, results land in `<outdir>/<run_name>/`. This guide walks through the key output files in the order you should read them.

---

## Start with group_detail.tsv

`group_detail.tsv` gives the clearest picture first: one row per logical entity (organism or distractor group), showing whether it was detected.

```
group_name                        category    expected  detected  classification  member_count  total_reads
Influenza_A_H3N2                  sample      true      true      TP              1             738
SARS_CoV_2                        sample      true      true      TP              1             787
Dengue_virus_2                    target      false     false     TN_target       1             0
distractors                       distractor  false     true      FP_distractor   2             475
```

**Columns:**

| Column | Meaning |
|--------|---------|
| `category` | `sample` (expected), `target` (non-sample, in panel), `distractor` (off-target) |
| `expected` | Was this entity expected to be detected? (`true` = sample) |
| `detected` | Was at least one read assigned to this entity? |
| `classification` | TP, FN, FP_target, FP_distractor, TN_target, or TN_distractor |
| `member_count` | Number of sequences in this group |
| `detected_member_count` | How many member sequences individually had reads |
| `total_reads` | Total reads assigned across all members |

If you did not provide `--groups`, distractor contigs from the same file are still automatically grouped by file stem. Target sequences each form their own singleton group.

---

## Then check detected_detail.tsv

`detected_detail.tsv` has one row per individual sequence and adds coverage statistics — useful for understanding how uniformly the probes cover each target.

```
reference_id              classification  reads_assigned  avg_coverage  pct_covered_5x  pct_covered_20x
SARS_CoV_2                TP              787             47.2×         97.5%           95.2%
Influenza_A_H3N2          TP              738             44.3×         97.2%           96.0%
Human_chr1_frag           FP_distractor   250             15.0×         93.5%           14.1%
Human_chr22_frag          FP_distractor   225             13.5×         92.1%           11.2%
```

**The coverage columns are the key diagnostic:**

- `pct_covered_5x`: fraction of bases with at least 5 reads — indicates whether the probe set tiles the target
- `pct_covered_20x`: fraction with at least 20 reads — indicates whether depth is sufficient for confident detection

**Genuine capture vs background noise:**

| Signal type | Typical pattern |
|-------------|-----------------|
| Genuine capture | High average depth, high 5× AND 20× breadth (e.g., 97% at 5×, 95% at 20×) |
| Background noise | Moderate average depth, high 5× but low 20× breadth (e.g., 93% at 5×, 14% at 20×) |

Background reads spread thinly across an entire sequence. Probe-biased reads concentrate at binding sites. Even a sequence with many background reads will show low 20× breadth because no individual position accumulates sufficient depth.

---

## Summary metrics (results.tsv)

`results.tsv` is a single-row TSV with every metric in one place. The key classification counts and rates:

### Classification counts

| Column | What it counts |
|--------|----------------|
| `tp_count` | Sample targets that were detected |
| `fn_count` | Sample targets that were NOT detected |
| `fp_target_count` | Non-sample targets (in panel, not in specimen) that were detected |
| `fp_distractor_count` | Distractor sequences that were detected |
| `tn_target_count` | Non-sample targets correctly not detected |
| `tn_distractor_count` | Distractors correctly not detected |

A genome (or group) is "detected" if at least one read maps to it after capture and mapping.

### Performance rates

| Metric | Formula | What it means |
|--------|---------|----------------|
| `sensitivity` | TP / (TP + FN) | Fraction of sample targets found. Low = probes miss targets. |
| `specificity` | TN / (TN + FP) | Fraction of non-sample entities correctly rejected. Low = false positives. |
| `precision` | TP / (TP + FP) | Of everything detected, fraction that is a real target. Low = detections are mostly noise. |
| `f1_score` | 2 × (prec × sens) / (prec + sens) | Harmonic mean of precision and sensitivity. Balanced overall score. |

### Fragment and read counts

| Column | What it tracks |
|--------|----------------|
| `sample_captured` | Fragments from sample targets that passed the capture filter |
| `distractor_captured` | Fragments from distractors that passed the capture filter |
| `reads_correctly_mapped` | Reads that mapped back to their source reference |
| `reads_incorrectly_mapped` | Reads that mapped to a different reference (cross-mapping) |

`reads_incorrectly_mapped` is a sensitive indicator of cross-reactivity. A panel can show perfect genome-level metrics (all TP, no FP) while still having substantial cross-mapping — reads from virus A map to virus B's reference. If you see a large `reads_incorrectly_mapped` count alongside good sensitivity, investigate the `detected_detail.tsv` to find which references are collecting misassigned reads.

---

## Interpreting common result patterns

### Perfect sensitivity, imperfect specificity

```
sensitivity: 1.0   specificity: 0.5   fp_distractor_count: 1
```

All targets detected; one distractor group also detected. Check `detected_detail.tsv` to see how many reads the distractor received and whether the coverage pattern looks like genuine capture (high 20× breadth) or background noise (low 20× breadth).

### Imperfect sensitivity

```
sensitivity: 0.8   fn_count: 1
```

One sample target was missed. Check the corresponding row in `detected_detail.tsv` — `reads_assigned` will be 0. Common causes: probes don't cover that target, or the target is very short relative to fragment length.

### FP_target detections

```
fp_target_count: 1   tn_target_count: 0
```

A non-sample target (in the panel, not in the specimen) was detected. This means probes designed for other panel members accidentally capture this target. Run `baitbench xreact` to find which probes are responsible.

### High `reads_incorrectly_mapped`

Reads generated from one source are mapping to a different reference. This does not necessarily affect genome-level TP/FP counts, but indicates cross-mapping between targets. Consider running `baitbench panel-qc` to assess target uniqueness.

---

## Reading the JSON output (results.json)

`results.json` contains the same data as `results.tsv` in a nested structure, suitable for downstream processing:

```json
{
  "metrics": {
    "sensitivity": 1.0,
    "specificity": 0.5,
    "tp": 2, "fn": 0, "fp_target": 0, "fp_distractor": 1
  },
  "details": [ ... ]
}
```

---

## Coverage depth (coverage.tsv)

`coverage.tsv` is run-length encoded: consecutive positions with the same depth are collapsed into one interval.

```
reference_id        start  end   depth
Influenza_A_H3N2    1      50    0
Influenza_A_H3N2    51     170   12
Influenza_A_H3N2    171    350   47
...
```

Use this to plot coverage profiles, identify uncovered gaps, or calculate any depth threshold you need. The HTML report does this automatically if R is available.
