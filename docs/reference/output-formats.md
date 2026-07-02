# Output Formats

Column-level specifications for all BaitBench output files.

---

## Output Directory Structure

Each `baitbench run` creates a timestamped subdirectory under `--outdir`:

```
results/run_20250101_120000/
├── combined_reference.fa       # All sequences merged for fragment generation
├── mapping_reference.fa        # Targets + distractors for mapping (genome mode only)
├── weights.txt                 # Per-sequence sampling weights (TSV: id<TAB>weight)
├── targets.txt                 # Target sequence IDs (one per line)
├── distractors.txt             # Distractor sequence IDs (one per line)
├── sample.txt                  # Sample sequence IDs (one per line)
├── genomes.txt                 # Genome IDs (genome mode only)
├── sample_target_map.txt       # Genome-to-target mapping (genome mode only)
├── target_groups.tsv           # Target group assignments (if --groups)
├── distractor_groups.tsv       # Distractor group assignments (always written)
├── fragments.fa                # Simulated DNA fragments (probe-biased + background)
├── reads.fa                    # Sequencing reads (trimmed to read length)
├── filtered.fa                 # Host-filtered reads (if --host-fasta)
├── mapped.sam                  # Read alignments to reference
├── detected.list               # Read counts per reference (TSV: id<TAB>count)
├── run_params.tsv              # Run configuration (used by report)
├── results.tsv                 # Summary metrics (one row)
├── detected_detail.tsv         # Per-reference detection and coverage detail
├── group_detail.tsv            # Per-group summary
├── results.json                # Machine-readable JSON metrics
├── coverage.tsv                # Run-length encoded read depth intervals
├── report.html                 # HTML report (--report full, requires R)
├── report.Rmd                  # Editable RMarkdown (--report rmd)
├── species_calls.tsv           # Species-level calls (if --identify)
├── species_calls.json          # Species calls JSON (if --identify)
├── target_similarity.tsv       # Target pairwise similarity (if --identify)
├── capture.log                 # Capture alignment log
├── mapping.log                 # Read mapping log
└── host_filter.log             # Host filtering log (if --host-fasta)
```

With `--cleanup`, intermediate files (fragments.fa, reads.fa, filtered.fa, mapped.sam, *.log) are deleted after the run.

---

## results.tsv

One-row TSV with all summary metrics. The first row is a header.

| Column | Description |
|--------|-------------|
| `run_name` | Run identifier |
| `timestamp` | Completion time (ISO 8601) |
| `num_fragments` | Fragments requested (`--num-fragments`) |
| `seed` | Random seed used (or `NA` if not set) |
| `fragments_generated` | Fragments actually generated |
| `fragments_captured` | Fragments passing the capture step |
| `capture_rate` | `fragments_captured / fragments_generated` |
| `sample_generated` | Fragments generated from sample targets |
| `nonsample_target_generated` | Fragments generated from non-sample targets |
| `distractor_generated` | Fragments generated from distractor sequences |
| `untargeted_generated` | Fragments generated from untargeted genomes (genome mode) |
| `sample_captured` | Captured fragments originating from sample targets |
| `nonsample_target_captured` | Captured fragments from non-sample targets (should be 0) |
| `distractor_captured` | Captured fragments from distractor sequences |
| `untargeted_captured` | Captured fragments from untargeted genomes (genome mode) |
| `reads_correctly_mapped` | Reads that mapped back to their source reference |
| `reads_incorrectly_mapped` | Reads that mapped to a different reference |
| `reads_sequenced` | Reads after the sequencing step (0 if not tracked) |
| `reads_after_filter` | Reads after host filtering (0 if filter not used) |
| `reads_mapped` | `reads_correctly_mapped + reads_incorrectly_mapped` |
| `reads_unmapped` | Reads that entered mapping but did not map to any reference |
| `sample_total` | Number of distinct sample targets |
| `nonsample_target_total` | Number of non-sample targets |
| `distractors_total` | Number of distractor sequences |
| `tp_count` | True Positives: sample targets detected |
| `fn_count` | False Negatives: sample targets not detected |
| `fp_target_count` | False Positives: non-sample targets detected |
| `fp_distractor_count` | False Positives: distractors detected |
| `fp_total` | `fp_target_count + fp_distractor_count` |
| `tn_target_count` | True Negatives: non-sample targets not detected |
| `tn_distractor_count` | True Negatives: distractors not detected |
| `tn_total` | `tn_target_count + tn_distractor_count` |
| `sensitivity` | `TP / (TP + FN)` |
| `specificity` | `TN_total / (TN_total + FP_total)` |
| `precision` | `TP / (TP + FP_total)` |
| `f1_score` | `2 × (precision × sensitivity) / (precision + sensitivity)` |

---

## detected_detail.tsv

One row per: (1) every detected reference sequence, (2) FN sequences (sample targets with zero reads), (3) untargeted genomes (genome mode). TN_target and TN_distractor sequences are NOT included.

| Column | Description |
|--------|-------------|
| `reference_id` | Sequence ID (as in the targets/distractors FASTA) |
| `group` | Group name this sequence belongs to (its own ID if no groups file) |
| `category` | `sample` (sample targets, TP or FN), `target` (non-sample targets, FP), `distractor` (FP), `untargeted` (genome mode), or `unknown` |
| `expected` | `"true"` if expected to be detected (sample target), `"false"` otherwise |
| `detected` | `"true"` if at least one read maps here, `"false"` otherwise |
| `fragments_generated` | Fragments generated from this sequence |
| `fragments_captured` | Fragments captured by probes from this sequence |
| `reads_assigned` | Reads mapped to this reference |
| `classification` | `TP`, `FN`, `FP_target`, `FP_distractor`, `TN_target`, `TN_distractor`, or `untargeted` |
| `ref_length` | Reference sequence length (bp) |
| `avg_coverage` | Average read depth across the full reference |
| `pct_covered_5x` | % of positions with ≥ 5× depth |
| `pct_covered_20x` | % of positions with ≥ 20× depth |

---

## group_detail.tsv

One row per logical group. Always written (distractor groups are auto-created from FASTA file stems even without `--groups`).

| Column | Description |
|--------|-------------|
| `group_name` | Group identifier |
| `category` | `sample` (sample groups), `target` (non-sample target groups), or `distractor` |
| `expected` | `"true"` if the group is expected to be detected (sample group), `"false"` otherwise |
| `detected` | `true` if at least one member sequence has ≥ 1 read |
| `classification` | `TP`, `FN`, `FP_target`, `FP_distractor`, `TN_target`, or `TN_distractor` |
| `member_count` | Number of sequences in this group |
| `detected_member_count` | Number of member sequences individually detected |
| `total_reads` | Sum of reads assigned to all members |

---

## results.json

Structured JSON with nested sections. Fields mirror results.tsv plus the per-sequence detail array.

```json
{
  "run_info": {
    "run_name": "run_20250101_120000",
    "timestamp": "2025-01-01T12:00:00Z",
    "num_fragments": 10000,
    "seed": "42"
  },
  "capture_stats": {
    "fragments_generated": 10000,
    "fragments_captured": 3500,
    "capture_rate": 0.35
  },
  "read_level": {
    "reads_correctly_mapped": 3400,
    "reads_incorrectly_mapped": 100,
    "reads_mapped": 3500,
    "reads_unmapped": 0,
    "reads_sequenced": 3500,
    "reads_after_filter": 0
  },
  "metrics": {
    "tp_count": 5, "fn_count": 0,
    "fp_target_count": 0, "fp_distractor_count": 1, "fp_total": 1,
    "tn_target_count": 10, "tn_distractor_count": 50, "tn_total": 60,
    "sensitivity": 1.0,
    "specificity": 0.98,
    "precision": 0.83,
    "f1_score": 0.91
  },
  "details": {
    "true_positives": ["SARS_CoV_2", "Influenza_A_H3N2"],
    "false_negatives": [],
    "fp_targets": [],
    "fp_distractors": ["distractors"],
    "tn_targets": [...],
    "tn_distractors": [...],
    "unknown_detected": [],
    "untargeted_genomes": [],
    "detail_rows": [ ... ]
  }
}
```

---

## coverage.tsv

Run-length encoded read depth. Consecutive positions at the same depth are merged into one interval. Coordinates are 1-based, inclusive.

```
reference_id	start	end	depth
dengue_1	1	50	0
dengue_1	51	100	3
dengue_1	101	200	5
```

| Column | Description |
|--------|-------------|
| `reference_id` | Sequence ID |
| `start` | Start position (1-based, inclusive) |
| `end` | End position (1-based, inclusive) |
| `depth` | Read depth for all positions in this interval |

This format is typically 100–1000× smaller than per-position output, making it practical for large target panels.

---

## Intermediate Files

These are written during a run and can be removed with `--cleanup`:

| File | Format | Written by |
|------|--------|-----------|
| `combined_reference.fa` | FASTA | `prepare` |
| `mapping_reference.fa` | FASTA | `prepare` (genome mode) |
| `weights.txt` | TSV: `id<TAB>weight` | `prepare` |
| `fragments.fa` | FASTA, IDs: `{seq_id}_fragment_{n} start={pos} length={len}` | `simulate` |
| `reads.fa` | FASTA | `sequence` |
| `reads_R2.fa` | FASTA | `sequence` (paired-end) |
| `filtered.fa` | FASTA | `filter` |
| `mapped.sam` | SAM | `map` |
| `detected.list` | TSV: `id<TAB>count` | `list` |
| `*.log` | Plain text | various steps |
