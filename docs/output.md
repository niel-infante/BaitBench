# Output Files

## Run Output Directory

Each `baitbench run` creates a timestamped subdirectory:

```
results/run_20250101_120000/
├── combined_reference.fa       # All sequences merged for fragment generation
├── mapping_reference.fa        # Targets + distractors for mapping (genome mode only)
├── weights.txt                 # Per-sequence sampling weights
├── targets.txt                 # Target sequence IDs
├── distractors.txt             # Distractor sequence IDs
├── sample.txt                  # Sample sequence IDs
├── genomes.txt                 # Genome IDs (genome mode only)
├── sample_target_map.txt       # Genome-to-target mapping (genome mode only)
├── fragments.fa                # Simulated DNA fragments (probe-biased + background)
├── reads.fa                    # Sequencing reads (trimmed to read length)
├── filtered.fa                 # Host-filtered reads (if --host-fasta)
├── mapped.sam                  # Read alignments to reference
├── detected.list               # Read counts per reference
├── run_params.tsv              # Run configuration (used by report)
├── target_groups.tsv           # Target group assignments (if --groups)
├── distractor_groups.tsv       # Distractor group assignments (always; auto or from --distractor-groups)
├── results.tsv                 # Summary metrics
├── detected_detail.tsv         # Per-reference detection and coverage detail
├── group_detail.tsv            # Per-group summary (if groups are present)
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

## results.tsv Columns

| Column | Description |
|--------|-------------|
| `run_name` | Run identifier |
| `timestamp` | Completion time |
| `num_fragments` | Fragments requested |
| `seed` | Random seed (or "NA") |
| `fragments_generated` | Fragments actually generated |
| `fragments_captured` | Fragments passing capture |
| `capture_rate` | fragments_captured / fragments_generated |
| `sample_captured` | Captured fragments from sample targets |
| `nonsample_target_captured` | Captured fragments from non-sample targets |
| `distractor_captured` | Captured fragments from distractors |
| `untargeted_captured` | Captured fragments from untargeted genomes (genome mode) |
| `reads_correctly_mapped` | Reads mapping to their source reference |
| `reads_incorrectly_mapped` | Reads mapping to a different reference |
| `sample_total` | Number of distinct sample targets |
| `nonsample_target_total` | Number of non-sample targets |
| `distractors_total` | Number of distractor sequences |
| `tp_count` | True Positives: sample targets detected |
| `fn_count` | False Negatives: sample targets not detected |
| `fp_target_count` | False Positives: non-sample targets detected |
| `fp_distractor_count` | False Positives: distractors detected |
| `fp_total` | fp_target_count + fp_distractor_count |
| `tn_target_count` | True Negatives: non-sample targets not detected |
| `tn_distractor_count` | True Negatives: distractors not detected |
| `tn_total` | tn_target_count + tn_distractor_count |
| `sensitivity` | TP / (TP + FN) |
| `specificity` | TN_total / (TN_total + FP_total) |
| `precision` | TP / (TP + FP_total) |
| `f1_score` | 2 * (precision * sensitivity) / (precision + sensitivity) |
| `reads_sequenced` | Number of reads after the sequencing step (0 if not tracked) |
| `reads_after_filter` | Number of reads after host filtering (0 if filter not applied) |
| `reads_mapped` | reads_correctly_mapped + reads_incorrectly_mapped |
| `reads_unmapped` | Reads that entered mapping but did not map to any reference |

## detected_detail.tsv Columns

One row per reference sequence:

| Column | Description |
|--------|-------------|
| `reference_id` | Sequence ID |
| `group` | Group name this sequence belongs to (sequence's own ID if no groups file provided) |
| `category` | `sample`, `nonsample_target`, `distractor`, or `untargeted` |
| `expected` | 1 if expected to be detected (sample target), 0 otherwise |
| `detected` | 1 if at least one read maps to this reference, 0 otherwise |
| `fragments_generated` | Number of fragments generated from this sequence |
| `fragments_captured` | Number of fragments captured by probes |
| `reads_assigned` | Number of reads mapped to this reference |
| `classification` | `TP`, `FN`, `FP_target`, `FP_distractor`, `TN_target`, `TN_distractor`, or `untargeted` |
| `ref_length` | Reference sequence length (bp) |
| `avg_coverage` | Average read depth across reference |
| `pct_covered_5x` | % positions with >= 5x depth |
| `pct_covered_20x` | % positions with >= 20x depth |

## group_detail.tsv Columns

Written when group files are present (`target_groups.tsv` or `distractor_groups.tsv`). One row per group:

| Column | Description |
|--------|-------------|
| `group_name` | Group identifier |
| `category` | `sample`, `nonsample_target`, or `distractor` |
| `expected` | `true` if the group is expected to be detected (sample group) |
| `detected` | `true` if at least one member sequence has reads mapped to it |
| `classification` | `TP`, `FN`, `FP_target`, `FP_distractor`, `TN_target`, or `TN_distractor` |
| `member_count` | Number of sequences in this group |
| `detected_member_count` | Number of member sequences individually detected |
| `total_reads` | Sum of reads assigned to all members of this group |

## results.json Structure

Structured JSON output with nested sections:

```json
{
  "run_info": {
    "run_name": "...",
    "timestamp": "...",
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
    "sensitivity": 1.0,
    "specificity": 0.95,
    "precision": 0.8,
    "f1_score": 0.89,
    "tp": 5, "fn": 0,
    "fp_target": 2, "fp_distractor": 1,
    "tn_target": 10, "tn_distractor": 50
  },
  "details": [ ... ]
}
```

## coverage.tsv Format

Run-length encoded read depth intervals. Consecutive positions with the same depth are collapsed into a single interval (1-based inclusive coordinates):

```
reference_id	start	end	depth
dengue_1	1	50	0
dengue_1	51	100	3
dengue_1	101	200	5
...
```

This format is typically 100-1000x smaller than per-position output, making it feasible for large target panels.
