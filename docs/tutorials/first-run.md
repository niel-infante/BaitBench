# Your First Simulation

This tutorial runs BaitBench on a small included dataset and walks through what every part of the output means. By the end you will understand the 7-step pipeline, what the metrics measure, and how to read the two key output files.

**Time:** about 5 minutes.  
**Prerequisite:** [Installation](installation.md) complete, conda environment active.

---

## The dataset

The tutorial dataset is in `examples/tutorial/` and represents a small respiratory virus panel:

| File | Contents |
|------|----------|
| `targets.fa` | 3 target sequences: Influenza A H3N2, SARS-CoV-2, Dengue virus 2 |
| `distractors.fa` | 2 human genomic sequences (background) |
| `probes.fa` | 96 probes tiling all 3 target sequences |
| `sample.tsv` | Sample manifest: only Influenza A and SARS-CoV-2 are present in this specimen |

Dengue virus 2 is in the probe panel but *not* in the sample manifest — it is a non-sample target. This lets us test whether BaitBench can distinguish a specimen containing two viruses from a panel designed for three.

---

## Run the simulation

From the BaitBench directory:

```bash
conda activate baitbench

./target/release/baitbench run \
  --targets examples/tutorial/targets.fa \
  --distractors examples/tutorial/distractors.fa \
  --probes examples/tutorial/probes.fa \
  --sample examples/tutorial/sample.tsv \
  --distractor-fraction 0.5 \
  --num-fragments 2000 \
  --seed 42 \
  --report none \
  --outdir results/tutorial
```

`--seed 42` makes the simulation reproducible. `--report none` skips HTML report generation so you do not need R installed for this tutorial. `--distractor-fraction 0.5` sets half the simulated fragments to come from the human background sequences.

!!! note
    If you added `baitbench` to your PATH, omit the `./target/release/` prefix.

---

## The pipeline output

You will see 7 steps in the log:

```
[INFO ] Step 1/7: Preparing reference...
[INFO ] Step 2/7: Simulating fragments (mode=thermodynamic, capture_fraction=0.50)...
[INFO ] Step 3/7: Sequencing fragments (simulator=Perfect)...
[INFO ] Step 4/7: Skipping host filtering (no host genome provided)
[INFO ] Step 5/7: Mapping reads to reference...
[INFO ] Step 6/7: Generating detection list...
[INFO ] Step 7/7: Calculating metrics and coverage...
```

Each step produces an intermediate file in the output directory. Here is what each one does:

| Step | What happens |
|------|-------------|
| **Prepare** | Combines targets and distractors into a single reference FASTA; assigns abundance weights to each sequence based on the sample manifest and distractor fraction |
| **Simulate** | Aligns probes to the reference; generates fragments biased toward probe binding sites using thermodynamic scoring (SantaLucia 1998 nearest-neighbor model) |
| **Sequence** | Trims fragments to read length (120 bp by default) |
| **Filter** | Removes host-matching reads — skipped here because we did not provide a host genome |
| **Map** | Aligns reads back to the combined reference |
| **List** | Counts how many reads mapped to each reference sequence |
| **Metrics** | Classifies each target and distractor as TP/FP/FN/TN and computes sensitivity, specificity, precision, and F1 |

Near the end of the log you will see the classification summary:

```
[INFO ]   True Positives (sample detected): 2
[INFO ]   False Negatives (sample missed): 0
[INFO ]   FP targets (non-sample target detected): 0
[INFO ]   FP distractors (distractor detected): 1
[INFO ]   TN targets (non-sample target not detected): 1
[INFO ]   TN distractors (distractor not detected): 0
[INFO ]   Sensitivity: 1.0000
[INFO ]   Specificity: 0.5000
[INFO ]   Precision: 0.6667
[INFO ]   F1 Score: 0.8000
```

---

## Reading the results

Results are written to `results/tutorial/<run_name>/`. The run name includes a timestamp, so it will differ from examples below.

### group_detail.tsv

This is the best file to read first. Each row is one logical entity (a target or distractor group), and shows whether it was detected:

| group_name | category | expected | detected | classification | total_reads |
|------------|----------|----------|----------|---------------|-------------|
| SARS_CoV_2 | sample | true | true | **TP** | 787 |
| Influenza_A_H3N2 | sample | true | true | **TP** | 738 |
| distractors | distractor | false | true | **FP_distractor** | 475 |
| Dengue_virus_2 | target | false | false | **TN_target** | 0 |

Both sample viruses are detected (TP). Dengue is correctly *not* detected (TN_target) — the probes are present in the panel, but no Dengue fragments were generated because it was absent from the sample manifest. The human distractor group was detected (FP_distractor) — more on this below.

### detected_detail.tsv

This shows the same information at the individual sequence level, with per-sequence coverage statistics:

| reference_id | classification | reads_assigned | avg_coverage | pct_covered_5x | pct_covered_20x |
|---|---|---|---|---|---|
| SARS_CoV_2 | TP | 787 | 47.2× | 97.5% | 95.2% |
| Influenza_A_H3N2 | TP | 738 | 44.3× | 97.2% | 96.0% |
| Human_chr1_frag | FP_distractor | 250 | 15.0× | 93.5% | 14.1% |
| Human_chr22_frag | FP_distractor | 225 | 13.5× | 92.1% | 11.2% |

Notice the coverage depth difference:

- **Sample targets**: ~45× average depth, >95% of bases covered at 20×
- **Human distractors**: ~14× average depth, only ~12% of bases covered at 20×

This is an important pattern. The human sequences received reads because probes occasionally align to background sequences by chance — every random 120 bp probe has some probability of a short fortuitous match. But those reads spread thinly across the whole sequence rather than concentrating in specific probe-binding sites. High breadth at low depth (93% at 5×, 12% at 20×) is a hallmark of background noise, not genuine capture.

In a real evaluation you would use this depth difference, along with cross-reactivity analysis (`baitbench xreact`), to distinguish genuine target detection from background.

---

## What the metrics mean

**Sensitivity** (1.0): the fraction of sample targets that were detected. Both Influenza A and SARS-CoV-2 were found — perfect recall.

**Specificity** (0.5): the fraction of non-sample entities that were correctly *not* detected. There are two non-sample entities: Dengue (correctly negative) and the human distractor group (incorrectly positive). 1 of 2 = 0.5.

**Precision** (0.67): of the three entities detected (Influenza, SARS-CoV-2, human), two were true positives. 2/3 = 0.67.

**F1** (0.80): the harmonic mean of sensitivity and precision.

---

## Next steps

- Try removing `--sample` from the command. All three viruses will be treated as sample targets, so Dengue will generate reads and appear as a TP. This shows how the sample manifest controls the discrimination test.
- Try `--distractor-fraction 0.9` (the default) to simulate a lower-abundance specimen and observe how coverage on the targets changes.
- Continue to the [How-To Guides](../how-to/index.md) to learn how to tune parameters, interpret results in depth, and run genome mode with the bacterial dataset in `examples/tutorial-genome/`.
