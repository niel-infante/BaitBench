# Genome Mode: Simulating Bacterial Capture

Standard mode works well for viruses, where the full genome is also the probe target. Bacteria are different: a typical metagenomics probe panel targets a short conserved region — the 16S rRNA gene (~1500 bp) — but the organism itself has a genome of several megabases. Probe hybridization occurs against whatever DNA is in the sample, not just the 16S gene.

**Genome mode** (`--genomes`) models this correctly. Fragments are generated from full bacterial genomes, but reads are mapped back to the 16S target sequences. This tutorial walks through a genome mode run using two bacteria included in the tutorial dataset.

**Time:** about 10 minutes.  
**Prerequisite:** Complete [Your First Simulation](first-run.md) first — this tutorial builds on what was introduced there.

---

## Why genome mode matters

In standard mode, every fragment comes from a target sequence. For a bacterium this would mean simulating reads only from the 16S gene, which misses an important reality: most of the DNA in a sample comes from the rest of the genome, not from the probe-targeted region. Genome mode generates fragments proportionally from the full genome, applies thermodynamic scoring to find which fragments overlap probe-binding sites, and then maps reads back to the 16S reference — just as a real capture-sequencing experiment would.

---

## The dataset

The tutorial bacterial dataset is in `examples/tutorial-genome/`:

| File | Contents |
|------|----------|
| `genomes.fa` | 2 full bacterial genome fragments (~8500 bp each): *M. tuberculosis* H37Rv and *S. aureus* MRSA252 |
| `targets.fa` | 2 16S rRNA sequences (~1500 bp each), one per bacterium |
| `distractors.fa` | 1 human mitochondrial DNA fragment (background) |
| `probes.fa` | 70 probes tiling the two 16S target sequences |
| `sample.tsv` | Sample manifest: both bacteria are present in this specimen |
| `sample_target_map.tsv` | Maps each genome ID to its corresponding 16S target ID |

### The two-file split

This is the key difference from standard mode:

| Role | Standard mode | Genome mode |
|------|--------------|-------------|
| Fragment generation | `--targets` | **`--genomes`** |
| Read mapping | `--targets` | `--targets` |

In genome mode `--targets` is still required, but it serves only as the mapping reference — the sequences you expect reads to map *to*. Fragment generation draws from `--genomes` instead.

### The sample-target-map

`--sample` accepts genome IDs (not target IDs) in genome mode. The sample-target-map tells BaitBench which target ID corresponds to each genome:

```
# sample_target_map.tsv
Mycobacterium_tuberculosis_H37Rv    Mycobacterium_tuberculosis_H37Rv_16S
Staphylococcus_aureus_MRSA252       Staphylococcus_aureus_MRSA252_16S
```

This matters for classification: a read originating from *M. tuberculosis* H37Rv that maps to its 16S target is correctly mapped; a read from *S. aureus* that accidentally maps to the MTB 16S target would be flagged as incorrectly mapped.

If you omit `--sample-target-map`, BaitBench tries to auto-link by name — a genome ID matches a target if the target ID starts with `{genome_id}|` or is identical. For this dataset the names differ, so we supply the map explicitly.

---

## Run the simulation

```bash
conda activate baitbench

./target/release/baitbench run \
  --targets   examples/tutorial-genome/targets.fa \
  --genomes   examples/tutorial-genome/genomes.fa \
  --distractors examples/tutorial-genome/distractors.fa \
  --probes    examples/tutorial-genome/probes.fa \
  --sample    examples/tutorial-genome/sample.tsv \
  --sample-target-map examples/tutorial-genome/sample_target_map.tsv \
  --num-fragments 3000 \
  --capture-fraction 0.6 \
  --seed 42 \
  --report none \
  --outdir results/tutorial-genome
```

---

## Pipeline differences in genome mode

The log will look familiar, but two things are different in the **Prepare** step:

```
[INFO ] Combining genomes + distractors to .../combined_reference.fa...
[INFO ] Combining targets + distractors to .../mapping_reference.fa...
[INFO ] Sample-target-map: 2 explicit, 0 auto-linked, 0 untargeted
```

Genome mode creates **two** reference files instead of one:

- `combined_reference.fa` — genomes + distractors, used for fragment simulation (step 2)
- `mapping_reference.fa` — targets + distractors, used for read mapping (step 5)

The sample-target-map summary confirms 2 explicit links were loaded, 0 were resolved by auto-naming, and 0 genomes are untargeted (a genome without a matching target would still generate fragments but be excluded from TP/FP/FN/TN accounting).

In the **Simulate** step, all 70 probes aligned successfully to the combined reference (the 16S sequences are embedded within the genome fragments):

```
[INFO ] Loaded 70 probe hits (0 skipped) from .../fragments.probe_hits.sam
[INFO ] Sampling 1800 capture fragments (capture_fraction=0.60) + 1200 background fragments
```

---

## Reading the results

### group_detail.tsv

| group_name | category | detected | classification | total_reads |
|---|---|---|---|---|
| Staphylococcus_aureus_MRSA252_16S | sample | true | **TP** | 969 |
| Mycobacterium_tuberculosis_H37Rv_16S | sample | true | **TP** | 897 |
| distractors | distractor | true | **FP_distractor** | 813 |

Both bacteria are detected via their 16S targets. The human mitochondrial distractor is also detected — human mitochondria carry their own 12S and 16S ribosomal RNA genes, which share enough similarity with bacterial 16S probes to produce cross-reactive reads. This is a real phenomenon in 16S metagenomic panels and one reason probe cross-reactivity analysis (`baitbench xreact`) is worth running.

### detected_detail.tsv

| reference_id | reads_assigned | avg_coverage | pct_covered_5x | pct_covered_20x |
|---|---|---|---|---|
| Staphylococcus_aureus_MRSA252_16S | 969 | 76.7× | 100.0% | 97.7% |
| Mycobacterium_tuberculosis_H37Rv_16S | 897 | 71.0× | 98.9% | 96.1% |
| Human_mtDNA_frag | 813 | 48.8× | 96.0% | 92.7% |

Coverage depth on the 16S targets is high and uniform: both bacteria reach >97% breadth at 20×. The human distractor also has high depth (92.7% at 20×) — unlike the thin background signal seen in the virus tutorial, this reflects genuine cross-reactive probe binding to mitochondrial ribosomal sequences rather than accidental k-mer matches.

### Correctly vs incorrectly mapped reads

The log reports a metric not present in standard mode:

```
[INFO ]   Reads correctly mapped: 1866
[INFO ]   Reads incorrectly mapped: 813
[INFO ]   Reads unmapped: 321
```

In genome mode, a read is **correctly mapped** if it came from genome G and mapped to a target T where G→T is in the sample-target-map. The 813 incorrectly mapped reads are the distractor reads: they originated from the human mitochondrial sequence and mapped to bacterial 16S targets (the genome-target link for human DNA is not in the map).

The 321 unmapped reads came from the non-16S flanking regions of the bacterial genomes. Those fragments were generated as background, but the flanking genomic sequence has no probes and little similarity to the 16S mapping reference, so they do not map — which is the expected behaviour.

---

## Key metrics

| Metric | Value | Meaning |
|---|---|---|
| Sensitivity | 1.0000 | Both bacteria detected |
| Specificity | 0.0000 | The only non-sample entity (human distractor) was detected |
| Precision | 0.6667 | 2 of 3 detected entities are true positives |
| F1 | 0.8000 | Harmonic mean of sensitivity and precision |

Specificity is 0 here because there is only one distractor group and it was detected. In a real panel with several distractor organisms, only a subset would typically cross-react, giving a non-zero specificity.

---

## Next steps

- Try omitting `--sample-target-map` and observe whether BaitBench auto-links names correctly (it will not for this dataset, since the names do not follow the `{genome_id}|{target_id}` convention).
- Try `--ct 25` instead of the default distractor fraction to simulate a clinical specimen at a specific CT value.
- Run `baitbench xreact --probes examples/tutorial-genome/probes.fa --against examples/tutorial-genome/distractors.fa` to quantify the probe cross-reactivity against human mitochondrial DNA.
- See [Analyze Cross-Reactivity](../how-to/analyze-cross-reactivity.md) for the full cross-reactivity workflow.
