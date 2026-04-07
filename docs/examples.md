# Usage Examples

## Basic Probe Evaluation

Test whether probes capture all targets and reject distractors:

```bash
baitbench run \
  --targets targets.fa \
  --distractors distractors.fa \
  --probes probes.fa \
  --num-fragments 10000 \
  --seed 42 \
  --outdir results
```

All targets are treated as "present" (no `--sample`). The default distractor fraction is 0.9 (90% background, 10% target).

## Sample Discrimination Testing

Test whether probes can detect specific targets while rejecting others in the panel:

```bash
# Inline sample IDs
baitbench run \
  --targets all_viruses.fa \
  --distractors bacteria.fa \
  --probes probes.fa \
  --sample dengue_1 zika_virus \
  --num-fragments 10000 \
  --outdir results

# With custom weights (dengue at 5x abundance)
baitbench run \
  --targets all_viruses.fa \
  --distractors bacteria.fa \
  --probes probes.fa \
  --sample dengue_1 5 zika_virus \
  --num-fragments 10000 \
  --outdir results

# Using a TSV manifest file
baitbench run \
  --targets all_viruses.fa \
  --distractors bacteria.fa \
  --probes probes.fa \
  --sample sample.tsv \
  --num-fragments 10000 \
  --outdir results
```

Non-sample targets will have FP_target classification if detected, testing cross-reactivity within the panel.

## Clinical Specimen Simulation with CT

Simulate specimens at different viral loads:

```bash
# High viral load (CT 20)
baitbench run \
  --targets targets.fa \
  --distractors distractors.fa \
  --probes probes.fa \
  --ct 20 \
  --num-fragments 10000 \
  --outdir results_ct20

# Low viral load (CT 30)
baitbench run \
  --targets targets.fa \
  --distractors distractors.fa \
  --probes probes.fa \
  --ct 30 \
  --num-fragments 10000 \
  --outdir results_ct30
```

## Genome Mode for Bacteria

When probe targets are sub-regions of large genomes:

```bash
# targets.fa: 16S gene sequences
# genomes.fa: full bacterial genomes
# mapping.tsv links genome IDs to target gene IDs

baitbench run \
  --targets 16S_targets.fa \
  --genomes bacteria_genomes.fa \
  --distractors human.fa \
  --probes probes.fa \
  --sample-target-map mapping.tsv \
  --sample e_coli s_aureus \
  --num-fragments 50000 \
  --outdir results
```

Use higher `--num-fragments` for bacteria since genomes are much larger than target regions, requiring more fragments to achieve adequate target coverage.

## Mixed Panels (Virus + Bacteria)

Genome mode handles mixed panels naturally. Virus genomes that match their target IDs auto-link:

```bash
# genomes.fa: influenza_a (13kb), e_coli (5Mb)
# targets.fa: influenza_a (same seq), e_coli_16S (1.5kb subsequence)
# mapping.tsv only needs the e_coli entry (influenza_a auto-links)

baitbench run \
  --targets targets.fa \
  --genomes genomes.fa \
  --distractors human.fa \
  --probes probes.fa \
  --sample-target-map mapping.tsv \
  --sample influenza_a e_coli \
  --num-fragments 50000 \
  --outdir results
```

## Multiple Distractor Sources

Provide multiple distractor FASTA files:

```bash
baitbench run \
  --targets targets.fa \
  --distractors bacteria.fa \
  --distractors fungi.fa \
  --distractors protozoa.fa \
  --probes probes.fa \
  --num-fragments 10000 \
  --outdir results
```

All distractor sequences are concatenated and share the same per-sequence weight.

## Capture Fraction Sweep

Control what fraction of simulated fragments come from probe binding sites. With thermodynamic mode (default), high-affinity probe-reference alignments receive higher weight.

```bash
# Default: 50% probe-biased, 50% background (thermodynamic mode)
baitbench run \
  --targets targets.fa \
  --distractors distractors.fa \
  --probes probes.fa \
  --num-fragments 10000 \
  --capture-fraction 0.5 \
  --outdir results_thermo

# High capture fraction with simple (uniform) weighting
baitbench run \
  --targets targets.fa \
  --distractors distractors.fa \
  --probes probes.fa \
  --simulate-mode simple \
  --capture-fraction 0.8 \
  --num-fragments 10000 \
  --outdir results_simple

# Sweep capture fractions with coverage-curve
baitbench coverage-curve \
  --targets targets.fa \
  --distractors distractors.fa \
  --probes probes.fa \
  --sample target_1 target_2 \
  --capture-fraction-values 0.2 0.4 0.6 0.8 \
  --ct 25 \
  --outdir cf_sweep
```

## Sequencing Depth Control

Control the number of reads output by the sequencing step:

```bash
# Sample 5000 reads with replacement (models limited sequencing)
baitbench run \
  --targets targets.fa \
  --distractors distractors.fa \
  --probes probes.fa \
  --num-fragments 50000 \
  --num-sequences 5000 \
  --seed 42 \
  --outdir results
```

## Host Filtering

Remove host reads before mapping:

```bash
baitbench run \
  --targets targets.fa \
  --distractors distractors.fa \
  --probes probes.fa \
  --host-fasta human_genome.fa \
  --num-fragments 10000 \
  --outdir results
```

## Coverage Curve Analysis

Sweep parameters to understand how conditions affect coverage:

```bash
# Sweep CT values only
baitbench coverage-curve \
  --targets targets.fa \
  --distractors distractors.fa \
  --probes probes.fa \
  --sample dengue_1 zika_virus \
  --ct-values 20 25 30 35 \
  --num-fragments 10000 \
  --seed 42 \
  --outdir coverage_ct

# Sweep CT and capture fraction (combinatorial)
baitbench coverage-curve \
  --targets targets.fa \
  --distractors distractors.fa \
  --probes probes.fa \
  --sample dengue_1 \
  --ct-values 20 25 30 \
  --capture-fraction-values 0.3 0.5 0.7 \
  --num-fragments 10000 \
  --outdir coverage_ct_cf

# Sweep all three parameters
baitbench coverage-curve \
  --targets targets.fa \
  --distractors distractors.fa \
  --probes probes.fa \
  --sample dengue_1 \
  --ct-values 20 25 30 \
  --capture-fraction-values 0.3 0.5 0.7 \
  --num-sequences-values 500 1000 5000 \
  --num-fragments 10000 \
  --outdir coverage_full

# Fixed CT with capture fraction sweep
baitbench coverage-curve \
  --targets targets.fa \
  --distractors distractors.fa \
  --probes probes.fa \
  --sample dengue_1 \
  --ct 25 \
  --capture-fraction-values 0.1 0.3 0.5 0.7 0.9 \
  --num-fragments 10000 \
  --outdir coverage_cf
```

## Probe Design QC

Evaluate probe tiling independently of the simulation:

```bash
baitbench probe-coverage \
  --targets targets.fa \
  --probes probes.fa \
  --proximity 100 \
  --outdir probe_qc
```

## Cross-Reactivity Analysis

Check whether probes have off-target homology to specific genomes:

```bash
# Probe-to-genome: which probes match the human genome?
baitbench xreact \
  --probes probes.fa \
  --against human_genome.fa \
  --threshold 80 \
  --outdir xreact_human

# Probe-to-genome: check against multiple references
baitbench xreact \
  --probes probes.fa \
  --against human_genome.fa mouse_genome.fa \
  --threshold 80 \
  --outdir xreact_hosts

# Probe-to-probe: find probes that are too similar to each other
baitbench xreact \
  --probes probes.fa \
  --self \
  --threshold 80 \
  --outdir xreact_self

# Both modes together
baitbench xreact \
  --probes probes.fa \
  --against human_genome.fa \
  --self \
  --threshold 80 \
  --outdir xreact_full
```

## Target Panel QC

Assess whether a target panel can distinguish between species before running simulations:

```bash
# Basic panel QC
baitbench panel-qc \
  --targets gene_targets.fa \
  --sample-target-map mapping.tsv \
  --outdir panel_qc_results

# Stricter similarity threshold (95% instead of default 90%)
baitbench panel-qc \
  --targets gene_targets.fa \
  --sample-target-map mapping.tsv \
  --identity-threshold 95 \
  --outdir panel_qc_strict

# Skip HTML report (just produce TSV files)
baitbench panel-qc \
  --targets gene_targets.fa \
  --sample-target-map mapping.tsv \
  --report none \
  --outdir panel_qc_tsv
```

The HTML report includes a species discriminability chart, confusion matrix heatmap, and target composition breakdown.

## Species Identification

Call species from existing pipeline results or as part of `baitbench run`:

```bash
# Standalone: using pre-computed similarity from panel-qc
baitbench identify \
  --detected-detail results/run/detected_detail.tsv \
  --sample-target-map mapping.tsv \
  --target-similarity panel_qc_results/target_similarity.tsv \
  --outdir identify_results

# Standalone: compute similarity on-the-fly from target FASTA
baitbench identify \
  --detected-detail results/run/detected_detail.tsv \
  --sample-target-map mapping.tsv \
  --targets gene_targets.fa \
  --outdir identify_results

# Integrated into pipeline (genome mode)
baitbench run \
  --targets gene_targets.fa \
  --genomes full_genomes.fa \
  --distractors human.fa \
  --probes probes.fa \
  --sample-target-map mapping.tsv \
  --sample e_coli influenza_a \
  --identify \
  --num-fragments 50000 \
  --outdir results

# With stricter calling threshold (require 2 unique markers)
baitbench run \
  --targets gene_targets.fa \
  --genomes full_genomes.fa \
  --distractors human.fa \
  --probes probes.fa \
  --sample-target-map mapping.tsv \
  --sample e_coli influenza_a \
  --identify \
  --min-unique-targets 2 \
  --num-fragments 50000 \
  --outdir results
```

When `--identify` is used with `baitbench run`, the species calls are compared against the ground-truth `--sample` manifest and included in the HTML report.

## Probe Assessment

Run combined probe coverage + cross-reactivity analysis on an existing probe set:

```bash
# Basic assessment (probe coverage + self-homology)
baitbench assess-probes \
  --targets targets.fa \
  --probes probes.fa \
  --outdir assess_results

# With cross-reactivity against genomes
baitbench assess-probes \
  --targets targets.fa \
  --probes probes.fa \
  --genomes human_genome.fa other_genomes.fa \
  --threshold 80 \
  --outdir assess_results

# Skip HTML report (produce only TSV outputs)
baitbench assess-probes \
  --targets targets.fa \
  --probes probes.fa \
  --report none \
  --outdir assess_results

# Refinement: re-run coverage on low-coverage targets 3 times
baitbench assess-probes \
  --targets targets.fa \
  --probes probes.fa \
  --refine-iterations 3 \
  --refine-threshold 80 \
  --outdir assess_results

# Refinement: repeat until no targets remain below 80% 1X coverage
baitbench assess-probes \
  --targets targets.fa \
  --probes probes.fa \
  --refine-until-stable \
  --refine-threshold 80 \
  --outdir assess_results
```

Build probes and automatically assess them:

```bash
# Build + assess (default behavior)
baitbench build-probes \
  --targets targets.fa \
  --outdir probes_output

# Build + assess with cross-reactivity against genomes
baitbench build-probes \
  --targets targets.fa \
  --genomes human_genome.fa \
  --outdir probes_output

# Build only, skip assessment
baitbench build-probes \
  --targets targets.fa \
  --skip-assess \
  --outdir probes_output

# Build using CATCH method
baitbench build-probes \
  --targets targets.fa \
  --method catch \
  --probe-length 120 \
  --outdir probes_output

# Build using CATCH with custom parameters
baitbench build-probes \
  --targets targets.fa \
  --method catch \
  --catch-stride 30 \
  --catch-mismatches 3 \
  --catch-extension 10 \
  --outdir probes_output
```

## Running Individual Steps

Run pipeline steps independently for custom workflows:

```bash
# 1. Prepare
baitbench prepare \
  --targets targets.fa \
  --distractors distractors.fa \
  --distractor-fraction 0.95 \
  --outdir prep

# 2. Simulate (probe alignment + TNN scoring + multinomial sampling)
baitbench simulate \
  --reference prep/combined_reference.fa \
  --weights prep/weights.txt \
  --probes probes.fa \
  --num-fragments 50000 \
  --capture-fraction 0.5 \
  --seed 42 \
  --output prep/fragments.fa

# 3. Sequence
baitbench sequence \
  --input prep/fragments.fa \
  --read-length 150 \
  --output prep/reads.fa

# 4. Map
baitbench map \
  --reference prep/combined_reference.fa \
  --reads prep/reads.fa \
  --output prep/mapped.sam

# 5. List
baitbench list \
  --sam prep/mapped.sam \
  --output prep/detected.list

# 6. Metrics
baitbench metrics \
  --targets prep/targets.txt \
  --distractors prep/distractors.txt \
  --sample prep/sample.txt \
  --detected prep/detected.list \
  --fragments prep/fragments.fa \
  --captured prep/fragments.fa \
  --sam prep/mapped.sam \
  --run-name custom_run \
  --num-fragments 50000 \
  --output-summary prep/results.tsv \
  --output-detail prep/detected_detail.tsv
```

## Reproducible Runs

Use `--seed` for reproducibility. The same seed with the same inputs produces identical results:

```bash
baitbench run \
  --targets targets.fa \
  --distractors distractors.fa \
  --probes probes.fa \
  --num-fragments 10000 \
  --seed 42 \
  --outdir results
```

## Batch Comparisons

Run multiple configurations and compare results:

```bash
# Compare different CT values
for ct in 20 25 30 35; do
  baitbench run \
    --targets targets.fa \
    --distractors distractors.fa \
    --probes probes.fa \
    --ct $ct \
    --num-fragments 10000 \
    --seed 42 \
    --report none \
    --outdir "results_ct${ct}"
done

# Aggregate results
head -1 results_ct20/*/results.tsv > comparison.tsv
for ct in 20 25 30 35; do
  tail -1 "results_ct${ct}"/*/results.tsv >> comparison.tsv
done
```
