# Choose a Simulation Mode

BaitBench has several mode decisions to make before running. This guide helps you pick the right combination for your experiment.

---

## Standard mode vs genome mode

### Standard mode (default)

Fragments are generated from target sequences directly. Use this when your probe targets ARE the sequences you want to detect — as is typical for viruses.

```bash
baitbench run \
  --targets virus_targets.fa \
  --distractors human.fa \
  --probes probes.fa \
  --outdir results
```

### Genome mode (`--genomes`)

Fragments are generated from full genome sequences, but reads are mapped to a smaller set of target subsequences (e.g., 16S rRNA genes). Use this when your probes target a specific region of a much larger genome — typical for bacteria.

```bash
baitbench run \
  --targets 16S_sequences.fa \
  --genomes bacterial_genomes.fa \
  --distractors human.fa \
  --probes probes.fa \
  --sample-target-map mapping.tsv \
  --outdir results
```

**Choose genome mode when:**

- Your targets are a small region of a larger genome (16S, a single gene, a genomic island)
- The full genome is megabases but probes only target a few kilobases
- You want fragments to represent the full genomic context of the sample, not just the targeted region

See the [Genome Mode tutorial](../tutorials/genome-mode.md) for a complete walkthrough.

---

## CT score vs distractor fraction

Both flags control how abundant target DNA is relative to background. They are mutually exclusive.

### Distractor fraction (`--distractor-fraction`)

Sets the fraction of simulated fragments that come from distractor sequences directly.

```bash
--distractor-fraction 0.9   # 90% background, 10% target (default)
--distractor-fraction 0.5   # 50/50 split
--distractor-fraction 0.99  # ~1% target, challenging scenario
```

Use distractor fraction when you want to directly control the mixing ratio.

### CT score (`--ct`)

Converts a qPCR cycle-threshold value to a distractor fraction using a calibrated formula. More intuitive when working with clinical specimens where abundance is expressed as a CT value.

```bash
--ct 20   # ~1% target (default calibration)
--ct 25   # ~0.03% target
--ct 30   # ~0.001% target (near limit of detection)
```

**Default calibration:** CT 20 = 1% target. Each CT unit is a 2-fold change.

| CT | Approx. target fraction |
|----|------------------------|
| 15 | 32% |
| 20 | 1% (baseline) |
| 25 | 0.03% |
| 30 | 0.001% |
| 35 | 0.00003% |

Adjust the calibration if your experimental system differs:

```bash
# If your system shows CT 25 = 0.1% target
--ct 30 --ct-baseline 25 --ct-baseline-fraction 0.001
```

Use `--ct-calibration "CT1,FRAC1" "CT2,FRAC2"` for two-point calibration from reference standards.

**Use CT scores when:**

- You want results that map directly to clinical CT measurements
- You're sweeping a range of abundances using `coverage-curve --ct-values 20 25 30 35`
- You're comparing probes across CT ranges

---

## Thermodynamic mode vs simple mode

### Thermodynamic (default)

Probe binding sites are scored using the SantaLucia (1998) nearest-neighbor free energy model. Fragments are sampled proportionally to `exp(-ΔG / RT)` — stronger-binding probes generate more reads. The hybridization temperature (`--hybridization-temperature`, default 70°C) affects how scores are weighted.

```bash
--simulate-mode thermodynamic --hybridization-temperature 70
```

This is physically motivated and produces more realistic enrichment patterns. Use it for most evaluations.

### Simple mode

All probe binding sites are weighted equally regardless of binding affinity.

```bash
--simulate-mode simple
```

**Use simple mode when:**

- You want to test uniform coverage without thermodynamic bias
- You are debugging a pipeline or input format issue and want reproducible, easy-to-reason-about output
- You don't have a good estimate of hybridization conditions

---

## When to use coverage curves

`baitbench coverage-curve` runs the full pipeline at multiple parameter values in one pass and plots coverage depth across conditions. Use it when you want to answer questions like:

- "At what CT value does coverage drop below acceptable depth?"
- "How does capture fraction affect uniformity?"
- "What sequencing depth do I need?"

```bash
baitbench coverage-curve \
  --targets targets.fa \
  --distractors distractors.fa \
  --probes probes.fa \
  --sample Influenza_A_H3N2 SARS_CoV_2 \
  --ct-values 20 25 30 35 \
  --outdir coverage_results
```

`--sample` is required for `coverage-curve`. The `--ct-values` flag is mutually exclusive with `--distractor-fraction` and `--ct`.

---

## Decision guide

| Question | Answer |
|----------|--------|
| Is your target organism a virus or small genome? | **Standard mode** |
| Are probes targeting one region of a large bacterial genome? | **Genome mode** |
| Do you have CT values from qPCR? | Use `--ct` |
| Do you want to set the mixing ratio directly? | Use `--distractor-fraction` |
| Do you want thermodynamically realistic enrichment? | **Thermodynamic** (default) |
| Do you want simple, reproducible results for debugging? | **Simple mode** |
| Do you need to sweep multiple conditions at once? | **coverage-curve** |
