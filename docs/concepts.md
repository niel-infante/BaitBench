# Concepts

## Pipeline Overview

BaitBench runs a multi-step simulation pipeline:

1. **Prepare** -- Combine target and distractor sequences; generate sampling weights
2. **Simulate** -- Align probes to reference; score binding sites by TNN thermodynamics; generate probe-biased fragments + background (controlled by `--capture-fraction`)
3. **Sequence** -- Trim fragments to read length; optionally sample to model sequencing depth
4. **Filter** (optional) -- Remove reads mapping to a host genome
5. **Map** -- Align reads back to reference sequences
6. **List** -- Count reads per reference
7. **Metrics** -- Classify each reference as TP/FP/FN/TN; compute summary statistics
8. **Report** (optional) -- Generate HTML report with figures

The `baitbench run` command chains all steps automatically. Each step is also available as a standalone subcommand for custom workflows.

## Standard Mode vs Genome Mode

**Standard mode** (default): Fragments are generated from target sequences and distractors. Reads are mapped back to the same sequences. Use this for viruses and other small genomes where the probe target IS the genome.

**Genome mode** (`--genomes`): Fragments are generated from full genome sequences, but reads are mapped back to probe target subsequences. Use this for bacteria and other large pathogens where probes target specific gene regions (e.g., 16S rRNA) rather than the whole genome. A `--sample-target-map` links genome IDs to their target IDs.

## Sample Manifest

The `--sample` flag specifies which targets (or genomes) are "present" in the simulated specimen. Only sample entries generate fragments; non-sample targets become negatives that should NOT be detected.

Without `--sample`, all targets (or genomes) are treated as present with equal weight. This tests basic capture efficiency. With `--sample`, the tool tests discrimination -- can the probes detect sample targets while rejecting non-sample targets within the same panel?

See [Input File Formats](reference.md#sample-manifest-format) for syntax details.

## 3-Way Classification

BaitBench classifies each reference sequence into one of three categories, then evaluates detection:

| Category | Detected | Classification | Meaning |
|----------|----------|----------------|---------|
| Sample target | Yes | **TP** | Correctly detected |
| Sample target | No | **FN** | Missed detection |
| Non-sample target | Yes | **FP_target** | Cross-reactive within panel |
| Non-sample target | No | **TN_target** | Correctly rejected |
| Distractor | Yes | **FP_distractor** | Off-target capture |
| Distractor | No | **TN_distractor** | Correctly rejected |
| Untargeted genome | -- | **untargeted** | No expected target (genome mode only) |

This distinguishes two types of false positives:
- **FP_target**: Cross-reactivity within the target panel (e.g., probe for virus A captures virus B)
- **FP_distractor**: True off-target capture (e.g., probe captures bacterial DNA)

Without `--sample`, all targets are in the sample, reducing to a 2-way classification (TP/FP/FN/TN with no FP_target).

## CT Scores

CT (cycle threshold) scores from qPCR provide a natural way to express target abundance. BaitBench converts CT values to distractor fractions using a calibrated exponential formula. Lower CT = more target DNA = easier to detect.

See [CT Score Calculation](parameters.md#ct-score-calculation) for the formula, default calibration, and how to customize it.

## Capture Fraction and Thermodynamic Simulation

`--capture-fraction` (default 0.5) controls what fraction of simulated fragments come from probe binding sites. The remaining fraction are background fragments drawn uniformly by sequence weight × length.

Probe binding sites are scored using the SantaLucia (1998) nearest-neighbor thermodynamic model: ΔG is computed from consecutive Watson-Crick stacking interactions along each probe-reference alignment, and the Boltzmann factor `exp(-ΔG / RT)` weights sampling toward high-affinity sites. Use `--simulate-mode simple` to skip TNN scoring and use uniform weights instead (no temperature required).

Target enrichment is emergent from TNN affinity × sequence weights rather than being imposed post-hoc — sequences with weight 0.0 (non-sample targets) never generate probe-biased fragments. Fold enrichment is no longer a parameter.

### Thermodynamic model details

The Boltzmann score `exp(-ΔG / RT)` is used as a relative sampling weight — only the ratio between scores matters, not the absolute value. A score of **1** (ΔG = 0) is the neutral baseline: no thermodynamic tendency to bind. Scores above 1 represent favorable hybridization; the model never produces scores below 1.

ΔG has three contributions:

**1. Nearest-neighbor stacking.** The dominant term. Each pair of consecutive Watson-Crick base pairs contributes a stacking free energy from the SantaLucia (1998) Table 2 parameters. Stacking is always stabilising (negative ΔG). A mismatch breaks the stacking chain — BaitBench uses the SkipStacking strategy, meaning mismatches interrupt but do not reverse accumulated stacking energy.

**2. Initiation penalty.** Forming a duplex from two separate strands carries a one-time nucleation cost: the strands lose translational and rotational freedom on association, and the terminal base pairs are less constrained than interior ones. This is captured by the SantaLucia initiation parameters, applied once for each end of the duplex. AT termini are more costly (+2.3 kcal/mol enthalpy) than GC termini (+0.1 kcal/mol) because A-T pairs have two hydrogen bonds versus three for G-C and are weaker duplex anchors.

The initiation penalty and salt correction are only applied when at least one stacking step exists — i.e., when there are at least two consecutive Watson-Crick pairs. A single isolated complementary base pair flanked by mismatches cannot sustain a stable duplex, so there is no nucleation event to penalise and no duplex over which to apply a salt correction. Without this guard, isolated complementary positions at low salt would produce scores below 1 (implying the probe is repelled from that site), which is unphysical.

**3. Salt correction.** At salt concentrations below 1 M Na+, the negatively charged DNA backbone is less shielded, making hybridization less favorable. BaitBench applies the Owczarzy et al. (1997) correction `ΔS += 0.368 × (n_wc − 1) × ln([Na+])`. At 1 M Na+ this term is zero; at 50 mM Na+ (a typical hybridization buffer, controlled by `--salt-concentration`) it meaningfully weakens binding, with the effect scaling with the number of Watson-Crick pairs. Set `--salt-concentration 1000` to disable the correction.

### Comparison with RAmpSim

RAmpSim (Rooney et al. 2025) uses the same SantaLucia stacking table and SkipStacking strategy but omits initiation and salt correction entirely. For relative scoring, the constant part of the initiation penalty cancels when normalising sampling weights across probes. BaitBench preserves the AT vs GC terminal distinction (~2 kcal/mol effect on ΔG) because it is physically real: GC-terminated duplexes are genuinely more stable and should be sampled more frequently. The salt correction has a larger practical impact — at 50 mM Na+ it can shift scores by tens of kcal/mol for longer probe hits — and meaningfully re-ranks long vs short binding sites relative to omitting it.

## Weight Calculation

Sampling weights determine how many fragments each sequence generates. The number of fragments from a sequence is proportional to `weight * sequence_length`.

**Standard mode:**
- Sample targets: weight from sample manifest (default 1.0)
- Non-sample targets: weight = 0 (no fragments)
- Distractors: calculated to achieve the requested distractor fraction

**Genome mode:**
- Sample genomes: weight from sample manifest (default 1.0)
- Non-sample genomes: weight = 0
- Distractors: same formula as standard mode

The distractor weight formula ensures the requested fraction of total fragments come from distractors:

```
distractor_weight = (distractor_fraction * total_sample_weight) / (n_distractors * (1 - distractor_fraction))
```

## Sequence ID Conventions

Sequence IDs are taken from the first whitespace-delimited word of each FASTA header (everything after `>` up to the first space). These IDs must be unique within each file and consistent across input files.

**Sequence names must not contain spaces.** Use underscores or other delimiters: `>Zika_virus` not `>Zika virus`.

Fragment names follow the pattern `{seq_id}_fragment_{n}`, using the last occurrence of `_fragment_` as the delimiter. This allows sequence IDs to contain the substring `_fragment_` without ambiguity.

---

## Pipeline Flowcharts

### Standard Mode

```
INPUT FILES                    STEP                          OUTPUT FILES
=============                  ====                          ============

targets.fa ──────────┐
                     │
distractors.fa ──────┤
                     ├──── 1. PREPARE ──────────────────── combined_reference.fa
sample (optional) ───┤         │                            weights.txt
                     │         │                            targets.txt
--distractor-fraction│         │                            distractors.txt
  or --ct ───────────┘         │                            sample.txt
                               │
                               ▼
combined_reference.fa ──┐
probes.fa ──────────────┤
weights.txt ────────────┤ 2. SIMULATE ────────────────── fragments.fa
--num-fragments ────────┤      │      (probe-biased + background)
--capture-fraction ─────┤      │
--simulate-mode ────────┤      │
--hybridization-temp ───┤      │
--fragment-length-* ────┤      │
--seed ─────────────────┘      │
                               │
                               ▼
fragments.fa ───────────┐
                        ├─ 3. SEQUENCE ──────────────────── reads.fa
--read-length ──────────┤      │
--num-sequences ────────┤      │
--seed ─────────────────┘      │
                             │
                    ┌────────┴────────┐
                    │  --host-fasta   │
                    │   specified?    │
                    └──┬──────────┬───┘
                   yes │          │ no
                       ▼          │
reads.fa ───────┐                 │
host.fa ────────┤ 4. FILTER       │
--host-minimap- ┤    │            │
  preset ───────┘    │            │
                     ▼            │
              filtered.fa         │
                     │            │
                     ▼            ▼
              (filtered or reads).fa
                     │
combined_            ├─ 5. MAP ────────────────────────── mapped.sam
  reference.fa ──────┤      │
--minimap-preset ────┘      │
                            │
                            ▼
mapped.sam ──────────── 6. LIST ───────────────────────── detected.list
                            │
                            ▼
targets.txt ─────────┐
distractors.txt ─────┤
sample.txt ──────────┤
detected.list ───────┤ 7. METRICS ────────────────────── results.tsv
fragments.fa ────────┤                                    detected_detail.tsv
fragments.fa ────────┤                                    results.json
mapped.sam ──────────┘                                    coverage.tsv
                            │
                            ▼
results.tsv ─────────┐
detected_detail.tsv ─┤ 8. REPORT (optional) ──────────── report.html
run_params.tsv ──────┤
coverage.tsv ────────┘
```

### Genome Mode

Genome mode adds a separate mapping reference and genome-aware classification:

```
INPUT FILES                    STEP                          OUTPUT FILES
=============                  ====                          ============

targets.fa ──────────┐
                     │
genomes.fa ──────────┤
                     │
distractors.fa ──────┤
                     ├──── 1. PREPARE ──────────────────── combined_reference.fa
sample ──────────────┤                                        (genomes + distractors)
                     │                                      mapping_reference.fa
sample-target-map ───┤                                        (targets + distractors)
                     │                                      weights.txt
--distractor-fraction│                                      targets.txt
  or --ct ───────────┘                                      distractors.txt
                                                            genomes.txt
                                                            sample.txt
                                                            sample_target_map.txt

    Steps 2-4 are identical to standard mode:
      - Simulate uses combined_reference.fa (genomes + distractors); probes align to genomes

                     ... (steps 2-4) ...

              (filtered or reads).fa
                     │
mapping_             ├─ 5. MAP ────────────────────────── mapped.sam
  reference.fa ──────┤       (targets + distractors)
                     │
                            │
                            ▼
                     ... (step 6 same) ...
                            │
                            ▼
targets.txt ─────────┐
distractors.txt ─────┤
sample.txt ──────────┤
sample_target_map ───┤ 7. METRICS ────────────────────── results.tsv
detected.list ───────┤   (genome-aware classification)    detected_detail.tsv
fragments.fa ────────┤                                    results.json
fragments.fa ────────┤                                    coverage.tsv
mapped.sam ──────────┘

                     ... (step 8 same) ...
```

Key differences in genome mode:
- **combined_reference.fa** = genomes + distractors (fragments generated from full genomes)
- **mapping_reference.fa** = targets + distractors (reads mapped to target regions)
- A read from genome G mapping to target T is correct if T is linked to G in the sample-target-map
- Untargeted genomes (no target mapping) are tracked separately and do not affect TP/FP/FN/TN
