# Parameter Reference

## Input Files

| Parameter | Flag | Default | Applies to | Description |
|-----------|------|---------|------------|-------------|
| Targets | `--targets`, `-t` | required | run, prepare, probe-coverage, coverage-curve | FASTA of target sequences the probes are designed to capture |
| Genomes | `--genomes`, `-g` | none | run, prepare, coverage-curve | FASTA of full genomes for fragment generation (genome mode) |
| Distractors | `--distractors`, `-d` | required | run, prepare, coverage-curve | FASTA of background sequences that should not be captured. Can be specified multiple times to provide multiple distractor files |
| Probes | `--probes`, `-p` | required | run, capture, probe-coverage, coverage-curve | FASTA of probe sequences |
| Sample | `--sample` | all targets | run, coverage-curve | Sample targets or genomes: TSV file path OR inline IDs with optional weights. See [Sample Manifest Format](reference.md#sample-manifest-format) |
| Sample-target map | `--sample-target-map` | none | run, prepare, coverage-curve | TSV mapping genome IDs to target IDs (genome mode). See [Sample-Target Map Format](reference.md#sample-target-map-format) |
| Host FASTA | `--host-fasta` | none | run, coverage-curve | Host genome for read filtering |

## Fragment Generation

| Parameter | Flag | Default | Description |
|-----------|------|---------|-------------|
| Num fragments | `--num-fragments` | 10000 | Number of fragments to generate. More fragments = better statistical power but slower |
| Fragment length mean | `--fragment-length-mean` | 175 | Mean fragment length in bp. Center of the normal distribution |
| Fragment length min | `--fragment-length-min` | 150 | Minimum fragment length in bp. Fragments shorter than this are discarded |
| Fragment length max | `--fragment-length-max` | 200 | Maximum fragment length in bp. Fragments longer than this are truncated |

## Target Abundance

| Parameter | Flag | Default | Description |
|-----------|------|---------|-------------|
| Distractor fraction | `--distractor-fraction`, `-f` | 0.9 | Fraction of fragments from distractor sequences (0-1). Higher = lower target abundance. **Mutually exclusive with `--ct`** |
| CT score | `--ct` | none | qPCR CT (cycle threshold) score. Converted to distractor fraction via calibration formula. Lower CT = more target. **Mutually exclusive with `--distractor-fraction`** |

If neither `--distractor-fraction` nor `--ct` is specified, defaults to a distractor fraction of 0.9 (10% target).

## CT Score Parameters

These parameters calibrate the CT-to-fraction conversion. Only relevant when using `--ct`.

| Parameter | Flag | Default | Description |
|-----------|------|---------|-------------|
| CT baseline | `--ct-baseline` | 20.0 | The CT value at which the target fraction equals the baseline fraction |
| CT baseline fraction | `--ct-baseline-fraction` | 0.01 | The target fraction at the baseline CT value |

See [CT Score Calculation](#ct-score-calculation) for details.

## Simulation Parameters

| Parameter | Flag | Default | Description |
|-----------|------|---------|-------------|
| Simulate mode | `--simulate-mode` | thermodynamic | `thermodynamic` (TNN Boltzmann weighting) or `simple` (uniform probe-site weights) |
| Hybridization temperature | `--hybridization-temperature` | 70.0 | Hybridization temperature in °C; only used in thermodynamic mode |
| Capture fraction | `--capture-fraction` | 0.5 | Fraction of fragments from probe binding sites (0.0–1.0); remainder are background |

## Sequencing Parameters

| Parameter | Flag | Default | Description |
|-----------|------|---------|-------------|
| Read length | `--read-length` | 120 | Trim captured fragments to this length (bp). Fragments shorter than this are kept as-is |
| Num sequences | `--num-sequences` | all | Number of reads to sample with replacement. If not set, all captured fragments become reads. Models sequencing depth control |

## Execution Parameters

| Parameter | Flag | Default | Description |
|-----------|------|---------|-------------|
| Threads | `--threads` | 1 | Number of threads for external tools (minimap2, BLAST) |
| Output dir | `--outdir`, `-o` | ./results | Output directory. A timestamped subdirectory is created for each run |
| Output prefix | `--output-prefix` | (empty) | String prepended to every auto-generated output filename. Available on `run`, `prepare`, `probe-coverage`, `coverage-curve`, `xreact`, `panel-qc`, `identify`. E.g., `--output-prefix myrun_` produces `myrun_results.tsv` instead of `results.tsv` |
| Run name | `--run-name` | auto | Custom name for the run. Default: `run_YYYYMMDD_HHMMSS` |
| Report mode | `--report` | full | Report output: `full` (render HTML), `none` (skip), `rmd` (editable RMarkdown file) |
| Seed | `--seed`, `-s` | random | Random seed for reproducibility. If not set, results vary between runs |
| Verbose | `--verbose` | false | Enable debug logging (global flag) |
| Minimap preset | `--minimap-preset` | sr | Minimap2 preset for read mapping |
| Host minimap preset | `--host-minimap-preset` | sr | Minimap2 preset for host read filtering |
| Cleanup | `--cleanup` | false | Delete intermediate files after completion, keeping only report inputs and final outputs. Available on `run`, `coverage-curve`, `probe-coverage`, and `xreact` |
| Identify | `--identify` | false | Enable species-level identification after metrics (genome mode only, requires `--sample-target-map`). Available on `run` |
| Identity threshold | `--identity-threshold` | 90.0 | Minimum sequence identity % to consider targets "similar" for species identification. Available on `run`, `panel-qc`, `identify` |
| Min unique targets | `--min-unique-targets` | 1 | Minimum unique target detections required to call a species PRESENT. Available on `run`, `identify` |

---

## CT Score Calculation

CT (cycle threshold) scores from qPCR provide an intuitive way to express target abundance. In qPCR, each cycle doubles the DNA, so a CT difference of 1 corresponds to a 2-fold change in DNA quantity. Lower CT = more target DNA.

### The Formula

```
target_fraction = ct_baseline_fraction * 2^(ct_baseline - ct)
distractor_fraction = 1 - target_fraction
```

Where:
- `ct_baseline` is a known CT value (default: 20.0)
- `ct_baseline_fraction` is the target fraction at that CT (default: 0.01)
- `ct` is the CT value you want to simulate

### Default Calibration

With defaults (`--ct-baseline 20.0`, `--ct-baseline-fraction 0.01`), the interpretation is: "at CT 20, 1% of DNA is from targets."

### CT Reference Table

| CT | Target fraction | Distractor fraction | Interpretation |
|----|-----------------|---------------------|----------------|
| 10 | 100%* | 0% | Pure target (capped at 100%) |
| 15 | 32% | 68% | Very high abundance |
| 18 | 4% | 96% | High abundance |
| 20 | 1% | 99% | Moderate (baseline) |
| 22 | 0.25% | 99.75% | Low-moderate |
| 25 | 0.031% | 99.97% | Low abundance |
| 28 | 0.004% | 99.996% | Very low |
| 30 | 0.001% | 99.999% | Near limit of detection |
| 35 | 0.00003% | ~100% | Extremely low |
| 40 | 0.000001% | ~100% | At qPCR detection limit |

*Target fractions above 100% are capped at 100% (distractor fraction = 0).

### Custom Calibration

The default calibration assumes CT 20 = 1% target. If your experimental system has different characteristics, use `--ct-baseline` and `--ct-baseline-fraction` to calibrate:

**Example: Your lab data shows CT 25 = 0.1% target reads:**

```bash
baitbench run \
  --targets targets.fa \
  --distractors distractors.fa \
  --probes probes.fa \
  --ct 30 \
  --ct-baseline 25 \
  --ct-baseline-fraction 0.001 \
  --num-fragments 10000 \
  --outdir results
```

This shifts the entire curve:
- CT 25 = 0.1% target (your calibration point)
- CT 30 = 0.003% target (5 CT higher = 32x less)
- CT 20 = 3.2% target (5 CT lower = 32x more)

**Example: Calibrate with a strong-positive sample (CT 15 = 50% target):**

```bash
baitbench run \
  --targets targets.fa \
  --distractors distractors.fa \
  --probes probes.fa \
  --ct 25 \
  --ct-baseline 15 \
  --ct-baseline-fraction 0.5 \
  --num-fragments 10000 \
  --outdir results
```

### Tips for Using CT Scores

- **Match your experimental system.** If you have empirical data linking CT to target fraction, use `--ct-baseline` and `--ct-baseline-fraction` to match your curve.
- **Use coverage-curve to sweep.** The `coverage-curve` command with `--ct-values` lets you visualize performance across a range of CT values in a single analysis.
- **Remember the log scale.** Each CT unit represents a 2-fold change. A 10-CT range spans ~1000-fold differences in abundance.
