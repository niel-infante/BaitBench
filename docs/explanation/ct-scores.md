# CT Scores

BaitBench can express sample abundance as a CT (cycle threshold) value from a qPCR assay, rather than a raw distractor fraction. This page explains the conversion, the underlying assumptions, and when to use CT vs distractor fraction directly.

---

## What Is a CT Score?

CT (cycle threshold) is the number of PCR amplification cycles required for a fluorescence signal to cross a detection threshold in a quantitative PCR assay. It is inversely proportional to the initial amount of target DNA:

- **Low CT** (e.g., CT 15–20): abundant target — the assay detects the signal quickly
- **High CT** (e.g., CT 30–35): scarce target — many cycles are needed to amplify enough signal

CT is the standard clinical metric for reporting viral or bacterial load in diagnostic assays. BaitBench's CT support lets you answer the question: "At a patient CT of 28, can my probe panel still detect the pathogen?"

---

## CT to Distractor Fraction Conversion

BaitBench uses a one-point calibration model. You provide a single reference point (a CT value and the corresponding target fraction at that CT), and BaitBench interpolates to any other CT value using the PCR doubling model:

```
target_fraction(CT) = f_baseline × (1 + E)^(CT_baseline - CT)
distractor_fraction = 1 - target_fraction
```

Where:
- `CT_baseline` is the reference CT (default: 20)
- `f_baseline` is the target fraction at `CT_baseline` (default: 0.01, i.e., 1%)
- `E` is the PCR efficiency (default: 1.0 = 100%)

With 100% efficiency, each cycle doubles the template (×2 per cycle). With efficiency 0.9, each cycle multiplies by 1.9 instead of 2.

### Reference CT table (default calibration: CT 20 = 1% target, 100% efficiency)

| CT | Target fraction | Distractor fraction |
|----|----------------|---------------------|
| 10 | 100% (capped at 1.0) | 0% |
| 15 | ~3.2% | ~96.8% |
| 20 | 1.0% | 99.0% |
| 25 | ~0.031% | ~99.97% |
| 30 | ~0.001% | ~99.999% |
| 35 | ~3×10⁻⁵% | — |

Each 5-CT increase reduces target fraction by ~30×; each 10-CT increase reduces it by ~1000×.

---

## PCR Efficiency

Real PCR assays don't double perfectly each cycle. Efficiency E is typically 0.90–0.98 for well-optimised assays. Lower efficiency means slower amplification — the same CT at lower efficiency corresponds to a *higher* initial template concentration.

```bash
--ct 25 --ct-efficiency 0.95   # use 95% PCR efficiency
```

If you don't know the efficiency of the assay that generated your CT values, use the default (1.0). This gives a conservative (worst-case) target fraction estimate.

---

## One-Point Calibration

The default calibration (CT 20 = 1% target, 100% efficiency) may not match your assay. Override it with:

```bash
--ct 28 --ct-baseline 20 --ct-baseline-fraction 0.01 --ct-efficiency 1.0
```

Or change the reference point:

```bash
# If CT 25 corresponds to 0.1% target in your assay:
--ct 30 --ct-baseline 25 --ct-baseline-fraction 0.001
```

---

## Two-Point Calibration

If you have two (CT, target-fraction) reference points from your assay, BaitBench can derive the efficiency automatically:

```bash
--ct 30 --ct-calibration "20,0.01" "25,0.0003"
```

BaitBench computes the efficiency from the two points:

```
E = (f1 / f2)^(1 / (CT2 - CT1)) - 1
```

And uses the first point as the baseline. The derived efficiency is logged so you can verify it matches expectations.

`--ct-calibration` is mutually exclusive with `--ct-baseline`, `--ct-baseline-fraction`, and `--ct-efficiency` — all three are blocked when calibration is active.

---

## Practical CT Values and Their Meanings

For respiratory viruses in clinical nasopharyngeal swabs, typical CT values and their clinical interpretations:

| CT range | Approximate meaning | BaitBench distractor fraction |
|----------|--------------------|-----------------------------|
| < 20 | Very high viral load | < 99% |
| 20–25 | High viral load | 99–99.97% |
| 25–30 | Moderate viral load | 99.97–99.999% |
| 30–35 | Low viral load, near limit of detection | > 99.999% |
| > 35 | Near or below detection limit | — |

BaitBench is most informative for CT values in the 25–33 range, where the target is present but rare. Very low CT values (< 20) give near-100% sensitivity regardless of probe design; very high CT values (> 35) will typically show FN regardless of probe quality.

---

## CT vs Distractor Fraction: Which to Use?

Use **`--ct`** when:
- Your experimental context uses CT values (clinical diagnostics, qPCR data)
- You want to express abundance in terms of a specific assay's calibration
- You are running a `coverage-curve` sweep across CT values to find the limit of detection

Use **`--distractor-fraction`** when:
- You want to directly control the ratio of target to background without thinking about qPCR
- You are running exploratory simulations where the exact clinical correspondence doesn't matter
- You are debugging or reproducing an earlier run

The two flags are mutually exclusive: you can specify one or the other, not both. If neither is specified, BaitBench defaults to `--distractor-fraction 0.9` (90% background, 10% target).

---

## Coverage-Curve CT Sweeps

`baitbench coverage-curve` can sweep across multiple CT values to produce a sensitivity-vs-CT curve:

```bash
baitbench coverage-curve \
  --targets targets.fa \
  --distractors distractors.fa \
  --probes probes.fa \
  --ct-values 15 20 25 28 30 32 \
  --outdir ct_sweep
```

This runs the full pipeline at each CT and produces a curve showing how sensitivity degrades as the sample becomes more dilute. The resulting plot (in `ct_sweep_report.html`) shows the panel's limit of detection.
