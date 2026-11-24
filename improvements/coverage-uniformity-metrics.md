# Coverage uniformity: fold-80 and friends

Raised 2026-11-24. **Not started.** Written up and parked.

## What fold-80 is

`FOLD_80_BASE_PENALTY` from Picard's `CollectHsMetrics` is the hybrid-capture field's
standard uniformity statistic: *the fold over-coverage needed to raise 80% of target
bases to the mean coverage of those targets*. Computed as mean depth divided by the
20th-percentile depth over target bases, restricted to targets with non-zero coverage.
1.0 is perfect uniformity; a well-behaved exome lands around 1.2–1.8.

**Check Picard's `TargetMetricsCollector` source before implementing.** The mean/P20
shape and the 1.0-is-perfect convention are solid, but the exact percentile handling and
the rules for restricting to non-zero-coverage targets are worth reading rather than
reconstructing from the one-line doc string.

## Why BaitBench could report it today

The input already exists and is already being iterated over:

| piece | where |
|---|---|
| per-base depth vector per reference | `coverage::compute_coverage` (`src/alignment/coverage.rs`) |
| a function already walking that vector | `coverage::calculate_stats`, for `pct_covered_5x` / `20x` |
| run-length-encoded profile on disk | `coverage.tsv`, written by `metrics` |

Adding fold-80 is a percentile over a vector the code already holds. Nothing resembling
a uniformity metric exists anywhere in `src/`, `R/` or `docs/` today.

## The part that actually needs thought

Fold-80 is only as meaningful as the non-uniformity the simulator reproduces.

**Modelled:**
- probe tiling density and gaps — `p_retain` is a function of the probes overlapping a
  fragment, so sparse tiling genuinely produces lower coverage
- probe affinity, including GC dependence through the nearest-neighbour ΔG
- target-edge fall-off, via `--min-overlap-frac` and fragment geometry
- mappability and repeat behaviour, since reads are really mapped
- Poisson sampling noise

**Not modelled:**
- **PCR amplification GC bias.** Hybridization GC bias is in the model via ΔG;
  *amplification* GC bias is a different mechanism and is absent.
- **PCR duplicates and jackpotting** — plan Stage 6, not started.
- Fragment-length realism — plan Stage 5, not started: lengths are a clamped Normal with
  a hardwired sd, plus the contig-end mass bug found 2026-11-24.

Those first two are among the largest real contributors. **A simulated fold-80 will
therefore read systematically optimistic — too close to 1.0 — against a Picard number
from a real BAM.** Reported without that label, someone will make exactly that
comparison and conclude the tool is broken.

### The depth trap

Fold-80 has a pure-noise floor set by sequencing depth, because the 20th percentile of a
Poisson is well below its mean at low λ:

| mean depth | approximate fold-80 from sampling noise alone |
|---|---|
| 10x | ~1.4 |
| 30x | ~1.2 |
| 100x | ~1.1 |

Below reasonable depth the metric measures `--num-sequences`, not probe design. Real
Picard has the same property, which is why it is quoted at sensible depth — but
BaitBench users set depth freely, so this needs a guard that warns (or refuses) below
some mean target depth rather than silently reporting noise.

## Proposed shape: two metrics, two different claims

They have different epistemic status and should not be presented as one number.

**1. Static tiling uniformity, in `assess-probes`.** Probes-per-base variance across
target sequences. No simulation, no reads, no claim about a wet-lab number — purely a
property of the panel against its targets. Fits beside the existing gap classification
and `--gap-min-length` work, and is the honest pre-experiment answer.

**2. Simulated fold-80, in the pipeline metrics.** From `coverage.tsv`, labelled as
capture-only, with the exclusions above stated in the column documentation
(`docs/reference/output-formats.md`). Pair it with coefficient of variation or a Gini
coefficient, which are better conditioned at the depths people actually simulate —
fold-80 earns its place by being the statistic the field recognises, not by being the
best-behaved one.

**Frame #2 as comparative, not absolute.** Ranking two probesets against the same
targets is defensible, because the missing terms are approximately common to both and
largely cancel. Predicting the number a sequencing core will report is not. That
comparative framing is exactly the `build-probes` → `assess-probes` use case, so it fits
the tool's existing purpose instead of stretching it.

## The reason this is more than a nice-to-have

**It is a second validation axis for the paper.** Everything in the current revision
scores *composition* — which organisms got captured. Fold-80 scores *uniformity* — how
evenly. Real capture data exists for both validation datasets, so Picard fold-80 on the
real BAMs against simulated fold-80 would be an independent test speaking to reviewer
item 3, and **the gap between them would quantify how much of real non-uniformity is PCR
rather than probes** — a number worth having on its own account.

**Sequencing: do not attempt that validation before plan Stages 5 and 6.** Without
fragment-length realism and a duplication model, the measured gap is mostly "features
not yet implemented", which teaches nothing. The metric itself can land any time; the
validation should wait.

## Work estimate

- Metric #1 (static tiling uniformity): small. New stats in the probe-tiling half of
  `coverage.rs`, a column in the assess-probes tables, a section in
  `R/assess_probes.Rmd`.
- Metric #2 (simulated fold-80 + CV): small. Extend `CoverageStats` and
  `calculate_stats`, add columns to `detected_detail.tsv` and the `results.json` block,
  add the depth guard, document in `output-formats.md` and
  `docs/explanation/metrics-and-classification.md`.
- The real-vs-simulated validation: medium-to-large, and blocked on Stages 5–6.

## Related

- Plan Stage 6 (PCR duplicates) and Stage 7 (configurable detection and coverage
  thresholds — the hardcoded 5x/20x live in the same function fold-80 would).
- [`partial-duplex-capture.md`](partial-duplex-capture.md) — edge behaviour at fragment
  boundaries is one of the modelled contributors to simulated non-uniformity, so
  changing it moves this metric.
