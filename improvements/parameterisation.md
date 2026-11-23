# Parameterisation of the capture model

**Status: the next substantial piece of simulation work.** Everything below came out of
getting the per-locus gate to pass; the fix that worked exposed that the model's
absolute scale is not interpretable.

## Why this is the blocker

Capture is currently driven by equilibrium Langmuir occupancy,
`θ = K/(1+K)` with `K = [probe]·exp(−ΔG/RT)`, and `[probe]` is a fitted constant
(`DEFAULT_PROBE_CONCENTRATION`, e^-2).

That constant is **not scale-free**. `ln K = ln[probe] − ΔG(T, salt)/RT`, so the same
number lands in a different place on the saturation curve for every temperature, salt
concentration and probe length. Measured: at the Borrelia conditions (58 °C, 1 M Na+)
it reproduces per-locus capture well (log-RMSE 2.71 -> 1.02); at the minimal example's
conditions (70 °C, 50 mM) both targets saturate to identical efficiency, where the
previous model separated them 9-fold.

The deeper problem is what the constant is compensating for. A perfect 120-mer has
ΔG ~ -116 kcal/mol, so at any physically real probe concentration `K ~ e^170` and
**every well-matched probe saturates**. The fitted value is small precisely because it
is absorbing that. Equilibrium occupancy of a long duplex is close to a constant 1, so
it cannot be what discriminates loci — which is consistent with everything else
measured here: capture discrimination is kinetic (nucleation, off-rate under finite
wash), not thermodynamic.

**So this is not "pick a better constant".** It is: decide what quantity drives capture,
given that equilibrium occupancy is the wrong one.

### The trap to avoid

Anchoring the scale to the panel's own ΔG distribution — "set `[probe]` so the median
placement sits at θ = 0.5" — self-scales across conditions, and is therefore tempting.
It also re-centres whenever the ΔG basis changes, which is the **same absorption
failure** that made the original solver useless:

- Solving `[probe]` to hit `--capture-fraction` made the model a rank-threshold device.
  A 10 bp melting window changed *nothing*; `[probe]` moved 22 orders of magnitude and
  the output was bit-identical.
- Any rule that re-derives the scale from the thing being changed will do this again.

Whatever replaces it must be anchored to something that does **not** move when the
affinity basis moves.

## A. Make the occupancy scale intrinsic

Options, roughly in order of how principled they are:

1. **Replace equilibrium occupancy with a kinetic capture propensity.** Nucleation
   requires a contiguous match; off-rate under a finite wash is set by local stability.
   Both discriminate where equilibrium ΔG does not. The 10 bp melting window is already
   a crude stand-in for exactly this — it works because 10 bp is the nucleation length,
   not because the duplex melts in 10 bp units.
2. **Keep occupancy but express the scale in physical units** (nM probe, kcal/mol) and
   accept that an explicit efficiency term carries the rest. Honest about what is
   fitted, but concedes that the thermodynamics alone does not predict capture.
3. **Dimensionless reparameterisation** — carry `ln K` relative to a fixed reference
   duplex (say a perfect 60-mer at the run's conditions) so the number means the same
   thing at any temperature and salt.

**Acceptance test, all three required:** reproduce the LOD per-locus fit, do *not*
saturate the minimal example, and hold SimVal.

## A2. `--capture-fraction` and `--background-retention` are not separately identifiable

To first order `p_retain ~ p_ns + f * p_bound`, so output **composition depends only on
the ratio `f / p_ns`**. The absolute scale moves yield, and yield is never validated
against real data because the real library yield is unknown.

Measured on the D6331 mock community, with everything else held fixed:

| `--capture-fraction` | Spearman | Lin's CCC |
|---|---|---|
| 1.0 | 0.9500 | 0.8027 |
| 0.207 | 0.9098 | 0.7932 |

Those differ only by the multiplier; the model is identical. With `p_ns = 1e-3`, the
better fit corresponds to `f / p_ns ~ 1000`.

**Measured response (MOXTAA, p_ns held at 1e-3, 2026-11-23):**

| `--capture-fraction` | `f / p_ns` | Spearman | Lin's CCC |
|---|---|---|---|
| 0.05 | 50 | 0.7278 | 0.7574 |
| 0.1 | 100 | 0.8105 | 0.7806 |
| 0.2 | 200 | 0.9188 | 0.7872 |
| 0.4 | 400 | 0.9489 | 0.8009 |
| **0.7** | **700** | **0.9522** | **0.8116** |
| 1.0 | 1000 | 0.9500 | 0.8027 |

A genuine interior optimum at `f / p_ns ~ 700`, not a boundary — so `p_ns = 1e-3` is not
simply too high, and `f = 0.7` (70% of bound fragments surviving wash and elution) is
physically plausible.

But the curve is **flat**: CCC moves 0.801 -> 0.812 -> 0.803 across a 2.5x range of the
parameter. Enriched composition constrains the ratio only weakly. Fitting against the
per-organism *enrichment factor* (enriched / unenriched, section A3) should sharpen it,
since the input composition cancels.

**Do not set the harness value from this sweep.** It is fitted on MOXTAA, which the gate
then reports — the circular form. Fit against enrichment using the non-enriched arm and
validate on held-out replicates.

**Consequences:**
- No amount of fitting against composition will pin both parameters. One must come from
  protocol knowledge, or they must be reparameterised into (ratio, scale) with only the
  ratio fitted.
- Every calibrated `--capture-fraction` in the validation harness was derived under the
  *old* meaning (fraction of fragments drawn from probe sites) and does not carry over
  to the new one (fraction of probe-bound fragments surviving wash and elution).
  `tests/validation/simval/run.sh` still passes 0.207 on that basis. Re-derive before
  treating SimVal as a pass.

## A3. There IS a non-enriched arm — use it (2026-11-23)

Roadmap decision **D6 is wrong and should be reopened.** It was resolved as "TELSeq has
no non-enriched arm for the mock community, so the manufacturer's theoretical D6310
composition is the only independent input available."

`$SIMVAL_REAL` holds 24 real count files, not the 3 the gate uses. The `*NEG*` samples
are non-enriched libraries of the same mock community — 143,777 reads over 1,088
references, and their composition tracks the manufacturer's theoretical table almost
exactly:

| organism | theoretical | NEG | enriched |
|---|---|---|---|
| F. prausnitzii | 0.140 | 0.132 | 0.038 |
| V. rogosae | 0.140 | 0.142 | 0.037 |
| B. fragilis | 0.140 | 0.184 | 0.038 |

NEG vs theoretical: Pearson r(log10) = **+0.989** (n=15, E. coli strains excluded — the
theoretical table lists them separately while `groups.tsv` collapses them). The enriched
arm gives +0.931, distorted by capture as expected.

`SimVal/RESULTS.md` was already drawing fragment lengths from "the non-enriched runs",
so the data was in use without its significance being noticed.

**Why this changes the plan:**

1. **Reviewer item 3 (circularity) gets a far stronger answer.** Input weights become
   the *measured* unenriched composition, and the prediction target is the *measured*
   enriched composition. A genuine held-out prediction rather than "we used the
   manufacturer's numbers".
2. **Better than the theoretical table**, which carries none of the real extraction, GC
   or fragment-length bias of the actual prep. Those biases are not the simulation's to
   reproduce, and scoring against theory charges it for them.
3. **Gives a per-organism empirical enrichment factor** (enriched / unenriched) — the
   exact quantity the capture model predicts, and a sharper target than absolute
   composition.
4. **Breaks the A2 identifiability deadlock.** With input pinned independently,
   `f / p_ns` can be fitted against enrichment instead of against the composition that
   is then reported as validation. Fitting the ratio on enriched composition and
   reporting that same composition would have been circular — structurally the same
   criticism as reviewer item 3.

**Sample naming, confirmed from SRA (PRJNA1081843, fetched 2026-11-23).** The
BioProject holds 96 runs, 24 of them synthetic mock community — exactly the 24 local
count files, in 8 groups of 3 replicates:

| SRA experiment title | runs | local prefix |
|---|---|---|
| Non-enrichment of synthetic mock community: Negative control | SRR28150698, 699, 700, 732, 733, 735 | `*NEG*` |
| Target-enrichment: Probe type ARGs: Enrichment XT | SRR28150694, 695, 696 | `MOXTA*` |
| Target-enrichment: Probe type MGEs: Enrichment XT | SRR28150688, 689, 690 | `MOXTM*` |
| Target-enrichment: Probe type ARGs and MGEs: Enrichment XT | SRR28150691, 692, 693 | `MOXTAM*` |
| Target-enrichment: Probe type ARGs: Enrichment XTHS2 | SRR28150707, 709, 710 | `MOV2A*` |
| Target-enrichment: Probe type MGEs: Enrichment XTHS2 | SRR28150701, 702, 703 | `MOV2M*` |
| Target-enrichment: Probe type ARGs and MGEs: Enrichment XTHS2 | SRR28150704, 705, 706 | `MOV2AM*` |

So `NEG` is literally "Non-enrichment ... Negative control", `A` = ARGs (the RES panel),
`M` = MGEs (MOB), `AM` = both (combo) — consistent with `MOXTAA_res` and
`MOV2AMA_combo` in `SimVal/logs/`. `MOXT`/`MOV2` map to the XT and XTHS2 enrichment
chemistries; that last pairing is inferred from the clean 12/12 split rather than
stated outright.

The composition evidence and the metadata agree: `MOXTNEGA` is an unenriched library.

**Still to confirm:** whether a `NEG` run and an enriched run share the *same* input
library, or are separate preps of the same standard. The per-organism enrichment factor
(enriched / unenriched) is sharper if they are matched pairs. Also check whether the
Borrelia LOD series has an equivalent non-enriched arm.

## B0. The target should be enrichment factor, not composition (2026-11-23)

Pooling the six non-enriched runs gives a well-determined input composition, and
dividing the enriched arm by it gives the quantity the capture model actually predicts.

**Replicate agreement sets the noise floor.** All six non-enriched runs agree at
r(log10) >= 0.989, most >= 0.996 — and *across* chemistry (XT vs XTHS2) as closely as
within it, so all six pool into one estimate. This also confirms the local
`NEGA`/`NEGM`/`NEGAM` suffixes carry no panel meaning, as physics requires.

**Measured enrichment, MOXTAA (ARGs panel, XT) over the pooled non-enriched input:**

| organism | input | enriched | enrichment |
|---|---|---|---|
| Escherichia coli | 0.16427 | 0.73646 | 4.48 |
| Salmonella enterica | 0.00014 | 0.00024 | 1.65 |
| Lactobacillus fermentum | 0.04745 | 0.05947 | 1.25 |
| Clostridioides difficile | 0.01100 | 0.01072 | 0.97 |
| Enterococcus faecalis | 0.00002 | 0.00001 | 0.72 |
| Bifidobacterium adolescentis | 0.08136 | 0.02765 | 0.34 |
| Faecalibacterium prausnitzii | 0.12960 | 0.03828 | 0.30 |
| Veillonella rogosae | 0.13588 | 0.03691 | 0.27 |
| Roseburia hominis | 0.10524 | 0.02635 | 0.25 |
| Bacteroides fragilis | 0.17805 | 0.03764 | 0.21 |
| Prevotella corporis | 0.07755 | 0.01561 | 0.20 |
| Fusobacterium nucleatum | 0.02619 | 0.00516 | 0.20 |
| Candida albicans | 0.01188 | 0.00161 | 0.14 |
| Saccharomyces cerevisiae | 0.01127 | 0.00143 | 0.13 |
| Akkermansia muciniphila | 0.01924 | 0.00242 | 0.13 |
| Methanobrevibacter smithii | 0.00086 | 0.00006 | 0.08 |

3 enriched, 13 depleted, **60x dynamic range**. Mechanistically sensible: an ARG panel
should enrich ARG-rich lab *E. coli* and little else. The depletion of everything else
is compositional — if *E. coli* takes 74% of reads, the rest must shrink.

**Why this is the right objective:**
- It isolates the capture effect from input composition, which is what the model
  predicts; scoring absolute composition also scores how well the input weights were
  guessed.
- Non-circular: input comes from the non-enriched arm, not from the data being
  predicted. Fitting `f / p_ns` to enriched composition and then reporting that same
  composition was the circular form (section A2).
- Sharp: 60x signal against an r ~ 0.996 noise floor, versus a flat CCC response that
  barely moved over a 2.5x parameter range.
- It answers reviewer item 3 directly.

**Held-out split available:** 3 replicates x 3 panels x 2 chemistries = 18 enriched
runs. Fit on one panel/chemistry, validate on the rest. The panels differ (ARGs, MGEs,
both), so generalising across them is a real test rather than a replicate check.

## B. Decide what we optimise against

The current objective is log-RMSE over eight Borrelia loci. Two parameters have already
been fitted to it. That is thin, and the harness has twice produced results that did not
survive contact with the pipeline.

Needed:

- **A joint objective** over LOD titration, LOD per-locus and SimVal composition, so
  parameters are fitted rather than hand-swept.
- **Held-out data.** The two unused SimVal replicates and the eight unused LOD
  conditions are the natural split.
- **A weighting decision.** Titration spans ~5 orders of magnitude, per-locus ~1,
  composition ~5. An unweighted log-RMSE silently optimises whichever has most spread.
- Extend `tests/locus_calibration.rs`, which currently scores one dataset. Keep its
  caveat in view: its first version averaged over all starts spanning a locus while
  production averages over covered positions only, and that mismatch produced a
  confident 29x improvement that the pipeline did not show.

B is tempting to defer as bookkeeping. It is what makes A's answer believable.

## C. Molecule counts instead of `--num-fragments`

`--num-fragments` is really sequencing depth, and depth is a real choice. The deeper
fiction is that there is **no molecule count at all**: weights are continuous densities,
so a reference is an infinitely divisible source of fragments.

At NTC the real sample holds tens of *B. burgdorferi* genome copies and yielded 37
reads; the model draws 20,000 fragments from a smooth profile. Poisson sampling noise,
locus dropout and jackpotting are all absent — in precisely the regime the LOD gate is
meant to validate.

Replace input weights with **copy numbers**, keep depth separate. Independent of A and
B, and the change most likely to improve the low end. Overlaps Stage 6 (PCR duplicates)
of the simulation-rework plan.

## D. Prerequisites for trusting any of it

- [x] **Re-run the gates after the multiplier landed** (2026-11-21, fast subsets,
      5 threads). Titration 0.00010 / 0.00205 / 0.83475 against real 0.0000416 /
      0.00109 / 0.79565; per-locus log-RMSE 2.71 -> 1.00; SimVal 0.9098 / 0.7932.
- [ ] **Re-derive the harness `--capture-fraction` values.** Added by the re-run above.
      `simval/run.sh` still passes 0.207, calibrated under the old meaning (fraction of
      fragments drawn from probe sites). Under the new meaning the same model scores
      0.9500 / 0.8027 at `cf=1.0` versus 0.9098 / 0.7932 at 0.207, so the SimVal gate
      is currently measuring the wrong scenario rather than the model. Note A2: only
      `f / p_ns` is identifiable from composition, so this is a choice about the ratio.
- [ ] **Full 11-condition LOD and all three SimVal replicates**, not the fast subsets.
      Everything measured so far is 3 LOD conditions and 1 replicate.
- [ ] **Decide whether `--capture-fraction 0.5` remains a sensible default** now that it
      means protocol efficiency (wash stringency, hybridisation time, elution losses)
      rather than a solved target.

## Suggested order

**D -> A -> B -> C.** D is cheap now that probe placements are cached and gives a
trustworthy baseline. Nothing else is interpretable until A fixes the scale. B stops the
hand-fitting before more parameters accumulate. C is independent and can be taken
whenever.

## What is already settled, so it is not re-litigated

- `--capture-fraction` is a **uniform multiplier** on specific capture, not a solved
  target. It is the scenario knob — how good or bad you expect the run to be — paired
  with `--num-sequences` for depth. It provably cannot reshape relative capture between
  loci.
- Non-specific carryover (`--background-retention`) is a genuinely separate physical
  path, combined as a union of independent events rather than a second sampling arm.
  Splitting it into its own arm is what capped the realized distractor fraction at
  `1 − capture_fraction` (reviewer item 6).
- Tested and rejected, with evidence, in `~/.claude/plans/`: the two-hit seed heuristic,
  longer seeds, a dangling-arm/flap penalty, probe density alone, GC bias, lower
  calibration operating points, and swapping the probe-combination rule for a Poisson
  sum.
