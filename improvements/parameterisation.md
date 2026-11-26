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

### A1. Partly answered: the scale is already behaving intrinsically (2026-11-24)

The worry was that `e^-2` is conditions-specific — a per-dataset fudge absorbing
temperature rather than a concentration. Swept it against the B0b enrichment objective,
on a dataset it was never fitted to: different organisms, different panel (ARGs vs
Borrelia MLST), 65 °C instead of 58 °C, 5.8 kb fragments instead of 225 bp, and a
different objective function.

| ln [probe] | Spearman rho | r(log10) | log10 RMSE | centred RMSE |
|---|---|---|---|---|
| -8 | +0.919 | 0.860 | 0.639 | 0.464 |
| -6 | +0.943 | 0.869 | 0.588 | 0.413 |
| -4 | +0.936 | 0.900 | 0.411 | 0.281 |
| -3 | +0.936 | 0.914 | 0.297 | 0.217 |
| -2.5 | +0.936 | 0.917 | 0.256 | 0.200 |
| **-2 (shipped)** | **+0.936** | **0.914** | **0.233** | **0.195** |
| -1.5 | +0.936 | 0.907 | 0.228 | 0.200 |
| -1 | +0.936 | 0.895 | 0.235 | 0.213 |
| 0 | +0.925 | 0.858 | 0.268 | 0.247 |
| +2 | +0.890 | 0.791 | 0.332 | 0.300 |
| +4 | +0.890 | 0.764 | 0.365 | 0.319 |

**`e^-2` sits at the minimum of the shape error on a dataset that did not fit it.** The
basin is shallow locally — anything from `e^-2.5` to `e^-1.5` is within 0.005 of the
best — and rises clearly outside it, by 44% at `e^-4` and 27% at `e^0`. Spearman rho is
*exactly* +0.936 across the whole `e^-3 … e^-1` band: within the basin the concentration
moves the scale, not the ranking.

The same sweep over `--melting-window`, again never fitted to this dataset:

| window | Spearman rho | r(log10) | log10 RMSE | centred RMSE |
|---|---|---|---|---|
| 0 (whole duplex) | +0.886 | 0.748 | 0.389 | 0.326 |
| **10 (shipped)** | **+0.936** | **0.914** | **0.233** | **0.195** |
| 15 | +0.930 | 0.831 | 0.285 | 0.272 |
| 20 | +0.886 | 0.787 | 0.328 | 0.305 |
| 30 | +0.916 | 0.749 | 0.371 | 0.334 |
| 60 | +0.881 | 0.749 | 0.383 | 0.325 |

Both constants, fitted on eight Borrelia loci, come out independently optimal here. The
window is a *sharp* optimum — 15 bp is already 40% worse on shape — and it lands on the
canonical nucleation length.

This also **reverses the earlier reading that the window hurt SimVal** (CCC 0.791 ->
0.586). That was measured under the composition objective with `[probe]` solved per run,
which section A2 shows is the circular configuration; on the non-circular enrichment
objective with `[probe]` fixed, the window is clearly right.

**What this does and does not settle.** It is real cross-dataset validation of both
fitted constants, which is what "still provisional" was waiting on, and it is evidence
the scale behaves like a concentration rather than a temperature fudge, since one value
holds at 58 °C and at 65 °C. It does *not* make the parameterisation intrinsic in the
sense option 3 means: the number is still in arbitrary units absorbed into `K`, and two
datasets at two temperatures do not establish that it transfers to a third. Option 3
(carry `ln K` relative to a fixed reference duplex) remains worth doing, but as
tidying-up rather than as a fix for a measured failure.

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

3 enriched, 13 depleted, **60x dynamic range** (63x as the harness scripts compute it — see B0b). Mechanistically sensible: an ARG panel
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

## B0b. First predicted-vs-measured enrichment: rho = +0.94 (2026-11-24)

Ran the model against the B0 target. Input weights from the pooled non-enriched arm
(section A3), ARGs panel (`res_probes.fasta`), MOXTAA as the enriched comparison,
40 000 fragments, `--capture-fraction 0.7`, 5 threads.

Because the weights encode the non-enriched composition, the simulation's own input
share reproduces B0's input column exactly (ratio 1.000 for all 17 organisms, verified
from `weights.txt` x reference lengths). So any departure of the output from 1.0 is
capture, and predicted enrichment is directly comparable to measured.

| organism | input | predicted | measured | pred/meas |
|---|---|---|---|---|
| Escherichia coli | 0.16329 | 3.778 | 4.510 | 0.84 |
| Salmonella enterica | 0.00016 | 1.745 | 1.520 | 1.15 |
| Lactobacillus fermentum | 0.04576 | 0.800 | 1.299 | 0.62 |
| Clostridioides difficile | 0.01152 | 0.820 | 0.930 | 0.88 |
| Bifidobacterium adolescentis | 0.08054 | 1.063 | 0.343 | 3.10 |
| Faecalibacterium prausnitzii | 0.12984 | 0.580 | 0.295 | 1.97 |
| Veillonella rogosae | 0.13696 | 0.479 | 0.269 | 1.78 |
| Roseburia hominis | 0.10362 | 0.404 | 0.254 | 1.59 |
| Bacteroides fragilis | 0.17968 | 0.203 | 0.209 | 0.97 |
| Prevotella corporis | 0.07741 | 0.218 | 0.202 | 1.08 |
| Fusobacterium nucleatum | 0.02733 | 0.287 | 0.189 | 1.52 |
| Akkermansia muciniphila | 0.01922 | 0.283 | 0.126 | 2.25 |
| Saccharomyces cerevisiae | 0.01142 | 0.112 | 0.125 | 0.90 |
| Methanobrevibacter smithii | 0.00089 | 0.175 | 0.072 | 2.42 |
| Enterococcus faecalis | 0.00001 | 1.168 | 1.063 | 1.10 |
| *Candida albicans* | 0.01235 | 0.029 | 0.130 | excluded |

**On the 15 scorable organisms: Spearman rho = +0.939, Pearson r(log10) = +0.918,
log10 RMSE = 0.230 (1.70x), centred 0.193, bias +0.126.** Against a 63x measured
dynamic range and an r >= 0.989 replicate noise floor, that is the first time the model
has been scored on a non-circular quantity and come out meaningfully better than chance.

*E. faecalis* is worth singling out: at an input share of 1.1e-5 — four orders of
magnitude below the dominant organism — it is predicted at 1.17 against a measured 1.06.
Inflating rare organisms was one of the original failures (*E. faecalis* 1360x under the
quota model), so this is the specific regression that is now gone.

Scoring the **sampled reads** instead of the analytic prediction gives +0.917 / +0.900 /
0.250 / 0.212 on 14 organisms, and there *E. faecalis* must be dropped: at this input
share 40 000 fragments give 0.44 expected fragments, so the run cannot sample it. That
is a depth limit of the fast harness, not of the model — raise `--num-fragments` and it
comes back. The remaining analytic-vs-sampled gap is Poisson noise plus a flat 8%
mapping loss.

Reproduced by `tests/validation/simval/input_weights.py` (pools the non-enriched arm,
derives weights and input shares) and `enrichment.py` (scores predicted against
measured). The numbers above come from those scripts, so they differ in the last digit
from the first hand-computed pass.

For scale, the composition-based comparison this replaces gave Lin's CCC 0.80 while
simultaneously under-predicting *E. coli* 3x and inflating *S. enterica* 135x — the two
facts could coexist because composition is dominated by whether the input weights were
guessed well.

**One organism excluded, for a diagnosed reason, not for being inconvenient:**

- *Candida albicans* — **the reference assembly is 13 062 contigs over 12.58 Mb**, mean
  963 bp, against a 5.8 kb fragment length. Only **25.4%** of its bases sit on a contig
  long enough to start a mean-length fragment, versus 99.7-99.8% for every other
  organism in the panel. Predicted/floor = 0.028/0.115 = **0.244** against an available-
  starts fraction of **0.254** — the shortfall is the assembly, to two significant
  figures, with nothing left for the capture model to explain.

The exclusion was checked rather than assumed. Per-organism mapping rate is flat at
88-94% across all 15 detected organisms, so the loss is not a mapping artifact; and the
realised input composition matches the intended one at ratio 1.000 for all 17 organisms
(computed from `weights.txt` x reference lengths), so it is not the weights.

### What the residual says

The first guess was compositional closure: if *E. coli* is under-captured (3.78 vs 4.48)
then, because enrichments must satisfy `sum(input x enrichment) = 1` (measured 0.9996,
predicted 1.0000), every other organism has to be over-predicted — by
`(1-0.620)/(1-0.736) = 1.44x`. **That explains the centre and not the spread.** The
geometric mean of the 13 non-*E. coli* ratios is 1.21, close to the prediction, but they
run from 0.62 to 3.10, and centring the log residuals (the CLR transform, which removes
closure exactly) barely moves the error: 0.230 -> 0.193. So most of the error is real
shape, not an artifact of the constraint.

The shape error has a direction: **the model does not deplete the untargeted bacteria
enough.** Bifidobacterium 3.10x, Methanobrevibacter 2.42x, Akkermansia 2.25x,
Faecalibacterium 1.97x. In the real data these four sit at 0.08-0.34, which is close to
where an organism with *no* probe binding at all should sit — the carryover floor is
`p_ns / yield = 1e-3 / 8.699e-3 = 0.115`.

### A fast way to test this, built first

Predicted enrichment does not need the pipeline. The sampler draws a reference with
probability proportional to `weight x retention_mass`, so an organism's expected output
share is its pooled `w x mass`, and the whole table is computable from
`fragments.retention.tsv`. `simulate` alone with a warm `--probe-hits` cache therefore
replaces a ~25 minute run with ~90 seconds, and because it is analytic it carries no
sampling noise, and reaches organisms the sampled run is too shallow to see. Verified
against the sampled run on the same 14 organisms: rho +0.939 vs +0.917, log10 RMSE 0.238
vs 0.250, per-organism ratios within ~5% (the gap is Poisson noise plus the flat 8%
mapping loss). Script: `tests/validation/simval/enrichment.py`, with `--reads` to score
the sampled side instead.

One caveat it caught immediately: `simulate` run without
`--fragment-length-mean/-min/-max` does not reproduce a `run` that set them, and the
retention profile changes enough to invert conclusions (yield 0.13% vs 0.87%, log RMSE
0.564 vs 0.233). Any sweep like this must carry an unfiltered baseline row and check it
against the reference run before the other rows are read.

### The predicted enrichment does not depend on the input composition

Found while moving the input-weights derivation into the harness. Writing `w_o = s_o /
L_o` (share over group length, which is what the weights must be — see
`input_weights.py`) and substituting into the sampler's draw probability:

```
predicted_share_o  = w_o M_o / sum_j w_j M_j
predicted_enrich_o = predicted_share_o / s_o = (M_o / L_o) / sum_j (s_j M_j / L_j)
                   = efficiency_o / mean_efficiency
```

The share `s_o` cancels. **Predicted enrichment is retention mass per base, divided by
one composition-weighted scalar common to every organism.** So all relative predicted
enrichments — the ranking, and every ratio between organisms — are completely
independent of the input composition, and only the common scale depends on it.

Measured, as a check rather than a derivation: pooling the six non-enriched runs by
count instead of averaging their shares moves individual input shares by up to 5% (and
much more for the two organisms near one read), and the score does not move — rho
+0.9361 both ways, r(log10) +0.9143 both ways, centred RMSE 0.195 both ways, log10 RMSE
0.233 vs 0.234. Only the bias shifts, by 0.002.

Three consequences:

1. **The objective is maximally non-circular.** The input composition enters only as the
   denominator of the *measured* enrichment. It plays no part in the prediction, so
   there is no route by which guessing the input well could flatter the model.
2. **It is robust to how the non-enriched arm is pooled**, which was the one judgement
   call in deriving the input. Pooling by count is still the better estimator — it is
   the ML estimate under a shared multinomial — but nothing hinges on it.
3. **Shape can be scored straight from `fragments.retention.tsv`**: the `efficiency`
   column *is* the predicted enrichment up to that one scalar. Verified against the
   table above — *E. coli* 0.03231/0.008699 = 3.71 vs 3.76 predicted, *Akkermansia*
   0.00243/0.008699 = 0.279 vs 0.282, *M. smithii* 0.00157/0.008699 = 0.180 vs 0.183.

### Two candidate causes, both rejected by measurement

> The sweep tables below, and those in A1 and D, were scored before the input
> derivation moved into `input_weights.py`, so their absolute values sit ~0.005 from
> B0b's headline. Every one of them is a *relative* comparison, and the
> composition-independence result above shows the ranking cannot depend on that
> derivation — rho was bit-identical across both versions. The baseline row in each
> table is the shipped configuration, so read the rows against it, not against B0b.

**1. Spurious cross-reactive placements.** The residual does not track probe signal:
`resid` vs per-organism specific fraction r = +0.026 (+0.277 in log), vs covered-start
fraction r = +0.104. Direct inspection then made the placement story look compelling
anyway — the three organisms the real data puts at or below the carryover floor
(Saccharomyces 0.13, Akkermansia 0.13, Methanobrevibacter 0.08) are *exactly* the three
with zero placements at >= 110 bp and >= 95% identity, while every enriched organism has
them; and median aligned length is 78-96 bp on untargeted organisms against 118 bp on
*E. coli*. A hand-built density statistic confirmed it: measured enrichment correlates
with `hits/Mb` at r(log) = +0.595, rising monotonically to +0.784 when restricted to
hits >= 110 bp at >= 95% identity. The perfect-identity hits on untargeted organisms are
not low-complexity artifacts either — they are 16S rRNA, which a conserved region makes
genuinely universal.

So: filter the cache and re-simulate. **It does not work.**

| min aligned frac | min identity | hits kept | Spearman rho | log10 RMSE | centred | bias |
|---|---|---|---|---|---|---|
| — | — | 100.0% | +0.936 | 0.233 | 0.195 | +0.128 |
| 0.50 | 0.90 | 42.1% | +0.899 | 0.215 | 0.195 | +0.091 |
| 0.75 | 0.90 | 37.6% | +0.870 | 0.208 | 0.203 | +0.046 |
| 0.90 | 0.90 | 34.4% | +0.901 | 0.210 | 0.208 | +0.030 |
| 0.90 | 0.95 | 27.6% | +0.901 | 0.215 | 0.215 | +0.015 |
| 0.90 | 0.98 | 23.4% | +0.903 | 0.206 | 0.205 | +0.012 |
| 0.95 | 0.95 | 27.2% | +0.822 | 0.219 | 0.219 | -0.001 |

Raw RMSE improves ~10%, but **every bit of it is the bias going to zero** — centred
RMSE is flat to worse (0.195 -> 0.219) and Spearman rho *degrades* (0.936 -> 0.822).
Discarding three quarters of the placements barely moves composition, because with
5.8 kb fragments and `--min-overlap-frac 0.5` a single hit covers ~2.9 kb of start
positions: roughly one hit per 4 kb saturates coverage, so hit density is not the
model's discriminating variable even though it is the hand-built statistic's. This is
the lesson of correction #8 again — a statistic that correlates with the target is not
evidence about the mechanism until the model is run with it changed.

**2. The 10 bp melting window as a compressor.** A probe at 85% identity almost always
contains some perfectly matched 10-mer, so a window criterion is length-blind by
construction and should make marginal placements look capture-competent. Swept it
(table in section A1): window 10 is not the compressor, it is the discriminator —
**the best value on this dataset too**, by 40% in shape over 15 bp and 67% over the
whole duplex.

**The residual's cause is therefore still open**, and it is now known not to be
placement stringency, not the melting window, not probe concentration (A1: identical
ranks across the basin), not `--capture-fraction` (D: identical ranks across a 10x
range), not compositional closure (centring barely moves it), not mapping, and not the
input weights. What remains untested is the combination rule itself at fragment scale:
the product over probes overlapping a 5.8 kb fragment saturates for any organism with a
few hundred scattered hits, which would compress the top and lift the bottom — the
observed pattern. That is the same "product form is inherently flattening" concern the
plan recorded for Borrelia per-locus capture, reappearing at organism scale, and it is
the one hypothesis the sweeps above cannot address because every one of them leaves the
product rule intact.

### A real defect found on the way: fragmented assemblies lose mass silently

The Candida diagnosis is not only about this dataset. Retention mass counts fragment
**start positions** (`n_starts = L - f + 1`, with `f` clamped to `L` for short contigs,
giving exactly 1 start), while the input weights encode DNA **mass** (proportional to
`L`). The two agree to 0.1% when `L >> f` and diverge badly when they are comparable:

| | bases | starts@5.8kb | realised share |
|---|---|---|---|
| closed bacterial chromosome | 2-5 Mb | 99.7-99.9% | as intended |
| Candida albicans draft | 12.58 Mb / 13 062 contigs | 25.4% | 0.25x intended |

`sample_retention_fragments` bails only when *every* sequence is unsampleable; losing
three quarters of one organism is silent. Draft assemblies are the norm in real
reference panels, so this will hit users. Two things to do, recorded in the plan's
Stage 5 (it is the same contig-end arithmetic as D10):

1. **Warn.** Report the fraction of weighted input mass that cannot yield a
   mean-length fragment, per reference and in aggregate, whenever it is non-trivial.
2. **Fix the arithmetic.** Make mass proportional to `L` for every contig rather than to
   `L - f + 1`, which also removes the inward bias at contig ends that D10 already
   flags. The uniform factor `f` cancels in composition, so only the end effect matters.

## B0e. The mapper confound is real but negligible (2026-11-25)

Every SimVal comparison scores BaitBench-mapped simulated reads against count tables
produced by *external* minimap2 (`map-hifi`, `--secondary=no`, per
`SimVal/scripts/map_real_reads.sh`), while BaitBench mapped with rammap's `sr` default —
a **short-read preset applied to 5.8 kb reads**. Three differences stacked: aligner,
preset, and read error model. All of it sat inside the residual, unattributed.

With the raw reads in hand this is directly measurable. Mapping MOXTAA's **real** CCS
reads against the same reference, varying only the mapper:

| organism | external mm2 `map-hifi` | bb `sr` | bb `map-hifi` |
|---|---|---|---|
| E. coli | 0.73646 | 0.73404 | 0.73645 |
| S. cerevisiae | 0.00143 | 0.00267 (**1.87x**) | 0.00143 |
| E. faecalis | 0.00001 | 0.00002 (1.40x) | 0.00001 |
| S. enterica | 0.00024 | 0.00028 (1.20x) | 0.00024 |
| mapped reads | 342,020 | 365,407 | 341,960 |

**rammap reproduces minimap2 exactly at the same preset** — every organism at ratio
1.00, mapped-read counts within 0.02%. The embedded aligner is not a source of
discrepancy, which is worth having measured given reviewer item 4 asks about it.

**The preset does matter on real reads**: `sr` over-maps by 6.8% and inflates the rare
organisms, up to 1.87x on *S. cerevisiae*. That looked like it might explain the model's
over-prediction of low-abundance untargeted organisms.

**It does not, because the simulated side is insensitive to it.** Re-running the gate's
simulated arm under both presets:

| | Spearman rho | r(log10) | Lin's CCC |
|---|---|---|---|
| `sr` | 0.9353 | 0.9821 | 0.9691 |
| `map-hifi` | 0.9324 | 0.9811 | 0.9688 |

The reason is the read model: simulated reads are **error-free**, so they map
unambiguously under either preset, while `sr`'s inflation on real reads comes from its
short-read scoring interacting with CCS's ~1% error rate. The asymmetry is visible in
the totals — on real reads `sr` maps 6.8% *more* than `map-hifi`, on simulated reads it
maps 8.7% *fewer* — but in both directions the losses are near-proportional across
organisms, and composition is a ratio.

**So: a real confound, correctly identified, that turns out not to explain the
residual.** The scale offset stays where B0d left it.

The gate now passes `--minimap-preset map-hifi` regardless — matching the real side is
free, `sr` drops 8.7% of simulated reads, it is ~35% slower, and "we mapped 5.8 kb reads
with a short-read preset" is not a sentence worth defending to a reviewer.

## B0f. Capture shortens the real library 5x; simulating that is wrong (2026-11-25)

Measured from the raw reads, for the first time. Read length is **completely
systematic** across the BioProject:

| arm | runs | mean | median |
|---|---|---|---|
| non-enriched | 6 | 5,508–6,109 | 5,312–5,911 |
| enriched | 18 | 925–1,909 | 614–1,667 |

Every enriched library is ~5x shorter than every non-enriched one, across both
chemistries and all three panels. The gate simulates at 5,800 bp and compares against
enriched data that is ~1,200 bp, so the simulated output is 5x too long — and the model
has no mechanism that would shorten it, since retention does not depend on length.

That looked like a strong candidate for the residual: fragment length sets how many
probes overlap a fragment, and at 5.8 kb a single hit covers ~2.9 kb of start positions
while at 1.2 kb it covers ~600 bp, so probe density would matter far more at the real
length.

**Tested, and it is wrong.** Same panel and seed, only the length changed:

| `--fragment-length-mean` | comp rho | comp r | comp CCC | enrich rho | enrich centred |
|---|---|---|---|---|---|
| 5800 (current) | 0.9324 | 0.9811 | **0.9688** | +0.8179 | **0.227** |
| 1200 (measured enriched) | 0.9265 | 0.9582 | 0.9221 | +0.7036 | 0.320 |

Worse on every metric. (Absolute values here are below the gate's because this run used
the theoretical `--sample` rather than the non-enriched weights; only the 5800-vs-1200
contrast is meaningful.)

**Why, in hindsight:** the input library *is* ~6 kb — that is exactly what the
non-enriched arm measures, and it is what the probes see. The shortening happens during
or after capture (post-capture PCR strongly favours short templates, and there may be a
size-selection step). Feeding the *output* length back in as the *input* length
misrepresents the physical library, and the data says so.

**What this does leave:** a real, unmodelled, systematic effect — capture shortens the
library 5x and BaitBench preserves input length exactly. It belongs with the Stage 6
amplification work, as a length-dependent amplification bias, not as a fragment-length
parameter. It is not the composition residual's cause.

Fourth candidate eliminated by measurement, after placement stringency, the melting
window, and the mapper.

## B0g. Read realism does not move the metrics either (2026-11-26)

The `pacbio-hifi` badread profile exists now (`pacbio2021` + `--read-identity`), so the
question of whether error-free reads were costing us anything is directly testable
rather than argued from first principles.

**The profile works.** Identity measured from NM over aligned length, same method on
both sides:

| reads | mean | median | p5 |
|---|---|---|---|
| `perfect` | 100.00% | 100.00% | 100.00% |
| `badread pacbio-hifi` | 99.24% | **99.88%** | 96.40% |
| real CCS | 98.53% | **99.83%** | 91.53% |

The median is near-exact. The simulated left tail is lighter than reality — the profile
default is `99.5,100,1.5` against a measured sd of 2.82 — so `--read-identity
98.5,100,2.8` would match more closely. Good enough to answer the question.

**It changes nothing measurable:**

| | Spearman rho | r(log10) | Lin's CCC |
|---|---|---|---|
| `perfect` | 0.9324 | 0.9811 | 0.9688 |
| `badread pacbio-hifi` | 0.9324 | 0.9781 | 0.9675 |

Spearman is identical to four decimals; CCC moves by 0.0013 — in the *worse* direction,
within noise.

**So `--read-simulator perfect` stays, but now for a defensible reason.** The old
justification ("badread's pacbio profile models CLR, further from the data than
error-free reads") was wrong about badread's capability — it was a limitation of this
wrapper, which never passed `--identity` and so ran at 95% whatever the error model. The
correct justification is that the composition metrics cannot tell the difference: real
CCS is 99.83% median identity, and the ~0.2% of bases that differ do not change which
reference a 5.8 kb read maps to.

Fifth candidate eliminated by measurement, after placement stringency, the melting
window, the mapper, and fragment length. The pattern is consistent and worth stating
plainly: **everything BaitBench reports from reads is compositional, and compositional
metrics are insensitive to read-level detail.** Near-proportional effects cancel in a
ratio. Future hypotheses of this shape should be expected to fail, and the scale offset
should be looked for somewhere other than the read pipeline.

## B0d. Gates re-run after Stage 5 — all pass (2026-11-25)

Stage 5 changed retention mass on every reference, so the gates were re-run. **All
pass.** Composition improved on all three replicates; LOD and enrichment held.

| | pre-Stage-5 | post-Stage-5 |
|---|---|---|
| MOXTAA composition r / CCC | 0.974 / 0.961 | **0.982 / 0.969** |
| MOXTAB composition r / CCC | 0.965 / 0.958 | **0.988 / 0.982** |
| MOXTAC composition r / CCC | 0.975 / 0.966 | **0.989 / 0.982** |
| LOD Pearson r (log10) | +0.9993 | +0.9993 |
| LOD sim/real, dilute | 1.92x | 2.02x |

Enrichment now scores **15** organisms rather than 14: *Candida albicans* no longer
needs excluding, because the mass fix put it at ratio 1.000 of its intended input share
(from 0.244) and it predicts 0.156 against a measured 0.130. Per-replicate enrichment
rho +0.929 / +0.889 / +0.771.

Per-locus efficiency remains stable across the whole dilution series — spread
(max/min over 11 conditions) 1.0–2.4x — which is the separability of efficiency from
abundance that the quota model could not produce.

**The global scale offset is unchanged and remains the one open item:** 2.02x in the
LOD dilute regime, ~1.3x as the enrichment bias. Stage 5 did not touch it, which is
expected — it is a property of `f / p_ns`, not of fragment geometry.

**Harness footprint:** a full gate run now costs **3.5 GB**, down from ~50 GB, after
`combined_reference.fa` and `mapping_reference.fa` were added to `--cleanup`. They are
4.2 GB per LOD condition on this dataset and had been excluded from the cleanup list,
which is what filled the disk mid-run on 2026-11-24.

## B0c. Full gates, all 11 LOD conditions and 3 SimVal replicates (2026-11-24)

First time the full gates have ever actually run: `run_full.sh` exported
`VALIDATION_MODE=full` and then exec'd `run_fast.sh`, which began with an unconditional
`export VALIDATION_MODE=fast`. Every previous "full" run was the fast subsets. Three
harness defects fixed before these numbers mean anything (see the plan's harness list).

### LOD titration — 11 conditions

| | |
|---|---|
| real dynamic range | **27,035x** (4.16e-5 at NTC to 0.796 at Pos1) |
| Pearson r (log10) | **+0.9993** |
| log10 RMSE | 0.243 (1.75x) |
| sim/real, 8 dilute conditions | **1.92x** (range 1.82–2.06) |
| sim/real, 3 saturated conditions | 1.09x (range 1.05–1.12) |

The shape is essentially exact over four and a half orders of magnitude. Nearly all the
RMSE is a **constant multiplicative offset** in the dilute regime that compresses as
capture saturates — a scale error, not a shape error.

### Per-locus efficiency at Pos1, relative to rplB

r(log10) = **+0.964**, Spearman rho = +0.786, log10 RMSE 0.139 (**1.38x**), worst locus
clpX 1.70x. Per-locus efficiency is stable across the whole dilution series — spread
(max/min over 11 conditions) is 1.0–2.1x for seven of eight loci — which is the
separability of efficiency from abundance that the real data shows and the quota model
could not produce.

*Not comparable to the previously recorded "per-locus log-RMSE 1.00"*: that came from
`locus_efficiency.py`'s own statistic, and these gates also ran with the shipped
10 bp window, which the old `config.sh` was overriding. Re-derive that figure the same
way before claiming an improvement.

### SimVal, 3 replicates

Composition (E. coli strains collapsed), and enrichment with `--min-real-reads 10`:

| replicate | comp. rho | comp. r | comp. CCC | enrich rho | enrich r | enrich centred |
|---|---|---|---|---|---|---|
| MOXTAA | 0.929 | 0.974 | 0.961 | +0.939 | +0.917 | 0.204 |
| MOXTAB | 0.929 | 0.965 | 0.958 | +0.917 | +0.916 | 0.186 |
| MOXTAC | 0.949 | 0.975 | 0.966 | +0.780 | +0.863 | 0.203 |

**MOXTAC's composition figures were an artifact, and are fixed.** They first read
r = 0.753 / CCC = 0.751 against ~0.97 / ~0.96 for the other two. Cause: `compare.py`
floored a zero proportion at `eps = 1e-9`. Real proportions here span ~1e-5 to 0.8, so a
floored zero landed four orders of magnitude below the smallest real observation and
acted as an enormous leverage point in both log-scale statistics. MOXTAC is the only
replicate with zero-read organisms (*E. faecalis* and *C. perfringens*), so only it was
affected. Zeros are now floored at that arm's **detection limit** — the proportion one
read represents — which is what a zero is actually bounded by, and the censored count is
reported. MOXTAC becomes r = 0.975 / CCC = 0.966, in line with the others; MOXTAA and
MOXTAB are unchanged, having nothing censored. Spearman was never affected: it runs on
raw proportions.

**The enrichment gap is real, though.** The first suspicion there was
rare-organism noise, and that was *part* of it: *E. faecalis* has 4 / 1 / **0** reads
across A / B / C and *C. perfringens* 0 / 0 / 1, so the scored set differed between
replicates and single-read organisms carried ~100% measurement error. Hence the new
`--min-real-reads` floor (default 10), which puts all three on the same 14 organisms.
It improves MOXTAC (rho 0.768 -> 0.780, r 0.822 -> 0.863, centred 0.239 -> 0.203) but
does not close the gap. MOXTAC's measured dynamic range is 34x against 63x and 43x, and
*C. difficile* sits at 0.323 against 0.930 and 0.985 — the experiment itself enriched
less that time. So the two MOXTAC anomalies have different causes: the composition one
was ours, the enrichment one is the data's.

### The number that matters: model error against the replicate noise floor

Three replicates make it possible to ask whether the residual is model error or
experimental noise. Per organism, across replicates:

| quantity | log10 | fold |
|---|---|---|
| replicate noise, pooled within-organism SD | 0.125 | 1.33x |
| SD of the 3-replicate **mean** (what the prediction is scored against) | 0.072 | 1.18x |
| prediction vs replicate mean, RMSE | 0.208 | 1.61x |
| — of which bias | +0.113 | 1.30x |
| — centred (shape) | 0.175 | 1.50x |
| **excess model error after removing replicate noise** | **0.159** | **1.44x** |

**The model is not yet at the noise floor**: its shape error is about twice the
irreducible noise. An earlier, looser comparison — median replicate spread (max/min,
1.59x) against median |prediction / replicate mean| (1.41x) — appeared to show the
model matching replicate variability, but those are two different statistics and the
comparison flattered the model. The variance decomposition above is the honest version.

### Both datasets say the same thing

| | shape | scale offset |
|---|---|---|
| LOD titration, 11 conditions | r(log10) = +0.9993 | 1.92x (dilute regime) |
| SimVal enrichment, 3 replicates | rho = +0.78 to +0.94 | 1.30x (bias term) |

Two independent datasets, different organisms, panels, chemistries, fragment lengths and
temperatures, agree that **relative capture is predicted well and there is one global
scale offset of roughly 1.3–2x.** That is a single scalar, and section A2 already
identified which one: only `f / p_ns` is identifiable from composition, and section D
showed `--capture-fraction` moves bias without touching shape (Spearman identical to
four decimals across a 10x range). So the offset is correctable in principle and
**cannot be fitted from this data** — it needs an absolute anchor: a spike-in, a
measured library yield, or an on-target fraction against a distractor arm.

That is now the single clearest open item for the model, and it is a parameterisation
question rather than a physics one.

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
- **Build it on the gate's analytic scorers, not a separate harness.**
  `tests/locus_calibration.rs` was deleted (2026-11-25). It existed because "running
  the whole pipeline to test a change to the capture law takes hours" — no longer true:
  the full gate, all 11 LOD conditions and 3 SimVal replicates, now runs in ~32 minutes
  since read mapping was parallelised. Its header also still described the per-locus
  defect as open, which it is not. And it had already diverged from production once:
  its first version averaged over all starts spanning a locus while production averages
  over covered positions only, producing a confident 29x improvement the pipeline did
  not show. A second implementation of the sweep is a second thing to keep in step.
  `simulate` with a warm `--probe-hits` cache plus `simval/enrichment.py` gives the same
  seconds-per-evaluation loop while running the production code path.

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
- [x] **Re-derive the harness `--capture-fraction` values** (2026-11-24) — and the
      answer is that **the enrichment data cannot pin it, by construction.** Swept
      against B0b with `[probe]` and the window fixed:

      | `--capture-fraction` | Spearman rho | r(log10) | log10 RMSE | centred | bias |
      |---|---|---|---|---|---|
      | 0.1 | +0.9361 | 0.921 | 0.445 | 0.301 | +0.328 |
      | 0.3 | +0.9361 | 0.925 | 0.308 | 0.219 | +0.217 |
      | 0.5 | +0.9361 | 0.920 | 0.256 | 0.198 | +0.162 |
      | 0.7 | +0.9361 | 0.914 | 0.233 | 0.195 | +0.128 |
      | 0.9 | +0.9361 | 0.909 | 0.224 | 0.198 | +0.104 |
      | 0.99 | +0.9361 | 0.907 | 0.222 | 0.201 | +0.095 |

      Spearman rho is *identical to four decimals* at every value across a 10x range,
      and centred RMSE varies by 0.006 over 0.5–0.99. All of the movement in raw RMSE
      is the bias term. This is a direct measurement of what the flag's help text
      asserts by design: it is a uniform multiplier on specific capture and does not
      reshape relative capture. Enrichment is scale-free by compositional closure, so a
      pure scale knob is invisible to it — exactly A2's identifiability problem,
      now confirmed on the other axis.

      **Consequence: stop trying to fit it.** Pinning `--capture-fraction` needs an
      absolute quantity the compositional data does not contain — a spike-in, a measured
      library yield, or an on-target fraction against a distractor arm. The harness
      should therefore pass a *documented scenario value*, not a fitted one, which is
      also what the flag is for. `simval/run.sh`'s 0.207 was calibrated under the old
      meaning and is simply stale; replace it with the default rather than refit.
- [ ] **Full 11-condition LOD and all three SimVal replicates**, not the fast subsets.
      Everything measured so far is 3 LOD conditions and 1 replicate.
- [x] **Decide whether `--capture-fraction 0.5` remains a sensible default** — keep it.
      It is mid-range for a knob that now means protocol efficiency, the sweep above
      shows the choice does not affect any shape metric, and no available dataset
      constrains it. Changing it would imply a precision the data does not support.

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
