# BaitBench — Paper Revision Roadmap

> Governing document for the manuscript revision. Unlike the other files in this
> directory, which are design sketches, this one tracks committed work. Phases run in
> order where dependencies demand it; Phases 3 and 4 can be picked up at any time.
> Nothing here has been started yet.

## Context

The BaitBench Application Note (submitted to *Bioinformatics*, `MTEC-25-04-MT-025`)
was rejected. The AI-disclosure failure was clerical — the PI checked "no AI used"
on the wrong submission. The reviewer's other five comments are substantive, and two
are worse than the letter suggests:

- **Item 6 is a real modelling defect**, not a clamp. The realized distractor read
  fraction saturates at `1 − capture_fraction` (0.5 at the default). Root cause is
  already recorded in [capture-probe-competition.md](capture-probe-competition.md).
- **Items 2 and 3 are not a mis-citation.** The validation ran against an in-house
  ONT run of the ZymoBIOMICS **Gut** Standard (D6331, 21 strains), not TELSeq's
  public PacBio data on the CSII standard (D6310, 10 organisms) — and the input
  weights were derived from the very reads being predicted. Reviewer comment 3 is
  correct.

**Items 6 and 2–3 are the same problem.** The capture model is abundance-blind, so
the simulation largely reproduces whatever weights it is given; ρ = 0.916 measured
weight propagation, not prediction. Fixing the model is what makes the validation a
real test — and is also a prerequisite for it, since CSII spans ~6 orders of
magnitude of input abundance, exactly the regime the current model cannot represent.

| # | Reviewer comment | Phase |
|---|---|---|
| 1 | AI use not disclosed | 6 |
| 2 | Mock community discrepancies (20 vs 10 organisms, *M. smithii*) | 5 |
| 3 | Weights may have come from the output (circularity) | 5 |
| 4 | minimap2 cited, rammap used | 3 |
| 5 | Reimplementations unvalidated | 4 |
| 6 | `--distractor-fraction` > 0.5 silently ignored | 1, 2 |

---

## Decisions

| ID | Decision | Blocks | Status |
|---|---|---|---|
| **D1** | How far to take the capture-model fix (see analysis below) | Phase 2 | **Resolved**: full per-fragment retention model, quota removed outright (no compatibility flag). Shipped 2026-11-19; see the simulation-rework plan. |
| **D2** | TELSeq probe panel: request from authors, or reconstruct from MEGARes v1.0 + BacMet + VFDB and declare as a limitation | Phase 5 | **Open — largest risk** |
| **D3** | Target venue: appeal to *Bioinformatics* vs fresh submission. Sets the AI-disclosure heading and wording | Phase 6 | **Open** |
| **D4** | Whether to keep the in-house D6331 ONT run as a deposited secondary dataset (legitimate once circularity is removed; needs an SRA deposit) | Phase 5 | **Open** |
| **D5** | Which run is canonical for any interim citation of the old result | — | **Open** (see Phase 5 note) |
| **D6** | Independent input weights for the validation | Phase 5 | **REOPENED 2026-11-23**: a non-enriched arm exists. SRA PRJNA1081843 lists 6 runs titled "Non-enrichment of synthetic mock community: Negative control" (SRR28150698/699/700/732/733/735), matching the local `*NEG*` count files. Composition tracks the manufacturer's theoretical table at r(log10) = +0.989 vs +0.931 for the enriched arm. Reviewer item 3 can therefore be answered by predicting measured-enriched from measured-unenriched, rather than from theory. See `parameterisation.md` A3. |
| **D7** | Whether the Borrelia LOD titration becomes a second manuscript figure (see "Folding in" below). Needs an SRA deposit, like D4 | Phase 5 | **Open — recommend yes** |

---

## D1 background: what is actually wrong with the capture model

`--distractor-fraction` is never clamped and never compared against a threshold. The
ceiling is structural, and three problems compound.

**(a) `capture_fraction` fixes output composition, not efficiency.**
`src/commands/simulate.rs:88-97` computes
`n_capture = round(num_fragments × capture_fraction)` and fills each of the two
quotas to size with replacement. Neither quota consults `weights.txt`. The split
between the probe-mediated arm and the background arm is a user-supplied constant,
independent of how much target DNA exists.

**(b) The capture arm is abundance-blind.**
`src/sampling/thermo_sim.rs:296` is
`WeightedIndex::new(vec![1.0f64; probe_names.len()])` — level-1 probe selection is
uniform. Abundance enters only at level 2, *within* one probe's hit list
(`thermo_sim.rs:315`). A probe whose hits are all on targets yields a target
fragment no matter how much distractor DNA is present, so the distractor share of
the capture arm `q_c` is bounded by the fraction of probes with a cross-reactive
distractor hit — a property of the panel, not of the flag. For a clean panel
`q_c = 0` exactly.

**(c) The background arm is mass-weighted, not copy-weighted.**
`src/sampling/thermo_sim.rs:405` samples ∝ `weight × seq.len()`, so
`q_b = df·L_d / (df·L_d + (1−df)·L_t)`, which equals `df` only when target and
distractor mean lengths match. This contradicts
`docs/explanation/simulation-modes.md:109` and `CLAUDE.md`, which claim the weight
formula "ensures" a `df : (1−df)` *fragment* ratio. It ensures a *weight* ratio only.

Net: `D_realized = cf·q_c + (1−cf)·q_b ≤ 1 − cf`. At the default `cf = 0.5`:

| `--distractor-fraction` | realized distractor reads | realized target reads |
|---|---|---|
| 0.5 | 0.250 | 0.750 |
| 0.9 (default) | 0.450 | 0.550 |
| 0.99 | 0.495 | 0.505 |
| 0.999 | 0.4995 | 0.5005 |

Above `df = 0.5` everything compresses into a 0.05-wide band, indistinguishable from
sampling noise. The target fraction has a hard floor at `capture_fraction`. The
reviewer's "0.5" is numerically `1 − capture_fraction`.

**Why it matters beyond one flag.** The advertised use case is limit-of-detection
work — rare pathogen in a large background — which is the regime the quota model
cannot represent. `docs/how-to/set-parameters.md:34` advertises
`--distractor-fraction 0.999  # 0.1% target — near limit of detection`; the pipeline
delivers ~50% target. Manuscript Figure 1 sweeps `capture_fraction` rather than
abundance for exactly this reason.

### D1 options

**Option 1 — Physical retention model** *(recommended)*

Replace the two quotas with per-fragment retention: draw a pre-capture library from
all sequences ∝ `weight × length`; for each fragment take its best overlapping probe
hit, compute ΔG via the existing `thermodynamics::delta_g`, convert to a Langmuir
occupancy `θ = K/(1+K)` with `K = [probe]·exp(−ΔG/RT)`; retain with probability `θ`,
or with `--background-retention` (wash carryover, default ~1e-3) if no probe
overlaps; draw until `--num-fragments` are retained. Keep `--capture-fraction` as
the calibration knob setting `[probe]`, and keep the old behaviour as
`--capture-model quota`.

- *Pros.* Output composition emerges instead of being imposed, so `df` is monotonic
  across its whole range and limit-of-detection simulation becomes meaningful.
  On-target fraction becomes a **prediction**, which answers comments 3 and 6
  together. Langmuir occupancy is citable physics and removes the need for the
  `boltzmann_score` clamp at `src/thermodynamics.rs:277-279`, which currently makes
  it impossible for any probe to deplete. Resolves the manuscript's own stated
  limitation rather than defending it.
- *Cons.* Largest blast radius — every simulation output changes, so Figures 1 and 2
  and all `coverage-curve` results must be regenerated. Two new parameters for users
  to understand. Yield becomes variable, shifting `--num-fragments` semantics
  (mitigated by retry-to-quota). Most new testing.
- *Note.* This is the Bernoulli-filter idea from
  [capture-probabilistic-model.md](capture-probabilistic-model.md) (Option C), but
  deriving the retention probability from ΔG rather than from a sigmoid. That removes
  Option C's main stated drawback — its sigmoid had to be re-calibrated by binary
  search on every run, because the alignment-score distribution changes with the
  probes, targets and distractors. A Langmuir occupancy needs no per-run
  calibration. **Option 1 supersedes Option C**, which can be retired.

**Option 2 — Abundance-weighted probe selection + guards**

The "D1" fix from [capture-probe-competition.md](capture-probe-competition.md):
weight level-1 by `Σ_hits(boltzmann × seq_weight)`, fix (c) in
`src/sampling/weights.rs:57-67`, fix the background retry/warning, tighten the CLI
range, hard-error when `distractor_fraction > 1 − capture_fraction`.

- *Pros.* Small, surgical, reviewable. Fixes the unambiguous bugs. Keeps the engine
  and published numbers broadly intact. Ships quickly.
- *Cons.* Does **not** remove the ceiling — the capture arm is still a fixed quota,
  so the floor at `capture_fraction` survives and the honest response becomes "we
  now error instead of silently ignoring," inviting "then the model cannot do
  limit-of-detection." Also re-introduces the problem the manuscript already names:
  abundance-weighted probe selection over-selects common species, and with nothing
  counteracting it, rare organisms drop out. On log-distributed CSII that likely
  makes the new validation *worse*.
- *Decisive against.* [capture-probe-competition.md](capture-probe-competition.md)
  (line 13) already records that **D1 cannot fix the limit-of-detection case**: a
  single target diluted against a distractor has no competing probes, so "there's no
  alternative probe for D1's reweighting to shift mass toward — BB's probes remain
  the only eligible pool at every dilution step, and the quota-filling resampling
  persists regardless of weighting." The doc concludes a capacity mechanism or the
  Bernoulli filter is needed. That is Option 1.

**Option 3 — `--enrichment-factor` reparameterization**

Keep the engine; derive `capture_fraction` from `(1−df, E)` via
`on_target_out = f_in·E / (f_in·E + (1−f_in))`.

- *Pros.* Cheapest physically-interpretable option; ~one function; makes `df`
  meaningful end-to-end at first order; gives users a parameter measurable from a
  pilot run.
- *Cons.* Cosmetic — `q_c` is still panel-determined and within-arm composition is
  still abundance-blind. A reviewer reading the code sees the same structure. Does
  nothing for rare organisms. Best used *alongside* Option 1 as a convenience alias.

**Option 4 — Document the ceiling only**

- *Pros.* Hours, not weeks. Nothing can regress.
- *Cons.* Concedes the headline scenario is unsupported. The reviewer has already
  read the code; "we documented it" is unlikely to move a rejection.

**Recommendation.** Option 2's bug fixes **and** Option 1, with Option 3 as a cheap
alias. Option 2 is a strict prerequisite subset — fixes (b), (c), the shortfall and
the CLI range are needed under either model — so land those first as independently
reviewable commits (Phase 1), then build the retention model on top (Phase 2).
Default to `retention`, keep `quota` for v1.0 reproducibility.

---

## Folding in the existing improvement docs

Reviewed all seven. Verdicts:

| Doc | Verdict |
|---|---|
| [capture-probe-competition.md](capture-probe-competition.md) (Option D) | **Already the core of Phases 1–2.** Its LOD finding is what rules out Option 2 as a standalone fix. Retire on completion. |
| [capture-probabilistic-model.md](capture-probabilistic-model.md) (Option C) | **Superseded by Option 1** — same Bernoulli filter, ΔG instead of a calibrated sigmoid. Retire. |
| [capture-quality-weighted-resampling.md](capture-quality-weighted-resampling.md) (Option A) | **Retire without implementing.** Its own recommendation says the realism "likely doesn't justify the complexity" for a tool where the user *specifies* enrichment rather than predicts it — and moving from specifying to predicting is the whole point of this revision. Leaving three competing capture designs in the repo is itself a bad look for a code-reading reviewer. |
| [mismatch-stacking-thermo.md](mismatch-stacking-thermo.md) | **Fold into Phase 4 — highest-value add.** See below. |
| [rna-dna-thermodynamics.md](rna-dna-thermodynamics.md) | **Defer the implementation; state the limitation now.** See below. |
| [mutation-distribution-validation.md](mutation-distribution-validation.md) | **Keep deferred.** Unscoped (no panel, method, or divergence dataset chosen) and not a reviewer comment. |
| [INDEX.md](INDEX.md) | Three entries point at files that do not exist. Phase 0. |

### Fold in: mismatch stacking (high value, shared cost)

`src/thermodynamics.rs` currently uses "SkipStacking" — any nearest-neighbour step
with a mismatch at either end contributes **zero** stacking energy, so a probe with
a few scattered mismatches scores like a completely non-binding one. Four reasons to
do this now rather than later:

1. **The data-entry cost is already being paid.** Phase 4 requires building a
   reference ΔG table from SantaLucia & Hicks (2004) cross-checked against
   biopython's `MeltingTemp.DNA_IMM1` — the *same* source and the *same* dictionary
   the mismatch table comes from. Doing them separately means opening Table 4 twice.
2. It converts a reviewer-visible deviation from the cited model ("we skip stacking")
   into an implementation of it ("SantaLucia & Hicks 2004, Table 4").
3. `thoughts.md` already flags SkipStacking as a pre-submission **[MAJOR]** concern.
4. Mismatch tolerance is what the panel-breadth claims rest on, and Phase 2's
   retention model makes ΔG matter *more* — it now sets a retention probability per
   fragment rather than a relative weight within one probe's hit list.

Scope per the doc: ~15 lines of logic in `delta_g()`, a ~256-entry `MM_TABLE` as
`Option<NnParams>` (~48 independent values), 2–3 existing tests revised, ~5 new.
Keep the `has_stacking` guard — mismatch stacking alone does not nucleate a duplex.

### Fold in: run the existing RAmpSim thermodynamics comparison (nearly free)

`~/Projects/CrossLab/RAmpSim/compare/thermo_compare.py` (25 KB) **already exists and
already does this** — it reimplements both RAmpSim's and BaitBench's scoring in pure
Python, parses real probe alignments, and reports how the scores compare. It has
never appeared in the paper. Running it and reporting the result costs almost
nothing and defuses a cluster of exposure the reviewer has not yet joined up:
BaitBench reimplements RAmpSim's thermodynamic approach, depends on `rammap-core`
from a related lab, was validated on the same TELSeq dataset RAmpSim used, and cites
neither. `thoughts.md` lists this as **[CRITICAL]** + **[MAJOR]**. Add to Phase 4;
port the script into `validation/` alongside the others.

### State now, implement later: RNA chemistry

There is no `DuplexChemistry` support in `src/thermodynamics.rs` — everything uses
SantaLucia (1998) **DNA-DNA** parameters. Manuscript Figure 1 is a **Zika virus**
example, i.e. an RNA target scored with DNA-DNA parameters. A reviewer already
reading the code may well notice. Adding Sugimoto (1995) RNA:DNA and Xia (1998)
RNA:RNA tables plus the Nakano (1999) salt correction is real work and is not among
the six comments, so: **state it as an explicit limitation in the revised
manuscript** (cheap, pre-empts the hit) and keep the implementation deferred. Do not
add a one-variant enum speculatively.

### Promote: the Borrelia LOD titration is the best validation asset on disk

[capture-probe-competition.md](capture-probe-competition.md) (lines 9–15) points at
a second dataset that is more valuable than the composition correlation currently in
the paper:

- Real data: `~/Projects/CrossLab/MTEC_probes/borrelia_LOD/4Path/` — an 11-condition
  *B. burgdorferi* LOD titration into real tick matrix, with quantified BB-mapped
  read counts spanning **~20,000-fold** (NTC = 37 → Pos1 = 748,104), plus
  `target_analysis/ontarget_summary.tsv` reporting measured **target enrichment of
  37–51×** per sample.
- Simulation side: `~/Projects/CrossLab/MTEC_runs/BB_LOD_deep_dive/` with a thorough
  README. Current result: a flat ~78.5% plateau across NTC–s1E5 instead of the real
  drop-off.

Why this matters more than its current status suggests:

1. **It is the acceptance test for Phase 2.** Pass/fail, unambiguous: does the
   simulated titration curve reproduce the real one? Make it a required gate.
2. **It cannot be circular.** The input is known spike copy number; the output is
   observed reads. No weights are derived from the thing being predicted — which is
   precisely the criticism in reviewer comment 3.
3. **The measured 37–51× enrichment is an empirical anchor** for the quantity
   Option 1's retention model predicts. Phase 2 can be calibrated against it.
4. **It is the stronger second figure.** Reproducing a real LOD curve across four
   orders of magnitude answers "insufficient to convince this reader that it works as
   described" far better than a composition correlation does. Hence **D7**.
5. **Bonus:** the README documents constant *Rickettsia* background tracking with
   tick-matrix presence rather than with BB dilution — probe cross-reactivity with
   the host genome. That is a ready-made real-world test of `xreact` and of the
   FP_distractor machinery, which nothing currently validates.

Cost: an SRA deposit (same call as D4) and porting `BB_LOD_deep_dive/scripts/` into
`validation/borrelia_lod/`.

---

## Phase 0 — Housekeeping (no dependencies, S)

- [x] Write this roadmap to `improvements/paper-revision-roadmap.md` and add the
      `improvements/INDEX.md` line.
- [ ] Fix the three `INDEX.md` entries pointing at files that do not exist
      (`sweep-read-length.md`, `paired-end-support.md`, `long-read-testing.md`).
- [ ] Mark `capture-quality-weighted-resampling.md` (Option A) and
      `capture-probabilistic-model.md` (Option C) as superseded, with a pointer to
      this roadmap, so the repo no longer presents three competing capture designs.
- [ ] Add `.github/workflows/test.yml` running `cargo test`,
      `cargo clippy -- -D warnings`, `cargo fmt --check`. Neither existing workflow
      runs tests; `RELEASE_CHECKLIST.md` never mentions `cargo test`.
- [ ] Create a `tests/` integration directory. There is none today — all 70 tests are
      inline `#[cfg(test)]` across 7 files, and nothing tests the simulation core.
- [ ] Add `cargo test` to `RELEASE_CHECKLIST.md` §1.

*Done when:* CI runs tests on push and `tests/` holds at least one end-to-end run.

---

## Phase 1 — Capture-model bug fixes — **DONE** (except one item, reopened below)

- [x] Tightened the `df` bound check to `[0, 1)` in `src/main.rs`, with a message
      naming the flag and suggesting 0.9999. `1.0` previously passed the early check
      and failed inside `prepare` after the run directory and combined reference had
      been written.
- [x] Short-contig shortfall. Resolved differently than planned: there is no longer a
      background arm with a retry loop. `build_profiles` now builds each reference's
      profile at `min(frag_len, seq_len)`, so a contig shorter than the nominal
      fragment length stays in the library instead of getting zero mass and
      disappearing. Covered by `short_references_are_not_dropped`.
- [x] Realized composition is now reported. `simulate` computes yield, on-target
      fraction and enrichment analytically and writes `fragments.model_stats.tsv`;
      `metrics` carries them into `results.tsv` as `model_yield_fraction`,
      `model_on_target_fraction`, `model_enrichment`, and into the JSON as a
      `capture_model` block. Per-reference efficiency goes to
      `fragments.retention.tsv`. `capture_rate` is kept but is 1.0 by construction and
      documented as such — under per-fragment retention `fragments.fa` *is* the
      retained library, so counting it cannot recover what was filtered out.
- [x] The regression test exists as `tests/simulation.rs::distractor_fraction_is_honoured`,
      sweeping `df ∈ {0.1, 0.5, 0.9, 0.99, 0.999, 0.9999}` at a fixed seed. Its
      assertion changed: see the note below.
- [ ] **Reopened — `src/sampling/weights.rs:57-67`, length-normalizing the distractor
      weight.** The original goal, "so the realized *fragment* ratio matches the
      documented `df : (1−df)`", is no longer the right target: under retention the
      realized output ratio is *supposed* to differ from the input ratio, because that
      difference is the enrichment. What remains worth deciding is whether
      `--distractor-fraction` denotes a **molar** fraction (equal weight per contig, as
      now, so long contigs contribute proportionally more DNA) or a **mass** fraction
      (weight ∝ 1/length). This changes what the flag means, so it needs a decision
      rather than a patch. Not blocking: it affects only how distractor material is
      apportioned *among* distractors.

**Note on the regression test's assertion.** It originally asserted realized ≈
requested distractor fraction. That was wrong: `--distractor-fraction` sets the *input*
library, and capture enriches targets, so realized ≠ requested at moderate values is
correct physics. The real defect in reviewer item 6 was the **ceiling** at
`1 − capture_fraction` — the parameter went inert over the top half of its range
whatever the input. The test now asserts monotonicity, the absence of that ceiling, and
convergence at an input so extreme that finite enrichment cannot rescue the target.

*Done:* measured titration 0.00000 / 0.00125 / 0.01700 / 0.16700 / 0.66775 / 0.95500
across requested 0.1 … 0.9999, against the quota model's 0.043 / 0.242 / 0.450 / 0.495
/ 0.495 / 0.495.

---

## Phase 2 — Capture-model rework — **DONE** (gates met 2026-11-21)

Superseded in detail by the simulation-rework plan, which expanded this phase into
Stages 0–3 (harness, ΔG correctness, divergent-probe alignment, internal-mismatch
thermodynamics, then the capture model). What shipped differs from the sketch below in
two ways worth recording:

- **No `--capture-model {retention,quota}` flag.** Per D1 the quota model was removed
  outright, so there is no v1.0 bit-for-bit reproduction target. `boltzmann_score` was
  deleted with it — the quota sampler was its only consumer, and its clamp-to-1.0
  ("probes can only enrich") is wrong under an occupancy model, which is bounded in
  [0,1] by construction.
- **`--capture-fraction` was redefined, not retired.** It is now capture *efficiency*
  (how much of the material it has probes for the panel pulls down), back-solved to a
  probe concentration by bisection. An intermediate version targeted the realized
  *output* on-target fraction and reproduced the quota ceiling almost exactly
  (0.494 vs 0.495), because adding distractor made the solver compensate. Any parameter
  targeting an output composition does this.

- [x] `src/thermodynamics.rs` — `langmuir_occupancy`, `ln_k_unit_conc`,
      `occupancy_from_ln_k`; `boltzmann_score` removed.
- [x] `src/sampling/thermo_sim.rs` — both quota samplers deleted; `sweep_retention`,
      `build_retention_profile`, `build_profiles`, `capture_efficiency_at`,
      `solve_probe_concentration`, `retention_stats`, `sample_retention_fragments`,
      `write_retention_report`.
- [x] `src/cli.rs` — `--background-retention`; hidden `--probe-concentration` and
      `--min-overlap-frac`.
- [x] `src/commands/simulate.rs` — quota split deleted, concentration solver added.
- [x] `src/commands/metrics.rs` — `model_*` columns (yield, on-target, enrichment).
- [x] Diagrams 8, 9 and 12 regenerated for the retention model. **Note:** the
      `poster_*` variants have no generator in this repo and still show the two-arm
      quota model — regenerate wherever they are authored before reusing the poster.
- [ ] Regenerate manuscript Figure 1 (`coverage_curve.png`).
- [ ] Port `~/Projects/CrossLab/MTEC_runs/BB_LOD_deep_dive/scripts/` into
      `tests/validation/lod/` and use the titration as the gate below. Calibrate
      against the measured 37–51× target enrichment in
      `borrelia_LOD/4Path/target_analysis/ontarget_summary.tsv`.

*Done when:*

1. ✅ `df` is monotonic across its whole range and the `1 − capture_fraction` ceiling is
   gone (`distractor_fraction_is_honoured`, un-ignored and green).
2. ✅ Per-locus efficiency is emergent and separable — on `examples/minimal` the two
   targets come out at 0.896 and 0.097 from probe density and ΔG alone, with
   distractors exactly at the carryover floor. `fragments.retention.tsv` exports this
   per reference, so the Borrelia per-locus ratios can be compared directly without
   running reads through.
3. ✅ **the Borrelia LOD titration resolves from a flat ~78.5% plateau into a real
   titration curve.** Measured 2026-11-21 (fast subset, 3 conditions):

   | condition | real BB fraction | simulated |
   |---|---|---|
   | NTC | 0.0000416 | 0.00010 |
   | s1E3 | 0.00109 | 0.00205 |
   | Pos1 | 0.79565 | 0.83475 |

   Monotone across three orders of magnitude where the old model gave 78.50% for every
   condition. Pos1 within 1.05x.

**Phase 2 gates are met.** Two things from the sketch above remain unticked and are
tracked rather than blocking: manuscript Figure 1 has not been regenerated, and the
per-locus fit — beyond the stated gate — reached log-RMSE 1.00 only after a capture-law
change whose parameterisation is not yet interpretable across temperature and salt. See
`parameterisation.md`; it is the next substantial piece of simulation work, and Phase 5
should not be written up before it settles.

---

## Phase 3 — Reviewer item 4: cite and name rammap (parallel, M)

`Cargo.toml:19` depends on `rammap-core = "1.1.1"` (crates.io; repo
`github.com/jwanglab/rammap`, MIT). Everything funnels through
`src/external/minimap2.rs`, whose `check_available()` is a no-op stub — there is
**no** external minimap2 anywhere in `src/`.

- [ ] `paper/BaitBenchRef.json` — add a `rammap` entry; keep the minimap2 citation
      only where the manuscript refers to the algorithm rammap implements. While in
      the file: delete the junk Zotero duplicate `MethodsApplicationsComputational`
      and add Gnirke et al. 2009 (flagged in `thoughts.md`, never added).
- [ ] `BB_Paper/paper.md` — `### Assessing Probes` and `### Simulating Capture` both
      say minimap2. Switch to rammap, cite minimap2 as its basis, and add the
      version, preset, and flags the manuscript currently omits entirely.
- [ ] Code — rename the module to `src/external/aligner.rs`; fix log strings at
      `assess_probes.rs:74` and `xreact.rs:163`; fix the `run_params.tsv` value at
      `assess_probes.rs:270` / `xreact.rs:431` (this is what the HTML reports show);
      fix `R/assess_probes.Rmd:363`. Keep `--aligner minimap2` and `--minimap-preset`
      as flag *spellings* for compatibility, with help text saying rammap.
- [ ] Docs — `docs/reference/{commands,parameters,input-formats}.md`,
      `docs/how-to/{design-and-assess-probes,analyze-cross-reactivity,troubleshoot}.md`,
      `docs/diagrams/prepare_diagrams.R` (diagram boxes are rendered into the
      published site), `gui/src/lib/views/RunView.svelte`, `README.md:285`,
      `ARCHITECTURE.md:63,80,229,258-260,634`, `CLAUDE.md:107,112`.
- [ ] Bug: `docs/tutorials/installation.md:48` links rammap to `github.com/lh3/rammap`
      — wrong org, correct to `jwanglab/rammap`.
- [ ] Dead weight: `gui/src-tauri/src/commands/setup.rs:531-601` still downloads
      `minimap2.exe` from `lh3/minimap2` releases on Windows; nothing invokes it.
- [ ] `ARCHITECTURE.md` is stale on two other counts — it lists
      `catch.rs`/`syotti.rs`/`probetools.rs` at top level (they are in `src/probes/`)
      and lists minimap2 in the External Tool Wrappers table.

*Done when:* `grep -ri minimap2` returns only flag spellings, algorithm citations,
and the rammap-implements-minimap2 note.

---

## Phase 4 — Reviewer item 5: validate the reimplementations (parallel, M–L)

Scope: the comparisons needing no new external tooling, plus thermodynamics.
Upstream Syotti and ProbeTools head-to-heads deferred (no wrapper; needs C++/SDSL
and VSEARCH).

- [ ] **catch-lite vs CATCH** — `validation/catch_comparison/`. Feasible today with
      no new code: `catch` is in `environment.yml` and `external_catch::design_probes`
      has an identical signature to the native one, so `--method catch` vs
      `--method catch-lite` is a matched comparison. Report panel size, per-target
      breadth of coverage at matched mismatch tolerance, probe-sequence Jaccard, and
      runtime. **First verify** `src/external/catch.rs:54` — it passes `-ps` where
      the doc comment claims `-s`. Report the hardcoded internals so the result is
      interpretable (`SEED_LEN = 10`, `MINHASH_NUM_FUNCS = 100`, `MINHASH_KMER_LEN = 8`,
      `MAX_N_FRAC = 0.5`). Name the known divergence honestly: catch-lite implements
      the Hamming coverage model (equivalent to LCF-k with `lcf_threshold = probe_len`,
      the CATCH default), not general LCF-k.
- [ ] **sdust vs dustmasker** — masked-interval agreement (Jaccard, boundary offsets)
      at matched `T` and window. `dustmasker` ships with the already-required
      `blast>=2.12`. Cheapest credible "the reimplementation reproduces the original"
      result in the whole roadmap.
- [ ] **Thermodynamics vs published values** — a reference table of duplex ΔG°₃₇ / Tm
      from SantaLucia (1998) plus an independent calculator, asserted as unit tests.
      The only numeric test today (`src/thermodynamics.rs:400-414`) recomputes the
      formula from the same constants, which proves nothing. Document and bound the
      remaining deliberate deviations: the `has_stacking` guard suppressing
      initiation and salt correction below 2 consecutive WC pairs
      (`thermodynamics.rs:219-259`) and the `boltzmann_score` clamp to ≥ 1.0
      (`:277`), which Phase 2 removes from the retention path anyway.
- [ ] **Mismatch stacking** (folded in from
      [mismatch-stacking-thermo.md](mismatch-stacking-thermo.md)) — add `MM_TABLE`
      from SantaLucia & Hicks (2004) Table 4, cross-checked against biopython
      `MeltingTemp.DNA_IMM1`, and the third branch in `delta_g()`. Do this in the
      *same sitting* as the reference-value table above: same paper, same cross-check
      dictionary. Sanity check: a probe with one central mismatch must score strictly
      between a perfect match and a full mismatch, rather than equal to the latter as
      today.
- [ ] **BaitBench vs RAmpSim thermodynamics** — port and run
      `~/Projects/CrossLab/RAmpSim/compare/thermo_compare.py`, which already
      implements both scorers and compares them on real probe alignments. Report the
      result. Addresses the uncited-overlap exposure alongside Phase 3.
- [ ] **State the RNA limitation** in the manuscript: Figure 1's Zika example is an
      RNA target scored with DNA-DNA parameters; RNA:DNA (Sugimoto 1995) and RNA:RNA
      (Xia 1998) tables are future work per
      [rna-dna-thermodynamics.md](rna-dna-thermodynamics.md).
- [ ] **Unit tests for the untested reimplementations.** `src/probes/catch.rs`,
      `syotti.rs`, `probetools.rs` have **zero** tests. Priority: Syotti's
      no-false-negative seeding claim (`src/probes/syotti.rs:151-153`) and
      ProbeTools' pure helpers (`parse_clstr`, `percentile_coverage_for_targets`,
      `find_low_coverage_regions`, `expand_run`). De-duplicate `build_kmer_index` /
      `load_fasta`, copied verbatim between `catch.rs` and `syotti.rs`.

*Done when:* two agreement tables exist in-repo and the thermodynamics reference
values pass in CI.

---

## Phase 5 — Reviewer items 2–3: rebuild the validation (blocked by Phase 2 and **D2**, L)

### What exists now — do not reuse as-is

`~/Projects/CrossLab/RAmpSim/compare/`, **untracked**, with a hardcoded
`MM2="/Users/niel/miniconda3/bin/minimap2"`. `run_comparison.sh` steps 4–6 are the
circularity: map real reads → per-species counts → write those counts as `--sample`
weights → set `--capture-fraction` to the observed on-target fraction (0.0535) →
correlate against the same counts. Targets are `TestRAmpSim/D6331.refseq/ssrRNAs`
(Gut standard); reads are `TestRAmpSim/NN07_2025-07-16_FBA62613/barcode07` (in-house
ONT); probes are `Jared_MEGARESv3_Jan2024_1_Probes.txt`.

### Steps

- [ ] **Reads:** SRA BioProject **PRJNA751055**. Mock community = ZymoBIOMICS CSII /
      **D6310**, PacBio Sequel, **3 replicates** at 2 / 5 / 8 kb fragment lengths —
      the replicate structure doubles as a test of `--fragment-length-*`.
- [ ] **Input weights:** manufacturer's theoretical D6310 composition from the Zymo
      product sheet (resolves ~~D6~~). Independent of the reads being predicted.
- [ ] **Probe panel:** per **D2**. The TELSeq probes are not public (Agilent ELID
      0792071, "previously described"). Either obtain them from the authors, or
      reconstruct from MEGARes v1.0 + BacMet + VFDB with `baitbench build-probes` and
      declare the reconstruction as a limitation.
- [ ] **Capture parameters:** stop supplying the observed on-target fraction. Fit
      `--capture-fraction` / `--background-retention` **once** on replicate 1, report
      the fitted values, and apply them unchanged to replicates 2–3 as held-out
      predictions. State this explicitly in the manuscript.
- [ ] **Metrics:** Spearman ρ, Pearson on log₁₀, Lin's CCC, per-replicate and pooled,
      with an explicit statement that the input is theoretical composition.
- [ ] **Make the harness reproducible:** move it into the repo as `validation/telseq/`
      — SRA download, MEGARes fetch + `build-probes`, theoretical-composition TSV,
      the BaitBench invocation, a port of `real_data_comparison.py`, and a README.
      Drop the hardcoded minimap2 path; use the embedded aligner. This is what
      answers "insufficient to convince this reader."
- [ ] Regenerate Figure 2 and rewrite `BB_Paper/paper.md:64-72`.
- [ ] **D4:** decide whether the in-house D6331 ONT run becomes a deposited secondary
      dataset.
- [ ] **D7:** decide whether the Borrelia LOD titration becomes a second figure. It
      is the stronger of the two candidates — non-circular by construction, four
      orders of magnitude of dynamic range, and it demonstrates that the Phase 2 fix
      works rather than merely asserting it. Needs the same SRA deposit as D4, plus
      `validation/borrelia_lod/` from Phase 2.

### The existing numbers are broken regardless (**D5**)

- `paper.md` says ρ = 0.884 — matches **no** surviving results table.
- The submitted Oct 1 `.docx` contains `0.916884` and `p=38.4×10-67`: a mangled
  tracked-changes merge of `0.916` over `0.884`, and `3`+`8.4×10-6`+`7`. The same
  corruption repeats in the Figure 2 caption. This went out in the submitted PDF.
- `letter.docx` quotes the clean ρ = 0.916, CCC = 0.954 — the **badread** run
  (`results/badread_comparison_table.tsv`).
- `n = 20` is really 21 species with *C. albicans* at 0/0 silently dropped.

All superseded by the new numbers, but worth a line in the response letter as
evidence the figures were re-derived. **Work from `paper.md`, never the Oct 1 `.docx`.**

---

## Phase 6 — Reviewer item 1: AI disclosure and missing declarations (blocked by **D3**, S)

Full and specific. Written for the editor and reviewer.

What a reviewer can already see on the repo: `CLAUDE.md` (29 KB) and
`memory/MEMORY.md` are both committed; `improvements/*.md` are AI-authored design
notes; 2 of 148 commits carry `Co-Authored-By: Claude …`, which is what puts "claude"
in the GitHub contributor list. Development ran Feb–Oct 2026 across ~105 k lines of
Rust, ~12 k of R, ~3 k of GUI TypeScript/Svelte.

- [ ] **Manuscript — new "Use of generative AI" subsection.** State the tool (Claude
      Code with Claude models, Anthropic); what it was used for (the Rust CLI, the
      R/RMarkdown report templates, the Tauri GUI, documentation, and the
      reimplementations of published algorithms); that algorithm selection, model
      design, parameter choices, and validation were directed and verified by the
      authors; that the authors take full responsibility for correctness; and that no
      AI tool is an author. Match the venue's required heading per **D3**.
- [ ] **Cover letter.** A short, plain paragraph: the declaration was checked
      incorrectly during simultaneous submissions, this was an error rather than
      concealment, and the manuscript now carries a full disclosure. Lead with it.
- [ ] **Repo.** Same statement in `README.md` and `CITATION.cff`.
- [ ] **Response letter.** Note that the AI-authored planning documents
      (`improvements/`) are retained deliberately, as a record of which design
      questions were considered and rejected.
- [ ] **Add the missing declaration sections**, absent from the manuscript entirely
      and required by most venues: Data Availability, Code Availability (currently
      only the prose `## Distribution`), Author Contributions, Competing Interests
      (exists only in the cover letter), Acknowledgements (only "Funding
      Acknowledgement" exists).

---

## Phase 7 — Docs, manuscript rebuild, response letter (blocked by 1–6, M)

Per `CLAUDE.md`'s update table, the Phase 2 model change touches:

- [ ] `src/cli.rs` help for `--capture-fraction`, `--distractor-fraction`, and the
      new flags
- [ ] `docs/reference/{parameters,commands,output-formats}.md` (new metrics columns)
- [ ] `docs/explanation/simulation-modes.md:103-111` (the incorrect "ensures a
      fragment ratio" claim), `thermodynamic-scoring.md:104-120` (describes the
      uniform level-1 draw), `pipeline-overview.md:50`,
      `metrics-and-classification.md`
- [ ] `docs/how-to/set-parameters.md:27-40` and
      `docs/how-to/choose-simulation-mode.md:45-59` — both advertise unreachable `df`
- [ ] `CLAUDE.md` Metrics Definitions + Weight Generation; `ARCHITECTURE.md`
- [ ] `gui/src/lib/views/RunView.svelte` for the new flags, then `make copy-sidecar`
- [ ] `BB_Paper/paper.md`: rewrite `## Validation Using Public Data` (64–72), replace
      the limitation paragraph that defends the abundance-blind capture model, update
      `### Simulating Capture`, add the Phase 6 sections, rebuild via `build.sh`
- [ ] Point-by-point response letter covering all six comments

---

## Dependency graph

```
Phase 0 ──┬─> Phase 1 ──> Phase 2 (D1) ──> Phase 5 (D2) ──┐
          ├─> Phase 3                                      ├─> Phase 7
          ├─> Phase 4                                      │
          └─> Phase 6 (D3) ─────────────────────────────────┘
```

## Verification

```bash
conda activate baitbench
export PATH="$HOME/.cargo/bin:$PATH"
cargo build --release && cargo test && cargo clippy -- -D warnings
```

- **Item 6 regression:** sweep `--distractor-fraction` over
  `{0.1, 0.5, 0.9, 0.99, 0.999}` on `examples/minimal/` at fixed seed; assert
  `realized_distractor_fraction` tracks the request monotonically across the whole
  range under `--capture-model retention`, and that `quota` reproduces v1.0 numbers
  bit-for-bit at a fixed seed.
- **Backwards compatibility:** re-run the three `CLAUDE.md` smoke commands
  (standard, simple, sample-manifest) plus a genome-mode run and a `coverage-curve`
  run; confirm no panics and that reports render.
- **Borrelia LOD gate (Phase 2):** `validation/borrelia_lod/` reproduces the real
  titration shape across NTC → Pos1. This is the primary evidence that the capture
  rework works.
- **Phase 4:** `validation/catch_comparison/` and the dustmasker comparison produce
  agreement tables; thermodynamics reference values and the new mismatch-stacking
  tests pass as unit tests; the RAmpSim thermo comparison produces a reported number.
- **Phase 5:** `validation/telseq/run.sh` reproduces Figure 2 end-to-end from SRA
  accessions on a clean checkout, with no absolute paths.
- **Docs:** `mkdocs build --strict`.
