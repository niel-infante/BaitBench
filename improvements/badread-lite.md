# badread-lite: native multi-threaded long-read simulation

Scoped 2026-11-25. **Future work — parked 2026-11-26, design settled, not scheduled.**

No trigger to start. The question that would have forced it — whether read realism
affects the validation — was measured and answered no (`parameterisation.md` B0g), so
nothing is blocked on this. Revisit if a reviewer asks for error-model sensitivity
analysis, or if a future feature makes read simulation throughput matter.

**If it is built: build design C** (parametric). It is smaller, faster than a faithful
port, and avoids the licence entirely. Reasoning below.

## The measured problem

Badread is single-threaded Python. Measured on this machine, simulating PacBio HiFi
reads from 150 000 captured fragments:

| | |
|---|---|
| throughput | **36 reads/s** |
| one SimVal replicate (150k reads) | **~70 min** |
| full SimVal gate (3 replicates) | **~3.5 h**, in read simulation alone |
| CPU | 99.1% — one core, on a 10-core machine |

For comparison, the rest of that pipeline run takes ~6 minutes: probe placement is
cached, mapping is parallel (6.2x since it was fixed), and simulation is seconds.
**Read simulation would be ~90% of the gate's wall time.**

This is why the gate uses `--read-simulator perfect`. The cost is not that error-free
reads are indefensible — real CCS median identity is 99.83%, so they are close — but
that we *cannot afford to check*, which is a bad position to argue from in a paper.

Reads are independent, so this is embarrassingly parallel. A Rust implementation at
8 workers would plausibly turn 70 minutes into well under a minute.

## What Badread actually is

Small, and simpler than its output suggests. The core is ~1 300 lines:

| module | lines | role |
|---|---|---|
| `simulate.py` | 536 | main loop, fragment assembly, adapters, chimeras, glitches |
| `qscore_model.py` | 334 | context → quality-score distribution |
| `error_model.py` | 229 | k-mer → observed-variant distribution |
| `fragment_lengths.py` | 117 | gamma-distributed lengths |
| `identities.py` | 103 | beta-distributed target identity |

The algorithm is straightforward: draw a fragment length, draw a target identity, walk
the sequence in overlapping 7-mers substituting each for a sampled variant, tune the
substitution rate until the read hits its target identity, then assign quality scores
from the qscore model.

**The realism does not live in the algorithm — it lives in two empirical data files:**

- `error_models/<chem>.gz` — **16 384 lines**, one per 7-mer (4^7), each mapping that
  7-mer to a distribution over what it gets sequenced as. Example:
  `AAAAAAA,0.963678;AAAAAA,0.035852;AAAAA,0.000052;…` — note the second entry is a
  deletion, which is how homopolymer-length error is encoded.
- `qscore_models/<chem>.gz` — **494 lines**, keyed by alignment context (`overall`,
  `=`, `===`, …), each mapping to a quality-score distribution.

Five chemistries ship: `nanopore2018/2020/2023`, `pacbio2016`, `pacbio2021`.

## The licensing constraint — this decides the design

**Badread is GPL-3.0. BaitBench is MIT.** Reproduced from `simulate.py`:

> Copyright 2018 Ryan Wick … Badread is free software: you can redistribute it and/or
> modify it under the terms of the GNU General Public License … version 3

Two separable questions:

**The algorithm is fine to reimplement.** It is published — Wick (2019), *Badread:
simulation of error-prone long reads*, JOSS 4(36):1316 — and the project already has
this pattern twice over. `src/probes/catch.rs` implements Metsky et al. 2019 from the
Nature Biotechnology paper; `src/probes/syotti.rs` implements Alanko et al. 2022 from
the Bioinformatics paper and explicitly swaps the FM-index for a k-mer hash. Both cite
the paper, not the source. Same approach applies here.

**The model files are not.** They are GPL-licensed data, and they are the substance of
badread's realism. Shipping them inside an MIT binary would be a derivative-work
problem, not a grey area.

That single fact rules out the obvious "port it and vendor the models" plan, and forces
a choice between three real designs.

## Where the time actually goes — read the loop before porting it

It is tempting to assume 36 reads/s is Python interpreter overhead and that a port
therefore wins by a large constant factor. **It is mostly alignment.** `edlib` is
imported by three of the five core modules and called in two hot paths:

- `sequence_fragment` realigns every 25 accepted changes (`ALIGNMENT_INTERVAL = 25`)
  over a 1 kb window (`ALIGNMENT_SIZE = 1000`), to correct its own running estimate of
  the read's identity.
- `get_qscores` performs a **full read-versus-fragment alignment for every read**. The
  author's own comment sits directly above it: `# TODO: I fear this full sequence
  alignment will be slow for long and inaccurate sequences.`

A faithful port pays the same algorithmic cost. It would still win — a Rust edit-distance
implementation plus 8 workers is a large improvement — but the win comes from threading
and constant factors, not from removing the dominant cost.

Why does Badread need the alignments at all? Because it applies errors by **k-mer
substitution**: it picks a random position, replaces that 7-mer with a sampled variant
(which may insert or delete), and therefore loses track of the true identity as
coordinates shift. The targeting loop is correspondingly heuristic:

```python
errors += new_errors * (estimated_identity ** 1.5)
```

That exponent is an empirical fudge with no first-principles justification. Matching
Badread means porting the fudge verbatim; doing something defensible instead means
designing a different targeting scheme and validating that it hits the requested
identity.

**This reverses the ranking of the options below.** Under design C errors are decided
per base at a known rate, so identity is exact by construction and the per-base error
state is known directly — which means *both* edlib call sites disappear. No realignment
to estimate identity, no alignment to assign qscores. C is therefore not the compromise
option; it is the technically superior one, and by a wide margin, **unless the 7-mer
error spectrum turns out to matter**.

Remaining fiddliness, in any design: `new_fragment_bases` holds variable-width entries
(`''` for a deletion, `'AT'` for an insertion), and qscore assignment builds a
`seq_pos -> alignment_pos` map with margin pullback at read ends. Mechanical, but easy
to get subtly wrong, and worth property-testing.

## Three designs

### A. Reimplement the algorithm, load badread's models at runtime

Parse the installed `error_models/` and `qscore_models/` files at startup.

- **Gets:** the full speed and parallelism win, bit-comparable output to badread.
- **Does not get:** dependency removal. Badread must still be installed, so
  `environment.yml` is unchanged and so is the install burden.
- **Licensing:** reading files the user already has is fine; we ship nothing.
- **Effort:** small-to-medium. The file formats are simple (CSV-ish, gzipped).

Honest framing: this solves the *speed* problem and nothing else. That may be the whole
point — speed is the measured problem.

### B. Regenerate equivalent models from public data

Badread documents a `badread error_model` subcommand that builds a model from real reads
plus alignments. Running that on public data produces our own tables under our own
licence.

- **Gets:** everything — speed, parallelism, dependency removal, MIT-clean.
- **Costs:** a reference dataset per chemistry, a model-building pipeline, and the
  burden of validating that our models behave like badread's. Model files are also
  ~100 KB each compressed, which is fine to ship.
- **Effort:** large, and most of it is data curation rather than code.

### C. Parametric error model, no empirical tables

Drop 7-mer conditioning. Sample errors to hit a target identity with a realistic
indel:substitution ratio and an explicit homopolymer-length bias term.

- **Gets:** speed, parallelism, dependency removal, entirely original code and data —
  **and it removes both edlib call sites**, which is the dominant cost (see above).
  Identity is exact by construction rather than estimated and corrected.
- **Loses:** the 7-mer error spectrum. Reads would have the right *identity* and the
  right *kind* of errors in aggregate, but not the exact context-dependent spectrum.
  Homopolymer-length bias has to be modelled explicitly rather than falling out of the
  k-mer table.
- **Effort:** small. Probably 400–600 lines including tests, versus ~800–1200 plus an
  edit-distance dependency for a faithful port.

## The question that decides between them

**Does BaitBench need 7-mer-level error realism, or only the identity distribution?**

Everything BaitBench reports from reads is compositional — per-reference counts, ratios,
enrichment factors. The evidence so far says the error *spectrum* is unlikely to matter:

- Changing the mapping preset dropped **8.7%** of simulated reads and moved Lin's CCC by
  **0.0003** (`parameterisation.md` B0e). Composition is a ratio; near-proportional
  losses cancel.
- The same section found rammap and minimap2 agree to within 0.02% at the same preset.

If identity distribution is what matters, **C is sufficient and is also the cheapest and
cleanest option.** If the spectrum turns out to matter — most plausibly for assigning
reads between near-identical references, like the five E. coli strains — then B is the
only option that is both faithful and MIT-clean.

**This has now been measured** (`parameterisation.md` B0g). A SimVal replicate under
`pacbio-hifi` versus the same run under `perfect`:

| | Spearman rho | r(log10) | Lin's CCC |
|---|---|---|---|
| `perfect` | 0.9324 | 0.9811 | 0.9688 |
| `badread pacbio-hifi` | 0.9324 | 0.9781 | 0.9675 |

Spearman identical to four decimals, CCC moving 0.0013 in the worse direction. The
simulated reads were genuinely realistic — 99.88% median identity against a real 99.83%
— so this is not a case of the test failing to apply the treatment.

**The error spectrum does not matter at BaitBench's resolution, so C is justified on
evidence**, not merely on convenience. Combined with the fact that C removes both edlib
call sites, it is the design to build if this is ever built.

## Implementation notes if it proceeds

- **Determinism.** Reads are independent, so parallelise with the same pattern as
  `map_one_file`: render in parallel with an indexed `par_iter`, write serially in input
  order. The mapping fix learned this the hard way — the previous code iterated a
  `HashMap` and produced non-reproducible output at a fixed `--seed`.
- **Seeding.** Each read needs a deterministic per-read seed derived from the global
  seed and the read index, or output depends on worker scheduling.
- **Scope boundary.** Badread's adapters, chimeras, glitches and junk reads are already
  exposed as flags. A reimplementation must either cover them or declare them
  unsupported; silently ignoring flags that currently do something would be worse than
  not having the feature.
- **Naming.** Follow the existing convention: `badread-lite` as a `--read-simulator`
  value, with badread itself kept available, so the two can be compared directly. That
  comparison is also the acceptance test.

## Should we do this at all

Not yet, and probably not before the manuscript.

It is a performance and dependency improvement, not a correctness one. Nothing currently
blocked is blocked *by* badread's speed — the gate uses `perfect`, and the open
modelling question is the global scale offset, which read simulation does not touch.

The case strengthens considerably if either of these becomes true:

1. The in-flight measurement shows read realism **does** move the validation metrics, in
   which case the gate needs badread and 3.5 hours per run is untenable.
2. A reviewer asks for error-model sensitivity analysis, which would mean many runs.

Until then this is a well-understood piece of work sitting on the shelf, which is where
it should stay. Revisit when the measurement lands.

## Related

- `improvements/parameterisation.md` B0e — mapper and preset sensitivity, the evidence
  that composition is insensitive to read-level detail.
- `src/probes/catch.rs`, `src/probes/syotti.rs` — the established
  reimplement-from-the-paper pattern and how provenance is documented.
- Plan Stage 6 (PCR duplicates) — a duplication model would also multiply read counts,
  which makes simulation throughput matter more than it does today.
