# Future Improvements

**[Paper revision roadmap](paper-revision-roadmap.md)** — governing document for the
manuscript revision: the six reviewer comments, seven open decisions, and eight
phases of work. Supersedes several of the capture entries below; read it first.

- [Capture: quality-weighted resampling](capture-quality-weighted-resampling.md) — weight fragment retention by alignment quality score rather than pass/fail (Option A) **(superseded by the shipped per-fragment retention model, Stage 3 — kept for the rejected-alternatives record)**
- [Capture: probabilistic model](capture-probabilistic-model.md) — sigmoid-based continuous capture probability function calibrated to target fold enrichment (Option C) **(superseded by Stage 3; its per-run calibration is exactly what the shipped model avoids)**
- [Capture: probe binding competition](capture-probe-competition.md) — weight probe selection by expected yield to fix over-enrichment of rare-species-specific probes (Option D; D1 is a one-line change) **(superseded by Stage 3: probe competition is implicit in `1 - prod(1 - theta)`)**
- [Thermodynamics: RNA/DNA hybrid parameters](rna-dna-thermodynamics.md) — add Sugimoto (1995) RNA:DNA and Xia (1998) RNA:RNA NN tables for RNA-target capture
- [Thermodynamics: mismatch stacking](mismatch-stacking-thermo.md) — replace SkipStacking with SantaLucia & Hicks (2004) internal mismatch parameters (~4–6 hours, data-entry-heavy)
- [Sweep: read length](sweep-read-length.md) — add `--read-length-values` to `coverage-curve` to sweep across read lengths and measure organism discrimination vs read length
- [Sequencing: paired-end support](paired-end-support.md) — extend simulate/sequence/filter/map pipeline to emit R1+R2 FASTQ pairs; already partially threaded through (`SequenceArgs.output_r2`)
- [Sequencing: long read testing](long-read-testing.md) — validate and benchmark the Badread long-read simulator path; add long-read presets and a `coverage-curve` long-read mode
- [Validation: mutation distribution tolerance](mutation-distribution-validation.md) — test whether BaitBench shows TITAN-RNA's uniform-vs-natural mutation tolerance gap (10% vs 27%, PMID 40501929)
- [Aligner: reuse seed-extend](seed-extend-wider-use.md) — the local seed-and-extend written for `simulate` probably fits xreact, assess-probes, panel-qc; also consolidates three duplicate k-mer index implementations
- [Metrics: coverage uniformity / fold-80](coverage-uniformity-metrics.md) — report the field's standard hybrid-capture uniformity statistic; the depth vector already exists, but simulated fold-80 reads optimistic until PCR GC bias and duplication are modelled, so frame it comparatively
- [Sequencing: badread-lite](badread-lite.md) — **future work, design settled, not scheduled.** Native multi-threaded long-read simulation. Badread does 36 reads/s single-threaded (~70 min per SimVal replicate) and is alignment-bound, not interpreter-bound. Build the *parametric* design if ever built: it sidesteps badread's GPL error-model tables, and removes both edlib call sites rather than inheriting them. Measured to be sufficient — read realism moves Lin's CCC by 0.0013 (B0g). No trigger to start; revisit if a reviewer asks for error-model sensitivity.
- [Capture: partial-duplex occupancy](partial-duplex-capture.md) — scale θ by how much of a probe's site a fragment actually carries, replacing the min-overlap threshold shipped in Stage 3
- [Cache: per-genome probe placements](per-genome-probe-cache.md) — **before next release**: index each genome once per panel and compose, instead of one cache per experimental mix; note the global `max_hits_per_probe` composition trap
- [Capture: parameterisation](parameterisation.md) — **next major piece**: the occupancy scale is not interpretable across temperature/salt, equilibrium ΔG cannot discriminate loci, and the fitting objective needs held-out data
