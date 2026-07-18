# Future Improvements

- [Capture: quality-weighted resampling](capture-quality-weighted-resampling.md) — weight fragment retention by alignment quality score rather than pass/fail (Option A)
- [Capture: probabilistic model](capture-probabilistic-model.md) — sigmoid-based continuous capture probability function calibrated to target fold enrichment (Option C)
- [Capture: probe binding competition](capture-probe-competition.md) — weight probe selection by expected yield to fix over-enrichment of rare-species-specific probes (Option D; D1 is a one-line change)
- [Thermodynamics: RNA/DNA hybrid parameters](rna-dna-thermodynamics.md) — add Sugimoto (1995) RNA:DNA and Xia (1998) RNA:RNA NN tables for RNA-target capture
- [Thermodynamics: mismatch stacking](mismatch-stacking-thermo.md) — replace SkipStacking with SantaLucia & Hicks (2004) internal mismatch parameters (~4–6 hours, data-entry-heavy)
- [Sweep: read length](sweep-read-length.md) — add `--read-length-values` to `coverage-curve` to sweep across read lengths and measure organism discrimination vs read length
- [Sequencing: paired-end support](paired-end-support.md) — extend simulate/sequence/filter/map pipeline to emit R1+R2 FASTQ pairs; already partially threaded through (`SequenceArgs.output_r2`)
- [Sequencing: long read testing](long-read-testing.md) — validate and benchmark the Badread long-read simulator path; add long-read presets and a `coverage-curve` long-read mode
