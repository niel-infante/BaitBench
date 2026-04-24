# BaitBench Paper — Outline Thoughts

---

## Introduction

### Background: Target capture sequencing
- Hybridization-based target capture (bait capture) is a widely used method for enriching specific genomic regions prior to sequencing
- Essential for detecting low-abundance targets in complex backgrounds: clinical metagenomics, pathogen surveillance, antimicrobial resistance (AMR) gene detection, viral genomics, ancient DNA, whole-exome sequencing (WES)
- Probes (baits) bind to complementary sequences in fragmented DNA; unbound material is washed away; bound fragments are sequenced
- Reference review for the field: Bravo et al. (2025) — covers applications, limitations, and open computational problems

### Applications / why it matters
- Clinical diagnostics: detecting low-titer pathogens (Rickettsia, SARS-CoV-2, arboviruses) directly from clinical specimens [Paskey 2024, Nagy-Szakal 2021, Kamaraj 2019]
- Pathogen surveillance panels: broad-range viral or bacterial detection from environmental and clinical metagenomes
- Epidemiological genomics: full-genome recovery from complex samples (Lassa outbreak, ZIKV/DENV co-infections [Metsky 2019, Kamaraj 2019])
- AMR surveillance from microbiomes
- Ancient DNA recovery from degraded samples

### The design–experiment gap
- Probe panel design is relatively mature (CATCH, Syotti, tiling, ProbeTools), but evaluating performance before wet-lab experiments remains largely manual or ad hoc
- No standardized way to predict: sensitivity for a given target abundance, specificity against distractors, within-panel cross-reactivity, coverage uniformity
- Wet-lab iteration is expensive and slow; in-silico evaluation is underutilized
- Existing simulators (if any are cited) typically assume uniform coverage — don't model binding affinity or thermodynamics
- RAmpSim [Zhang 2025] introduced thermodynamic nearest-neighbor (TNN) Boltzmann-weighted simulation for metagenomics but is scope-limited to the simulation step; does not include probe design, probe QC, or structured performance metrics

### What BaitBench does
- Unified, end-to-end tool: design probes → assess probes → simulate a capture experiment → quantitative performance metrics → HTML report
- Thermodynamic simulation: binding affinities derived from SantaLucia (1998) nearest-neighbor model, Boltzmann-weighted fragment sampling — more realistic than uniform-coverage approaches
- Two operational modes: standard (small genomes such as viruses) and genome mode (large pathogens such as bacteria where probes target a gene region within a full genome)
- Direct integration of clinical context: qPCR CT scores translate to target abundances, making simulations interpretable alongside real diagnostic data
- Written in Rust for performance; optional desktop GUI for accessibility

---

## Methods

### Overview
- BaitBench is a command-line tool written in Rust; interactive GUI via Tauri desktop application
- Modular subcommand design: each pipeline step is independently callable; `baitbench run` chains all steps automatically
- External dependencies: minimap2 [Li 2018] for alignment, cd-hit [Li & Godzik 2006] for sequence clustering, R/ggplot2/RMarkdown for report generation; BLAST+ optionally for cross-reactivity

### Probe building (`baitbench build-probes`)
- Input: target FASTA(s); output: filtered, deduplicated probe FASTA
- Pre-processing steps: remove sequences with >5% ambiguous bases (N-filter), remove sequences shorter than probe length, collapse near-identical sequences with cd-hit-est
- Probe design methods (user-selectable via `--method`):
  - **Tiling**: sliding window across each target with configurable stride; simple and exhaustive
  - **CATCH-lite**: native Rust reimplementation of the CATCH optimization algorithm [Metsky 2019] — tile candidates → MinHash LSH deduplication → greedy set-cover to minimize probe count while guaranteeing coverage
  - **Syotti-lite**: native Rust reimplementation of the Syotti greedy set-cover design [Alanko 2022] — k-mer hash index, seed-and-extend coverage tracking; linear scaling
  - **CATCH** (external): calls the original Python CATCH tool if installed
- Post-processing filters: GC content bounds, low-complexity masking via sDUST [Morgulis 2006], final deduplication with cd-hit-est
- Automatically chains into `assess-probes` unless `--skip-assess` is specified

### Probe assessment (`baitbench assess-probes`)
- **Probe coverage analysis**: align probes to all targets with minimap2; compute per-target coverage depth, percent covered at various thresholds (5×, 20×), gap identification, multimapping probes
- **Cross-reactivity**: probes vs. genomes (host or other organisms), probes vs. themselves (internal competition); highlights probes with off-target or self-complementary binding
- **Panel QC** (`baitbench panel-qc`): all-vs-all target similarity; per-species discriminability scores; confusion matrix of species pairs sharing targets — identifies which organisms may be indistinguishable before any experiment is run
- Refinement options for panels with many similar strains: rerun mapping on low-coverage targets, or `--all-individual-targets` (map probes to each target in isolation)
- Output: single self-contained HTML report with interactive tables (DT) and ggplot2 figures

### Simulation pipeline (`baitbench run`)

#### Step 1 — Prepare
- Combine target and distractor FASTAs into a single reference
- Generate per-sequence sampling weights:
  - Sample targets: weight from user manifest (default 1.0); non-sample targets: weight = 0
  - Distractors: `distractor_weight = (distractor_fraction × total_sample_weight) / (n_distractors × (1 − distractor_fraction))`
- CT score conversion (alternative to `--distractor-fraction`): `target_fraction = ct_baseline_fraction × 2^(ct_baseline − ct)`; default calibration: CT 20 → 1% target DNA
- Genome mode: generates two references — a combined genome+distractor reference for fragment generation, and a target+distractor reference for read mapping; sample-target-map links genome IDs to target IDs

#### Step 2 — Simulate (thermodynamic fragment generation)
- Align probes to combined reference with minimap2; parse CIGAR + MD tags to reconstruct per-position (probe_base, ref_base) pairs for each alignment
- **Thermodynamic scoring**: compute ΔG (Gibbs free energy) for each probe-reference alignment using the SantaLucia (1998) nearest-neighbor model via a `ThermoModel` struct (temperature + salt concentration)
  - *NN stacking*: accumulate stacking energy over consecutive Watson-Crick pairs; mismatches break the stacking chain (SkipStacking strategy)
  - *Initiation terms*: add AT (+2.3 kcal/mol ΔH, +4.1 cal/mol/K ΔS) or GC (+0.1, −2.8) initiation penalty for the first and last WC terminal of each alignment (SantaLucia 1998 Table 2)
  - *Salt correction*: adjust ΔS for actual Na+ concentration via Owczarzy et al. (1997): `ΔS += 0.368 × (n_wc−1) × ln([Na+])`; user-specified via `--salt-concentration` (mM, default 50 mM); at 1 M the correction is exactly zero
  - Convert to Boltzmann binding score: `score = exp(−ΔG / RT)` at user-specified hybridization temperature
- **Two-level multinomial fragment sampling** for captured reads:
  1. Sample a probe uniformly from probes with ≥1 alignment hit
  2. Sample an alignment hit for that probe, weighted by `Boltzmann_score × sequence_weight`
  3. Fragment center: alignment center ± uniform jitter (±fragment_length/4)
  4. Fragment length: sampled from truncated normal distribution (user-specified mean, SD, min, max)
- Background fragments (fraction `1 − capture_fraction`): sampled uniformly weighted by `sequence_weight × sequence_length`
- Capture fraction (`--capture-fraction`): controls ratio of probe-biased to background fragments; models incomplete capture efficiency in real experiments
- Target enrichment is emergent — no imposed fold-enrichment parameter

#### Step 3 — Sequence
- Trim fragments to read length (current implementation); architecture is designed for drop-in read simulators (ART-modern, PBSIM2) for realistic error profiles
- Optional: subsample to target sequencing depth

#### Step 4 — Filter (optional)
- Map reads against host genome(s); discard mapping reads; models host depletion step in real workflows

#### Step 5 — Map
- Align reads to combined reference with minimap2; configurable preset and secondary alignment settings

#### Step 6 — List
- Parse SAM; count reads per reference sequence

#### Step 7 — Metrics (3-way classification)
- Classification at genome/group level (was each target detected?):
  - **TP**: sample target detected
  - **FN**: sample target not detected
  - **FP_target**: non-sample target within panel detected (within-panel cross-reactivity)
  - **TN_target**: non-sample target within panel not detected
  - **FP_distractor**: distractor detected (off-target capture)
  - **TN_distractor**: distractor not detected
  - **Untargeted**: genome-mode genomes with no target mapping (tracked separately)
- Summary metrics: sensitivity, specificity, precision, F1-score
- Coverage statistics: per-reference depth, pct_covered_5x, pct_covered_20x
- Read-level tracking: correctly mapped, incorrectly mapped, source vs. mapping destination

#### Step 8 — Report
- Self-contained HTML generated via RMarkdown/ggplot2
- Sankey diagram of fragment flow (generation → capture → sequencing → filtering → mapping)
- Performance metrics bar charts, detection detail lollipop chart, coverage depth plots
- Parameters section with reconstructed CLI command for reproducibility

### Additional modules
- **Coverage curves** (`baitbench coverage-curve`): parameter sweeps over CT × hybridization temperature × capture fraction × sequencing depth; produces coverage depth curve plots to identify detection limits and experimental requirements
- **Species identification** (`baitbench identify`): ordered-explanation algorithm — call PRESENT if unique marker targets detected, ABSENT when all hits explained by cross-reactivity, AMBIGUOUS when indeterminate; uses cross-reactivity knowledge from xreact module
- **Cross-reactivity** (`baitbench xreact`): standalone probe cross-reactivity check against genomes and/or other probes

### Group-level metrics
- Optional `--groups` TSV: maps sequence IDs to logical group names (e.g., multiple strain variants → one species group)
- All classification (TP/FP/FN/TN) computed at group level; a group is detected if any member has mapped reads
- Distractor contigs auto-grouped by source FASTA file; `--distractor-groups` overrides with explicit mapping

### Sample manifest and discrimination testing
- TSV of `id<TAB>weight` specifying which targets are "present" in the simulated sample (and at what relative abundance)
- Enables discrimination testing: can probe set distinguish dengue serotype 1 from dengue 2 within a panel?
- Without manifest: all targets are in-sample, reducing to traditional 2-way sensitivity/specificity

### Genome mode for large pathogens
- For bacteria and other organisms where probes target a gene region (e.g., 16S, ompB) rather than the full genome
- Fragments generated from full genomes; reads mapped back to target regions only
- `--sample-target-map` TSV links genome IDs to their target IDs (supports 1-to-1, 1-to-many, many-to-1)
- Auto-linking by exact name match or prefix convention (`genome_id|target_name`)

### Use cases
- **Basic probe panel evaluation**: test capture efficiency on a set of viral targets with host as distractor
- **Sample discrimination**: test whether a panel can specifically detect a focal virus within a multi-virus panel
- **Clinical specimen modeling**: use CT score to set realistic target abundance; compare predicted coverage depth to assay detection thresholds
- **Bacterial targeted capture**: genome mode for 16S or other gene targets within full genome sequences
- **Panel optimization**: iterate probe design and assess-probes to close coverage gaps before synthesis
- **Host depletion validation**: test specificity after simulated host filtering
- **Multi-species surveillance panels**: group-level metrics for panels covering many closely related strains
- **Detection limit estimation**: coverage curve sweeps to determine minimum sequencing depth for reliable detection
- **Pre-experiment QC**: panel-qc to identify pairs of organisms that cannot be discriminated by the current probe set

---

## Discussion

### Strengths and innovations
- **End-to-end, unified workflow**: no other tool covers probe design → QC → thermodynamic simulation → classification metrics → reports in one package
- **Thermodynamic realism**: TNN nearest-neighbor scoring (SantaLucia 1998) produces probe-site-specific binding affinities; fragment enrichment near high-affinity sites emerges naturally, without an imposed fold-enrichment parameter — more faithful to real hybridization physics than uniform-coverage models
- **3-way classification**: separately quantifying FP_target (within-panel cross-reactivity) vs. FP_distractor (true off-target capture) provides actionable diagnostic information that a simple positive/negative call cannot
- **CT score integration**: bridging qPCR-derived abundance estimates to simulation inputs makes predictions directly comparable to real clinical data; enables "what CT score is detectable?" queries
- **Genome mode**: extends applicability to bacteria and other large pathogens; practically important for clinical surveillance panels that mix viral and bacterial targets
- **Rust implementation**: memory-safe, fast; streaming I/O avoids loading large files into memory; practical for large reference databases
- **Accessibility**: desktop GUI (Tauri) lowers barrier for non-command-line users; precompiled binaries planned

### Comparison to related tools
- **RAmpSim** [Zhang 2025]: shares TNN thermodynamic core and Boltzmann-weighted sampling; BaitBench differs by providing the complete design/QC/metrics pipeline, 3-way classification, clinical CT integration, genome mode, and interactive reports. RAmpSim hands off to external read simulators for platform-specific error modeling; BaitBench's sequencing step is currently trimming-only but is architecture-ready for external simulator wrapping.
- **CATCH** [Metsky 2019] and **Syotti** [Alanko 2022]: probe design tools; BaitBench reimplements both natively in Rust and integrates them into a quality-controlled build pipeline followed by performance assessment
- **ProbeTools**: probe design for diverse/hypervariable viral taxa; focused on design, not simulation or performance evaluation
- Uniform-coverage simulation approaches: do not model binding affinity differences across probe sites; overestimate coverage uniformity

### Limitations
- Read sequencing is currently simulated as simple trimming (no PCR amplification bias, adapter artifacts, GC-dependent coverage bias, platform error profiles) — a drop-in interface for ART-modern (Illumina) or PBSIM2 (long reads) is planned
- Library preparation steps (end repair, A-tailing, adapter ligation, size selection) are not modeled
- Hybridization kinetics model is equilibrium-based (Boltzmann); does not model probe concentration, hybridization time, or wash stringency dynamics explicitly (capture fraction parameter is a pragmatic proxy)
- Syotti-lite does not use an FM-index (unlike the original Syotti); performance on very large datasets (>1 GB) is reduced
- CATCH-lite reimplementation may differ subtly from the original Python CATCH in probe selection when probe sets are equivalent by the set-cover objective
- Pan-genome and structural variation not currently modeled

### Future directions (from paper.md notes)
- Wrap external read simulators (ART-modern for Illumina short reads, PBSIM2 for long reads) to produce realistic error profiles
- Native long-read support throughout the pipeline
- Probe editing: tools to remove redundant or low-utility probes, and to fill coverage gaps with targeted probes
- GUI enhancements
- Documentation: example CLI recipes (e.g., generating a groups file from a FASTA)
- Benchmarking against real capture datasets to validate thermodynamic predictions

---

## Key Innovations Summary

1. **Thermodynamic probe-site scoring** via SantaLucia (1998) nearest-neighbor free energy model, converted to Boltzmann binding affinity — applied at base-pair resolution to minimap2 alignments
2. **Two-level multinomial sampling**: probe selection → alignment site selection → fragment generation; enrichment emerges from binding physics, not from an imposed parameter
3. **3-way genome-level classification**: distinguishes within-panel cross-reactivity from true off-target capture — clinically meaningful distinction
4. **CT score → abundance conversion**: directly links simulation parameters to qPCR-measurable clinical abundances
5. **Genome mode**: handles the biologically important case where probes target a genetic locus within a full genome (bacteria, large DNA viruses)
6. **Sample manifest + group-level metrics**: enables discrimination testing and handles multi-strain panels as logical taxonomic units
7. **Native reimplementations of CATCH and Syotti** in Rust within an integrated, quality-controlled probe build pipeline
8. **Species identification module**: ordered-explanation algorithm using cross-reactivity knowledge to call PRESENT/ABSENT/AMBIGUOUS from multi-target detection patterns
9. **Coverage curve parameter sweeps**: systematic exploration of experimental parameter space to find detection limits
10. **End-to-end pipeline with interactive HTML reports**: actionable visualization from a single command

---

## Suggested Citations

### Already in BaitBenchRef.json (verify all before submitting)
- **Alanko et al. 2022** — Syotti bait design
- **Bravo et al. 2025** — Review of bait capture enrichment methods and computational challenges (excellent for framing the introduction)
- **Kamaraj et al. 2019** — Targeted enrichment for Dengue/Zika/Chikungunya
- **Li & Godzik 2006** — cd-hit
- **Li 2018** — minimap2
- **Metsky et al. 2019** — CATCH
- **Nagy-Szakal et al. 2021** — Hybrid capture for SARS-CoV-2
- **Paskey et al. 2024** — Targeted enrichment for rickettsial pathogens
- **ProbeTools** (Springer Nature link in refs — confirm full citation)
- **Zhang et al. 2025** — RAmpSim (thermodynamic simulator)

### Need to add to BaitBenchRef.json
- **SantaLucia JR. 1998** — "A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor thermodynamics." *PNAS* 95(4):1460–1465. doi:10.1073/pnas.95.4.1460 — **critical**: provides all ΔH/ΔS stacking parameters used in the TNN model
- **Morgulis et al. 2006** — sDUST low-complexity sequence masking — "A Fast and Symmetric DUST Implementation to Mask Low-Complexity DNA Sequences." *J. Computational Biology* 13(5):1028–1040. doi:10.1089/cmb.2006.13.1028
- **Gnirke et al. 2009** — "Solution hybrid selection with ultra-long oligonucleotides for massively parallel targeted sequencing." *Nature Biotechnology* 27:182–189 — **seminal** bait capture method paper
- **Albert et al. 2007** or **Hodges et al. 2007** — original exome capture papers (context for WES applications)
- Read simulator citations if/when integrated: **Huang et al. 2012** (ART) or ART-modern, **Ono et al. 2021** (PBSIM2/PBSIM3)
- Any paper describing the original tiling/sliding window probe design approach
- Possibly **Bolger et al. 2014** (Trimmomatic) or similar for context on sequencing artifact modeling

### Possible additional citations depending on use cases presented in results
- Papers using the specific probe panels BaitBench is validated against
- Clinical metagenomics review papers for introduction context
- AMR surveillance papers if genome-mode bacteria use case is featured

---

## What Else Should Be Added to the Paper

### Abstract (currently empty)
Suggested elements:
- One sentence on target capture sequencing importance
- One sentence stating the gap (no unified in-silico evaluation tool)
- One sentence describing BaitBench (what it does, thermodynamic approach)
- One sentence on key capabilities (probe design, assessment, simulation, metrics)
- One sentence on availability

### Availability / Implementation section
- Source code (GitHub URL)
- Precompiled binaries
- Conda environment installation (`conda activate baitbench`)
- Operating systems supported
- License
- Desktop GUI availability

### Figure plan (not yet outlined)
Consider:
1. Pipeline overview schematic (design → assess → simulate → report)
2. Thermodynamic scoring illustration (probe aligned to reference, ΔG calculation, Boltzmann weighting)
3. Fragment sampling schematic (two-level multinomial)
4. Example output figures from the HTML report (Sankey, lollipop chart, coverage plot)
5. Coverage curve example (parameter sweep)
6. Probe assessment report example (coverage heatmap or gap visualization)

### Methods: implementation details worth mentioning
- FASTA I/O is native Rust (no seqtk dependency); streaming reduces memory footprint
- SAM parsing extracts CIGAR + MD tags natively to reconstruct base-pair alignments for TNN scoring
- All intermediate files are documented; pipeline can be re-entered at any step

### Possibly a "Quick Start" or example walkthrough
- A minimal example (the `examples/minimal/` data) showing the full workflow from command line
- Could be a supplementary figure or box

### Limitations section or paragraph
- Currently in Discussion bullets above, but worth making it a named subsection for transparency

### Note on RAmpSim relationship
- May be worth a brief "Note" or paragraph explicitly clarifying the relationship with RAmpSim [Zhang 2025] since the thermodynamic core is conceptually very similar — the distinction (full pipeline vs. simulation-only, different scope and audience) should be stated clearly and positively

---

## Weaknesses, Reviewer Concerns, and Pre-submission Checklist

Labels: **[CRITICAL]** = reviewer may reject without this; **[MAJOR]** = must be addressed or explicitly justified; **[MINOR]** = should be discussed or noted; **[NIT]** = polish, unlikely to block acceptance.

---

### Thermodynamic Model Limitations

- ~~**[CRITICAL] No initiation penalty.**~~ **FIXED.** SantaLucia (1998) initiation terms (AT: +2.3/+4.1 kcal/mol; GC: +0.1/−2.8 kcal/mol) are now applied to the first and last WC-paired terminal of each alignment. All 8 unit tests pass including two new tests verifying exact values.

- **[MAJOR] SkipStacking is a significant simplification.** The current strategy breaks the stacking chain completely on any mismatch — consecutive WC pairs separated by even a single mismatch contribute zero energy for that step. The published literature has more nuanced treatments (e.g., unified nearest-neighbor models that include mismatch stacking penalties). This may underestimate binding affinity for probes with scattered single mismatches and overestimate the penalty of mismatches in otherwise high-affinity alignments. This choice is defensible but needs explicit justification and ideally a brief sensitivity analysis.

- ~~**[MAJOR] No salt concentration correction.**~~ **FIXED.** The Owczarzy et al. (1997) entropy correction `ΔS += 0.368 × (n_wc−1) × ln([Na+])` is now applied. A new `--salt-concentration` flag (mM, default 50 mM) was added to `baitbench run`, `baitbench simulate`, and `baitbench coverage-curve`. The `ThermoModel` struct encapsulates both temperature and salt concentration. At 1 M Na+ the correction is exactly zero (ln(1)=0), preserving backward compatibility. Salt concentration is logged and written to `run_params.tsv`.

- **[MINOR] No initiation for duplex ends (dangling end terms).** Dangling end stacking contributions from unpaired terminal bases are ignored. These are small but non-negligible for short probes.

- **[MINOR] Single-stranded secondary structure of probes not modeled.** Probe hairpins and self-dimers reduce effective hybridization affinity; this is not accounted for. For short (80-mer) probes this is usually small, but can matter for GC-rich probes.

- **[MINOR] Probe concentration is implicit.** The Boltzmann score captures binding affinity but not probe concentration. When multiple probes compete for the same site, relative concentrations affect the effective enrichment. BaitBench treats all probes as equally concentrated; unequal probe concentrations in synthesized panels are not modeled.

- **[MINOR] Equilibrium assumption.** The model assumes thermodynamic equilibrium (Boltzmann weighting). Real hybridization is kinetically controlled; probes reaching equilibrium depends on hybridization time, temperature ramp, and probe concentration. A single equilibrium temperature is an approximation of a temperature ramp protocol.

---

### Simulation Realism

- **[CRITICAL] No experimental validation.** This is almost certainly the first thing reviewers will ask. The thermodynamic model predicts coverage distributions and enrichment patterns, but without a benchmark against real capture sequencing data these remain theoretical claims. RAmpSim (Zhang 2025) validated against empirical coverage distributions using earth mover's distance. BaitBench will need at least one comparable validation showing predicted vs. observed coverage profiles, enrichment fold, or sensitivity/specificity, on a real dataset. Without this, the paper is a methods description with unverified assumptions.

- **[CRITICAL] Read sequencing is just trimming.** The current "sequencing" step trims fragments to read length with no error model. Real Illumina reads have: base-quality score degradation toward 3' ends, GC-dependent coverage bias, PCR duplicate patterns, adapter contamination, and index hopping. The claimed sensitivity/specificity metrics are computed on error-free reads mapped back to the exact reference they were generated from — this is highly optimistic and may not reflect real-world performance. This limitation should be prominently stated and ideally partially addressed (e.g., by integrating ART-modern or running a sensitivity analysis with ART-generated reads).

- **[MAJOR] Detection threshold is hardcoded at ≥1 read.** A reference is classified as "detected" if it has at least one mapped read. This is very permissive — a single spuriously mapped read calls a target positive. There is no configurable minimum coverage threshold, no statistical test, no FDR control. This inflates sensitivity and deflates FP counts. At minimum this should be configurable; ideally a brief analysis showing how metrics change across thresholds (e.g., ≥1, ≥5, ≥10 reads) would support the choice.

- **[MINOR, well-motivated] Capture fraction is a single global scalar.** The `--capture-fraction` parameter applies uniformly to all sequences. Reviewers may ask whether this is redundant with thermodynamic scoring — it is not. Boltzmann weighting models *differential enrichment at the site level* (which probe-binding positions are preferentially captured); `--capture-fraction` models *overall pull-down efficiency* (what fraction of total DNA enters the captured pool at all, regardless of probe affinity). These are orthogonal parameters: the former emerges from ΔG; the latter depends on probe concentration, hybridization time, and wash stringency, which are outside the thermodynamic model. The paper should make this distinction explicit to pre-empt reviewer confusion. The remaining limitation (per-probe or per-GC-content variation in capture efficiency within the captured fraction) is already substantially handled by the per-site Boltzmann weighting.

- **[MAJOR] No PCR duplicate modeling.** Real capture libraries undergo PCR amplification. Duplicate reads from the same original fragment introduce coverage non-uniformity (jackpotting), particularly for low-input libraries. This is one of the largest sources of coverage non-uniformity in real data and is completely absent from the simulation.

- **[MAJOR] Fragment size distribution.** Fragments are drawn from a truncated normal distribution. Real DNA fragmentation (sonication, enzymatic) produces distributions that are closer to log-normal or empirical, with a pronounced mode and heavier right tail. This affects coverage uniformity at target boundaries.

- **[MINOR] No library prep biases.** End repair, A-tailing, and adapter ligation have sequence-context biases that are not modeled. These create systematic coverage dips at certain sequence motifs (e.g., the ENCODE "low-mapability" regions often have library prep artifacts).

- **[MINOR] No index hopping / cross-contamination.** When multiplexing multiple samples, a small fraction of reads are assigned to the wrong sample index. This can create phantom low-level signals — relevant for the 3-way classification if the threshold is ≥1 read.

- **[MINOR] Circular genomes not handled at boundaries.** Fragments sampled near the start or end of a circular genome (most bacteria) will be truncated rather than wrapping around. This creates artificial coverage gaps at genome edges. Not a large effect for long genomes, but worth noting.

- **[NIT] Only single-end reads.** Paired-end sequencing is standard for most target capture protocols. Paired reads allow duplicate detection, improved mapping at repetitive regions, and realistic insert size distributions. The current single-end model is a recognized simplification.

---

### RAmpSim Overlap — The Hardest Reviewer Question

- **[CRITICAL] The thermodynamic core is nearly identical to RAmpSim.** Both BaitBench and RAmpSim use SantaLucia (1998) nearest-neighbor ΔG, Boltzmann-weighted fragment sampling, and a two-level multinomial procedure. A reviewer familiar with RAmpSim (Zhang 2025) will immediately ask: "What does BaitBench add beyond RAmpSim, and why should the community adopt it instead of RAmpSim?" The paper must have a clear, confident answer — not just scope (pipeline vs. simulation) but also whether the implementations differ in detail (e.g., initiation terms, stacking strategy, how MD tags are parsed). If BaitBench's TNN implementation differs from RAmpSim's in consequential ways, that should be stated.

- **[MAJOR] No head-to-head comparison with RAmpSim.** A benchmark showing that BaitBench's simulation output is comparable (or better) to RAmpSim's on a shared dataset would greatly strengthen the paper. If RAmpSim is not comparable because they serve different purposes, the paper should explain what RAmpSim produces that BaitBench does not, and vice versa.

---

### Probe Design Tool Validation

- **[MAJOR] CATCH-lite and Syotti-lite not benchmarked against originals.** The paper describes native Rust reimplementations of both tools. Reviewers will ask: are they equivalent? Under what conditions (dataset size, diversity level) do they diverge from the originals in probe count, coverage, or set composition? Without a validation experiment, these are unverified reimplementations. At minimum, a comparison on the example dataset or a public benchmark dataset is needed.

- **[MINOR] Syotti-lite lacks the FM-index used in original Syotti.** The original Syotti achieves linear scaling via succinct data structures. Without the FM-index, the Rust reimplementation's scaling properties on large (>1 GB) datasets are untested. This is acknowledged in paper.md but not formally characterized.

- **[MINOR] No runtime/memory benchmarks for probe design methods.** For a tool paper, providing timing on realistic datasets (a few thousand viral targets, a large bacterial collection) is expected.

---

### Metrics and Classification Design

- **[MAJOR] The 3-way classification rationale needs justification.** Why is separating FP_target from FP_distractor clinically meaningful? The paper should provide a concrete worked example where they diverge in implication (e.g., a false positive from a non-sample panel member is different from a false positive from a host sequence because the former implies the panel cross-reacts while the latter implies contamination). Reviewers may not find this obvious.

- **[MAJOR] Coverage thresholds (5×, 20×) are arbitrary.** The metrics pct_covered_5x and pct_covered_20x use thresholds that are standard in some fields (WES) but not universally applicable. For viral metagenomics (the main use case), these may not be the right thresholds. Justification or configurability is needed.

- **[MINOR] Species identification algorithm is novel but unvalidated.** The ordered-explanation algorithm for PRESENT/ABSENT/AMBIGUOUS calling is an original contribution but is not benchmarked against any ground truth or compared to existing taxonomic classifiers (Kraken2, Bracken, etc.). For a paper, at least a simulated ground-truth test (put known species in at known abundances, measure recall and precision of calls) would be expected.

- **[MINOR] Reads correctly mapped is defined by source == mapping destination, but minimap2 may map reads ambiguously.** When two targets are highly similar, a read generated from target A may correctly map to target B (because B is equally similar). The current metric would count this as incorrectly mapped, but it is biologically ambiguous. The paper should clarify how this case is handled.

- **[NIT] No statistical uncertainty on metrics.** Sensitivity/specificity are point estimates with no confidence intervals or bootstrap error bars. For small panels (5–10 targets), a single FN changes sensitivity by 20%. Reporting uncertainty would be more rigorous.

---

### Usability and Reproducibility

- **[MAJOR] No Docker/Singularity container.** The tool requires a Conda environment with minimap2, BLAST+, cd-hit, R, and multiple R packages. This is a non-trivial installation. Many journals now expect or recommend containers for reproducibility. A Dockerfile or Bioconda package would significantly lower the barrier.

- **[MINOR] No runtime or memory benchmarks.** For a Rust tool emphasizing performance, users expect to know how long a typical run takes and how much memory it needs. A benchmark table on representative dataset sizes (small virus panel, large surveillance panel, bacteria in genome mode) would be standard for a tool paper.

- **[MINOR] GUI not formally described.** The Tauri desktop GUI is mentioned but not documented in the paper. If it's a key feature for accessibility, it deserves a brief description and screenshot. If it's not ready for prime time, it may be better to omit to avoid reviewer questions about it.

- **[NIT] Seed reproducibility across platforms.** BaitBench uses a user-specified `--seed` for the Rust RNG. However, minimap2 has its own internal randomness (e.g., in tie-breaking for secondary alignments). Full bit-for-bit reproducibility across platforms or minimap2 versions is not guaranteed, which should be noted.

- **[NIT] Intermediate file formats not described in the paper.** The pipeline produces many intermediate files (fragments.fa, captured.fa, coverage.tsv, etc.). A supplementary table describing each file's format and purpose would be useful for users building custom workflows around individual steps.

---

### What a Reviewer Will Specifically Ask

1. "Please validate your thermodynamic predictions against at least one real capture sequencing dataset. What metric do you use to compare predicted vs. observed distributions?"

2. "RAmpSim (Zhang et al. 2025) uses the same SantaLucia nearest-neighbor thermodynamic model with Boltzmann-weighted sampling. Please clearly differentiate your approach, and provide a quantitative comparison."

3. "Your implementation omits initiation parameters from the SantaLucia (1998) model. Please justify this choice or correct it and re-evaluate."

4. "Detection is called at ≥1 mapped read. This threshold seems very permissive. Please provide a sensitivity analysis showing how sensitivity, specificity, and precision change as you increase the detection threshold."

5. "Have the CATCH-lite and Syotti-lite reimplementations been validated against the original CATCH and Syotti tools? Please provide a comparison."

6. "The sequencing step is simulated as simple trimming with no error model. How does the absence of sequencing errors affect your benchmarking results? Please discuss or demonstrate with a read-error-model comparison."

7. "Please provide runtime and memory benchmarks on typical use-case dataset sizes."

8. "The capture fraction parameter is a single global value. Please discuss the scenarios where this approximation breaks down and how sensitive your metrics are to its value."

---

### What Must Be Done Before Submitting

| Priority | Item |
|----------|------|
| Must | Validate thermodynamic simulation predictions against at least one real capture dataset |
| Must | Head-to-head comparison with RAmpSim (or clear argument for why comparison is inappropriate) |
| Must | Add SantaLucia (1998) to BaitBenchRef.json |
| ~~Must~~ | ~~Address or explicitly justify omission of initiation terms in the TNN model~~ — **DONE: initiation terms implemented** |
| ~~Must (was MAJOR)~~ | ~~Add salt correction to ΔG~~ — **DONE: Owczarzy et al. (1997) correction implemented via `--salt-concentration`** |
| Must | Discuss the ≥1-read detection threshold and its implications; ideally make it configurable |
| Should | Validate CATCH-lite and Syotti-lite against the original tools on a benchmark dataset |
| Should | Provide runtime/memory benchmarks |
| Should | Provide sensitivity analysis on capture fraction and hybridization temperature |
| Should | Characterize how metrics change with detection threshold |
| Should | Docker/Conda package for easy installation |
| Should | Make coverage thresholds (5×, 20×) configurable or justify them |
| Consider | Run a subset of results through ART-modern to show the effect of read errors |
| Consider | Validate species identification module against simulated ground truth |
| Consider | Formally describe (or exclude) the GUI |
| Must | Find a reference / justification of how CT is calculated |

---

### Extended Future Directions

*(In addition to those already in paper.md: read simulator integration, long reads, probe editing, GUI)*

- **Configurable detection threshold**: allow users to set minimum reads (or minimum coverage depth) for calling a target "detected"; report a precision-recall curve or ROC curve across thresholds
- **Confidence intervals on metrics**: bootstrap or Poisson-based uncertainty on sensitivity/specificity for small panels
- **Salt and buffer correction**: implement Owczarzy et al. (2004) or SantaLucia & Hicks (2004) salt correction to make temperature predictions more accurate under real buffer conditions
- **Initiation term and terminal AT/GC penalty**: complete the SantaLucia (1998) model; small code change, would improve accuracy for short or AT-rich probes
- **PCR duplicate modeling**: simulate amplification of captured fragments to produce realistic duplicate distributions
- **Paired-end read output**: generate paired FASTQ files with realistic insert size distributions
- **Probe concentration modeling**: allow unequal probe concentrations (e.g., from synthesis yield variation or deliberate spiking) to affect Boltzmann sampling weights
- **Circular genome support**: wrap fragment sampling at genome boundaries for circular chromosomes
- **Integration with sequencing simulators**: formal plugin interface for ART-modern (Illumina), PBSIM3 (PacBio/ONT), Badread (long read with error models)
- **FM-index for Syotti-lite**: would enable scaling to large bacterial datasets comparable to the original Syotti
- **Real-time parameter optimization**: given a target coverage depth at a given CT score, back-calculate the required number of sequences or capture fraction
- **Multi-sample simulation**: simulate pooled clinical samples with different organisms at different abundances; compare to multiplexed sequencing data
- **FDR control for species identification**: provide p-values or FDR-corrected calls for the PRESENT/ABSENT/AMBIGUOUS classifier
- **Benchmarking against Kraken2/Bracken**: compare BaitBench's species ID module to standard k-mer classifiers on capture-enriched data
- **Support for degenerate/IUPAC probe bases**: probes with IUPAC ambiguous bases (used for diverse target panels) are not currently modeled in the TNN scoring
- **Export to common probe synthesis formats**: output probe sequences in formats accepted by probe synthesis vendors (Agilent SureSelect, IDT xGen, Twist Bioscience)
