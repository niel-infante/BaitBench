# BaitBench: an easy tool for building and assessing probes, and predicting outcomes for target sequence capture


## Abstract

## Introduction

Target capture good. 
Target capture important.

not much to assess probes

virus/whole genome different than targeted site/wes

Installing tools can be hard


Simulating thermodynamics gives insight into realistic expectations, and an make predictions for needed sequencing depth.

Can simulate all steps in target sequence capture workflow.


 Hybridization-based target capture (bait capture) is a widely used method for enriching specific genomic regions prior to sequencing. It is essential for detecting low-abundance targets in complex backgrounds such as clinical metagenomics, pathogen surveillance, antimicrobial resistance gene detection, viral genomics, ancient DNA, or whole-exome sequencing.

 
- Probes (baits) bind to complementary sequences in fragmented DNA; unbound material is washed away; bound fragments are sequenced
- Reference review for the field: Bravo et al. (2025) — covers applications, limitations, and open computational problems

[@bravoMethodsApplicationsComputational2025]


### Applications / why it matters
- Clinical diagnostics: detecting low-titer pathogens (Rickettsia, SARS-CoV-2, arboviruses) directly from clinical specimens [Paskey 2024, Nagy-Szakal 2021, Kamaraj 2019]
- Pathogen surveillance panels: broad-range viral or bacterial detection from environmental and clinical metagenomes
- Epidemiological genomics: full-genome recovery from complex samples (Lassa outbreak, ZIKV/DENV co-infections [Metsky 2019, Kamaraj 2019])
- AMR surveillance from microbiomes
- Ancient DNA recovery from degraded samples

### The design–experiment gap

Probe panel design is relatively mature with tools such as CATCH [@metskyCapturingSequenceDiversity2019], Syotti [@alankoSyottiScalableBait2022], tiling, ProbeTools [@kuchinskiProbeToolsDesigningHybridization2022], but evaluating performance before wet-lab experiments remains largely manual or ad hoc. Further, there is no standardized way to predict sensitivity for a given target abundance, specificity against distractors, within-panel cross-reactivity, or coverage uniformity. Wet-lab iteration is expensive and slow; in-silico evaluation is underutilized and complicated by multiple tools with diverse dependencies and conventions. BaitBench aims to fill this gap with an easy to install, simple to use tool.



RAmpSim [Zhang 2025] introduced thermodynamic nearest-neighbor (TNN) Boltzmann-weighted simulation for metagenomics but is scope-limited to the simulation step; does not include probe design, probe QC, or structured performance metrics


## The Tool

BaitBench is a command line tool written primarily in Rust. It is available as source code or precompiled binaries. (Not yet)
Conda?


BaitBench is single tool in a single environment providing unified, end-to-end solution: design probes → assess probes → simulate a capture experiment → quantitative performance metrics → HTML report. BaitBench includes thermodynamic simulation of binding affinities derived from [@santaluciaUnifiedViewPolymer1998] as described by RAmpSim [@zhangRAmpSimThermodynamicSimulator2025], utilizing nearest-neighbor model, Boltzmann-weighted fragment sampling and more realistic than uniform-coverage approaches. BaitBench accommodates small genomes such as viruses and large pathogens such as bacteria where probes target gene regions within a set of genomes, or whole exome regions targeting a single genome. Direct integration of clinical context qPCR CT scores translate to target abundances, making simulations interpretable alongside real diagnostic data. BaitBench is written in Rust for performance, with an optional desktop GUI for accessibility.

BaitBench can be broadly split into three functionalities: building probes, assessing probes, and simulating capture.


### Building Probes

By default, the build pipeline will do QC and simplifications steps of cd-hit concatenation of sequences, removing sequences of more than 5% Ns, and removing short sequences. See FIG_BUILD The actual probe building is done with either catch [@metskyCapturingSequenceDiversity2019], a built-in tiling algorithm, or built in versions of catch  or syotti [@alankoSyottiScalableBait2022]. After building, probes are filtered for complexity via an internal sDust [@morgulisFastSymmetricDUST2006'] implementaion, and GC content. Assess-probes is automatically run, and in addition to the information that provides (see next section), it gives a summary of sequences input, filtered at each step, number of probes built, and number filtered at each step. 


FIG_BUILD:
 ![Build pipeline](../docs/diagrams/paper_build_probes.png)




### Assessing Probes

Assess-probes is automatically run after building probes, but it can also be run on probes built with other tools, or rerun on probesets using different parameters FIG_ASSESS. Host genomes can be added to test probe specificity. BaitBench first aligns all probes to all targets using minimap2 [@liMinimap2PairwiseAlignment2018]. This gives a wealth of information, presented in the report first with a small summary of target coverage coverage and multimapping probes. Then a full searchable table of all targets is giving listing coverage statistics. Graphs and tables show coverage and gaps for targets in multiple manners, to give insight into probe performance. Probes are then mapped to themselves to assess cross reactivity and competition between probes. Graphic summaries are given, and problematic probes identified.
If host genomes are given, probes are mapped against them, and any potential danger areas are identified.
When working with multiple highly similar targets, coverage can become problematic, as minimap2, even with high secondary alignment settings, will not always find all alignments for a probe. In order to address this issue, BaitBench has various refinement options where it will rerun the mapping with only low coverage targets, even supporting multiple rounds of refinement. In extreme circumstances the --all-individual-targets flag can be used, which will map probes to each target in isolation. The choice of which mode to use is driven by the Biological goal of the project, and how many of the targets are expected to be in any single sample. 



FIG_ASSESS
![Probe assessment diagram](../docs/diagrams/paper_assess_probes.png)

### Run Simulation

Assessing probes assumes a (near) perfect world. Here we try to simulate a more realistic capture experiment incorporating all steps in the process. BaitBench’s workflow is split into eight steps each of which can be run separately, and since all intermediate files are documented and retained (unless --cleanup is called) the pipeline can be re-entered at any step.



#### Step 1 — Prepare

The prepare step function is to create a single fasta file with all sequences, along with a weights file that specifies what is in the simulated sample, and in what proportion. The total amount of distractor sequence can specified by either a simple percentage, or a hypothetical CT (cycle threshold) value from qPCR quantifies pathogen abundance — higher CT means lower concentration, with each unit representing a two-fold dilution. BaitBench converts CT to a target DNA fraction using $\text{target\_fraction} = \text{ct\_baseline\_fraction} \times 2^{\text{ct\_baseline} - \text{ct}}$ defaulting to a calibration point of CT 20 = 1% target DNA, so CT 25 yields ~0.03% and CT 30 ~0.001%; the remainder becomes the distractor fraction, directly linking real clinical sample measurements to simulation parameters. in genome mode two references are generated, a combined genome+distractor reference for fragment generation, and a target+distractor reference for read mapping; sample-target-map links genome IDs to target IDs.


FIG_PREPARE_1

![Simple prepare diagram](../docs/diagrams/prepare_mode1_standard_nosample.png)



FIG_PREPARE_3

![genome prepare diagram](../docs/diagrams/prepare_mode3_genomes_nosample.png)





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





The coverage report is intended for assessing samples with a single, or small number of targets present. BaitBench simulates all stages of the target capture workflow. The user specifies what targets and distractors are in the sample with relative weights. BaitBench then randomly generates fragments of the given DNA selected weighted by the user provided weights and the genome sizes. Fragment size follows a normal distribution bound by user specified min, max,and mean. Capturing is simulated incorporating the thermodynamic properties of the probe - fragment binding using the RAmpSim [@zhangRAmpSimThermodynamicSimulator2025] algorithm. A tune-able fraction of all the DNA is randomly selected regardless of capture to represent the incomplete nature of the capture process. Sequencing is then simulated (currently just trimming all fragments to a read length, but drop-in ready for a read simulator such as ART-modern or PBSIM2), followed by an optional distractor filtering step, and finally mapping to the targets and assessment of the sensitivity, specificity, and precision of the experiment. We track every fragment from generation through to the mapping, so we can give a very fine tuned view 


### Coverage report
The coverage report is intended for assessing samples with a single, or small number of targets present. BaitBench simulates all stages of the target capture workflow. The user specifies what targets and distractors are in the sample with relative weights. BaitBench then randomly generates fragments of the given DNA. Capturing is simulated incorporating the thermodynamic properties of the probe - fragment binding using the RAmSim (ref) algorithm. A tune-able fraction of all the DNA is randomly selected regardless of capture to represent the incomplete nature of the capture process. Sequencing is then simulated (currently just trimming all fragments to a read length, but drop-in ready for a read simulator such as .… to give realistic Illumina or long read sequences), followed by an optional distractor filtering step, and finally mapping to the targets and assessment of the sensitivity, specificity, and precision of the experiment. 













### Coverage Curve

This module allows user to do a parameter sweep over some key parameters: capture fraction, temperature, number of sequences generated, and initial fraction of desired sample present. The resulting coverage curve give users insights into the effort needed to reach coverage sufficient for their downstream analyses. Capture fraction is a measure of how clean the capture is, and how much bleed through there is. Generally speaking, longer capture times and a cleaner wash will lead to a higher capture fraction. 



## Results
Testing on some real datasets
Use cases

### Differences with rust implementation
Speed up,
Not everything implemented
    Syotti - No FM-index, so not great on large data (> 1GB)
    Catch - Different, but how?

## Discussion


### Future Directions
wrap sequence simulator(s)
support long reads
GUI
editing probests 
    Remove redundant or useless probes
    Add probes to targeted, low coverage areas

In documentation, include some helpful cli commands, such as how to create a sample group file from a target.fa



## References


