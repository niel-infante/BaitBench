# BaitBench: an easy tool for building, assessing, and predicting outcomes for target sequence capture probes

Aniello M Infante^1^, Shaun T Cross^1*^ 

1. Separtment of Environmental, Agricultural, and Occupational Health, University of Nebraska Medical Center, Omaha, Nebraska, USA

## Abstract

## Introduction

Hybridization-based target capture, sometimes called bait capture or hybrid-capture enrichment, has become an essential technique for selectively sequencing genomic targets of interest from complex nucleic acid mixtures. In the protocol, a set of biotinylated oligonucleotide probes hybridizes to complementary nucleic acids is a library, streptavidin-coated beads then pull down probe-bound fragments enriching for select targets with the unbound nucleic acid washed away, and the enriched library is sequenced. The result is orders-of-magnitude increase in on-target read depth enabling applications that are otherwise prohibitively expensive or impossible with whole sample sequencing alone. Hybrid target capture is essential for detecting low-abundance targets in complex backgrounds such as clinical metagenomics, pathogen surveillance, antimicrobial resistance gene detection, viral genomics, ancient DNA, and whole-exome sequencing [@bravoMethodsApplicationsComputational2025]. For many of these settings the target represents as little as 0.001–1% of the total DNA in the sample, making the design and performance of the probe set the decisive factor between a successful assay and a failed one.

The computational side of probe design has matured considerably. Tools such as CATCH [@metskyCapturingSequenceDiversity2019], Syotti [@alankoSyottiScalableBait2022], and ProbeTools [@kuchinskiProbeToolsDesigningHybridization2022] can efficiently construct minimal probe sets with guaranteed coverage across diverse, rapidly evolving viral targets. However, a systematic gap remains between designing a probe panel and knowing how it will perform. Prior to synthesis and wet-lab validation, there is no broadly adopted approach for predicting sensitivity as a function of target abundance, specificity against a background of host and distractor sequences, within-panel cross-reactivity, or coverage uniformity across target sequences. 

Wet lab iteration is slow and expensive; insilico evaluation remains largely ad-hoc and beyond the abilities of some bench scientists. When simulations are used, they commonly assume uniform per-base capture probability, an approximation that ignores the sequence dependent binding affinities that govern real hybridization.

Recent work on RAmpSim @zhangRAmpSimThermodynamicSimulator2025 demonstrated that thermodynamic simulation using the SantaLucia (1998) nearest-neighbor (TNN) model   [@santaluciaUnifiedViewPolymer1998] can produce substantially more realistic coverage distributions than uniform models. By computing the Gibbs free energy ΔG for each probe-reference alignment and converting it to a Boltzmann-weighted binding probability, fragment enrichment near high-affinity sites emerges naturally from the physics of hybridization rather than from an imposed enrichment parameter. However, RAmpSim addresses only the simulation step; it does not include probe design, structured probe quality control, or quantitative performance metrics, leaving users to assemble these capabilities from separate tools with incompatible formats and dependencies.

Here we present BaitBench, an end-to-end computational suite for designing, assessing, and benchmarking probe panels for target capture sequencing. BaitBench integrates native Rust reimplementations of the CATCH and Syotti probe design algorithms with a quality-controlled build pipeline, a comprehensive probe assessment module, and a full thermodynamic simulation pipeline based on the SantaLucia (1998) nearest-neighbor model, including initiation terms, Boltzmann-weighted fragment sampling, and a sodium concentration correction following Owczarzy et al. [@owczarzyPredictingSequencedependentMelting1997]. The simulation pipeline supports both standard mode for small-genome targets such as viruses and a genome mode for bacteria and other large pathogens where probes target a genetic locus within a full genome. To directly connect simulations to clinical practice, BaitBench accepts qPCR cycle-threshold (CT) scores as input, automatically converting them to target DNA fractions and enabling the question "at what CT is my assay expected to work?" to be answered insilico. Performance is evaluated with a three-way classification scheme that separately quantifies within-panel cross-reactivity and true off-target capture, a clinically meaningful distinction that a binary detected/not-detected call obscures. BaitBench also offers the ability to simulate how different temperatures, capture efficiencies, and sequencing depths will effect target coverage, giving some direction for experiment design. Written in Rust for performance and distributed with an optional desktop GUI, BaitBench is designed to lower the barrier to rigorous probe panel evaluation for both command-line and non-specialist users.

Probe panel design is relatively mature with tools such as CATCH [@metskyCapturingSequenceDiversity2019], Syotti [@alankoSyottiScalableBait2022], tiling, ProbeTools [@kuchinskiProbeToolsDesigningHybridization2022], but evaluating performance before wet-lab experiments remains largely manual or ad hoc. Further, there is no standardized way to predict sensitivity for a given target abundance, specificity against distractors, within-panel cross-reactivity, or coverage uniformity. Wet-lab iteration is expensive and slow; in-silico evaluation is underutilized and complicated by multiple tools with diverse dependencies and conventions. BaitBench aims to fill this gap with an easy to install, simple to use tool.


## The Tool

BaitBench is a command line tool written primarily in Rust. It is available as source code or precompiled binaries. (Not yet)
Conda? Docker?


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

### Simulating Capture

Assessing probes assumes a (near) perfect world. Here we try to simulate a more realistic capture experiment incorporating all steps in the process. BaitBench’s workflow is split into eight steps each of which can be run separately, and since all intermediate files are documented and retained (unless --cleanup is called) the pipeline can be re-entered at any step.



#### Step 1 — Prepare

The prepare step function is to create a single fasta file with all sequences, along with a weights file that specifies what is in the simulated sample, and in what proportion. The total amount of distractor sequence can specified by either a simple percentage, or a hypothetical CT (cycle threshold) value from qPCR quantifies pathogen abundance — higher CT means lower concentration, with each unit representing a two-fold dilution. BaitBench converts CT to a target DNA fraction using $\text{target\_fraction} = \text{ct\_baseline\_fraction} \times 2^{\text{ct\_baseline} - \text{ct}}$ defaulting to a calibration point of CT 20 = 1% target DNA, so CT 25 yields ~0.03% and CT 30 ~0.001%; the remainder becomes the distractor fraction, directly linking real clinical sample measurements to simulation parameters. in genome mode two references are generated, a combined genome+distractor reference for fragment generation, and a target+distractor reference for read mapping; sample-target-map links genome IDs to target IDs.


FIG_PREPARE_1

![Simple prepare diagram](../docs/diagrams/prepare_mode1_standard_nosample.png)



FIG_PREPARE_3

![genome prepare diagram](../docs/diagrams/prepare_mode3_genomes_nosample.png)





#### Step 2 — Simulate (thermodynamic fragment generation)

The simulate step is modeled directly on RAmpSim [@zhangRAmpSimThermodynamicSimulator2025]. Probes are aligned to the combined reference with minimap2 [@liMinimap2PairwiseAlignment2018], CIGAR and MD tags are parsed via an internal tool to reconstruct per-position (probe_base, ref_base) pairs for each alignment.  Gibbs free energy (ΔG) is calculated for each probe-reference alignment using the SantaLucia (1998) nearest-neighbor model via a `ThermoModel` struct (temperature and salt concentration).  NN stacking accumulates stacking energy over consecutive Watson-Crick pairs, mismatches break the stacking chain (SkipStacking strategy) Initiation terms add AT (+2.3 kcal/mol ΔH, +4.1 cal/mol/K ΔS) or GC (+0.1, −2.8) initiation penalty for the first and last WC terminal of each alignment (SantaLucia 1998 Table 2) Salt correction adjusts ΔS for actual Na+ concentration via Owczarzy et al. [@owczarzyPredictingSequencedependentMelting1997]: `ΔS += 0.368 × (n_wc−1) × ln([Na+])`; user-specified via `--salt-concentration` (mM, default 50 mM). At 1 M the correction is exactly zero. Convert to Boltzmann binding score: `score = exp(−ΔG / RT)` at user-specified hybridization temperature. Now we can use a Two-level multinomial fragment sampling for captured reads:
  1. Sample a probe uniformly from probes with ≥1 alignment hit
  2. Sample an alignment hit for that probe, weighted by Boltzmann_score × sequence_weight
  3. Fragment center: alignment center ± uniform jitter (±fragment_length/4)
  4. Fragment length: sampled from truncated normal distribution (user-specified mean, SD, min, max)
- Background fragments (fraction `1 − capture_fraction`): sampled uniformly weighted by sequence_weight × sequence_length. To model incomplete capture efficiency in real experiments we use the single parameter.
  Target enrichment is and emergent property of the thermodynamic sampling method. 


FIG_THERMO

![Thermodynamics algorithm](../docs/diagrams/paper_thermodynamic_scoring.png)


#### Step 3 — Sequence
We need to fix this. Right now we simply trim to sequence length. It wont be hard to add some read simulators. Probably want to add short read, long read, paired end.

#### Step 4 — Filter (optional)
Map reads against host genome(s); discard mapping reads; models host depletion step in real workflows

#### Step 5 — Map
Align reads to combined reference with minimap2; configurable preset and secondary alignment settings

#### Step 6 — List
Parse SAM; count reads per reference sequence

#### Step 7 — Metrics (3-way classification)
- Classification at genome/group level (was each target detected?):
  - **TP**: sample target detected
  - **FN**: sample target not detected
  - **FP_target**: non-sample target within panel detected (within-panel cross-reactivity)
  - **TN_target**: non-sample target within panel not detected
  - **FP_distractor**: distractor detected (off-target capture)
  - **TN_distractor**: distractor not detected
  - **Untargeted**: genome-mode genomes with no target mapping (tracked separately)

From this we are able to calculate summary metrics, coverage statistics: per-reference depth, pct_covered_5x, pct_covered_20x we implemented read-level tracking so we can report correctly mapped, incorrectly mapped, source vs. mapping destination.

#### Step 8 — Report
BaitBench produces a self-contained HTML report generated via RMarkdown (or an .Rmd file the user can alter for custom graphics). The report contains a sankey diagram of fragment flow (generation → capture → sequencing → filtering → mapping), performance metrics bar charts, detection detail lollipop chart, a confusion matrix, coverage depth plots, and interactive tables of useful metrics. Every report BaitBench gernerates includes a parameters section with a reconstructed CLI command for reproducibility.

### Coverage Curve

This module allows user to do a parameter sweep over some key parameters: capture fraction, temperature, number of sequences generated, and initial fraction of desired sample present. The resulting coverage curve gives users insights into the effort needed to reach coverage sufficient for their downstream analyses. 

FIG_COVERAGE_CURVE   - Need a better one, this one still has fold enrichment which is no longer a parameter.

![Coverage Curve](FIG_CovCurve.png)

### Other Modules
**Species identification** (`baitbench identify`)   When working with similar species, and targeting the same genes in each, there is the concern that even with perfect capture you may not be able to tell species apart. This tool will look at all of the targets of every species, and consider the homology between them. Species are then called PRESENT if unique marker targets detected, ABSENT when all hits explained by cross-reactivity, AMBIGUOUS when indeterminate.
 **Cross-reactivity** (`baitbench xreact`): Standalone probe cross-reactivity check against genomes and/or other probes based on homology. Also useful to check if your targets are close.
## 



## Results
Testing on some real datasets
Use cases

### Differences with rust implementation
Speed up,
Not everything implemented
    Syotti - No FM-index, so not great on large data (> 1GB)
    Catch - Different, but how?

## Discussion

## Limitations

## Distribution


### Future Directions
wrap sequence simulator(s)
support long reads
~~GUI~~
editing probests 
    Remove redundant or useless probes
    Add probes to targeted, low coverage areas

In documentation, include some helpful cli commands, such as how to create a sample group file from a target.fa


## References


