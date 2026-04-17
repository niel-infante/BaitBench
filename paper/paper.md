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

Probe panel design is relatively mature with tools such as CATCH [@metskyCapturingSequenceDiversity2019], Syotti [@alankoSyottiScalableBait2022], tiling, ProbeTools [@ProbeToolsDesigningHybridization], but evaluating performance before wet-lab experiments remains largely manual or ad hoc
- No standardized way to predict: sensitivity for a given target abundance, specificity against distractors, within-panel cross-reactivity, coverage uniformity
- Wet-lab iteration is expensive and slow; in-silico evaluation is underutilized
- Existing simulators (if any are cited) typically assume uniform coverage — don't model binding affinity or thermodynamics
- RAmpSim [Zhang 2025] introduced thermodynamic nearest-neighbor (TNN) Boltzmann-weighted simulation for metagenomics but is scope-limited to the simulation step; does not include probe design, probe QC, or structured performance metrics






## The Tool

BaitBench is a command line tool written primarily in Rust. It is available as source code or precompiled binaries. (Not yet)
Conda?


### Building Probes

By default, BB build pipeline will do QC and simplifications steps of cd-hit concatenation of sequences, removing sequences of more than 5% Ns, and removing short sequences. The actual probe building is done with either catch [@metskyCapturingSequenceDiversity2019], a built-in tiling algorithm, or built in versions of catch  or syotti [@alankoSyottiScalableBait2022]. After building, probes are filtered for complexity (internal sDust(ref)) and GC content. Assess-probes is automatically run, and in addition to the information that provides (see next section), it gives a summary of sequences input, filtered at each step, number of probes built, and number filtered at each step. 

### Assessing Probes

Assess-probes is automatically run after building probes, but it can also be run on probes built with other tools, or rerun on probesets using different parameters. Host genomes can be added to test probe specificity. BaitBench first aligns all probes to all targets using minimap2 [@liMinimap2PairwiseAlignment2018]. This gives a wealth of information, presented in the report first with a small summary of target coverage coverage and multimapping probes. Then a full searchable table of all targets is giving listing coverage statistics. Graphs and tables show coverage and gaps for targets in multiple manners, to give insight into probe performance. Probes are then mapped to themselves to assess cross reactivity and competition between probes. Graphic summaries are given, and problematic probes identified.
If host genomes are given, probes are mapped against them, and any potential danger areas are identified.
When working with multiple highly similar targets, coverage can become problematic, as minimap2, even with high secondary alignment settings, will not always find all alignments for a probe. In order to address this issue, BaitBench has various refinement options where it will rerun the mapping with only low coverage targets, even supporting multiple rounds of refinement. In extreme circumstances the --all-individual-targets flag can be used, which will map probes to each target in isolation. The choice of which mode to use is driven by the Biological goal of the project, and how many of the targets are expected to be in any single sample. 



### Run Simulation

Assessing probes assumes a (near) perfect world. Here we try to simulate a more realistic capture experiment incorporating all steps in the process. 

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


