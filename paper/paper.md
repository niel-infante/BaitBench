# BaitBench: an easy tool for building and assessing probes, and predicting outcomes for target sequence capture


## Abstract

## Introduction

Target capture good. 
Target capture important.

not much to assess probes

virus/whole genome different than targeted site/wes

Installing tools can be hard





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

Assessing probes assumes a (near) perfect world. Here we try to simulate a more realistic capture experiment. 

The coverage report is intended for assessing samples with a single, or small number of targets present. BaitBench simulates all stages of the target capture workflow. The user specifies what targets and distractors are in the sample with relative weights. BaitBench then randomly generates fragments of the given DNA. Capturing is simulated incorporating the thermodynamic properties of the probe - fragment binding using the RAmpSim [@zhangRAmpSimThermodynamicSimulator2025] algorithm. A tune-able fraction of all the DNA is randomly selected regardless of capture to represent the incomplete nature of the capture process. Sequencing is then simulated (currently just trimming all fragments to a read length, but drop-in ready for a read simulator such as .… to give realistic Illumina or long read sequences), followed by an optional distractor filtering step, and finally mapping to the targets and assessment of the sensitivity, specificity, and precision of the experiment. 







### Coverage Curve





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


