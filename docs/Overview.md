# Overview

BaitBench is a tool to help construct and assess target sequence capture probesets. 

Target sequence capture, also called probe capture or hybrid capture, is a molecular technique used to selectively enrich specific DNA or RNA sequences from a complex mixture — such as a clinical sample or environmental metagenome — before sequencing, so that sequencing effort concentrates on regions of interest rather than the entire sample. It works by hybridizing biotinylated oligonucleotide probes (baits) to complementary target sequences in the sample, then using streptavidin-coated beads to pull down the probe-target hybrids and wash away non-binding background. The key inputs are the probe sequences (baits), the target sequences the probes are designed to capture, and background (distractor) sequences representing everything else in the sample that the probes should ideally ignore.



The tools are broadly grouped into three sections.  

## Construction of a probeset from target sequence

BaitBench assumes you know what you want to capture and have sequence for it. Typically this will consist of whole genomes for targets such as viruses, or specific genes for organisms with larger genomes. From this sequence will construct a probeset consisting of k-mers that hopefully cover the entire target sequence without cross-reactivity between probes or non-specific sequence. 

Currently BaitBench offers three different algorithms for probe construction: Simple tiling, Syotti, or Catch. See construction or the manual for more details.

## Assessment of probesets





## Simulation of Capture











