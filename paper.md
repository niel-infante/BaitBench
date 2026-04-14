# BaitBench: an easy tool for building and assessing probes, and predicting outcomes for target sequence capture


## Abstract

## Introduction

Target capture good. 
Target capture important.

not much to assess probes

virus/whole genome different than targeted site/wes

Installing tools can be hard





## The Tool

BaitBench is a command line tool written primarily in Rust. It is available as source code or precompiled binaries. 
Conda?


### Building Probes

By default, BB build pipeline will do QC and simplifications steps of cd-hit concatenation of sequences, removing sequences of more than 5% Ns, and removing short sequences. The actual probe building is done with either catch (ref), a built-in tiling algorithm, or built in versions of catch or syotti(ref). After building, probes are filtered for complexity (internal sDust(ref)) and GC content. Assess-probes is automatically run, and in addition to the information that provides (see next section), it ...

### Assessing Probes

Assess-probes is automatically run after building probes, but it can also be run on probes built with other tools, or rerun on probesets using different parameters. BaitBench first aligns all probes to all targets using minimap2 (ref). This gives a wealth of information, presented in the report first with a small summary of target coverage coverage and multimapping probes. Then a full searchable table of all targets is giving listing 


### Coverage Report
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





## References


