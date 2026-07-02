# What Is Probe Capture?

## Target sequence capture

Target sequence capture (also called hybridization capture or probe capture enrichment) is a molecular technique for selectively isolating specific DNA sequences from a complex mixture. A collection of short synthetic oligonucleotides called **probes** (or baits) is designed to be complementary to the sequences of interest. When probes are introduced into a solution of fragmented DNA, they hybridize to matching sequences by Watson-Crick base pairing. The probe-target duplexes are then physically pulled out of solution — typically using streptavidin beads that bind to biotin labels on the probes — and the captured DNA is eluted and sequenced.

The technique is widely used in clinical and research genomics:

- **Viral metagenomics**: enriching for viral sequences in samples that are predominantly host DNA
- **Bacterial detection**: targeting conserved gene regions (16S rRNA, *gyrB*, *rpoB*) in clinical specimens
- **Whole-exome sequencing**: capturing all protein-coding regions of a genome
- **Panel sequencing**: targeted capture of a defined gene panel for clinical diagnostics

## Why evaluate a probe panel in silico?

The efficiency of probe capture depends on many factors: probe GC content, hybridization temperature, probe coverage of the target, sequence divergence between probe and target, and the abundance of target DNA relative to background. Empirically evaluating a panel against all combinations of conditions is expensive and time-consuming.

In-silico simulation lets you:

- **Predict sensitivity** before wet-lab work: will the probes capture all intended targets?
- **Predict specificity**: will the probes accidentally capture off-target sequences (host DNA, non-target organisms)?
- **Test discrimination**: can the panel distinguish between two closely related species when only one is present in a specimen?
- **Find the limit of detection**: at what target abundance (CT value) does the panel fail?
- **Compare design strategies**: is a tiling probe set better than a set-cover design for this target?

## What BaitBench models

BaitBench simulates the following experimental steps:

1. **Fragmentation**: the DNA in a sample is broken into fragments of approximately 150–200 bp (as by sonication or enzymatic fragmentation)
2. **Probe hybridization**: probes are introduced; binding affinity is scored using nearest-neighbor thermodynamics
3. **Capture**: fragments with high-affinity probe binding sites are preferentially retained; background fragments are also captured at a lower rate
4. **Sequencing**: captured fragments are trimmed to read length and sequenced
5. **Mapping**: reads are aligned back to reference sequences to determine which organisms generated them
6. **Detection**: references with at least one mapped read are called as detected; TP/FP/FN/TN are computed

## What BaitBench does not model

BaitBench is a simulation, not a full physical model. Several real-world factors are simplified or omitted:

- **PCR amplification**: BaitBench does not model the amplification step that typically follows capture. Amplification can introduce GC bias and duplicate reads.
- **Probe concentration effects**: all probes are treated as equally abundant. In practice, probes in a pool may be unequally synthesized.
- **Secondary structure**: RNA or DNA hairpin structures in targets can block probe access. BaitBench does not model target secondary structure.
- **Probe competition**: multiple probes can compete for the same target region. BaitBench's thermodynamic model scores each probe-reference alignment independently.
- **Ligation and repair artefacts**: library preparation introduces end-repair and adaptor-ligation artefacts not modelled here.

For most panel evaluation purposes these simplifications are acceptable. The thermodynamic scoring model captures the dominant factors (GC content, mismatch tolerance, temperature) that determine whether a probe will capture its target.
