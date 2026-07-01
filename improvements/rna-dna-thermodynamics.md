# Option E: RNA-aware nearest-neighbor thermodynamics

## Overview

The current TNN model (`src/thermodynamics.rs`) uses SantaLucia (1998) DNA-DNA parameters. Two well-characterized parameter sets exist for RNA chemistry:

- **RNA-RNA**: Xia et al. (1998) *Biochemistry* 37(47):14719–14735 — same 10-parameter NN framework, U substitutes T, generally more stable than DNA-DNA.
- **RNA:DNA hybrid** (DNA probe + RNA target, or RNA probe + DNA target): Sugimoto et al. (1995) *Biochemistry* 34(35):11211–11216 — 8 asymmetric dinucleotide steps; stability intermediate between DNA-DNA and RNA-RNA. This is the practically relevant case for capture probes hybridizing to RNA virus targets or transcripts.

Note: probe and target sequences are conventionally supplied with T even when RNA chemistry is intended, so no U-handling changes are needed in sequence parsing.

## Implementation sketch

- Add `DuplexChemistry` enum (`DnaDna`, `RnaRna`, `RnaDna`) to `ThermoModel`
- Add Xia and Sugimoto NN tables alongside the existing SantaLucia table
- Dispatch in `delta_g()` based on the chemistry field
- The salt correction (Owczarzy 1997) applies to DNA-DNA; RNA duplexes use Nakano et al. (1999) — this would need updating for RNA modes
- Expose via `--duplex-chemistry` CLI flag
