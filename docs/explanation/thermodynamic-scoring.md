# Thermodynamic Scoring

BaitBench's thermodynamic simulation mode uses the SantaLucia (1998) nearest-neighbor (NN) model to score probe-target hybridization and weight the probability of fragment capture. This page explains the model and how it translates to capture efficiency predictions.

---

## Why Thermodynamics?

A probe captures its target if the probe-target duplex is stable enough under the hybridization conditions. Stability is determined by:

1. **Sequence composition**: GC base pairs form three hydrogen bonds (stronger) vs AT/AU two bonds (weaker)
2. **Context effects**: the stability of each base pair depends on its neighbours (nearest-neighbor effects)
3. **Temperature**: higher temperature destabilises duplexes; mismatches are more disruptive at high temperature
4. **Salt concentration**: Na⁺ ions screen the negative phosphate backbone charges, stabilising duplexes

A simple GC-content model misses all the context effects. The SantaLucia (1998) model captures them by summing stacking energies for each adjacent base-pair step — the dominant contribution to duplex stability.

---

## The SantaLucia Nearest-Neighbor Model

The nearest-neighbor model expresses ΔG of duplex formation as:

```
ΔG°(T) = ΔH° - T × ΔS°
```

Where ΔH° and ΔS° are the enthalpy and entropy of hybridization, each computed by summing over all adjacent dinucleotide steps in the sequence:

```
ΔH° = Σ ΔH°(stacking) + ΔH°(initiation)
ΔS° = Σ ΔS°(stacking) + ΔS°(initiation)
```

SantaLucia (1998) provides a unified table of 10 unique nearest-neighbor parameters (the 16 possible dinucleotide steps reduce to 10 by symmetry). BaitBench uses these published values directly.

---

## Computing ΔG: Three Contributions

### 1. Stacking terms

For each consecutive base-pair dinucleotide in the probe-target alignment (5'→3'), BaitBench looks up the stacking ΔH° and ΔS° from the SantaLucia (1998) table. Mismatched positions contribute zero stacking energy (they break the consecutive run).

This means a probe with 3 mismatches in a run of otherwise matched positions loses the stacking contributions of all bases adjacent to those mismatches, not just the mismatches themselves. Clustered mismatches are disproportionately destabilising.

### 2. Initiation terms

Two initiation penalties are added:
- One for the terminal AT (or AU) base pair at the 5' end (if applicable)
- One for the terminal AT base pair at the 3' end (if applicable)

GC terminal pairs incur no initiation penalty; AT terminal pairs add a small enthalpic and entropic cost (weaker terminal stacking).

### 3. Salt correction (Owczarzy et al.)

The ΔG calculation applies at 1 M NaCl. For realistic hybridization conditions (typically 0.1–1.0 M), the melting temperature is corrected using the Owczarzy et al. formula:

```
1/Tm(corrected) = 1/Tm(1M NaCl)
                + (4.29 × fGC - 3.95) × 10⁻⁵ × ln([Na⁺])
                + 9.40 × 10⁻⁶ × (ln([Na⁺]))²
```

Where `fGC` is the fraction of GC base pairs in the duplex. This correction shifts the melting temperature upward at physiological Na⁺ concentrations compared to the standard 1 M condition.

BaitBench applies this correction to the melting temperature, then uses the corrected Tm to evaluate ΔG at the actual hybridization temperature:

```
ΔG(T_hyb) = ΔH° × (1 - T_hyb / Tm_corrected)
```

---

## Boltzmann Weighting

After computing ΔG for each probe alignment at each position, BaitBench converts it to a capture score:

```
score = exp(-ΔG / (R × T))
```

Where:
- R = 1.987 cal/(mol·K) (gas constant)
- T = hybridization temperature in Kelvin
- ΔG is in cal/mol

A favourable (negative) ΔG gives a score > 1; an unfavourable ΔG gives a score < 1. Positions with no probe alignment have score 0.

These scores are used as sampling weights: positions with higher Boltzmann scores generate more captured fragments per unit of simulation time.

---

## From Scores to Fragment Sampling

The simulation divides fragments into two categories:

**Captured fragments** (governed by `--capture-fraction`): drawn from positions with probe alignments, weighted by Boltzmann score. A position covered by two overlapping probes accumulates the scores of both alignments — higher effective capture probability.

**Background fragments** (governed by `1 - capture-fraction`): drawn uniformly from all positions in a sequence, regardless of probe coverage. These represent non-specific co-capture: fragments that end up in the capture library without being bound by a probe (e.g., because they were in solution adjacent to a captured molecule, or due to non-specific binding).

At `--capture-fraction 0.5` (the default), 50% of a sequence's fragments come from probe-covered positions (thermodynamically weighted) and 50% come from uniform background sampling. Background reads spread thinly across the entire sequence, producing low depth per position. Probe-captured reads concentrate at binding sites, producing high local depth — this is what creates the characteristic "spike" pattern in captured vs non-captured sequences.

---

## Simple Mode vs Thermodynamic Mode

In `--simulate-mode simple`, the Boltzmann weighting step is skipped. All alignments above the minimum match threshold are treated as equally likely to capture. This removes temperature and sequence composition effects.

| Feature | Thermodynamic | Simple |
|---------|---------------|--------|
| Mismatch penalty | Yes — via ΔG | No |
| Temperature effect | Yes — `--hybridization-temperature` affects capture | No |
| GC content effect | Yes — high-GC probes capture more efficiently | No |
| Probe-position weighting | Boltzmann-weighted | Uniform over covered positions |
| Use case | Realistic efficiency prediction | Coverage testing, speed |

Thermodynamic mode is more accurate; simple mode is faster and useful for understanding which positions probes cover at all, independent of efficiency.
