# Possible Enhancements to Capture Enrichment

The current fold enrichment implementation (Option B) uses differential subsampling: run binary capture, then adjust the target:distractor ratio by subsampling captured distractors or adding back uncaptured distractors. This is simple and deterministic but doesn't model alignment-quality-dependent binding strength.

Below are two more sophisticated approaches that could be layered on in the future.

---

## Option A: Post-capture resampling with alignment-quality weighting

### How it works

1. Run alignment (minimap2/BLAST) as normal but **retain full alignment records** instead of just pass/fail IDs
2. Parse alignment quality for every fragment that has any hit (matching bases, mismatches, etc.)
3. Classify each fragment as target-origin or distractor-origin via `extract_source_id` + the targets/distractors ID lists
4. Assign a "binding score" to each fragment:
   - Aligned fragments: score based on matching bases / mismatches (higher = better binding)
   - Non-aligned fragments: score = small epsilon (very weak non-specific binding)
5. Calculate the target post-capture composition needed to hit the requested fold enrichment
6. Sample **without replacement** from all fragments, weighted by binding score, tuned to achieve the correct target:distractor ratio in the output
7. Output a pool of the same total size as the current capture would produce

### Pros

- Biologically nuanced: better probe matches = more likely to be retained
- Models non-specific binding (distractor fragments with partial homology get through at higher rates)
- Can achieve any fold enrichment from 1x (no enrichment) to very high
- Reuses existing alignment infrastructure
- Can bias selection toward fragments with stronger binding / fewer mismatches

### Cons

- Requires retaining full alignment data (modifying PAF/BLAST parsers to return scores, not just IDs)
- Need to pass target/distractor ID lists into the capture step (currently it doesn't know which is which)
- More complex implementation: new data structures for per-fragment scores, weighted sampling logic

### Implementation sketch

- Modify `paf::filter_paf` and `blastn::filter_blast_results` to return a `HashMap<String, f64>` of fragment_id -> binding_score instead of (or in addition to) `HashSet<String>`
- The binding score could be: `matching_bases / fragment_length` (0.0 to 1.0), penalized by mismatches
- In the enrichment step, use `rand::distributions::WeightedIndex` with these scores for selection
- Fragments with no alignment record get a minimal score (e.g. 0.001)

---

## Option C: Score-based probabilistic capture (full model)

### How it works

1. Run alignment on ALL fragments but **don't filter** — keep all alignment records with scores
2. For each fragment, compute a normalized binding score (0.0 to 1.0) from alignment quality
3. Apply a capture probability function parameterized by fold enrichment:
   - `P(capture | target, score) = sigmoid(a * score + b)`
   - `P(capture | distractor, score) = sigmoid(a * score + b - c)`
   - Where `c` is derived from the desired fold enrichment
4. For each fragment, draw a random number; if < probability, include in output
5. Calibrate parameters to match the target fold enrichment

### Pros

- Most biologically realistic — capture probability is a continuous function of binding affinity
- No duplicates, natural output size
- Fold enrichment emerges from the model rather than being forced
- Could be extended with additional parameters (wash stringency, probe concentration) in the future

### Cons

- Most complex to implement and test
- Stochastic — exact fold enrichment varies run-to-run
- Calibrating the sigmoid parameters requires iterative search (binary search on the offset parameter `c`) — must be done each run because the alignment score distribution depends on the specific probes, targets, and distractors
- Fragments with zero alignment (no PAF/BLAST record at all) need special handling
- Hardest to explain to users

### Calibration details

The sigmoid parameters can't be pre-computed because score distributions change every run (depending on probes, targets, distractors, fragment lengths/positions, and number of alignments). Each run requires:

1. Collect all alignment scores
2. Pick initial sigmoid parameters
3. Compute the *expected* fold enrichment under those parameters
4. Adjust via binary search on the offset parameter until expected FE matches target

This is a one-dimensional root-finding problem (~10-20 iterations), so it's not expensive computationally, but adds implementation complexity. The *realized* FE from random draws will also differ from the *expected* FE.

---

---

## Option D: Probe binding competition (concentration-aware enrichment)

### Background / motivation

Observed during real-data validation against ZymoBIOMICS D6331 (21-species gut mock community, 139,144 Nanopore reads). BaitBench's current Level-1 capture sampler selects probes *uniformly* — every probe with ≥1 hit is equally likely to be chosen regardless of how abundant its target species is in the sample. The sequence abundance weight only enters at Level 2 (choosing among a probe's hit locations), and only has effect when a probe hits multiple species.

The result: rare species (M. smithii, A. muciniphila) whose 16S probes are largely specific to them see disproportionate enrichment relative to their input weight. In a real bait-capture reaction, those probes would encounter far fewer target molecules and capture far fewer fragments. BaitBench doesn't model this because it doesn't account for probe saturation or competition between species for probe binding sites.

### Sketched approaches

**D1 — Expected-yield-weighted Level-1 probe selection**

Instead of uniform probe selection at Level 1, weight each probe by its *expected capture yield*: `Σ(boltzmann_score_i × seq_weight_i)` across all of its hits. This is a one-line change to the `WeightedIndex` construction in `sample_capture_fragments()`.

Weighting by `seq_weight` alone (abundance sum) is insufficient for two reasons:
- A probe with strong affinity for a rare species plus weak hits across many abundant species would be selected proportionally to hit count, then Level 2 might still choose the rare species because its Boltzmann score dominates.
- A promiscuous weakly-binding probe accumulates many low-weight hits and gets selected too often relative to a tight single-location probe. The tight probe's single high-affinity hit contributes far more expected yield but is penalized by having only one hit.

Weighting by `Σ(boltzmann × seq_weight)` fixes both: promiscuous weak probes have many small products; tight probes have one large product. Level 1 now reflects the total expected physical output of each probe, and Level 2 still correctly distributes that output across hit locations by the same per-hit scores.

Limitation: treats each probe independently. Doesn't model probe titration or competition between probes for the same target molecules.

**D2 — Per-probe molecule count budget**

Give each probe a "molecules available" count proportional to `sum(seq_weight × seq_count)` across its hits. Run a multinomial draw: each probe captures fragments in proportion to its molecule budget × its Boltzmann affinity, up to a capacity limit (probe concentration). Fragments beyond capacity are lost. This naturally suppresses rare-species over-enrichment when probes are titrated by abundant-species molecules.

More complex: requires choosing a probe concentration parameter and deciding what "titration" means when the same probe hits many species at different affinities.

**D3 — Competitive binding equilibrium**

Model the capture reaction as an equilibrium: for each probe, compute the fraction of probe molecules occupied by each target species at equilibrium, using concentrations (derived from `seq_weight`) and binding affinities (Boltzmann scores). Sample fragments from the occupied-probe distribution. This is the most physically grounded approach — it's essentially a multi-ligand binding model — but requires solving a system of equations per probe and introduces new parameters (probe concentration, reaction volume).

### Recommendation

D1 is a five-minute change and would likely eliminate most of the observed over-enrichment for rare-specific probes. Try it first before considering D2 or D3.

---

## Recommendation

Option A is the natural next step from the current Option B implementation. It adds alignment-quality weighting to the existing differential subsampling approach: instead of randomly selecting which distractors to add back (or which captured distractors to keep), weight the selection by alignment quality. This gives a biologically plausible model without the calibration complexity of Option C.

Option C is worth considering if BaitBench needs to model capture as a continuous physical process rather than a discrete enrichment step. For a simulation tool where the user *specifies* the enrichment rather than *predicts* it, the added realism likely doesn't justify the complexity.

---

## Option E: RNA-aware nearest-neighbor thermodynamics

The current TNN model (`src/thermodynamics.rs`) uses SantaLucia (1998) DNA-DNA parameters. Two well-characterized parameter sets exist for RNA chemistry:

- **RNA-RNA**: Xia et al. (1998) *Biochemistry* 37(47):14719–14735 — same 10-parameter NN framework, U substitutes T, generally more stable than DNA-DNA.
- **RNA:DNA hybrid** (DNA probe + RNA target, or RNA probe + DNA target): Sugimoto et al. (1995) *Biochemistry* 34(35):11211–11216 — 8 asymmetric dinucleotide steps; stability intermediate between DNA-DNA and RNA-RNA. This is the practically relevant case for capture probes hybridizing to RNA virus targets or transcripts.

Note: probe and target sequences are conventionally supplied with T even when RNA chemistry is intended, so no U-handling changes are needed in sequence parsing.

### Implementation sketch

- Add `DuplexChemistry` enum (`DnaDna`, `RnaRna`, `RnaDna`) to `ThermoModel`
- Add Xia and Sugimoto NN tables alongside the existing SantaLucia table
- Dispatch in `delta_g()` based on the chemistry field
- The salt correction (Owczarzy 1997) applies to DNA-DNA; RNA duplexes use Nakano et al. (1999) — this would need updating for RNA modes
- Expose via `--duplex-chemistry` CLI flag
