# Option D: Probe binding competition (concentration-aware enrichment)

## Background / motivation

Observed during real-data validation against ZymoBIOMICS D6331 (21-species gut mock community, 139,144 Nanopore reads). BaitBench's current Level-1 capture sampler selects probes *uniformly* — every probe with ≥1 hit is equally likely to be chosen regardless of how abundant its target species is in the sample. The sequence abundance weight only enters at Level 2 (choosing among a probe's hit locations), and only has effect when a probe hits multiple species.

The result: rare species (M. smithii, A. muciniphila) whose 16S probes are largely specific to them see disproportionate enrichment relative to their input weight. In a real bait-capture reaction, those probes would encounter far fewer target molecules and capture far fewer fragments. BaitBench doesn't model this because it doesn't account for probe saturation or competition between species for probe binding sites.

## Sketched approaches

### D1 — Expected-yield-weighted Level-1 probe selection

Instead of uniform probe selection at Level 1, weight each probe by its *expected capture yield*: `Σ(boltzmann_score_i × seq_weight_i)` across all of its hits. This is a one-line change to the `WeightedIndex` construction in `sample_capture_fragments()`.

Weighting by `seq_weight` alone (abundance sum) is insufficient for two reasons:
- A probe with strong affinity for a rare species plus weak hits across many abundant species would be selected proportionally to hit count, then Level 2 might still choose the rare species because its Boltzmann score dominates.
- A promiscuous weakly-binding probe accumulates many low-weight hits and gets selected too often relative to a tight single-location probe. The tight probe's single high-affinity hit contributes far more expected yield but is penalized by having only one hit.

Weighting by `Σ(boltzmann × seq_weight)` fixes both: promiscuous weak probes have many small products; tight probes have one large product. Level 1 now reflects the total expected physical output of each probe, and Level 2 still correctly distributes that output across hit locations by the same per-hit scores.

Limitation: treats each probe independently. Doesn't model probe titration or competition between probes for the same target molecules.

**D1 is a ~5-minute change and would likely eliminate most of the observed over-enrichment for rare-specific probes. Try it first before considering D2 or D3.**

### D2 — Per-probe molecule count budget

Give each probe a "molecules available" count proportional to `sum(seq_weight × seq_count)` across its hits. Run a multinomial draw: each probe captures fragments in proportion to its molecule budget × its Boltzmann affinity, up to a capacity limit (probe concentration). Fragments beyond capacity are lost. This naturally suppresses rare-species over-enrichment when probes are titrated by abundant-species molecules.

More complex: requires choosing a probe concentration parameter and deciding what "titration" means when the same probe hits many species at different affinities.

### D3 — Competitive binding equilibrium

Model the capture reaction as an equilibrium: for each probe, compute the fraction of probe molecules occupied by each target species at equilibrium, using concentrations (derived from `seq_weight`) and binding affinities (Boltzmann scores). Sample fragments from the occupied-probe distribution. This is the most physically grounded approach — it's essentially a multi-ligand binding model — but requires solving a system of equations per probe and introduces new parameters (probe concentration, reaction volume).
