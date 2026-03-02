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

## Recommendation

Option A is the natural next step from the current Option B implementation. It adds alignment-quality weighting to the existing differential subsampling approach: instead of randomly selecting which distractors to add back (or which captured distractors to keep), weight the selection by alignment quality. This gives a biologically plausible model without the calibration complexity of Option C.

Option C is worth considering if BaitBench needs to model capture as a continuous physical process rather than a discrete enrichment step. For a simulation tool where the user *specifies* the enrichment rather than *predicts* it, the added realism likely doesn't justify the complexity.
