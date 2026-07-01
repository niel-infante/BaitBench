# Option A: Post-capture resampling with alignment-quality weighting

## How it works

1. Run alignment (minimap2/BLAST) as normal but **retain full alignment records** instead of just pass/fail IDs
2. Parse alignment quality for every fragment that has any hit (matching bases, mismatches, etc.)
3. Classify each fragment as target-origin or distractor-origin via `extract_source_id` + the targets/distractors ID lists
4. Assign a "binding score" to each fragment:
   - Aligned fragments: score based on matching bases / mismatches (higher = better binding)
   - Non-aligned fragments: score = small epsilon (very weak non-specific binding)
5. Calculate the target post-capture composition needed to hit the requested fold enrichment
6. Sample **without replacement** from all fragments, weighted by binding score, tuned to achieve the correct target:distractor ratio in the output
7. Output a pool of the same total size as the current capture would produce

## Pros

- Biologically nuanced: better probe matches = more likely to be retained
- Models non-specific binding (distractor fragments with partial homology get through at higher rates)
- Can achieve any fold enrichment from 1x (no enrichment) to very high
- Reuses existing alignment infrastructure
- Can bias selection toward fragments with stronger binding / fewer mismatches

## Cons

- Requires retaining full alignment data (modifying PAF/BLAST parsers to return scores, not just IDs)
- Need to pass target/distractor ID lists into the capture step (currently it doesn't know which is which)
- More complex implementation: new data structures for per-fragment scores, weighted sampling logic

## Implementation sketch

- Modify `paf::filter_paf` and `blastn::filter_blast_results` to return a `HashMap<String, f64>` of fragment_id -> binding_score instead of (or in addition to) `HashSet<String>`
- The binding score could be: `matching_bases / fragment_length` (0.0 to 1.0), penalized by mismatches
- In the enrichment step, use `rand::distributions::WeightedIndex` with these scores for selection
- Fragments with no alignment record get a minimal score (e.g. 0.001)

## Recommendation vs Option C

Option A is the natural next step from the current implementation. It adds alignment-quality weighting to the existing differential subsampling approach: instead of randomly selecting which distractors to add back (or which captured distractors to keep), weight the selection by alignment quality. This gives a biologically plausible model without the calibration complexity of Option C (probabilistic model).

Option C is worth considering if BaitBench needs to model capture as a continuous physical process rather than a discrete enrichment step. For a simulation tool where the user *specifies* the enrichment rather than *predicts* it, the added realism likely doesn't justify the complexity.
