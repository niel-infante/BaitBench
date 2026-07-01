# Option C: Score-based probabilistic capture (full model)

## How it works

1. Run alignment on ALL fragments but **don't filter** — keep all alignment records with scores
2. For each fragment, compute a normalized binding score (0.0 to 1.0) from alignment quality
3. Apply a capture probability function parameterized by fold enrichment:
   - `P(capture | target, score) = sigmoid(a * score + b)`
   - `P(capture | distractor, score) = sigmoid(a * score + b - c)`
   - Where `c` is derived from the desired fold enrichment
4. For each fragment, draw a random number; if < probability, include in output
5. Calibrate parameters to match the target fold enrichment

## Pros

- Most biologically realistic — capture probability is a continuous function of binding affinity
- No duplicates, natural output size
- Fold enrichment emerges from the model rather than being forced
- Could be extended with additional parameters (wash stringency, probe concentration) in the future

## Cons

- Most complex to implement and test
- Stochastic — exact fold enrichment varies run-to-run
- Calibrating the sigmoid parameters requires iterative search (binary search on the offset parameter `c`) — must be done each run because the alignment score distribution depends on the specific probes, targets, and distractors
- Fragments with zero alignment (no PAF/BLAST record at all) need special handling
- Hardest to explain to users

## Calibration details

The sigmoid parameters can't be pre-computed because score distributions change every run (depending on probes, targets, distractors, fragment lengths/positions, and number of alignments). Each run requires:

1. Collect all alignment scores
2. Pick initial sigmoid parameters
3. Compute the *expected* fold enrichment under those parameters
4. Adjust via binary search on the offset parameter until expected FE matches target

This is a one-dimensional root-finding problem (~10-20 iterations), so it's not expensive computationally, but adds implementation complexity. The *realized* FE from random draws will also differ from the *expected* FE.
