# Partial-duplex capture: overlap-scaled occupancy

## The problem

A fragment is captured by a probe only if it carries the probe's target site. But a
fragment at a site boundary carries only *part* of it, and the retention model currently
credits that fragment the probe's **full** θ — which was computed from the whole 120 bp
duplex. A 20 bp sliver forms nothing stable at hybridization temperature, so capture is
overstated at every fragment boundary.

Stage 3 ships a **threshold** instead: `min_overlap_frac`, the fraction of the site a
fragment must carry before the probe pulls at all. That removes the unphysical sliver
case for one parameter, but it is still a step function — 59% of the site captures
nothing, 61% captures as if the whole site were present.

The size of the error scales with `hit_len / frag_len`. At the Borrelia LOD's 225 bp
fragments against 120 bp probes the partial-overlap zone is most of the window; at
SimVal's 5.8 kb fragments it is minor.

## The physical treatment

Recompute ΔG over the portion of the site the fragment actually carries, and take θ from
that:

```
θ_effective(s) = langmuir(ΔG(pairs ∩ fragment), [probe])
```

ΔG is not linear in duplex length — it has initiation terms and the stacking sum runs
over interior steps — so this genuinely requires re-evaluating `delta_g` on the
sub-range, not scaling.

## Why it is not a one-line change

The retention profile is piecewise constant *because* θ is constant per hit. Under
overlap-scaling, θ varies with start position throughout the partial-overlap zone, so
the profile gains a segment per base: for 28,933 hits at 5.8 kb fragments that is
~335M segments. Infeasible as written.

### Workable design

1. **Store aligned pairs per hit**, not just θ. 28,933 hits x 120 pairs x 2 bytes is
   about 7 MB — acceptable.
2. **Two truncation families.** For fragment `[s, s+L)` and hit `[h0, h1)`:
   - `s > h0` — the fragment carries the site's right portion (left-truncated hit)
   - `s + L < h1` — it carries the left portion (right-truncated hit)
   Both can only occur at once when `L < hit_len`.
3. **Quantise the overlap** into k levels (8–16) per family, precomputing θ at each.
   That is 2k floats per hit — under 2 MB at k=16 — and gives ~2k profile segments per
   hit rather than one per base. At k=16 that is ~460k segments total, which is the same
   order as the current profile.
4. The sweep in `build_retention_profile` then emits the quantised ramp on each side of
   a hit instead of a single rectangle.

### Validation

The Borrelia LOD per-locus efficiencies are the natural check, since that dataset has
short fragments relative to probe length and is where the threshold is crudest:
`tests/validation/expected/lod_locus_efficiency.tsv` (nifS 0.070 ... rplB 1.000).
Compare threshold vs overlap-scaled against those ratios.

## Relationship to other work

Shares machinery with `mismatch-stacking-thermo.md` — both are about evaluating ΔG over
a sub-range of aligned pairs rather than the whole alignment. If both are done, do them
together.

Not scheduled. Revisit once the capture model (Stage 3 of the simulation rework) has
landed and its gates are measured with the threshold in place, so the improvement can be
quantified against a working baseline.
