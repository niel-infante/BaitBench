# Per-genome probe-placement caches

**Status: do before shipping a new version.** The current cache is correct but reuses
far less than it could.

## Now

`--probe-hits PATH` stores one cache per *experimental mix*, keyed on a content hash of
the probes and of the whole weighted reference. Simulating tick+mosquito, then tick
alone, then tick+mosquito+host is three full scans of the tick genome — and on a 2.1 GB
genome that scan is the dominant cost of the run.

## Proposed

Cache per **(probe panel, single genome)** and compose at run time. Index a tick once, a
mosquito once; any mix of them is then assembled from the parts.

Placement already supports this: `find_local_hits` scans each reference sequence
independently, so a genome's placements do not depend on what else is in the mix.
Sequence IDs also survive `prepare` unchanged, so composition is concatenation by
`seq_id` with no ID remapping.

### Why it is worth it

- N genomes give 2^N mixes but only N indexes.
- Adding one distractor costs one scan, not a re-scan of everything.
- Each genome can be indexed on a different machine — the natural way to push the big
  references onto an HPC and ship small result files back.
- Editing one distractor stops invalidating every other genome's work.
- It subsumes the weight-filtering fix below: per-genome placements are independent of
  weights, so the key is just (panel, genome).

### The trap: `max_hits_per_probe` is global

Hits are merged **per probe across all references**, ranked by identity then length,
then truncated to 1000 (`seed_extend.rs`, end of `find_local_hits`). Capping per genome
and concatenating would give up to 1000 *per genome*, so results would depend on how the
mix was sliced — unacceptable for a cache.

Exactly fixable: store the top-K per genome and **re-rank the union at compose time**.
Ranking does not depend on which genome a hit came from, and the global top-1000 can
draw at most 1000 from any one genome, so re-ranking the union of per-genome top-1000s
reproduces the global top-1000 exactly. Be deliberate about this rather than finding it
later.

### Other costs

- Weight filtering must move out of alignment into compose/scoring, so a genome is
  indexed in full even for a run that excludes it. Amortised, but a real first-time cost.
- More storage (per-genome caps can exceed the combined cap) and a cache *directory*
  convention instead of a single path — probably `--probe-cache-dir`, keyed
  panel-hash x genome-hash.
- More staleness surface. One bug has already shipped here (below), so the per-genome
  version needs tests for composition, not just for individual key fields.

## Related fix already made

The cache key hashed the *whole* reference while the aligners skip zero-weight
sequences, so two runs differing only in `--sample` shared a key despite having
different placements — the second silently reused hits missing every sequence the first
had excluded. `CacheKey::new` now hashes only the positive-weight subset
(`zero_weight_sequences_are_outside_the_key`). Per-genome caching removes the class of
bug rather than patching this instance.
