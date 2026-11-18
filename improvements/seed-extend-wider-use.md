# Reuse seed-extend elsewhere in BaitBench

`src/alignment/seed_extend.rs` was written for probe placement in `simulate`, where it
replaced a read mapper that was discarding partial matches. It is ungapped **local**
seed-and-extend with a 2-bit CSR k-mer index over the probes, parallel across reference
sequences.

On a divergent panel against *B. burgdorferi* B31 it placed 77 of 228 probes against
blastn-short's 74 and the embedded mapper's 21, at ~85x blast's speed. That tradeoff —
blast-level recall for short divergent queries at mapper-level speed — is not specific
to `simulate`. Several other call sites make the same compromise today and would
likely benefit.

## Candidate call sites

| Command | Today | Why seed-extend might fit |
|---|---|---|
| `xreact` | `--aligner minimap2\|blast` | Identical tradeoff: cross-reactivity is about short divergent homology, and blast is the documented escape hatch precisely because minimizer seeding misses it. A third option could give blast's recall without the runtime. |
| `assess-probes` | `minimap2::probe_align` | Probe coverage of targets is a local-placement problem. Gap classification (`true_gap` vs `multimapper_gap`) depends on which probes are placed, so recall changes the reported gaps. |
| `probe-coverage` | `minimap2::probe_align` | Same as above. |
| `panel-qc` / `target_similarity.rs` | `xreact_align` | Target-vs-target similarity for discriminability scoring; divergent homology is the signal being measured. |
| `probes/probetools.rs` | `probe_align` (coverage step) | Already noted in ARCHITECTURE as a substitution away from upstream ProbeTools' BLAST; seed-extend is closer to BLAST than a read mapper is. |

## Consolidation opportunity

There are now **three** seed-and-extend implementations in the tree:

- `src/alignment/seed_extend.rs` — local, ungapped, CSR index over queries, parallel
- `src/probes/catch.rs` — `build_kmer_index` + `find_matches_with_model`, Hamming
- `src/probes/syotti.rs` — `build_kmer_index` + `find_matches`, Hamming

The latter two carry a **verbatim duplicate** of `build_kmer_index` (and of
`load_fasta`). Both index the *reference* rather than the queries, which is the wrong
direction for large references — the reason the new module inverts it.

Unifying them on the shared module would remove the duplication and give catch-lite and
syotti-lite the same index. Note the semantics differ: those two need full-length
Hamming coverage decisions for set cover, whereas seed-extend reports best local
segments, so the shared piece is the index and seeding, not the acceptance rule.

## Before adopting it anywhere else

Each of these sites has its own accuracy expectations, so none should switch on
assertion alone:

- Benchmark recall and runtime against the current backend on that site's own data, as
  was done for `simulate` (`--probe-aligner` comparison on the MTEC panel).
- Check what downstream metrics move. For `assess-probes`, more placed probes means
  different gap calls; for `panel-qc`, different discriminability scores.
- `max_seed_occ` is inactive on typical probe panels but was added to bound
  low-complexity inputs; target-vs-target comparisons are a different regime and may
  need a different value.

Not scheduled. Revisit after the capture-model rework lands, since `simulate` is the
only site whose correctness currently depends on it.
