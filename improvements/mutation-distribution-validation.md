# Validation: evenly-distributed vs. naturally occurring mutation tolerance

## Background

TITAN-RNA (PMID 40501929) reported a large split in probe tolerance depending on how target mutations are distributed: ~10% tolerance to evenly/uniformly distributed mutations vs. ~27% tolerance to naturally occurring divergence. If probe capture tools (including BaitBench) behave the same way, it implies that simulation benchmarks built on uniformly-distributed mutations understate a panel's real-world breadth against naturally divergent targets — a validity concern for how BaitBench itself gets evaluated.

## Proposed experiment

1. Read TITAN-RNA's methodology (PMID 40501929) to understand exactly how they generated and measured the two mutation distributions, so the BaitBench analog is comparable.
2. Build a probe panel (`baitbench build-probes`) against a target set.
3. Construct two divergent-target variants from the same starting targets:
   - evenly/uniformly distributed mutations at a given overall divergence level
   - naturally occurring mutation patterns (e.g. real strain divergence / phylogenetically structured mutations) at the same overall divergence level
4. Run both variants through the BaitBench pipeline (`run`, or swept via `coverage-curve`) and compare capture efficiency / detection tolerance between the two.
5. Check whether BaitBench shows a similar gap (uniform mutations tolerated less well than natural divergence at the same overall % divergence).

## Open questions

- No target panel, mutation-generation method, or natural-divergence source dataset chosen yet.
- Whether to implement mutation generation as a new `baitbench` tool/utility or as an external script feeding into existing inputs (`targets.fa` variants).
- Not scheduled — revisit when pursuing BaitBench validation studies.
