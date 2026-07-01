# Mismatch Stacking Thermodynamics

## Problem

The current `delta_g()` function in `src/thermodynamics.rs` uses "SkipStacking": whenever either end of a nearest-neighbor dinucleotide step is a mismatch, the step contributes **zero** stacking energy. This overestimates the thermodynamic penalty for probes with scattered single mismatches — treating a broken stack rather than a weakened one. The practical effect: probes with a few scattered mismatches get the same Boltzmann score as a completely non-binding probe for those positions, underestimating their real capture affinity.

## Difficulty: Moderate (~4–6 hours)

- Code change: Small and isolated to `src/thermodynamics.rs` (~40 lines)
- Data entry: Significant — mismatch stacking table (~96 entries from literature)
- Tests: 2–3 tests need updating; ~5 new tests needed
- Docs: 3 files need minor updates

## Approach

### Logic change in `delta_g()`

The stacking loop (lines 186–211) currently has two cases:
- `prev_wc && curr_wc` → use `NN_TABLE` (unchanged)
- anything else → skip (zero energy)

Add a third case:
- `prev_wc && !curr_wc` OR `!prev_wc && curr_wc` → look up `MM_TABLE`
- `!prev_wc && !curr_wc` → still skip zero (no published parameters for consecutive mismatches; rare in real probe alignments)

The `has_stacking` flag continues to require at least one WC|WC step before applying initiation terms and salt correction. Mismatch stacking alone doesn't nucleate a duplex, so this guard stays.

### New parameter table

```rust
// Indexed by [prev_top][prev_bot][curr_top][curr_bot], all 0..=3.
// Only WC|mismatch and mismatch|WC entries are populated (Some).
static MM_TABLE: [[[[Option<NnParams>; 4]; 4]; 4]; 4] = ...;
```

A 4^4 = 256-entry table as `Option<NnParams>`. Most entries are `None`.

### Parameter data source

**SantaLucia & Hicks (2004)** "The Thermodynamics of DNA Structural Motifs", *Annu. Rev. Biophys. Biomol. Struct.* 33:415–440, Table 4 (internal single mismatch parameters).

Cross-reference with **biopython's** `Bio/SeqUtils/MeltingTemp.py` (`DNA_IMM1` dictionary) for machine-readable values and to catch transcription errors.

The 12 mismatch types (top/bottom): `A/A`, `A/C`, `A/G`, `C/A`, `C/C`, `C/T`, `G/A`, `G/G`, `G/T`, `T/C`, `T/G`, `T/T`. Each in ~8 flanking contexts (RC symmetry halves independent parameters) → ~96 entries total, ~48 independent.

### Files to change

**`src/thermodynamics.rs`** only:
1. Add `MM_TABLE` constant (~100 lines of data)
2. Modify the stacking loop in `delta_g()` to add the mismatch stacking branch (~15 lines)
3. Update module doc comment (cite SantaLucia & Hicks 2004, remove SkipStacking framing)
4. Update tests: revise `all_mismatches_gives_zero_energy` and `isolated_wc_pairs_no_stacking_score_neutral_at_low_salt`; add new mismatch-stacking tests

**Documentation** (minor):
- `ARCHITECTURE.md` — update SkipStacking description
- `docs/concepts.md` (lines 62–82) — update user-facing mismatch handling description

### What stays the same

- `boltzmann_score()`, `ThermoModel`, `aligned_pairs` format — all unchanged
- `has_stacking` gating for initiation and salt correction — unchanged
- `thermo_sim.rs:265` call site — unchanged; gets more accurate scores automatically

## Verification

```bash
cargo test --lib thermodynamics
cargo build --release
./target/release/baitbench run \
  --targets examples/minimal/targets.fa \
  --distractors examples/minimal/distractors.fa \
  --probes examples/minimal/probes.fa \
  --num-fragments 1000 --capture-fraction 0.5 \
  --seed 42 --report none --outdir test_results_mm
cat test_results_mm/*/results.tsv
```

Sanity check: a probe with one central mismatch should score strictly between a perfect-match probe and a fully-mismatched probe (rather than equal to the latter as before).
