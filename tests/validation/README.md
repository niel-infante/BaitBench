# Validation harness

Three real datasets that gate the simulation rework (see the plan's Stage 0 and the
per-stage gates). Two of them currently **fail by design** — they are the regressions
the rework must turn green.

## Nothing large is vendored

Every dataset is referenced **in place** by path. The inputs total several gigabytes
(the *Ixodes* distractor genome alone is 2.25 GB, the D6331 references 93 Mb, the mock
FASTQs 2.5 GB), so copying them into the repo is not an option.

Only small derived reference tables live here, in `expected/` (~4 KB total).

Paths are declared in `config.sh` and every one can be overridden by exporting it:

```bash
bash tests/validation/check_data.sh                      # what is missing, if anything
CROSSLAB=/data/crosslab bash tests/validation/run_fast.sh # relocate the whole tree
LOD_COUNTS=/elsewhere/all_counts.tsv bash tests/validation/run_fast.sh  # one path
```

Symlinking the sibling projects into place works equally well.

## Running

```bash
bash tests/validation/run_fast.sh   # minutes: 3 LOD conditions, 1 mock run, 20k fragments
bash tests/validation/run_full.sh   # hours:  11 LOD conditions, 3 mock runs, 150k fragments
```

Use fast for the inner loop and full before declaring a stage complete. Results are
cached under `$VALIDATION_OUT` (default `$TMPDIR/baitbench_validation`); delete a
subdirectory to force a re-run.

## The three gates

### 1. `--distractor-fraction` monotonicity (`cargo test --test simulation -- --ignored`)

Runs on `examples/minimal`, needs no external data. Requesting more distractor must
deliver more distractor across the whole range. Pre-rework the realized fraction
saturates at `1 - capture_fraction`:

```
requested 0.100 -> realized 0.043
requested 0.500 -> realized 0.242
requested 0.900 -> realized 0.453
requested 0.990 -> realized 0.495
```

This is reviewer item 6. Fixed in Stage 3.

### 2. Borrelia LOD titration (`lod/run.sh`)

An 11-condition *B. burgdorferi* dilution into real tick matrix, spanning ~20,000-fold
(NTC = 37 BB reads, Pos1 = 748,104). Pre-rework the simulated BB fraction is a flat
**78.50%** for every condition from NTC to s1E5, because realized composition was
pinned by `--capture-fraction` (0.8183) and blind to abundance.

A second, independent gate on the same data: per-locus relative capture efficiency
(`expected/lod_locus_efficiency.tsv`, nifS 0.070 … rplB 1.000) must be reproduced *and*
stay constant across the dilution series, as it does in reality. The weak loci are weak
because the panel contains no *B. burgdorferi* probe for them, so their best probe has
6–10 mismatches — which is why Stages 2 and 2.5 (alignment sensitivity, internal
mismatch thermodynamics) come first. Background: `BB_Paper/LODLoci/FINDINGS.md`.

Per-locus comparison is not wired up yet; the titration column is the Stage 3 gate.

### 3. D6331 mock community composition (`simval/run.sh`)

Per-organism read composition against the real PacBio CCS reads of Slizovskiy et al.
2024 (PRJNA1081843), predicted from the **manufacturer's theoretical composition** —
not from the reads being predicted, which is what made the previous harness circular
(reviewer item 3).

Pre-rework the model under-predicts *E. coli* 3× (real 73.6%) while inflating
*S. enterica* 135× and *E. faecalis* 1360×, and agreement is *best with the capture
model switched off* (CCC 0.795 at cf=0 vs 0.362 at the physically-anchored cf=0.207).

Read model is `perfect`: the real reads are CCS/HiFi at mean Q20 (~99% accuracy), and
badread's `pacbio` profile models CLR (~85–90%), which is further from the data than
error-free reads are.

## Reference tables in `expected/`

| File | What it is |
|---|---|
| `lod_real_counts.tsv` | Real BB-mapped reads per condition, from `LOD_4Path/all_counts.tsv` |
| `lod_locus_efficiency.tsv` | Per-locus efficiency relative to rplB, from mosdepth coverage |
| `simval_real_composition.tsv` | Real per-organism proportions, minimap2 `map-hifi` primary alignments |

Regenerate them only if the underlying real data changes — they are measurements, not
expectations to be tuned.

## Provenance of the external data

| Variable | What it is | How it was produced |
|---|---|---|
| `LOD_*`, `REF_4PATH`, `REF_ISCAP`, `MTEC_PROBES` | Borrelia LOD run and the MTEC bacterial panel | In-house capture sequencing; `MTEC_probes/README.md` |
| `LOD_SIM_INPUTS` | Per-condition calibration for the LOD simulation | `MTEC_runs/BB_LOD_deep_dive/scripts/` |
| `SIMVAL_*` | D6331 genomes, groups, weights, real mapped counts | `BB_Paper/SimVal/scripts/` |
| `PANELS` | RES/MOB/Combo panels reconstructed from the published CATCH recipe | `BB_Paper/Comp/` — reproduces 12 published numbers exactly |
| `D6331_COMPOSITION` | Manufacturer's theoretical composition, 21 strains | Zymo D6331 instruction manual, Table 1 |
