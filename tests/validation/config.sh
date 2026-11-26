#!/usr/bin/env bash
# Paths to external validation data.
#
# NOTHING LARGE IS VENDORED. These datasets are referenced in place (the Ixodes genome
# alone is 2.25 GB, the D6331 references 93 Mb, the mock FASTQs 2.5 GB). Only small
# derived reference tables live in the repo, under expected/.
#
# Every path can be overridden by exporting the matching variable first, so the harness
# works on another machine without editing this file:
#
#   CROSSLAB=/data/crosslab bash tests/validation/run_fast.sh
#
# `check_data.sh` reports exactly which of these are missing.

# Root holding the sibling project directories.
: "${CROSSLAB:=$HOME/Projects/CrossLab}"

# --- Borrelia burgdorferi limit-of-detection titration -----------------------------
: "${LOD_DIR:=$CROSSLAB/MTEC_probes/borrelia_LOD}"
: "${LOD_COUNTS:=$LOD_DIR/LOD_4Path/all_counts.tsv}"
: "${LOD_COVERAGE:=$LOD_DIR/LOD_4Path/coverage}"
: "${LOD_TARGETS_BED:=$LOD_DIR/4Path/target_analysis/targets_merged.bed}"
: "${LOD_PROBES_BLAST:=$LOD_DIR/4Path/target_analysis/probes_blast6.tsv}"
: "${MTEC_PROBES:=$CROSSLAB/MTEC_probes/Targets_and_Probes/bacteria_probes.fa}"
: "${MTEC_TARGETS:=$CROSSLAB/MTEC_probes/Targets_and_Probes/bacteria_targets.fa}"
: "${REF_4PATH:=$CROSSLAB/MTEC_probes/Data/Refseq/4Path.fa}"
: "${REF_ISCAP:=$CROSSLAB/MTEC_probes/Data/Refseq/Iscapularis.fa}"
: "${LOD_SIM_INPUTS:=$CROSSLAB/MTEC_runs/BB_LOD_deep_dive/inputs}"

# --- ZymoBIOMICS D6331 mock community (Slizovskiy et al. 2024) ---------------------
: "${SIMVAL_DIR:=$CROSSLAB/BB_Paper/SimVal}"
: "${SIMVAL_GENOMES:=$SIMVAL_DIR/inputs/genomes.fa}"
: "${SIMVAL_GROUPS:=$SIMVAL_DIR/inputs/groups.tsv}"
: "${SIMVAL_SAMPLE:=$SIMVAL_DIR/inputs/sample.tsv}"
: "${SIMVAL_REAL:=$SIMVAL_DIR/real}"
# Raw CCS reads for the 24 mock-community runs, named <sample>_<accession>.fastq.gz.
# These already exist here; tests/validation/simval/fetch_reads.sh only re-fetches
# what is missing.
: "${SIMVAL_READS:=$SIMVAL_DIR/../Comp/sra/fastq}"
: "${PANELS:=$CROSSLAB/BB_Paper/Comp/panels}"
: "${D6331_COMPOSITION:=$CROSSLAB/BB_Paper/Comp/mock/d6331_theoretical_composition.tsv}"

# --- Harness settings ---------------------------------------------------------------
BB="${BB:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)/target/release/baitbench}"
VALIDATION_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EXPECTED_DIR="$VALIDATION_DIR/expected"
: "${VALIDATION_OUT:=${TMPDIR:-/tmp}/baitbench_validation}"

# Probe placement dominates runtime and scales close to linearly, but baitbench
# defaults to --threads 1. Without this the full gate runs ~6x slower than the machine
# allows (measured 215s -> 35s at 10 threads on the MTEC panel).
_ncpu=$( (command -v nproc >/dev/null && nproc) || sysctl -n hw.ncpu 2>/dev/null || echo 1 )
# Leave 2 cores free: these gates run for hours and the machine is in use meanwhile.
# Override by exporting THREADS explicitly.
: "${THREADS:=$(( _ncpu > 3 ? _ncpu - 2 : 1 ))}"
unset _ncpu
export THREADS

# Probe placements depend only on the panel and the reference, so all 11 LOD conditions
# and all 3 mock replicates share one each. These live outside VALIDATION_OUT so that
# wiping results does not throw the expensive part away.
: "${PROBE_CACHE_DIR:=${TMPDIR:-/tmp}/baitbench_probe_cache}"
mkdir -p "$PROBE_CACHE_DIR" 2>/dev/null || true
export PROBE_CACHE_DIR

# Model parameters. Both of these are EMPTY by default so the gates exercise the
# binary's own defaults (melting window 10, [probe] = e^-2) rather than pinning a
# configuration here that drifts out of step with what ships. Set either one to sweep
# it; the gate scripts pass the flag only when the variable is non-empty.
#
# This was previously `MELTING_WINDOW:=0`, which silently forced the whole-duplex
# occupancy on every gate run — not the shipped default, and the worst of the six
# values measured against the enrichment target (parameterisation.md A1).
: "${MELTING_WINDOW:=}"
export MELTING_WINDOW
: "${PROBE_CONC:=}"
export PROBE_CONC

# Capture efficiency for the mock-community gate. The flag is a scenario multiplier, and
# the enrichment objective provably cannot pin it — Spearman rho is identical to four
# decimals across a 10x range (parameterisation.md D) — so the gate passes the shipped
# default rather than a fitted value. It was 0.207, calibrated under the pre-rework
# meaning (fraction of fragments drawn from probe sites), which no longer exists.
: "${SIMVAL_CAPTURE_FRACTION:=0.5}"
export SIMVAL_CAPTURE_FRACTION

# Fast mode keeps the inner loop to minutes; full mode is the release gate.
if [[ "${VALIDATION_MODE:-fast}" == "fast" ]]; then
    LOD_CONDITIONS="NTC s1E3 Pos1"
    SIMVAL_RUNS="MOXTAA"
    NUM_FRAGMENTS=20000
else
    LOD_CONDITIONS="NTC s1E0 s1E001 s1E1 s1E2 s1E3 s1E4 s1E5 s1E6 Pos1 Pos2"
    SIMVAL_RUNS="MOXTAA MOXTAB MOXTAC"
    NUM_FRAGMENTS=150000
fi
export LOD_CONDITIONS SIMVAL_RUNS NUM_FRAGMENTS

# compare.py and locus_efficiency.py read these from the environment, so they must be
# exported, not just set: a child process does not inherit plain shell variables.
export CROSSLAB LOD_DIR LOD_COUNTS LOD_COVERAGE LOD_TARGETS_BED LOD_PROBES_BLAST
export MTEC_PROBES MTEC_TARGETS REF_4PATH REF_ISCAP LOD_SIM_INPUTS
export SIMVAL_DIR SIMVAL_GENOMES SIMVAL_GROUPS SIMVAL_SAMPLE SIMVAL_REAL SIMVAL_READS
export PANELS D6331_COMPOSITION
export BB VALIDATION_DIR EXPECTED_DIR VALIDATION_OUT
