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
: "${PANELS:=$CROSSLAB/BB_Paper/Comp/panels}"
: "${D6331_COMPOSITION:=$CROSSLAB/BB_Paper/Comp/mock/d6331_theoretical_composition.tsv}"

# --- Harness settings ---------------------------------------------------------------
BB="${BB:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)/target/release/baitbench}"
VALIDATION_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EXPECTED_DIR="$VALIDATION_DIR/expected"
: "${VALIDATION_OUT:=${TMPDIR:-/tmp}/baitbench_validation}"

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
