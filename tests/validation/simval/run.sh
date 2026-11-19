#!/usr/bin/env bash
# D6331 mock-community composition gate.
#
# Predicts per-organism read composition from the manufacturer's theoretical
# composition plus the reconstructed CATCH panel, and compares against the real
# PacBio CCS reads of Slizovskiy et al. 2024 (PRJNA1081843).
#
# Pass condition (Stage 3): E. coli (5 strains) rises toward its real ~73.6% while
# S. enterica stays near 0.02% and E. faecalis near 0.001%. Pre-rework the model
# under-predicted E. coli 3x and inflated the two rare organisms 135x / 1360x.
#
# Read model is `perfect`: the real data are CCS/HiFi at ~99% accuracy, and badread's
# pacbio profile models CLR (~85-90%), which is further from the data than error-free
# reads are.
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/../config.sh"

mkdir -p "$VALIDATION_OUT"
status=0

for run in $SIMVAL_RUNS; do
    real="$SIMVAL_REAL/${run}_counts.tsv"
    [[ -s "$real" ]] || { echo "missing real counts: $real" >&2; exit 1; }

    tag="${run}_res"
    outdir="$VALIDATION_OUT/$tag"
    if [[ -d "$outdir" ]]; then
        echo "[$(date +%H:%M:%S)] $tag: reusing $outdir"
    else
        echo "[$(date +%H:%M:%S)] $tag: simulating $NUM_FRAGMENTS fragments"
        "$BB" run \
            --targets "$SIMVAL_GENOMES" \
            --probes  "$PANELS/res_probes.fasta" \
            --sample  "$SIMVAL_SAMPLE" \
            --groups  "$SIMVAL_GROUPS" \
            --num-fragments "$NUM_FRAGMENTS" \
            --fragment-length-mean 5800 \
            --fragment-length-min 1000 \
            --fragment-length-max 20000 \
            --capture-fraction 0.207 \
            --read-simulator perfect \
            --read-length 20000 \
            --hybridization-temperature 65 \
            --salt-concentration 50 \
            --threads "$THREADS" \
            --probe-hits "$PROBE_CACHE_DIR/simval_res_panel.tsv" \
            --seed 42 \
            --report none \
            --outdir "$outdir" \
            > "$VALIDATION_OUT/${tag}.log" 2>&1 \
            || { echo "  simulate FAILED, see $VALIDATION_OUT/${tag}.log" >&2; status=1; continue; }
    fi

    echo
    VALIDATION_OUT="$VALIDATION_OUT" SIMVAL_REAL="$SIMVAL_REAL" SIMVAL_GROUPS="$SIMVAL_GROUPS" \
        python3 "$VALIDATION_DIR/simval/compare.py" "$tag" "$run" || status=1
    echo
done

echo "Reference composition: $EXPECTED_DIR/simval_real_composition.tsv"
exit $status
