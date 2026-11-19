#!/usr/bin/env bash
# Borrelia burgdorferi limit-of-detection gate.
#
# Two things must hold, and they are independent:
#
#   1. TITRATION. The simulated BB read fraction must track the real ~20,000-fold
#      range (NTC = 37 reads -> Pos1 = 748,104). The pre-rework model produced a flat
#      78.50% for every condition from NTC to s1E5, because realized composition was
#      pinned by --capture-fraction and blind to abundance.
#
#   2. PER-LOCUS EFFICIENCY. Relative capture efficiency per locus (nifS 0.070 ...
#      rplB 1.000) must be reproduced AND stay constant across the dilution series,
#      as it does in the real data. The weak loci are weak because the panel carries
#      no B. burgdorferi probe for them, so their best probe has 6-10 mismatches.
#      See BB_Paper/LODLoci/FINDINGS.md.
#
# Inputs are referenced in place. The per-condition calibration (num-fragments,
# distractor-fraction) is reused from the existing BB_LOD_deep_dive run rather than
# re-derived, so this harness measures the model, not the setup.
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/../config.sh"

CAL="$LOD_SIM_INPUTS/calibration.tsv"
[[ -s "$CAL" ]] || { echo "missing calibration: $CAL" >&2; exit 1; }

mkdir -p "$VALIDATION_OUT/lod"
echo "conditions: $LOD_CONDITIONS"
echo

printf '%-10s %12s %12s %14s\n' condition real_bb_frac sim_bb_frac note
printf '%s\n' "------------------------------------------------------------"

status=0
for cond in $LOD_CONDITIONS; do
    line=$(awk -F'\t' -v c="$cond" '$1==c' "$CAL") || true
    if [[ -z "$line" ]]; then
        printf '%-10s %12s\n' "$cond" "NO CALIBRATION"; status=1; continue
    fi
    nfrag=$(cut -f2 <<<"$line")
    df=$(cut -f3 <<<"$line")
    use_distractor=$(cut -f5 <<<"$line")

    real_frac=$(awk -F'\t' -v c="$cond" '$1==c {printf "%.8f", $4}' \
        "$EXPECTED_DIR/lod_real_counts.tsv")

    outdir="$VALIDATION_OUT/lod/$cond"
    if [[ ! -d "$outdir" ]]; then
        args=(run
            --genomes "$REF_4PATH"
            --targets "$LOD_SIM_INPUTS/targets_subset.fa"
            --sample-target-map "$LOD_SIM_INPUTS/sample_target_map.tsv"
            --sample "BB:AE000783.1"
            --probes "$MTEC_PROBES"
            --num-fragments "$([[ "${VALIDATION_MODE:-fast}" == fast ]] && echo "$NUM_FRAGMENTS" || echo "$nfrag")"
            --fragment-length-mean 225 --fragment-length-min 150 --fragment-length-max 350
            --hybridization-temperature 58 --salt-concentration 1000
            --read-simulator perfect --read-length 101
            --threads "$THREADS"
            --seed 42 --report none --outdir "$outdir")
        if [[ "$use_distractor" == "yes" ]]; then
            args+=(--distractors "$REF_ISCAP" --distractor-fraction "$df")
        fi
        "$BB" "${args[@]}" > "$VALIDATION_OUT/lod/${cond}.log" 2>&1 || {
            printf '%-10s %12s %12s   see %s\n' "$cond" "$real_frac" "FAILED" \
                "$VALIDATION_OUT/lod/${cond}.log"; status=1; continue; }
    fi

    detail=$(find "$outdir" -name group_detail.tsv -o -name detected_detail.tsv | head -1)
    if [[ -z "$detail" ]]; then
        printf '%-10s %12s %12s\n' "$cond" "$real_frac" "NO OUTPUT"; status=1; continue
    fi
    sim_frac=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++){if($i=="reads_assigned")r=i; if($i=="total_reads")r=i; if($i=="reference_id"||$i=="group_name")k=i}; next}
                           {t+=$r; if($k ~ /^BB:/) b+=$r} END{if(t>0) printf "%.8f", b/t; else print "0"}' "$detail")
    printf '%-10s %12s %12s\n' "$cond" "$real_frac" "$sim_frac"
done

echo
echo "Expected titration : $EXPECTED_DIR/lod_real_counts.tsv"
echo

echo "--- per-locus efficiency, relative to rplB ---"
# Targets are named {Species}|{gene}, so per-gene depth comes straight from the
# per-reference detail table; no chromosome re-mapping is needed. See the docstring in
# locus_efficiency.py for how this basis differs from the mosdepth reference values.
dirs=()
for cond in $LOD_CONDITIONS; do
    [[ -d "$VALIDATION_OUT/lod/$cond" ]] && dirs+=("$VALIDATION_OUT/lod/$cond")
done
if [[ ${#dirs[@]} -gt 0 ]]; then
    python3 "$(dirname "${BASH_SOURCE[0]}")/locus_efficiency.py" "${dirs[@]}" || status=1
else
    echo "no completed conditions to compare"
fi
exit $status
