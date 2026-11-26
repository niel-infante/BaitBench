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
#
# Fragment length sd is measured, not defaulted. The default (max - min) / 6 gives 3,167
# and a realised sd of 2,793, against 2,608 measured on the non-enriched arm — too wide
# in both tails. With the measured value the simulated distribution matches the real one
# closely: median 5,933 vs 5,948, mean 6,001 vs 6,141.
#
# Normal, not lognormal, despite the plan calling for lognormal: the real distribution is
# barely skewed (mean/median = 1.033), and a lognormal over-skews it — median 5,315
# against a real 5,948.
#
# Mapping preset is `map-hifi`, matching what produced the real counts
# (scripts/map_real_reads.sh). It was BaitBench's `sr` default — a *short-read* preset
# applied to 5.8 kb reads. Measured effect on the gate: negligible (CCC 0.9691 vs
# 0.9688), because composition is a ratio and the preset's losses are near-proportional.
# Changed anyway: `sr` drops 8.7% of simulated reads that `map-hifi` maps, it is
# indefensible to describe in a paper, and it is ~35% slower.
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
        # Expanded with the ${x+...} guard: macOS ships bash 3.2, where "${arr[@]}" on
        # an empty array is an unbound-variable error under `set -u`.
        conc_args=()
        [[ -n "$PROBE_CONC" ]] && conc_args=(--probe-concentration "$PROBE_CONC")
        [[ -n "$MELTING_WINDOW" ]] && conc_args+=(--melting-window "$MELTING_WINDOW")
        "$BB" run ${conc_args[@]+"${conc_args[@]}"} \
            --targets "$SIMVAL_GENOMES" \
            --probes  "$PANELS/res_probes.fasta" \
            --sample  "$SIMVAL_SAMPLE" \
            --groups  "$SIMVAL_GROUPS" \
            --num-fragments "$NUM_FRAGMENTS" \
            --fragment-length-mean 5800 \
            --fragment-length-min 1000 \
            --fragment-length-max 20000 \
            --fragment-length-sd 2608 \
            --capture-fraction "$SIMVAL_CAPTURE_FRACTION" \
            --read-simulator perfect \
            --read-length 20000 \
            --minimap-preset map-hifi \
            --hybridization-temperature 65 \
            --salt-concentration 50 \
            --threads "$THREADS" \
            --probe-hits "$PROBE_CACHE_DIR/simval_res_panel.tsv" \
            --seed 42 \
            --report none \
            --cleanup \
            --outdir "$outdir" \
            > "$VALIDATION_OUT/${tag}.log" 2>&1 \
            || { echo "  simulate FAILED, see $VALIDATION_OUT/${tag}.log" >&2; status=1; continue; }
    fi

    echo
    VALIDATION_OUT="$VALIDATION_OUT" SIMVAL_REAL="$SIMVAL_REAL" SIMVAL_GROUPS="$SIMVAL_GROUPS" \
    SIMVAL_TSV="$VALIDATION_OUT/${tag}_composition.tsv" \
        python3 "$VALIDATION_DIR/simval/compare.py" "$tag" "$run" || status=1

    # The scatter is the figure the manuscript needs, so generate it with the gate
    # rather than by hand afterwards. R is optional: a missing Rscript or package
    # should not fail a gate that has already produced its numbers.
    if command -v Rscript >/dev/null && [[ -s "$VALIDATION_OUT/${tag}_composition.tsv" ]]; then
        Rscript "$VALIDATION_DIR/simval/plot.R" \
            "$VALIDATION_OUT/${tag}_composition.tsv" \
            "$VALIDATION_OUT/${tag}_composition.png" \
            "D6331 mock community ($run): simulated vs real composition" \
            || echo "  composition plot failed (non-fatal)" >&2
    fi
    echo
done

# ---------------------------------------------------------------------------------
# Enrichment gate: the non-circular target (improvements/parameterisation.md B0).
#
# The composition gate above sets --sample from the manufacturer's theoretical table and
# scores absolute composition. This one takes input composition from the *non-enriched*
# arm of the same BioProject and scores the per-organism enrichment factor
# (enriched / non-enriched) — input and target are different libraries, so there is no
# route by which guessing the input well can flatter the model.
#
# It needs no reads: an organism's expected output share is its pooled
# `weight x retention_mass`, which `fragments.retention.tsv` already contains. So this
# runs `simulate` alone against the shared probe cache and scores analytically, which
# costs ~90s rather than the ~25 minutes a full pipeline run takes.
echo "=== enrichment gate (non-enriched arm as input) ==="
inputs="$VALIDATION_OUT/nonenriched_inputs"
if [[ ! -s "$inputs/input_share.json" ]]; then
    python3 "$VALIDATION_DIR/simval/input_weights.py" \
        --real "$SIMVAL_REAL" --groups "$SIMVAL_GROUPS" \
        --genomes "$SIMVAL_GENOMES" --outdir "$inputs" || status=1
fi

for run in $SIMVAL_RUNS; do
    real="$SIMVAL_REAL/${run}_counts.tsv"
    tag="${run}_enrich"
    outdir="$VALIDATION_OUT/$tag"
    if [[ ! -s "$outdir/fragments.retention.tsv" ]]; then
        echo "[$(date +%H:%M:%S)] $tag: simulating"
        mkdir -p "$outdir"
        # `prepare` is what normally writes weights.txt; simulate needs one, so derive
        # it here from the non-enriched sample file. Exponential notation is required:
        # these weights are ~1e-8 and a fixed six decimals rounds them all to zero.
        awk -F'\t' 'BEGIN{OFS="\t"} !/^#/ && NF>=2 {printf "%s\t%e\n", $1, $2}' \
            "$inputs/sample_nonenriched.tsv" > "$outdir/weights.txt"
        sim_args=(simulate
            --reference "$SIMVAL_GENOMES"
            --weights "$outdir/weights.txt"
            --probes "$PANELS/res_probes.fasta"
            --num-fragments "$NUM_FRAGMENTS"
            --fragment-length-mean 5800 --fragment-length-min 1000 --fragment-length-max 20000
            --fragment-length-sd 2608
            --capture-fraction "$SIMVAL_CAPTURE_FRACTION"
            --hybridization-temperature 65 --salt-concentration 50
            --threads "$THREADS"
            --probe-hits "$PROBE_CACHE_DIR/simval_res_panel.tsv"
            --seed 42 --output "$outdir/fragments.fa")
        [[ -n "$PROBE_CONC" ]] && sim_args+=(--probe-concentration "$PROBE_CONC")
        [[ -n "$MELTING_WINDOW" ]] && sim_args+=(--melting-window "$MELTING_WINDOW")
        "$BB" "${sim_args[@]}" > "$VALIDATION_OUT/${tag}.log" 2>&1 \
            || { echo "  simulate FAILED, see $VALIDATION_OUT/${tag}.log" >&2; status=1; continue; }
    else
        echo "[$(date +%H:%M:%S)] $tag: reusing $outdir"
    fi

    # Candida albicans used to be excluded here: its reference is 13,062 contigs over
    # 12.58 Mb against a 5.8 kb fragment length, and it came out at 0.244x its intended
    # input share. That was the simulator's bug, not the assembly's fault — retention
    # mass counted fragment starts rather than bases — and Stage 5 fixed it. Candida now
    # lands at ratio 1.000 of its input share and is predicted at 0.156 against a
    # measured 0.130, so it is scored like everything else.
    python3 "$VALIDATION_DIR/simval/enrichment.py" \
        --run "$outdir" --share "$inputs/input_share.json" \
        --enriched "$real" --groups "$SIMVAL_GROUPS" \
        --tsv "$VALIDATION_OUT/${tag}.tsv" || status=1

    if command -v Rscript >/dev/null && [[ -s "$VALIDATION_OUT/${tag}.tsv" ]]; then
        Rscript "$VALIDATION_DIR/simval/plot_enrichment.R" \
            "$VALIDATION_OUT/${tag}.tsv" \
            "$VALIDATION_OUT/${tag}.png" \
            "D6331 mock community ($run): predicted vs measured enrichment" \
            || echo "  enrichment plot failed (non-fatal)" >&2
    fi
    echo
done

echo "Reference composition: $EXPECTED_DIR/simval_real_composition.tsv"
exit $status
