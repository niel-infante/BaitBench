#!/usr/bin/env bash
# Verify that the externally-referenced validation datasets are present.
# Exits non-zero and names every missing path, so a failure is actionable.
set -uo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"

missing=0
check () {  # check <label> <path>
    if [[ -e "$2" ]]; then
        printf '  %-28s %s\n' "$1" "ok"
    else
        printf '  %-28s MISSING  %s\n' "$1" "$2"
        missing=1
    fi
}

echo "BaitBench validation data (CROSSLAB=$CROSSLAB)"
echo
echo "Borrelia LOD:"
check "real read counts"       "$LOD_COUNTS"
check "mosdepth coverage"      "$LOD_COVERAGE"
check "target intervals"       "$LOD_TARGETS_BED"
check "probe BLAST table"      "$LOD_PROBES_BLAST"
check "MTEC probes"            "$MTEC_PROBES"
check "4Path reference"        "$REF_4PATH"
check "Ixodes distractor"      "$REF_ISCAP"
echo
echo "D6331 mock community:"
check "genomes"                "$SIMVAL_GENOMES"
check "groups"                 "$SIMVAL_GROUPS"
check "sample weights"         "$SIMVAL_SAMPLE"
check "real mapped counts"     "$SIMVAL_REAL"
check "probe panels"           "$PANELS"
check "theoretical composition" "$D6331_COMPOSITION"
echo
echo "Binary:"
check "baitbench"              "$BB"

if [[ $missing -ne 0 ]]; then
    cat >&2 <<'MSG'

Some validation inputs are missing. They are referenced in place, never copied.
Point the harness at them by exporting CROSSLAB, or override individual paths:

    CROSSLAB=/data/crosslab bash tests/validation/run_fast.sh
    LOD_COUNTS=/elsewhere/all_counts.tsv bash tests/validation/run_fast.sh

See tests/validation/README.md for what each dataset is and how it was produced.
MSG
    exit 1
fi
echo
echo "All validation inputs present."
