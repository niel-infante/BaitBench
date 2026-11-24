#!/usr/bin/env bash
# Fast validation subsets — minutes, for the inner development loop.
# 3 LOD conditions, 1 mock run, 20k fragments each.
set -uo pipefail
# Default to fast, but do not clobber a mode already set by the caller. This was
# `export VALIDATION_MODE=fast`, which overwrote the `full` that run_full.sh exports
# immediately before exec'ing this script — so `run_full.sh` silently ran the fast
# subsets (3 LOD conditions, 1 mock run, 20k fragments) and reported them as the
# release gate.
: "${VALIDATION_MODE:=fast}"
export VALIDATION_MODE
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Refuse to run two gates against the same VALIDATION_OUT. Both would write the same
# result directories and the same logs, and a second run that wipes VALIDATION_OUT while
# the first is mid-flight produces an interleaved log and half-written outputs that still
# report "all gates passed". mkdir is atomic, so it is a sufficient lock here.
source "$HERE/config.sh"
LOCK="$VALIDATION_OUT.lock"
if ! mkdir "$LOCK" 2>/dev/null; then
    echo "A validation run already holds $LOCK." >&2
    echo "Wait for it, or set VALIDATION_OUT to a different path for this run." >&2
    echo "If no run is active, remove the lock: rmdir '$LOCK'" >&2
    exit 1
fi
trap 'rmdir "$LOCK" 2>/dev/null' EXIT

bash "$HERE/check_data.sh" || exit 1
echo
echo "=============================================================="
printf ' %s validation (VALIDATION_MODE=%s)\n' \
    "$([[ $VALIDATION_MODE == full ]] && echo FULL || echo FAST)" "$VALIDATION_MODE"
echo "=============================================================="

rc=0
echo; echo "--- distractor-fraction sweep (cargo) ---"
# Not `-- --ignored`: the distractor sweep was un-ignored in Stage 3, and that filter
# would now select zero tests and report "ok" while testing nothing.
( cd "$HERE/../.." && cargo test --test simulation -- --nocapture ) || rc=1

echo; echo "--- D6331 mock community ---"
bash "$HERE/simval/run.sh" || rc=1

echo; echo "--- Borrelia LOD titration ---"
bash "$HERE/lod/run.sh" || rc=1

echo
[[ $rc -eq 0 ]] && echo "$VALIDATION_MODE validation: all gates passed" \
               || echo "$VALIDATION_MODE validation: one or more gates failed"
exit $rc
