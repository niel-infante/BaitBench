#!/usr/bin/env bash
# Fast validation subsets — minutes, for the inner development loop.
# 3 LOD conditions, 1 mock run, 20k fragments each.
set -uo pipefail
export VALIDATION_MODE=fast
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

bash "$HERE/check_data.sh" || exit 1
echo
echo "=============================================================="
echo " FAST validation subsets (VALIDATION_MODE=fast)"
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
[[ $rc -eq 0 ]] && echo "fast validation: all gates passed" \
               || echo "fast validation: one or more gates failed (expected pre-Stage 3)"
exit $rc
