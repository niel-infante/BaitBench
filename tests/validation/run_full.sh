#!/usr/bin/env bash
# Full validation gates — hours. Run before declaring a stage complete.
# All 11 LOD conditions, 3 mock replicates, 150k fragments.
set -uo pipefail
export VALIDATION_MODE=full
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
exec bash "$HERE/run_fast.sh"
