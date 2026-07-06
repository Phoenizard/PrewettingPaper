#!/usr/bin/env bash
# Parallel small-scale verification over the top-N cases of result_cases.txt.
#
# Usage:  bash scripts/run_parallel.sh [N] [CORES]
#   N      number of cases from the top of result_cases.txt (0 or empty = all)
#   CORES  parallel workers (default: nproc, else 4)
#   PY     python launcher (env var, default: "conda run -n numenv python")
#
# Each case runs as an independent single-case verify.py invocation with
# --skip-existing (resume) and --no-summary (avoid concurrent SUMMARY writes).
# SUMMARY.csv is rebuilt once at the end from all pw_line.csv outputs.
# Parallelism lives here in the launch script, not in the Python code.
set -uo pipefail
cd "$(dirname "$0")/.."

N="${1:-0}"
CORES="${2:-$( (command -v nproc >/dev/null 2>&1 && nproc) || echo 4)}"
PY="${PY:-conda run -n numenv python}"
MANIFEST="result_cases.txt"

[ -f "$MANIFEST" ] || { echo "missing $MANIFEST" >&2; exit 1; }
if [ "${N}" -gt 0 ] 2>/dev/null; then
  WORK="$(head -n "$N" "$MANIFEST")"
else
  WORK="$(cat "$MANIFEST")"
fi
NLINES="$(printf '%s\n' "$WORK" | grep -c .)"

echo "cases=$NLINES  cores=$CORES  python='$PY'"
start=$(date +%s)
printf '%s\n' "$WORK" | xargs -P "$CORES" -L1 $PY scripts/verify.py --skip-existing --no-summary \
  || echo "warning: some cases exited non-zero"
elapsed=$(( $(date +%s) - start ))

echo "rebuilding SUMMARY.csv ..."
$PY scripts/verify.py --rebuild-summary

echo "done: $NLINES cases in ${elapsed}s on $CORES cores (~$(awk "BEGIN{printf \"%.1f\", $elapsed/$NLINES}")s/case wall)"
