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
# Compute writes DATA only (pw_line.csv) — no worker imports matplotlib, so
# there is no font-cache race and no compute/plot contention. After compute,
# SUMMARY.csv is rebuilt from all pw_line.csv, then figures are rendered once,
# single-process, by scripts/plot.py. Parallelism lives here, not in Python.
set -uo pipefail
cd "$(dirname "$0")/.."

# Pin each worker to a single BLAS/OpenMP thread. Without this, numpy/scipy spawn
# one thread per core PER process, so P parallel workers oversubscribe the machine
# (P*cores threads) and thrash. One thread each => P workers map cleanly onto P cores.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
       NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1

N="${1:-0}"
# nproc honors OMP_NUM_THREADS (pinned to 1 above), so use `nproc --all` for real cores.
CORES="${2:-$( (command -v nproc >/dev/null 2>&1 && nproc --all) || echo 4)}"
PY="${PY:-conda run -n numenv python}"
MANIFEST="${MANIFEST:-result_cases.txt}"   # override to run a case subset (e.g. the 10 baselines)

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

# Render figures once, single-process, from the CSVs. Decoupled from compute:
# if this step fails the data is already safe, and figures can be redrawn later.
echo "rendering figures ..."
$PY scripts/plot.py

echo "done: $NLINES cases in ${elapsed}s on $CORES cores (~$(awk "BEGIN{printf \"%.1f\", $elapsed/$NLINES}")s/case wall) [compute]"
