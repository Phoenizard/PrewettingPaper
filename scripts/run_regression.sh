#!/usr/bin/env bash
# Parallel regression gate: recompute every out/verify baseline in memory and
# assert its pre-wetting transition line is unchanged (<= TOL) vs the saved CSV.
#
# Per-case logs MIRROR the result tree:
#   result: out/verify/<chi>/<om>/<chibb>/pw_line.csv
#   log:    out/logs/regression/<chi>/<om>/<chibb>/run.log
# Parallelism lives here (xargs -P), not in Python. Read-only w.r.t. out/verify —
# writes only under out/logs/. Runs whatever optimization toggles the code ships
# with (equilibrium.USE_ANALYTIC_JAC / USE_WARM_START).
#
# Usage: bash scripts/run_regression.sh [CORES]
#   TOL env overrides the threshold (default 1e-3); PY overrides the launcher.
set -uo pipefail
cd "$(dirname "$0")/.."

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
       NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1

# nproc honors OMP_NUM_THREADS (which we pin to 1 above), so use `nproc --all` to
# get the real core count for the parallel fan-out.
CORES="${1:-$( (command -v nproc >/dev/null 2>&1 && nproc --all) || echo 4)}"
PY="${PY:-conda run -n numenv python}"
TOL="${TOL:-1e-3}"
LOGROOT="out/logs/regression"

mapfile -t RELS < <(find out/verify -name pw_line.csv | sed 's#^out/verify/##; s#/pw_line.csv$##')
echo "regression: ${#RELS[@]} baselines  cores=$CORES  tol=$TOL  python='$PY'"
[ "${#RELS[@]}" -gt 0 ] || { echo "no baselines under out/verify/"; exit 1; }

run_one() {
  rel="$1"
  a="${rel%%/*}"; rest="${rel#*/}"; b="${rest%%/*}"; c="${rest#*/}"
  d="$LOGROOT/$a/$b/$c"; mkdir -p "$d"
  $PY scripts/regression.py --tol "$TOL" "$a" "$b" "$c" > "$d/run.log" 2>&1
  grep -aE "^\[(PASS|FAIL)\]" "$d/run.log" | tail -1   # surface verdict to driver
}
export -f run_one
export PY TOL LOGROOT

start=$(date +%s)
printf '%s\n' "${RELS[@]}" | xargs -P "$CORES" -I{} bash -c 'run_one "$@"' _ {}
elapsed=$(( $(date +%s) - start ))

echo "=== summary (${elapsed}s) ==="
nfail=$(grep -rlaE "^\[FAIL\]" "$LOGROOT" 2>/dev/null | wc -l | tr -d ' ')
ntot="${#RELS[@]}"
echo "$(( ntot - nfail ))/$ntot passed (tol=$TOL)"
if [ "$nfail" -eq 0 ]; then echo "REGRESSION_PASS"; else echo "REGRESSION_FAIL ($nfail cases)"; fi