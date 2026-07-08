#!/usr/bin/env bash
# 12-case verification sweep. Parallel across cases (Python stays serial).
# Usage: bash scripts/run_verify.sh [PARALLEL]
# Env: PYTHON to override the interpreter (server: /root/miniconda3/envs/numenv/bin/python).
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
PY="${PYTHON:-python}"
PAR="${1:-12}"

# one BLAS thread per worker; workers provide the parallelism
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export MPLBACKEND=Agg

# warm the matplotlib font cache once before parallel workers import it
"$PY" -c "import matplotlib.pyplot"

CASES="$(grep -E '^[[:space:]]*-[[:space:]]' "$ROOT/config/verify_cases.yaml" | sed -E 's/^[[:space:]]*-[[:space:]]*//')"
N_CASES="$(echo "$CASES" | wc -l | tr -d ' ')"
echo "[sweep] $N_CASES cases, parallel=$PAR"

export ROOT PY
echo "$CASES" | xargs -P "$PAR" -I{} bash -c '
  rel="{}"
  out="$ROOT/out/verify/$rel"
  log_dir="$ROOT/log/verify/$rel"
  mkdir -p "$out" "$log_dir"
  echo "[start] $rel"
  if "$PY" "$ROOT/scripts/run_case.py" --case-rel "$rel" >"$log_dir/run.log" 2>&1 \
     && "$PY" "$ROOT/scripts/plot_case.py" --case-dir "$out" >>"$log_dir/run.log" 2>&1; then
    echo "[done] $rel"
  else
    echo "[FAIL] $rel (see $log_dir/run.log)"
  fi
'

"$PY" "$ROOT/scripts/build_summary.py"
