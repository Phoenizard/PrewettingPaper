#!/usr/bin/env bash
# Verification sweep. Parallel across cases AND within each case (line workers).
# Usage: bash scripts/run_verify.sh [PARALLEL] [CASES_FILE] [LINE_WORKERS]
# Env: PYTHON to override the interpreter (server: /root/miniconda3/envs/numenv/bin/python).
# Sizing: PARALLEL x LINE_WORKERS should fill the instance CPU quota
# (check cgroup cpu.max, not nproc — nproc shows the host, not the container).
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
PY="${PYTHON:-python}"
PAR="${1:-12}"
# one launch = one timestamp shared by the sweep log and all per-case logs;
# log/ is flat: log/<stamp>-verify[-<qualifier>].log
STAMP="${STAMP:-$(date +%Y%m%d-%H%M)}"

# one BLAS thread per worker; workers provide the parallelism
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export MPLBACKEND=Agg

# warm the matplotlib font cache once before parallel workers import it
"$PY" -c "import matplotlib.pyplot"

CASES_FILE="${2:-$ROOT/config/verify_cases.yaml}"
LINE_WORKERS="${3:-1}"
CASES="$(grep -E '^[[:space:]]*-[[:space:]]' "$CASES_FILE" | sed -E 's/^[[:space:]]*-[[:space:]]*//')"
N_CASES="$(echo "$CASES" | wc -l | tr -d ' ')"
echo "[sweep] $N_CASES cases, parallel=$PAR"

export ROOT PY STAMP LINE_WORKERS
mkdir -p "$ROOT/log"
echo "$CASES" | xargs -P "$PAR" -I{} bash -c '
  rel="{}"
  out="$ROOT/out/$rel"
  log_file="$ROOT/log/$STAMP-verify-${rel//\//+}.log"
  mkdir -p "$out"
  echo "[start] $rel"
  if "$PY" "$ROOT/scripts/run_case.py" --case-rel "$rel" --line-workers "$LINE_WORKERS" >"$log_file" 2>&1 \
     && "$PY" "$ROOT/scripts/plot_case.py" --case-dir "$out" >>"$log_file" 2>&1; then
    echo "[done] $rel"
  else
    echo "[FAIL] $rel (see $log_file)"
  fi
'

"$PY" "$ROOT/scripts/build_summary.py"
