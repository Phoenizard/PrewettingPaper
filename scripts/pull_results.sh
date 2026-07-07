#!/usr/bin/env bash
# Pull our results from the ssh compute server into the LOCAL results database.
# Run on the local Mac. Merge is by UNION (no --delete): different sessions or
# different ephemeral servers rsync into the same DB and combine cleanly, because
# each case path (chi/om/chibb) is deterministic and content is reproducible.
#
# Usage:
#   SSH="-p 32829 root@connect.cqa1.seetacloud.com" bash scripts/pull_results.sh [--logs]
#
# Env:
#   RESULTS_DB  database root (default: ./database; later e.g. /Volumes/<HDD>/PrewettingResults)
#   REMOTE_DIR  project dir on the server (default: ~/autodl-tmp/PrewettingPaper)
set -uo pipefail
cd "$(dirname "$0")/.."

: "${SSH:?set SSH, e.g. SSH=\"-p 32829 root@host\"}"
RESULTS_DB="${RESULTS_DB:-./database}"
REMOTE_DIR="${REMOTE_DIR:-~/autodl-tmp/PrewettingPaper}"   # ~ expands on the remote

PORT="$(printf '%s\n' "$SSH" | grep -oE '\-p[[:space:]]+[0-9]+' | grep -oE '[0-9]+' | head -1)"
TARGET="$(printf '%s\n' "$SSH" | grep -oE '[A-Za-z0-9._-]+@[A-Za-z0-9._.-]+' | head -1)"
[ -n "$TARGET" ] || { echo "could not parse user@host from SSH='$SSH'" >&2; exit 1; }
SSH_E="ssh -o ConnectTimeout=15${PORT:+ -p $PORT}"

mkdir -p "$RESULTS_DB/verify"
echo "pull: $TARGET:$REMOTE_DIR/out/verify/  ->  $RESULTS_DB/verify/"
rsync -a -e "$SSH_E" \
  "$TARGET:$REMOTE_DIR/out/verify/" "$RESULTS_DB/verify/"

if [ "${1:-}" = "--logs" ]; then
  mkdir -p "$RESULTS_DB/logs/verify"
  rsync -a -e "$SSH_E" \
    "$TARGET:$REMOTE_DIR/out/logs/verify/" "$RESULTS_DB/logs/verify/" || true
fi

n="$(find "$RESULTS_DB/verify" -name pw_line.csv 2>/dev/null | wc -l | tr -d ' ')"
echo "done. cases in DB: $n   (RESULTS_DB=$RESULTS_DB)"
