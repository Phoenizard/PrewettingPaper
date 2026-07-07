#!/usr/bin/env bash
# Serve the results gallery over local http, so the browser has no file:// image
# restriction (both ours and the result/ reference figures load). Server root is
# $RESULTS_DB; result/ resolves through the in-DB symlink created by build_gallery.py.
# Ctrl-C to stop.
#
#   RESULTS_DB=./database bash scripts/serve_gallery.sh [PORT]     # -> http://localhost:PORT
set -uo pipefail
cd "$(dirname "$0")/.."
RESULTS_DB="${RESULTS_DB:-./database}"
PORT="${1:-8000}"
[ -f "$RESULTS_DB/index.html" ] || {
  echo "no $RESULTS_DB/index.html — run build_gallery.py first" >&2; exit 1; }
echo "serving $RESULTS_DB at http://localhost:$PORT/  (Ctrl-C to stop)"
cd "$RESULTS_DB"
exec python3 -m http.server "$PORT"
