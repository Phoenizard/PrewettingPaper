#!/usr/bin/env python
"""Render phase-map figures from computed data (DATA -> figure).

Decoupled counterpart of scripts/verify.py: verify.py computes and writes
pw_line.csv only; this turns a case's inputs — experiment variables plus
computed results — into overlay.png, via plotting.render_phase_map. It is the
thin I/O adapter (locate case, reconstruct params from the dir name, load the
CSV); the reusable (variables + results) -> figure logic lives in
src/plotting.py and is not verify-specific.

Runs single-process and is cheap, so it never contends with the parallel
compute workers, and figures can be regenerated or restyled any time from the
CSVs without recomputing.

Usage:
  conda run -n numenv python scripts/plot.py                          # every out/verify case with pw_line.csv
  conda run -n numenv python scripts/plot.py <chi_dir> <om_dir> <chibb_dir>   # one case
  conda run -n numenv python scripts/plot.py --all [N]               # first N cases (sorted)
  --force   redraw overlay.png even if it already exists
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "src"))

import numpy as np  # noqa: E402
import plotting as P  # noqa: E402
import cases  # noqa: E402

# VERIFY_OUT overrides the results root (matches verify.py), so an optimized run
# in out/verify_opt gets its overlays plotted too.
VERIFY_ROOT = os.environ.get("VERIFY_OUT") or os.path.join(ROOT, "out", "verify")


def _load_pw(pw_csv):
    """Read a pw_line.csv back to an (n,2) array; empty file -> (0,2)."""
    raw = np.loadtxt(pw_csv, delimiter=",", skiprows=1, ndmin=2)
    return raw if raw.size and raw.shape[1] == 2 else np.empty((0, 2))


def _iter_done(root=VERIFY_ROOT):
    """Yield rel=(chi_dir, om_dir, chibb_dir) for every leaf that has pw_line.csv."""
    if not os.path.isdir(root):
        return
    for dirpath, _, files in os.walk(root):
        if "pw_line.csv" not in files:
            continue
        rel = os.path.relpath(dirpath, root).split(os.sep)
        if len(rel) == 3:
            yield tuple(rel)


def plot_case(rel, force=False):
    """Render <VERIFY_ROOT>/<rel>/overlay.png from its pw_line.csv. Returns a status."""
    out_dir = os.path.join(VERIFY_ROOT, *rel)
    overlay = os.path.join(out_dir, "overlay.png")
    pw_csv = os.path.join(out_dir, "pw_line.csv")
    if not os.path.exists(pw_csv):
        return "no_data"
    if os.path.exists(overlay) and not force:
        return "skip"
    chi, surf = cases.parse_case(*rel)
    P.render_phase_map(
        overlay, chi=chi, surf=surf, pw=_load_pw(pw_csv),
        title=f"Verify {rel[0]} | {rel[1]} | {rel[2]}",
        params_text="ours (DE solver) — compare result/.../overlay_omega1_omega2.png",
    )
    return "ok"


def main(argv):
    force = "--force" in argv
    argv = [a for a in argv if a != "--force"]

    if argv and argv[0] == "--all":
        limit = int(argv[1]) if len(argv) > 1 else None
        rels = sorted(_iter_done())
        if limit is not None:
            rels = rels[:limit]
    elif len(argv) == 3:
        rels = [tuple(argv)]
    else:
        rels = sorted(_iter_done())

    n = {"ok": 0, "skip": 0, "no_data": 0}
    for rel in rels:
        status = plot_case(rel, force)
        n[status] = n.get(status, 0) + 1
        if status != "skip":
            print(f"{'/'.join(rel)} -> {status}")
    print(f"done: {n['ok']} drawn, {n['skip']} skipped, {n['no_data']} missing-data")


if __name__ == "__main__":
    main(sys.argv[1:])
