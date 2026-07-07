#!/usr/bin/env python
"""Regression gate: recompute baselines in memory, assert the transition line is
unchanged within TOL. Read-only — NEVER writes under out/.

For each existing out/verify/**/pw_line.csv baseline, recompute prewetting_line
(with whatever optimization toggles equilibrium currently has) and compare the
transition phi1 per phi2 against the saved CSV. A case FAILS if any matched phi2
deviates by more than TOL, or if the row counts differ (branch-appearance shift).
Exit code is nonzero if any case fails.

Usage:
  conda run -n numenv python scripts/regression.py                     # all baselines
  conda run -n numenv python scripts/regression.py <chi_dir> <om_dir> <chibb_dir>  # one case
  --tol 1e-3     acceptance threshold on |dphi1| (default 1e-3)
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, HERE)

import numpy as np  # noqa: E402
import binodal as B  # noqa: E402
import cases  # noqa: E402
import verify  # noqa: E402

VERIFY_ROOT = os.path.join(ROOT, "out", "verify")


def _load_pw(pw_csv):
    raw = np.loadtxt(pw_csv, delimiter=",", skiprows=1, ndmin=2)
    return raw if raw.size and raw.shape[1] == 2 else np.empty((0, 2))


def _iter_baselines():
    for dirpath, _, files in os.walk(VERIFY_ROOT):
        if "pw_line.csv" in files:
            rel = os.path.relpath(dirpath, VERIFY_ROOT).split(os.sep)
            if len(rel) == 3:
                yield tuple(rel)


def _deviation(new, saved):
    """Max |dphi1| over phi2-matched rows, plus count of unmatched rows."""
    used = set()
    max_dev = 0.0
    for p1n, p2n in new:
        j = next((k for k, (_, p2s) in enumerate(saved)
                  if k not in used and abs(p2s - p2n) < 1e-6), None)
        if j is None:
            continue
        used.add(j)
        max_dev = max(max_dev, abs(p1n - saved[j][0]))
    unmatched = (len(new) - len(used)) + (len(saved) - len(used))
    return max_dev, unmatched


def check_case(rel, tol):
    pw_csv = os.path.join(VERIFY_ROOT, *rel, "pw_line.csv")
    saved = _load_pw(pw_csv)
    chi, surf = cases.parse_case(*rel)
    binodal = B.binodal_from_hull(chi)
    new = verify.prewetting_line(chi, surf, binodal,  # in memory; no writes
                                 progress="/".join(rel))
    max_dev, unmatched = _deviation(new, saved)
    ok = (max_dev <= tol) and (unmatched == 0)
    print(f"[{'PASS' if ok else 'FAIL'}] {'/'.join(rel)}  "
          f"n_saved={len(saved)} n_new={len(new)} "
          f"max_dev={max_dev:.2e} unmatched={unmatched}")
    return ok


def main(argv):
    tol = 1e-3
    if "--tol" in argv:
        i = argv.index("--tol")
        tol = float(argv[i + 1])
        argv = argv[:i] + argv[i + 2:]
    rels = [tuple(argv)] if len(argv) == 3 else sorted(_iter_baselines())
    if not rels:
        print("no baselines found under out/verify/")
        return 1
    results = [check_case(rel, tol) for rel in rels]
    n_fail = results.count(False)
    print(f"\n{len(results) - n_fail}/{len(results)} passed (tol={tol:g})")
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
