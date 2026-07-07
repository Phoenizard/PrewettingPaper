#!/usr/bin/env python
"""Compare two computed result sets — a plain diff of SAVED pw_line.csv, done
after the results exist (it reads result files, it does NOT recompute).

  conda run -n numenv python scripts/compare.py [pre_root] [post_root] [--tol 1e-3]

Defaults: pre=out/verify (baseline / 优化前结果), post=out/verify_opt (优化后结果).
For each case under pre_root that has a pw_line.csv, load both sides, match rows by
phi2 (within 1e-6), and report max |dphi1| with PASS/FAIL(<=tol). A missing post
side or a row-count mismatch is a hard fail. Exit code is nonzero if any case fails.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

import numpy as np  # noqa: E402


def _load_pw(pw_csv):
    """Read a pw_line.csv to an (n,2) array; missing file -> None; empty -> (0,2)."""
    if not os.path.exists(pw_csv):
        return None
    raw = np.loadtxt(pw_csv, delimiter=",", skiprows=1, ndmin=2)
    return raw if raw.size and raw.shape[1] == 2 else np.empty((0, 2))


def _iter_cases(root):
    for dirpath, _, files in os.walk(root):
        if "pw_line.csv" in files:
            rel = os.path.relpath(dirpath, root).split(os.sep)
            if len(rel) == 3:
                yield tuple(rel)


def _deviation(a, b):
    """Max |dphi1| over phi2-matched rows, plus unmatched counts on each side."""
    used = set()
    max_dev = 0.0
    for p1a, p2a in a:
        j = next((k for k, (_, p2b) in enumerate(b)
                  if k not in used and abs(p2b - p2a) < 1e-6), None)
        if j is None:
            continue
        used.add(j)
        max_dev = max(max_dev, abs(p1a - b[j][0]))
    unmatched_pre = len(a) - len(used)   # pre rows with no post match (regressions)
    unmatched_post = len(b) - len(used)  # post-only rows (recovered endpoints)
    return max_dev, unmatched_pre, unmatched_post


def main(argv):
    tol = 1e-3
    if "--tol" in argv:
        i = argv.index("--tol")
        tol = float(argv[i + 1])
        argv = argv[:i] + argv[i + 2:]
    # --allow-extend: pass a case when every PRE row is reproduced within tol, even
    # if POST added rows (the recovered prewetting-line endpoints). Without it, any
    # row-count difference is a hard fail (strict optimization-regression check).
    allow_extend = "--allow-extend" in argv
    argv = [x for x in argv if x != "--allow-extend"]
    pre = argv[0] if len(argv) >= 1 else os.path.join("out", "verify")
    post = argv[1] if len(argv) >= 2 else os.path.join("out", "verify_opt")

    rels = sorted(_iter_cases(pre))
    if not rels:
        print(f"no pw_line.csv under {pre}")
        return 1

    nfail = 0
    for rel in rels:
        a = _load_pw(os.path.join(pre, *rel, "pw_line.csv"))
        b = _load_pw(os.path.join(post, *rel, "pw_line.csv"))
        if b is None:
            print(f"[FAIL] {'/'.join(rel)}  post side missing")
            nfail += 1
            continue
        dev, um_pre, um_post = _deviation(a, b)
        ok = (dev <= tol) and (um_pre == 0) and (allow_extend or um_post == 0)
        nfail += 0 if ok else 1
        print(f"[{'PASS' if ok else 'FAIL'}] {'/'.join(rel)}  "
              f"n_pre={len(a)} n_post={len(b)} max_dev={dev:.2e} "
              f"unmatched_pre={um_pre} unmatched_post={um_post}")

    mode = "allow-extend" if allow_extend else "strict"
    print(f"\n{len(rels) - nfail}/{len(rels)} passed (tol={tol:g}, {mode})  "
          f"pre={pre} post={post}")
    return 1 if nfail else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
