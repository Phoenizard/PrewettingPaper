#!/usr/bin/env python
"""Verify result/ pre-wetting cases with our equilibrium-DE solver.

Output mirrors result/ 1:1 under out/verify/<chi_dir>/<om_dir>/<chibb_dir>/:
  overlay.png   binodal + our pre-wetting line
  pw_line.csv   computed (phi1_inf, phi2_inf) transition points
Top-level out/verify/SUMMARY.csv records every case run.

Usage:
  conda run -n numenv python scripts/verify.py                       # T-a template case
  conda run -n numenv python scripts/verify.py <chi_dir> <om_dir> <chibb_dir>
  conda run -n numenv python scripts/verify.py --all [N]             # first N (or all) cases
  conda run -n numenv python scripts/verify.py --rebuild-summary     # rebuild SUMMARY.csv from outputs
  --skip-existing   reuse cases already having overlay.png + pw_line.csv (resume)
  --no-summary      skip SUMMARY.csv write (for parallel workers; rebuild after)
"""
import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "src"))

import numpy as np  # noqa: E402
import thermo as T  # noqa: E402
import binodal as B  # noqa: E402
import equilibrium as E  # noqa: E402
import plotting as P  # noqa: E402
import cases  # noqa: E402

KAPPA = T.Kappa(1.0, 1.0)
SUMMARY = os.path.join(ROOT, "out", "verify", "SUMMARY.csv")
SUMMARY_COLS = ["chi_dir", "om_dir", "chibb_dir", "n_pw",
                "phi1_min", "phi1_max", "phi2_min", "phi2_max", "status"]


def _branches(binodal):
    """Dilute-left and dense-right phi1(phi2) interpolators + apex phi2."""
    pts = binodal.points
    apex_i = int(np.argmax(pts[:, 1]))
    apex_phi1, apex_phi2 = pts[apex_i, 0], pts[apex_i, 1]
    left = pts[pts[:, 0] <= apex_phi1]
    right = pts[pts[:, 0] >= apex_phi1]
    left = left[np.argsort(left[:, 1])]
    right = right[np.argsort(right[:, 1])]
    f_left = lambda p2: float(np.interp(p2, left[:, 1], left[:, 0]))
    f_right = lambda p2: float(np.interp(p2, right[:, 1], right[:, 0]))
    return f_left, f_right, apex_phi2


def prewetting_line(chi, surf, binodal):
    """(phi1*, phi2) points where gamma_thin = gamma_thick, scanning the dilute flank."""
    f_left, f_right, apex_phi2 = _branches(binodal)
    pw = []
    for phi2 in np.arange(0.01, 0.85 * apex_phi2, 0.01):
        bl = f_left(phi2)
        dense = [(f_right(phi2), phi2), (0.97 * f_right(phi2), phi2)]
        grid = bl + np.arange(-0.03, 0.012, 0.0025)
        grid = grid[grid > 1e-3]
        prev = None
        for phi1 in grid:
            st = E.find_states(chi, (phi1, phi2), surf, KAPPA, dense_seeds=dense)
            if len(st) >= 2:
                d = st[0].gamma - st[-1].gamma
                if prev is not None and prev[1] * d < 0:
                    p0, d0 = prev
                    pw.append((p0 + (phi1 - p0) * (0 - d0) / (d - d0), phi2))
                    break
                prev = (phi1, d)
            else:
                prev = None
    return np.array(pw)


def _row_from_pw(rel, pw):
    if len(pw):
        return [*rel, len(pw), f"{pw[:,0].min():.4f}", f"{pw[:,0].max():.4f}",
                f"{pw[:,1].min():.4f}", f"{pw[:,1].max():.4f}", "ok"]
    return [*rel, 0, "", "", "", "", "no_pw"]


def _load_pw(pw_csv):
    """Read a pw_line.csv back to an (n,2) array; empty file -> (0,2)."""
    raw = np.loadtxt(pw_csv, delimiter=",", skiprows=1, ndmin=2)
    return raw if raw.size and raw.shape[1] == 2 else np.empty((0, 2))


def verify_case(rel, chi, surf, skip_existing=False):
    out_dir = cases.verify_dir(rel, ROOT)
    overlay = os.path.join(out_dir, "overlay.png")
    pw_csv = os.path.join(out_dir, "pw_line.csv")
    if skip_existing and os.path.exists(overlay) and os.path.exists(pw_csv):
        return _row_from_pw(rel, _load_pw(pw_csv))

    os.makedirs(out_dir, exist_ok=True)
    binodal = B.binodal_from_hull(chi)
    pw = prewetting_line(chi, surf, binodal)

    np.savetxt(pw_csv, pw if len(pw) else np.empty((0, 2)),
               delimiter=",", header="phi1_inf,phi2_inf", comments="")
    branches = [{"points": pw, "state": "pw", "branch_id": 3}] if len(pw) else None
    P.plot_phase_map(
        overlay, binodal, pw_branches=branches,
        title=f"Verify {rel[0]} | {rel[1]} | {rel[2]}",
        params_text="ours (DE solver) — compare result/.../overlay_omega1_omega2.png",
    )
    return _row_from_pw(rel, pw)


def _write_summary(rows):
    os.makedirs(os.path.dirname(SUMMARY), exist_ok=True)
    existing = {}
    if os.path.exists(SUMMARY):
        with open(SUMMARY) as f:
            for r in csv.DictReader(f):
                existing[(r["chi_dir"], r["om_dir"], r["chibb_dir"])] = r
    for row in rows:
        existing[(row[0], row[1], row[2])] = dict(zip(SUMMARY_COLS, row))
    with open(SUMMARY, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=SUMMARY_COLS)
        w.writeheader()
        for k in sorted(existing):
            w.writerow(existing[k])


def _rebuild_summary():
    """Rebuild SUMMARY.csv from every out/verify/.../pw_line.csv (post-parallel merge)."""
    root = os.path.join(ROOT, "out", "verify")
    rows = []
    for chi_dir in sorted(os.listdir(root)) if os.path.isdir(root) else []:
        cp = os.path.join(root, chi_dir)
        if not os.path.isdir(cp):
            continue
        for om_dir in sorted(os.listdir(cp)):
            op = os.path.join(cp, om_dir)
            if not os.path.isdir(op):
                continue
            for chibb_dir in sorted(os.listdir(op)):
                pw_csv = os.path.join(op, chibb_dir, "pw_line.csv")
                if os.path.exists(pw_csv):
                    rows.append(_row_from_pw((chi_dir, om_dir, chibb_dir), _load_pw(pw_csv)))
    if os.path.exists(SUMMARY):
        os.remove(SUMMARY)
    _write_summary(rows)
    print(f"rebuilt {SUMMARY} from {len(rows)} cases")


def main(argv):
    if "--rebuild-summary" in argv:
        _rebuild_summary()
        return
    skip_existing = "--skip-existing" in argv
    no_summary = "--no-summary" in argv
    argv = [a for a in argv if a not in ("--skip-existing", "--no-summary")]
    if argv and argv[0] == "--all":
        limit = int(argv[1]) if len(argv) > 1 else None
        rows = []
        for i, (rel, chi, surf) in enumerate(cases.iter_cases(os.path.join(ROOT, "result"))):
            if limit is not None and i >= limit:
                break
            row = verify_case(rel, chi, surf, skip_existing)
            rows.append(row)
            print(f"[{i}] {'/'.join(rel)} -> {row[3]} pts ({row[8]})")
        if not no_summary:
            _write_summary(rows)
    else:
        rel = tuple(argv) if len(argv) == 3 else (
            "chi12_0p0__chi13_2p8__chi23_0p0",
            "om1_m0p30__om2_m0p30",
            "chibb11_0p00__chibb22_0p00__chibb12_0p00",
        )
        chi, surf = cases.parse_case(*rel)
        row = verify_case(rel, chi, surf, skip_existing)
        if not no_summary:
            _write_summary([row])
        print(f"{'/'.join(rel)} -> {row[3]} pw points ({row[8]})")
        print(f"wrote {cases.verify_dir(rel, ROOT)}/" + ("" if no_summary else "  and SUMMARY.csv"))


if __name__ == "__main__":
    main(sys.argv[1:])
