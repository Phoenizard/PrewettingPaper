#!/usr/bin/env python
"""Verify result/ pre-wetting cases with our equilibrium-DE solver (compute only).

Compute and plotting are decoupled: this script produces DATA only, never a
figure — so parallel workers never import matplotlib (no font-cache race) and
compute never contends with rendering. Turn the data into figures separately
with scripts/plot.py.

Output mirrors result/ 1:1 under out/<chi_dir>/<om_dir>/<chibb_dir>/:
  pw_line.csv   computed (phi1_inf, phi2_inf) transition points
Top-level out/SUMMARY.csv records every case run.

Usage:
  conda run -n numenv python scripts/verify.py                       # T-a template case
  conda run -n numenv python scripts/verify.py <chi_dir> <om_dir> <chibb_dir>
  conda run -n numenv python scripts/verify.py --all [N]             # first N (or all) cases
  conda run -n numenv python scripts/verify.py --rebuild-summary     # rebuild SUMMARY.csv from outputs
  --skip-existing   reuse cases already having pw_line.csv (resume)
  --no-summary      skip SUMMARY.csv write (for parallel workers; rebuild after)
"""
import csv
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "src"))

import numpy as np  # noqa: E402
import thermo as T  # noqa: E402
import binodal as B  # noqa: E402
import equilibrium as E  # noqa: E402
import cases  # noqa: E402

KAPPA = T.Kappa(1.0, 1.0)
SUMMARY = os.path.join(ROOT, "out", "SUMMARY.csv")
SUMMARY_COLS = ["chi_dir", "om_dir", "chibb_dir", "n_pw",
                "phi1_min", "phi1_max", "phi2_min", "phi2_max", "status"]


def _diag(**kw):
    """Emit one [PWDIAG] JSONL line to stderr when PW_DIAG is set (diagnosis only,
    no behaviour change, parallel-safe: stderr, no files)."""
    if os.environ.get("PW_DIAG"):
        print("[PWDIAG] " + json.dumps(kw), file=sys.stderr, flush=True)


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


# ---------------------------------------------------------------------------
# Pre-wetting line by DUAL-SCAN + CROSSING (replaces the old branch continuation).
#
# The continuation (seed multi-start + pw_point secant tracking two live branches)
# was fragile at the low-phi2 end: a third (middle) surface state appears, the thick
# branch is lost/mis-tracked, and the line falsely terminates. Instead we build the two
# gamma(phi_scan) curves INDEPENDENTLY — a thin branch (forward warm-start scan) and a
# thick branch (backward warm-start scan) — then read off the pre-wetting point as the
# gamma_thin == gamma_thick crossing, gated by a real film-thickness step (adsorption
# gap). A dropped (non-converged) point is just a gap; the crossing is still recovered
# from finite neighbours, so losing a branch point no longer kills the line.
#
# This mirrors the reference implementation's _run_hysteresis_scan +
# _extract_prewetting_crossings, re-implemented in our own code on solve_profile.
# ---------------------------------------------------------------------------

_trapz = getattr(np, "trapezoid", getattr(np, "trapz", None))  # numpy>=2 renamed trapz


def _adsorption(prof, res):
    """Total Gibbs adsorption cs = INT[(phi1-phi1_inf) + (phi2-phi2_inf)] dz (film size)."""
    return float(_trapz((prof.phi[0] - res[0]) + (prof.phi[1] - res[1]), prof.z))


def _branch_scan(chi, surf, fixed, scan_vals, axis, mode, f_right, progress=None):
    """Build ONE surface branch along a reservoir scan, warm-started point to point.

    axis="phi1": fixed=phi2, scan_vals are phi1_inf (reservoir phi1 swept).
    axis="phi2": fixed=phi1, scan_vals are phi2_inf.
    mode="thin"  scans scan_vals in the given (ascending) order; the cold seed is a small
                 wall enrichment (near-reservoir) — the thin basin.
    mode="thick" scans in REVERSE (descending); the cold seed is a phi1-rich wall plateau
                 (dense flank) — the thick basin. Reverse + high-phi1 start keeps the thick
                 attractor, exactly where independent multi-start loses it.

    Each point warm-starts from the previous converged profile (solve_profile(warm=...));
    a non-converged point is recorded as NaN and the previous warm seed is kept. Returns
    (gamma[], cs[]) aligned to scan_vals ASCENDING (thick is reversed back before return),
    NaN where the solve failed. Streams one line per scan point when progress is a tag."""
    order = list(range(len(scan_vals)))
    if mode == "thick":
        order = order[::-1]
    gamma = [float("nan")] * len(scan_vals)
    cs = [float("nan")] * len(scan_vals)
    warm = None
    for k in order:
        sv = float(scan_vals[k])
        res = (sv, fixed) if axis == "phi1" else (fixed, sv)
        if warm is None:  # cold seed for this branch's basin
            if mode == "thin":
                seed = (min(0.98, res[0] + 0.02), res[1])
                w = 2.0
            else:  # thick: phi1-rich wall plateau from the dense flank
                fr = f_right(res[1]) if axis == "phi1" else f_right(fixed)
                seed = (max(res[0] + 0.05, min(0.95, fr)), res[1])
                w = 4.0
            p = E.solve_profile(chi, res, surf, KAPPA, seed=seed, w=w)
        else:
            p = E.solve_profile(chi, res, surf, KAPPA, warm=warm)
        tag_hit = "None"
        if p is not None:
            gamma[k] = p.gamma
            cs[k] = _adsorption(p, res)
            warm = (p.sol_x, p.sol_y)
            tag_hit = f"g={p.gamma:+.5f} cs={cs[k]:.3f} wall={p.phi[0][0]:.3f}"
        if progress:
            print(f"[{progress}]   scan-{axis} {mode} fixed={fixed:.4f} "
                  f"{'phi1' if axis=='phi1' else 'phi2'}={sv:.4f} -> {tag_hit}",
                  file=sys.stderr, flush=True)
    return np.asarray(gamma), np.asarray(cs)


def _crossings(scan_vals, g_thin, g_thick, cs_thin, cs_thick,
               cs_threshold=0.1, min_pts=2, terminal_relax=0.67):
    """Pre-wetting crossings = sign changes of (g_thin - g_thick) where both branches are
    finite and the film-thickness gap |cs_thick - cs_thin| exceeds cs_threshold (a real
    thin->thick step, not spinodal noise). Linear-interpolate each sign change to the exact
    scan value. Falls back to a relaxed-gate boundary crossing (the line's endpoint) when no
    interior crossing is gated in. Returns a list of scan values (phi1* or phi2*)."""
    scan_vals = np.asarray(scan_vals, float)
    diff = g_thin - g_thick
    gap = np.abs(cs_thick - cs_thin)
    eligible = np.isfinite(diff) & np.isfinite(gap)
    valid = eligible & (gap > cs_threshold)
    out = []
    if int(valid.sum()) >= min_pts:
        idx = np.where(eligible)[0]
        splits = np.where(np.diff(idx) > 1)[0] + 1
        for region in np.split(idx, splits):
            rv = region[valid[region]]
            if rv.size < 2:
                continue
            sign_changes = np.where(np.diff(np.sign(diff[rv])) != 0)[0]
            for k in sign_changes.tolist():
                i, j = int(rv[k]), int(rv[k + 1])
                out.append(_interp_zero(scan_vals[i], scan_vals[j], diff[i], diff[j]))
    if not out and int(valid.sum()) > 0:  # terminal-boundary fallback (line endpoint)
        near = eligible & (gap > cs_threshold * terminal_relax)
        y1, y2 = diff[:-1], diff[1:]
        crossed = np.isfinite(y1) & np.isfinite(y2) & ((y1 * y2) <= 0.0)
        bnd = (valid[:-1] != valid[1:]) & near[:-1] & near[1:] & crossed
        for i in np.where(bnd)[0].tolist():
            out.append(_interp_zero(scan_vals[i], scan_vals[i + 1], y1[i], y2[i]))
    return out


def _interp_zero(x0, x1, y0, y1):
    """Linear-interpolate the zero crossing of y(x) between (x0,y0) and (x1,y1)."""
    if y1 == y0:
        return float(x0)
    return float(x0 - y0 * (x1 - x0) / (y1 - y0))


def prewetting_line(chi, surf, binodal, progress=None, max_lines=None,
                    cs_threshold=0.1, min_pts=2):
    """(phi1*, phi2) pre-wetting points (gamma_thin == gamma_thick) by DUAL-SCAN crossing.

    Two scan axes, unioned:
      fix_phi2 scan phi1 -> phi1*(phi2), the shallow parts of the line;
      fix_phi1 scan phi2 -> phi2*(phi1), the near-vertical low-phi2 tail toward phi2->0.
    Each axis builds a thin (forward) and thick (backward) branch via _branch_scan, then
    reads crossings via _crossings. progress truthy -> stream one line per scan point.
    max_lines caps the number of fixed values per axis (dry-run smoke test)."""
    f_left, f_right, apex = _branches(binodal)
    tag = progress if isinstance(progress, str) else "pw"
    t0 = time.perf_counter()

    def phi1_window(phi2):
        bl = f_left(phi2)
        lo = max(1e-3, bl - 0.03)
        hi = min(0.95, bl + 0.06)
        return np.arange(lo, hi + 1e-9, 0.0025)

    pts = []
    # --- axis A: fix phi2, scan phi1 ---
    phi2_band = np.arange(0.005, min(0.09, 0.95 * apex), 0.005)
    if max_lines:
        phi2_band = phi2_band[:max_lines]
    for i, phi2 in enumerate(phi2_band):
        phi2 = float(phi2)
        grid = phi1_window(phi2)
        if progress:
            print(f"[{tag}] axisA fix_phi2 {i+1}/{len(phi2_band)} phi2={phi2:.4f} "
                  f"scan {len(grid)} phi1 in [{grid[0]:.4f},{grid[-1]:.4f}]",
                  file=sys.stderr, flush=True)
        gt, ct = _branch_scan(chi, surf, phi2, grid, "phi1", "thin", f_right,
                              progress=tag if progress else None)
        gk, ck = _branch_scan(chi, surf, phi2, grid, "phi1", "thick", f_right,
                              progress=tag if progress else None)
        cr = _crossings(grid, gt, gk, ct, ck, cs_threshold, min_pts)
        for p1 in cr:
            pts.append((float(p1), phi2))
        if progress:
            print(f"[{tag}] axisA phi2={phi2:.4f} -> {len(cr)} crossing(s) "
                  f"phi1*={[round(x,4) for x in cr]}  elapsed={time.perf_counter()-t0:.0f}s",
                  file=sys.stderr, flush=True)

    # --- axis B: fix phi1, scan phi2 (recovers the near-vertical low-phi2 tail) ---
    phi1_lo = f_left(0.005)
    phi1_hi = f_left(min(0.08, 0.9 * apex))
    phi1_band = np.arange(min(phi1_lo, phi1_hi) - 0.005,
                          max(phi1_lo, phi1_hi) + 0.02, 0.005)
    if max_lines:
        phi1_band = phi1_band[:max_lines]
    for i, phi1 in enumerate(phi1_band):
        phi1 = float(phi1)
        grid = np.arange(0.002, min(0.09, 0.95 * apex) + 1e-9, 0.0025)
        if progress:
            print(f"[{tag}] axisB fix_phi1 {i+1}/{len(phi1_band)} phi1={phi1:.4f} "
                  f"scan {len(grid)} phi2 in [{grid[0]:.4f},{grid[-1]:.4f}]",
                  file=sys.stderr, flush=True)
        gt, ct = _branch_scan(chi, surf, phi1, grid, "phi2", "thin", f_right,
                              progress=tag if progress else None)
        gk, ck = _branch_scan(chi, surf, phi1, grid, "phi2", "thick", f_right,
                              progress=tag if progress else None)
        cr = _crossings(grid, gt, gk, ct, ck, cs_threshold, min_pts)
        for p2 in cr:
            pts.append((phi1, float(p2)))
        if progress:
            print(f"[{tag}] axisB phi1={phi1:.4f} -> {len(cr)} crossing(s) "
                  f"phi2*={[round(x,4) for x in cr]}  elapsed={time.perf_counter()-t0:.0f}s",
                  file=sys.stderr, flush=True)

    if not pts:
        if progress:
            print(f"[{tag}] no crossings on either axis", file=sys.stderr, flush=True)
        return np.array([])
    pw = _dedup_points(pts)
    if progress:
        el = time.perf_counter() - t0
        print(f"[{tag}] dual-scan done: {len(pw)} pw points  "
              f"phi2=[{pw[:,1].min():.4f},{pw[:,1].max():.4f}]  elapsed={el:.0f}s",
              file=sys.stderr, flush=True)
    return pw


def _dedup_points(pts, tol=2e-3):
    """Merge near-duplicate (phi1, phi2) points from the two scan axes; sort by phi2."""
    arr = np.asarray(sorted(pts, key=lambda p: (p[1], p[0])), float)
    keep = []
    for p in arr:
        if all(np.hypot(p[0] - q[0], p[1] - q[1]) > tol for q in keep):
            keep.append(p)
    return np.asarray(sorted(keep, key=lambda p: p[1]), float)


def _row_from_pw(rel, pw):
    if len(pw):
        return [*rel, len(pw), f"{pw[:,0].min():.4f}", f"{pw[:,0].max():.4f}",
                f"{pw[:,1].min():.4f}", f"{pw[:,1].max():.4f}", "ok"]
    return [*rel, 0, "", "", "", "", "no_pw"]


def _load_pw(pw_csv):
    """Read a pw_line.csv back to an (n,2) array; empty file -> (0,2)."""
    raw = np.loadtxt(pw_csv, delimiter=",", skiprows=1, ndmin=2)
    return raw if raw.size and raw.shape[1] == 2 else np.empty((0, 2))


def _out_dir(rel):
    # VERIFY_OUT overrides the output root (e.g. tmp/verify_opt) so a throwaway
    # side-by-side run does not clobber the baseline results under out/.
    base = os.environ.get("VERIFY_OUT")
    return os.path.join(base, *rel) if base else cases.verify_dir(rel, ROOT)


def verify_case(rel, chi, surf, skip_existing=False):
    out_dir = _out_dir(rel)
    pw_csv = os.path.join(out_dir, "pw_line.csv")
    # Resume keys on pw_line.csv (the only, and expensive, artifact here).
    if skip_existing and os.path.exists(pw_csv):
        return _row_from_pw(rel, _load_pw(pw_csv))

    os.makedirs(out_dir, exist_ok=True)
    binodal = B.binodal_from_hull(chi)
    # Progress streams by DEFAULT (per-step continuation lines to stderr): a long solve
    # must show live progress, not sit silent until done. VERIFY_PROGRESS=0 opts out
    # (e.g. a big parallel sweep where per-case interleaving would be noise).
    prog = None if os.environ.get("VERIFY_PROGRESS") == "0" else "/".join(rel)
    max_lines = int(os.environ.get("PW_MAX_PHI2", "0")) or None  # dry-run cap
    pw = prewetting_line(chi, surf, binodal, progress=prog, max_lines=max_lines)

    np.savetxt(pw_csv, pw if len(pw) else np.empty((0, 2)),
               delimiter=",", header="phi1_inf,phi2_inf", comments="")
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
    """Rebuild SUMMARY.csv from every out/<chi>/<om>/<chibb>/pw_line.csv (post-parallel merge)."""
    root = os.path.join(ROOT, "out")
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
        print(f"wrote {_out_dir(rel)}/" + ("" if no_summary else "  and SUMMARY.csv"))


if __name__ == "__main__":
    main(sys.argv[1:])
