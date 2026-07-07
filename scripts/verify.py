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


def _seed_scan(chi, surf, phi2, f_left, f_right):
    """Multi-start scan at fixed phi2; on the first thin/thick gamma sign change return
    (phi1_star, thin_warm, thick_warm, sep), else None. Only used to SEED the
    continuation with one well-separated point on the line (its converged thin/thick
    profiles), not to trace the whole line."""
    bl = f_left(phi2)
    dense = [(f_right(phi2), phi2), (0.97 * f_right(phi2), phi2)]
    grid = bl + np.arange(-0.03, 0.012, 0.0025)
    grid = grid[grid > 1e-3]
    prev = None
    prev_states = None
    for phi1 in grid:
        st = E.find_states(chi, (phi1, phi2), surf, KAPPA, dense_seeds=dense,
                           warm=prev_states)
        if len(st) >= 2:
            prev_states = ((st[0].sol_x, st[0].sol_y), (st[-1].sol_x, st[-1].sol_y))
            d = st[0].gamma - st[-1].gamma
            if prev is not None and prev[1] * d < 0:
                p0, d0 = prev
                phi1_star = p0 + (phi1 - p0) * (0 - d0) / (d - d0)
                sep = abs(st[0].phi[0][0] - st[-1].phi[0][0])
                return (phi1_star, (st[0].sol_x, st[0].sol_y),
                        (st[-1].sol_x, st[-1].sol_y), sep)
            prev = (phi1, d)
        else:
            prev = None
            prev_states = None
    return None


def _continue(chi, surf, seed, apex, direction, dphi2=0.005, phi2_floor=0.001,
              phi2_ceiling=0.15, max_steps=200, min_dphi2=6.25e-4, tag=None):
    """March the PW line from seed=(phi2, phi1_star, thin_warm, thick_warm) in
    `direction` (+1 up / -1 down), tracking the two branches with E.pw_point until they
    merge (terminus) or we leave the dilute flank. Returns [(phi1_star, phi2)].

    Adaptive step: a corrector failure does NOT immediately end the line — the phi2 step
    is halved and retried (phi2 not advanced), down to min_dphi2, and only a failure at
    that finest step is a real terminus. This is unconditional (every case uses it): on an
    easy line the first attempt succeeds and no halving happens, so behaviour is unchanged;
    on a line whose branches nearly coincide at the ends it recovers the points a coarse
    fixed step would skip over. After a success the step is restored toward dphi2 (doubled,
    capped) so the interior stays fast. Progress is streamed per accepted point / halving."""
    phi2_s, phi1_s, tw, kw = seed
    out = []
    p1_prev, p1_prev2 = phi1_s, None
    phi2 = phi2_s
    step = dphi2
    for _ in range(max_steps):
        phi2n = round(phi2 + direction * step, 10)
        if phi2n <= phi2_floor or phi2n >= min(apex, phi2_ceiling):
            break
        guess = p1_prev if p1_prev2 is None else (2.0 * p1_prev - p1_prev2)
        res = E.pw_point(chi, phi2n, surf, KAPPA, tw, kw, guess, max_it=20)
        if res is None:
            if step > min_dphi2 + 1e-12:
                step = max(min_dphi2, step * 0.5)  # shrink and retry from same phi2
                if tag:
                    print(f"[{tag}] {'down' if direction < 0 else 'up'} "
                          f"phi2~{phi2n:.4f} corrector miss -> step={step:.5f}",
                          file=sys.stderr, flush=True)
                continue
            break  # failed even at the finest step -> real line terminus
        phi1_star, thin, thick = res
        out.append((phi1_star, phi2n))
        if tag:
            print(f"[{tag}] {'down' if direction < 0 else 'up'} "
                  f"phi2={phi2n:.4f} phi1*={phi1_star:.4f} step={step:.5f} "
                  f"n={len(out)}", file=sys.stderr, flush=True)
        tw = (thin.sol_x, thin.sol_y)
        kw = (thick.sol_x, thick.sol_y)
        p1_prev2, p1_prev = p1_prev, phi1_star
        phi2 = phi2n
        if step < dphi2:
            step = min(dphi2, step * 2.0)  # recover toward the nominal step in the interior
    return out


def prewetting_line(chi, surf, binodal, progress=None, max_lines=None):
    """(phi1*, phi2) pre-wetting points (gamma_thin == gamma_thick) along the dilute
    flank, by branch CONTINUATION.

    A short multi-start seed pass over a middle phi2 band finds a well-separated
    thin/thick pair on the line; the continuation then tracks that pair outward in both
    phi2 directions (equilibrium.pw_point) down toward the binary limit and up to the
    surface-critical merge, so the line's endpoints — which per-line multi-start drops
    where the two branches nearly coincide — are recovered.

    progress truthy -> stream live per-step progress to stderr (each seed scan, each
    continuation point). max_lines caps the seed band (a dry-run smoke test; the
    continuation itself still traces the full line).
    """
    f_left, f_right, apex = _branches(binodal)
    band = np.arange(0.03, min(0.09, 0.85 * apex), 0.01)
    if max_lines:
        band = band[:max_lines]
    tag = progress if isinstance(progress, str) else "pw"
    t0 = time.perf_counter()
    seeds = []
    for i, phi2 in enumerate(band):
        r = _seed_scan(chi, surf, float(phi2), f_left, f_right)
        if r:
            phi1_star, tw, kw, sep = r
            seeds.append((float(phi2), phi1_star, tw, kw, sep))
        if progress:  # seed scan is slow (multi-start BVP per phi2) -> stream each one
            hit = f"phi1*={r[0]:.4f} sep={r[3]:.3f}" if r else "no pair"
            print(f"[{tag}] seed {i+1}/{len(band)} phi2={phi2:.3f} {hit} "
                  f"elapsed={time.perf_counter()-t0:.0f}s", file=sys.stderr, flush=True)
    if not seeds:
        if progress:
            print(f"[{tag}] no seed on band {band[0]:.2f}..{band[-1]:.2f}",
                  file=sys.stderr, flush=True)
        return np.array([])
    best = max(seeds, key=lambda s: s[4])  # widest branch separation = sturdiest seed
    seed = (best[0], best[1], best[2], best[3])
    down = _continue(chi, surf, seed, apex, -1, tag=tag if progress else None)
    up = _continue(chi, surf, seed, apex, +1, tag=tag if progress else None)
    pw = sorted([(best[1], best[0])] + down + up, key=lambda p: p[1])
    if progress:
        el = time.perf_counter() - t0
        print(f"[{tag}] seed@phi2={best[0]:.3f} sep={best[4]:.3f}  down+{len(down)} "
              f"up+{len(up)}  n={len(pw)}  phi2=[{pw[0][1]:.3f},{pw[-1][1]:.3f}]  "
              f"elapsed={el:.0f}s", file=sys.stderr, flush=True)
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
