"""T-f (chi12 = -8.5, chi13 = chi23 = 0) binodal critical points and a gapless loop.

The critical point of a bulk binodal is where the two coexisting compositions
merge and the tie line joining them shrinks to zero length. It is a property of
the bulk free energy alone -- omega does not enter f_b -- so it does not move
when the wall affinity changes.

Two independent routes, both run here so they can check each other.

Route 1, direct root find. A critical point satisfies two conditions on f_b:
the Hessian is singular,

    D = f11 f22 - f12^2 = 0,

and the third derivative along the Hessian's null direction v = (f12, -f11)
vanishes,

    C = f111 v1^3 + 3 f112 v1^2 v2 + 3 f122 v1 v2^2 + f222 v2^3 = 0.

Route 2, tie-line continuation. Two phases a, b coexist when both exchange
chemical potentials and the osmotic pressure match. That is three equations in
four unknowns, so the solutions form a one-parameter family -- the tie lines.
Continuing along the family from the symmetric tie line until its length goes to
zero lands on the same critical point, and the endpoints traced on the way out
are the binodal loop itself, with no detection heuristic and therefore no gaps.

The production binodal in src/bulk.py is left untouched: it builds a convex hull
over roughly five million grid points, which is minutes of work, and it accepts
a tie line only when consecutive hull vertex indices jump by more than a
threshold, which drops short tie lines. This script is analysis-only.

Usage:
  python scripts/tf_critical_point.py [--out-dir DIR]
"""

import argparse
import csv
import sys
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.optimize import root

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import model


@dataclass
class Bulk:
    """The fields of params.Physical that the bulk free energy actually uses.

    Declared here rather than imported so this script does not pull in the yaml
    config loader: nothing about a critical point depends on solver or scan
    settings, and the values below are the T-f row of the case table.
    """
    chi_12: float = -8.5
    chi_13: float = 0.0
    chi_23: float = 0.0
    n1: float = 1.0
    n2: float = 1.0
    n3: float = 1.0

SEED_GRID = 40          # per axis, for the route-1 seed sweep
ROOT_TOL = 1e-10        # max |residual| accepted for a route-1 root
DEDUPE_TOL = 1e-6       # two roots closer than this are the same root
FD_STEP = 1e-4          # finite-difference step for the self-check
TIE_MIN_LEN = 1e-4      # continuation stops here; the midpoint is the estimate
MAX_TIE_STEPS = 4000    # hard cap per branch, so a stalled step cannot spin
STALL_FRAC = 0.999      # a step that shortens the tie line by less than this
                        # much has stalled: the solver returns a length a hair
                        # above the target, so a bare "length > target" test
                        # never terminates.


# ------------------------------------------------------------ bulk derivatives

def third_derivatives(phi1, phi2, p):
    """f111, f112 (= f122), f222. The chi terms are quadratic and drop out."""
    p1, p2, p3 = model.clip_composition(phi1, phi2)
    f111 = -1.0 / (p.n1 * p1**2) + 1.0 / (p.n3 * p3**2)
    f112 = 1.0 / (p.n3 * p3**2)
    f222 = -1.0 / (p.n2 * p2**2) + 1.0 / (p.n3 * p3**2)
    return f111, f112, f222


def hessian_det(phi1, phi2, p):
    f11, f12, f22 = model.hessian(phi1, phi2, p)
    return f11 * f22 - f12**2


def null_direction(phi1, phi2, p):
    f11, f12, f22 = model.hessian(phi1, phi2, p)
    return f12, -f11


def third_along(phi1, phi2, p, v):
    """Third directional derivative of f_b along v = (v1, v2)."""
    f111, f112, f222 = third_derivatives(phi1, phi2, p)
    v1, v2 = v
    return (f111 * v1**3 + 3.0 * f112 * v1**2 * v2
            + 3.0 * f112 * v1 * v2**2 + f222 * v2**3)


def critical_residual(x, p):
    phi1, phi2 = x
    return np.array([hessian_det(phi1, phi2, p),
                     third_along(phi1, phi2, p, null_direction(phi1, phi2, p))])


# --------------------------------------------------------------- route 1

def find_critical_points(p):
    span = np.linspace(0.02, 0.96, SEED_GRID)
    seeds = [(a, b) for a in span for b in span if a + b < 0.98]

    print(f"  seeds: {len(seeds)}", flush=True)
    found = []
    for k, seed in enumerate(seeds):
        if (k + 1) % 200 == 0 or (k + 1) == len(seeds):
            print(f"  seed {k + 1}/{len(seeds)}, {len(found)} distinct root(s)",
                  flush=True)
        sol = root(critical_residual, np.array(seed), args=(p,), method="hybr")
        if not sol.success:
            continue
        phi1, phi2 = sol.x
        if not (phi1 > 1e-6 and phi2 > 1e-6 and phi1 + phi2 < 1.0 - 1e-6):
            continue
        res = critical_residual(sol.x, p)
        # C carries the cube of |v|, so scale it out before judging the residual.
        v = np.array(null_direction(phi1, phi2, p))
        scale = max(np.linalg.norm(v), 1e-30)
        if abs(res[0]) > ROOT_TOL or abs(res[1]) / scale**3 > ROOT_TOL:
            continue
        f11, _, f22 = model.hessian(phi1, phi2, p)
        if f11 <= 0.0 or f22 <= 0.0:      # not a stability boundary
            continue
        if not any(abs(phi1 - q[0]) < DEDUPE_TOL and abs(phi2 - q[1]) < DEDUPE_TOL
                   for q in found):
            found.append((float(phi1), float(phi2)))
    return sorted(found)


# --------------------------------------------------------------- route 2

def osmotic(phi1, phi2, p):
    """f - mu1 phi1 - mu2 phi2, the quantity that is equal across a tie line."""
    mu1, mu2 = model.chemical_potential(phi1, phi2, p)
    return model.free_energy(phi1, phi2, p) - mu1 * phi1 - mu2 * phi2


def coexistence_residual(x, p, constraint):
    a1, a2, b1, b2 = x
    mu1a, mu2a = model.chemical_potential(a1, a2, p)
    mu1b, mu2b = model.chemical_potential(b1, b2, p)
    return np.array([mu1a - mu1b,
                     mu2a - mu2b,
                     osmotic(a1, a2, p) - osmotic(b1, b2, p),
                     constraint(x)])


def symmetric_tie_line(p, n_grid=20001):
    """The tie line on the phi1 = phi2 diagonal, from a 1-D common tangent.

    On the diagonal both solutes are present in equal amount, so the mixture
    behaves as a pseudo-binary in the total solute fraction s. The two contact
    points of the lower convex envelope of f(s) are the tie-line ends.
    """
    s = np.linspace(1e-6, 1.0 - 1e-6, n_grid)
    f = model.free_energy(0.5 * s, 0.5 * s, p)

    # Lower convex hull of the 1-D curve by a monotone chain on the slopes.
    hull = []
    for i in range(len(s)):
        while len(hull) >= 2:
            j, k = hull[-2], hull[-1]
            if (f[k] - f[j]) * (s[i] - s[j]) >= (f[i] - f[j]) * (s[k] - s[j]):
                hull.pop()
            else:
                break
        hull.append(i)
    gaps = [(s[hull[i + 1]] - s[hull[i]], hull[i], hull[i + 1])
            for i in range(len(hull) - 1)]
    width, i, j = max(gaps)
    if width < 10.0 * (s[1] - s[0]):
        raise RuntimeError("no two-phase gap on the diagonal")

    def constraint(x):
        return (x[0] - x[1]) + (x[2] - x[3])      # midpoint offset = 0

    x0 = np.array([0.5 * s[i], 0.5 * s[i], 0.5 * s[j], 0.5 * s[j]])
    sol = root(coexistence_residual, x0, args=(p, constraint), method="hybr")
    if not sol.success:
        raise RuntimeError(f"symmetric tie line failed: {sol.message}")
    return sol.x


def trace_branch(p, x_sym, direction):
    """Continue the tie-line family away from the symmetric tie line.

    Two stages. First the midpoint offset from the diagonal is stepped, which
    breaks the mirror symmetry cleanly. Once the tie line has shortened enough
    the parameter switches to the tie-line length, which stays well conditioned
    as the two ends approach each other.
    """
    lines = []
    x = x_sym.copy()
    len0 = float(np.hypot(x[0] - x[2], x[1] - x[3]))

    def record(x):
        length = float(np.hypot(x[0] - x[2], x[1] - x[3]))
        lines.append((float(x[0]), float(x[1]), float(x[2]), float(x[3]), length))
        if len(lines) % 20 == 0:
            print(f"    branch {direction:+.0f}: {len(lines)} tie lines, "
                  f"length {length:.3e}", flush=True)
        return length

    record(x)

    # stage 1: step the midpoint offset t
    t, dt = 0.0, direction * 0.002
    while abs(dt) > 1e-6 and len(lines) < MAX_TIE_STEPS:
        t_new = t + dt

        def constraint(y, t_new=t_new):
            return (y[0] - y[1]) + (y[2] - y[3]) - 2.0 * t_new

        sol = root(coexistence_residual, x, args=(p, constraint), method="hybr")
        ok = sol.success and np.max(np.abs(
            coexistence_residual(sol.x, p, constraint))) < 1e-10
        if ok:
            x, t = sol.x, t_new
            if record(x) < 0.35 * len0:
                break
        else:
            dt *= 0.5

    # stage 2: step the tie-line length down
    length = float(np.hypot(x[0] - x[2], x[1] - x[3]))
    while length > TIE_MIN_LEN and len(lines) < MAX_TIE_STEPS:
        target = max(length * 0.7, TIE_MIN_LEN)

        def constraint(y, target=target):
            return (y[0] - y[2])**2 + (y[1] - y[3])**2 - target**2

        sol = root(coexistence_residual, x, args=(p, constraint), method="hybr")
        ok = sol.success and np.max(np.abs(
            coexistence_residual(sol.x, p, constraint))) < 1e-10
        if not ok:
            break
        x = sol.x
        new_length = record(x)
        if new_length > length * STALL_FRAC:
            break
        length = new_length

    tip = (0.5 * (x[0] + x[2]), 0.5 * (x[1] + x[3]), length)
    return lines, tip


# ----------------------------------------------------------------- self-checks

def fd_checks(phi1, phi2, p, h=FD_STEP):
    """Analytic D and C against central differences of model.free_energy."""
    def f(a, b):
        return float(model.free_energy(a, b, p))

    f11 = (f(phi1 + h, phi2) - 2 * f(phi1, phi2) + f(phi1 - h, phi2)) / h**2
    f22 = (f(phi1, phi2 + h) - 2 * f(phi1, phi2) + f(phi1, phi2 - h)) / h**2
    f12 = (f(phi1 + h, phi2 + h) - f(phi1 + h, phi2 - h)
           - f(phi1 - h, phi2 + h) + f(phi1 - h, phi2 - h)) / (4 * h**2)
    d_fd = f11 * f22 - f12**2

    v = np.array(null_direction(phi1, phi2, p))
    v = v / np.linalg.norm(v)
    hs = 10.0 * h

    def g(k):
        return f(phi1 + k * hs * v[0], phi2 + k * hs * v[1])

    c_fd = (g(2) - 2 * g(1) + 2 * g(-1) - g(-2)) / (2 * hs**3)
    c_an = third_along(phi1, phi2, p, v)
    return d_fd, hessian_det(phi1, phi2, p), c_fd, c_an


AUDIT_POINTS = [(0.20, 0.30), (0.35, 0.25), (0.15, 0.15), (0.40, 0.40)]


def fd_audit(p):
    """Check the algebra where it can actually be checked.

    At a critical point D and C are zero by construction, so comparing them
    against finite differences there tests nothing -- both sides are noise. The
    comparison is only meaningful at ordinary compositions, where both are large.
    """
    print("\n[check] analytic vs finite difference at ordinary compositions",
          flush=True)
    worst = 0.0
    for phi1, phi2 in AUDIT_POINTS:
        d_fd, d_an, c_fd, c_an = fd_checks(phi1, phi2, p)
        rel_d = abs(d_fd - d_an) / max(abs(d_an), 1e-30)
        rel_c = abs(c_fd - c_an) / max(abs(c_an), 1e-30)
        worst = max(worst, rel_d, rel_c)
        print(f"  ({phi1:.3f}, {phi2:.3f})  D fd={d_fd:+.6f} an={d_an:+.6f} "
              f"rel={rel_d:.2e} | C fd={c_fd:+.4f} an={c_an:+.4f} "
              f"rel={rel_c:.2e}", flush=True)
    print(f"  worst relative difference {worst:.2e}", flush=True)
    return worst


# ------------------------------------------------------------------- driver

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir",
                    default=str(ROOT / "out/analysis/Tf_omega_cross/critical"))
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    if not out_dir.is_dir():
        raise SystemExit(f"directory does not exist: {out_dir}")

    t_start = time.time()
    p = Bulk()
    print(f"case chi=({p.chi_12}, {p.chi_13}, {p.chi_23}) "
          f"n=({p.n1}, {p.n2}, {p.n3})", flush=True)

    # route 1
    t0 = time.time()
    roots = find_critical_points(p)
    print(f"\n[route 1] direct root find of D = 0, C = 0 "
          f"({time.time() - t0:.1f}s)", flush=True)
    for phi1, phi2 in roots:
        print(f"  critical point ({phi1:.8f}, {phi2:.8f})", flush=True)
    if not roots:
        raise SystemExit("route 1 found no critical point")

    # route 2
    t0 = time.time()
    x_sym = symmetric_tie_line(p)
    print(f"\n[route 2] tie-line continuation ({time.time() - t0:.1f}s)",
          flush=True)
    print(f"  symmetric tie line ({x_sym[0]:.6f}, {x_sym[1]:.6f}) -- "
          f"({x_sym[2]:.6f}, {x_sym[3]:.6f})", flush=True)
    branches, tips = [], []
    for direction in (+1.0, -1.0):
        lines, tip = trace_branch(p, x_sym, direction)
        branches.append(lines)
        tips.append(tip)
        print(f"  branch {direction:+.0f}: {len(lines)} tie lines, "
              f"tip ({tip[0]:.8f}, {tip[1]:.8f}) at length {tip[2]:.2e}",
              flush=True)

    # checks
    fd_audit(p)

    print("\n[check] residual of D and C at each route-1 root", flush=True)
    fd_rows = []
    for phi1, phi2 in roots:
        d_fd, d_an, c_fd, c_an = fd_checks(phi1, phi2, p)
        fd_rows.append((d_fd, d_an, c_fd, c_an))
        print(f"  ({phi1:.6f}, {phi2:.6f})  D={d_an:+.3e}  C={c_an:+.3e} "
              f"(finite difference there is pure noise: "
              f"D={d_fd:+.1e}, C={c_fd:+.1e})", flush=True)

    print("\n[check] route 1 vs route 2", flush=True)
    pairing = []
    for tip in tips:
        near = min(roots, key=lambda q: np.hypot(q[0] - tip[0], q[1] - tip[1]))
        gap = float(np.hypot(near[0] - tip[0], near[1] - tip[1]))
        pairing.append((near, tip, gap))
        print(f"  tip ({tip[0]:.8f}, {tip[1]:.8f}) vs root "
              f"({near[0]:.8f}, {near[1]:.8f})  distance {gap:.2e}", flush=True)

    print("\n[check] mirror symmetry of the route-1 root set", flush=True)
    for phi1, phi2 in roots:
        mate = min(roots, key=lambda q: np.hypot(q[0] - phi2, q[1] - phi1))
        print(f"  ({phi1:.8f}, {phi2:.8f}) mirrors ({mate[0]:.8f}, "
              f"{mate[1]:.8f})  residual "
              f"{np.hypot(mate[0] - phi2, mate[1] - phi1):.2e}", flush=True)

    # output
    with open(out_dir / "critical_points.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["phi1_c", "phi2_c", "D_analytic", "C_analytic",
                    "D_finite_diff", "C_finite_diff", "tie_tip_phi1",
                    "tie_tip_phi2", "tie_tip_length", "route_gap"])
        for (phi1, phi2), (d_fd, d_an, c_fd, c_an) in zip(roots, fd_rows):
            near, tip, gap = min(pairing,
                                 key=lambda r: np.hypot(r[0][0] - phi1,
                                                        r[0][1] - phi2))
            w.writerow([phi1, phi2, d_an, c_an, d_fd, c_fd,
                        tip[0], tip[1], tip[2], gap])

    with open(out_dir / "binodal_tie.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["branch", "phi1_a", "phi2_a", "phi1_b", "phi2_b", "length"])
        for k, lines in enumerate(branches):
            for a1, a2, b1, b2, length in lines:
                w.writerow([k, a1, a2, b1, b2, length])

    n_lines = sum(len(b) for b in branches)
    print(f"\nwrote {out_dir / 'critical_points.csv'} ({len(roots)} points)",
          flush=True)
    print(f"wrote {out_dir / 'binodal_tie.csv'} ({n_lines} tie lines)",
          flush=True)
    print(f"total elapsed {time.time() - t_start:.1f}s", flush=True)


if __name__ == "__main__":
    main()
