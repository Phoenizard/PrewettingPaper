#!/usr/bin/env python
"""Single-case demonstration of the pre-wetting line by branch CONTINUATION.

Hard case chi13=2.8, om1=om2=-0.30: the existing multi-start method only recovers
phi2=0.03-0.06; the reference spans phi2~0-0.09. This traces the thin/thick branch
pair from a well-separated middle seed outward (down toward phi2->0, up toward the
apex) with equilibrium.pw_point, stopping where the branches merge (the true
terminus). For cross-checking it also runs the old multi-start detection over the
middle band -- the two methods must agree there.

Outputs a point table + terminus to stdout, and out/pw_continue_test/overlay.png
(binodal + continuation line + old-method points) to eyeball against
result/.../overlay_omega1_omega2.png. Compute only; run on the server.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "src"))

import numpy as np  # noqa: E402
import thermo as T  # noqa: E402
import binodal as B  # noqa: E402
import equilibrium as E  # noqa: E402
import cases  # noqa: E402

KAPPA = T.Kappa(1.0, 1.0)
REL = ("chi12_0p0__chi13_2p8__chi23_0p0",
       "om1_m0p30__om2_m0p30",
       "chibb11_0p00__chibb22_0p00__chibb12_0p00")


def branches(binodal):
    pts = binodal.points
    ai = int(np.argmax(pts[:, 1]))
    a1, a2 = pts[ai, 0], pts[ai, 1]
    left = pts[pts[:, 0] <= a1]
    right = pts[pts[:, 0] >= a1]
    left = left[np.argsort(left[:, 1])]
    right = right[np.argsort(right[:, 1])]
    f_left = lambda p2: float(np.interp(p2, left[:, 1], left[:, 0]))
    f_right = lambda p2: float(np.interp(p2, right[:, 1], right[:, 0]))
    return f_left, f_right, a2


def seed_scan(chi, surf, phi2, f_left, f_right):
    """Old multi-start scan at fixed phi2; on the first thin/thick gamma sign change
    return (phi1_star, thin_warm, thick_warm, sep), else None."""
    bl = f_left(phi2)
    dense = [(f_right(phi2), phi2), (0.97 * f_right(phi2), phi2)]
    grid = bl + np.arange(-0.03, 0.012, 0.0025)
    grid = grid[grid > 1e-3]
    prev = None
    prev_states = None
    for phi1 in grid:
        st = E.find_states(chi, (phi1, phi2), surf, KAPPA,
                           dense_seeds=dense, warm=prev_states)
        if len(st) >= 2:
            prev_states = ((st[0].sol_x, st[0].sol_y), (st[-1].sol_x, st[-1].sol_y))
            d = st[0].gamma - st[-1].gamma
            if prev is not None and prev[1] * d < 0:
                p0, d0 = prev
                phi1_star = p0 + (phi1 - p0) * (0 - d0) / (d - d0)
                sep = abs(st[0].phi[0][0] - st[-1].phi[0][0])
                thin_w = (st[0].sol_x, st[0].sol_y)
                thick_w = (st[-1].sol_x, st[-1].sol_y)
                return phi1_star, thin_w, thick_w, sep
            prev = (phi1, d)
        else:
            prev = None
            prev_states = None
    return None


def continue_dir(chi, surf, seed, apex, direction, dphi2=0.005, phi2_floor=0.003):
    """March the PW line from `seed`=(phi2, phi1_star, thin_warm, thick_warm) in
    `direction` (+1 up / -1 down) until the branches merge or we leave the flank.
    Returns (points list [(phi1_star, phi2)], stop_reason, terminus_phi2)."""
    phi2_s, phi1_s, tw, kw = seed
    out = []
    p1_prev, p1_prev2 = phi1_s, None
    phi2 = phi2_s
    while True:
        phi2n = round(phi2 + direction * dphi2, 10)
        if phi2n <= phi2_floor:
            return out, "floor", phi2
        if phi2n >= apex:
            return out, "apex", phi2
        guess = p1_prev if p1_prev2 is None else (2.0 * p1_prev - p1_prev2)
        res = E.pw_point(chi, phi2n, surf, KAPPA, tw, kw, guess)
        if res is None:
            return out, "merge", phi2   # last surviving point was at phi2
        phi1_star, thin, thick = res
        out.append((phi1_star, phi2n))
        tw = (thin.sol_x, thin.sol_y)
        kw = (thick.sol_x, thick.sol_y)
        p1_prev2, p1_prev = p1_prev, phi1_star
        phi2 = phi2n


def main():
    chi, surf = cases.parse_case(*REL)
    bd = B.binodal_from_hull(chi)
    f_left, f_right, apex = branches(bd)
    print(f"case {'/'.join(REL)}  apex_phi2={apex:.4f}")

    # Seed over the middle band (also gives old-method points for cross-check).
    band = np.arange(0.02, 0.10, 0.01)
    old_pts, seeds = [], []
    for phi2 in band:
        r = seed_scan(chi, surf, float(phi2), f_left, f_right)
        if r:
            phi1_star, tw, kw, sep = r
            old_pts.append((phi1_star, float(phi2)))
            seeds.append((float(phi2), phi1_star, tw, kw, sep))
            print(f"  seed phi2={phi2:.3f}  phi1*={phi1_star:.4f}  sep={sep:.4f}")
    if not seeds:
        print("no seed found -- aborting")
        return
    best = max(seeds, key=lambda s: s[4])
    print(f"best seed: phi2={best[0]:.3f} phi1*={best[1]:.4f} sep={best[4]:.4f}")

    seed = (best[0], best[1], best[2], best[3])
    down, rd, td = continue_dir(chi, surf, seed, apex, -1)
    up, ru, tu = continue_dir(chi, surf, seed, apex, +1)
    line = sorted([(best[1], best[0])] + down + up, key=lambda p: p[1])
    print(f"\ncontinuation: {len(line)} pts, "
          f"phi2 in [{line[0][1]:.4f}, {line[-1][1]:.4f}]")
    print(f"  down stop={rd} at phi2={td:.4f}   up stop={ru} at phi2={tu:.4f}")
    print("  phi1_star  phi2")
    for p1, p2 in line:
        print(f"  {p1:.4f}    {p2:.4f}")

    outdir = os.path.join(ROOT, "out", "pw_continue_test")
    os.makedirs(outdir, exist_ok=True)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(bd.points[:, 0], bd.points[:, 1], s=6, c="#9db4c8",
               label="binodal", zorder=1)
    la = np.array(line)
    ax.plot(la[:, 0], la[:, 1], "-", color="#5a3d99", lw=2,
            label="PW continuation", zorder=3)
    if old_pts:
        oa = np.array(old_pts)
        ax.scatter(oa[:, 0], oa[:, 1], s=40, facecolors="none",
                   edgecolors="k", label="old multi-start", zorder=4)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.set_xlabel(r"$\phi_A$ ($\phi_1$)"); ax.set_ylabel(r"$\phi_R$ ($\phi_2$)")
    ax.set_title("PW continuation vs old method  (chi13=2.8, om=-0.30)")
    ax.legend(loc="upper right", fontsize=8)
    png = os.path.join(outdir, "overlay.png")
    fig.savefig(png, dpi=120, bbox_inches="tight")
    print(f"\nwrote {png}")


if __name__ == "__main__":
    main()
