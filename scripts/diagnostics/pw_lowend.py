#!/usr/bin/env python
"""Diagnose WHY the pre-wetting continuation stops at the low-phi2 end.

pw_point() returns a bare None on four very different failures; _continue() then just
stops, so we cannot tell a real physical terminus (the two branches genuinely merge,
consistent with continuous binary-limit wetting) from a numerical give-up. This probe
re-implements pw_point's inner logic with instrumentation: it marches phi2 DOWN from the
seed's low end toward the phi2=0 axis and, at each step, prints which exit is hit and the
quantities that decide it. It reads only; it changes no production code.

Exits (mirror src/equilibrium.pw_point):
  A  a branch failed to converge (solve_profile -> None) at the initial guess
  B  merge: |phi_thin(0) - phi_thick(0)| < merge_tol  (branches coincide -> terminus)
  C  a branch failed mid-secant
  D  secant used max_it without reaching |dphi1| < xtol
  OK a root was found (phi1*, gamma_thin==gamma_thick)

Run on the server, in a screen:
  python scripts/diagnostics/pw_lowend.py [<chi_dir> <om_dir> <chibb_dir>]  # default om2=-0.40
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.dirname(HERE)
ROOT = os.path.dirname(SCRIPTS)
sys.path.insert(0, SCRIPTS)                       # import verify
sys.path.insert(0, os.path.join(ROOT, "src"))

import numpy as np  # noqa: E402
import cases  # noqa: E402
import binodal as B  # noqa: E402
import equilibrium as E  # noqa: E402
import verify as V  # noqa: E402

KAPPA = E.T.Kappa(1.0, 1.0)


def _pw_point_diag(chi, phi2, surf, warm_thin, warm_thick, phi1_guess,
                   merge_tol=3e-3, xtol=1e-6, max_it=40, h=3e-4):
    """Instrumented copy of equilibrium.pw_point: returns (label, info, thin, thick).
    label in {A,B,C,D,OK}. info carries the deciding numbers for that step."""
    def g(phi1):
        thin, thick = E.track_branches(chi, (phi1, phi2), surf, KAPPA,
                                       warm_thin, warm_thick)
        if thin is None or thick is None:
            return None, thin, thick, ("branch", thin is None, thick is None)
        gap = abs(thin.phi[0][0] - thick.phi[0][0])
        if gap < merge_tol:
            return None, thin, thick, ("merge", gap)
        return thin.gamma - thick.gamma, thin, thick, ("ok", gap)

    x0 = phi1_guess
    g0, thin0, thick0, why0 = g(x0)
    if g0 is None:
        lab = "B" if why0[0] == "merge" else "A"
        return lab, why0, thin0, thick0
    x1 = phi1_guess + h
    g1, thin1, thick1, why1 = g(x1)
    if g1 is None:
        lab = "B" if why1[0] == "merge" else "A"
        return lab, why1, thin1, thick1
    for _ in range(max_it):
        if g1 == g0:
            break
        x2 = x1 - g1 * (x1 - x0) / (g1 - g0)
        step = x2 - x1
        if abs(step) > 0.02:
            x2 = x1 + (0.02 if step > 0 else -0.02)
        g2, thin2, thick2, why2 = g(x2)
        if g2 is None:
            lab = "B" if why2[0] == "merge" else "C"
            return lab, (why2, f"x2={x2:.4f}"), thin2, thick2
        x0, g0 = x1, g1
        x1, g1, thin1, thick1, why1 = x2, g2, thin2, thick2, why2
        if abs(x1 - x0) < xtol:
            return "OK", ("root", f"phi1*={x1:.5f}", why1), thin1, thick1
    return "D", ("no_converge", why1), thin1, thick1


def main(argv):
    if len(argv) >= 3:
        rel = tuple(argv[:3])
    else:
        rel = ("chi12_0p0__chi13_2p8__chi23_0p0", "om1_m0p30__om2_m0p40",
               "chibb11_0p00__chibb22_0p00__chibb12_0p00")
    chi, surf = cases.parse_case(*rel)
    bd = B.binodal_from_hull(chi)
    f_left, f_right, apex = V._branches(bd)

    # Seed exactly as production does: widest-separation pair on the middle band.
    seeds = []
    for phi2 in np.arange(0.03, min(0.09, 0.85 * apex), 0.01):
        r = V._seed_scan(chi, surf, float(phi2), f_left, f_right)
        if r:
            seeds.append((float(phi2), r[0], r[1], r[2], r[3]))
    if not seeds:
        print(f"{'/'.join(rel)}: no seed", flush=True)
        return
    best = max(seeds, key=lambda s: s[4])
    phi2_s, phi1_s, tw, kw, sep = best
    print(f"{'/'.join(rel[1:2])}  apex_phi2={apex:.4f}  seed phi2={phi2_s:.3f} "
          f"phi1*={phi1_s:.4f} sep={sep:.3f}", flush=True)

    # First walk DOWN with the production step to reach the low end (0.026-ish), warm
    # chaining as we go, then keep stepping past it toward the axis to see the failure.
    print("  phi2     exit  detail", flush=True)
    phi1_prev, phi1_prev2 = phi1_s, None
    phi2 = phi2_s
    dphi2 = 0.001            # fine, so we watch every step near the low end
    warm_t, warm_k = tw, kw
    while phi2 - dphi2 > 0.0:
        phi2n = round(phi2 - dphi2, 10)
        guess = phi1_prev if phi1_prev2 is None else (2.0 * phi1_prev - phi1_prev2)
        lab, info, thin, thick = _pw_point_diag(chi, phi2n, surf, warm_t, warm_k, guess)
        gapstr = ""
        if thin is not None and thick is not None:
            gapstr = f"  gap={abs(thin.phi[0][0]-thick.phi[0][0]):.4f}"
        print(f"  {phi2n:.4f}   {lab}    {info}{gapstr}", flush=True)
        if lab != "OK":
            print(f"  -> stops at phi2={phi2n:.4f} via exit {lab}", flush=True)
            break
        # accept: chain warm starts and predictor
        phi1_star = float(info[1].split("=")[1])
        warm_t = (thin.sol_x, thin.sol_y)
        warm_k = (thick.sol_x, thick.sol_y)
        phi1_prev2, phi1_prev = phi1_prev, phi1_star
        phi2 = phi2n
    print("DONE", flush=True)


if __name__ == "__main__":
    main(sys.argv[1:])
