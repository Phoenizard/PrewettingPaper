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


def _solve_profile_why(chi, res, surf, seed=None, w=2.0, L=12.0, n=300, warm=None,
                       tol=1e-8, max_nodes=30000):
    """Copy of equilibrium.solve_profile that, instead of a bare None, RETURNS WHY it
    failed so the autopsy can tell apart the failure modes:
      ('nonconv', sol)         solve_bvp did not converge (sol.success False)
      ('reject', reasons, y)   converged but the post-solve sanity check killed it
                               reasons subset of {phi1<0, phi2<0, sum>1}
      ('ok', profile, width)   converged & accepted; width = z where |phi-phi_res|
                               first falls below 1% of the wall excess (film width est.)
    Mirrors src line-for-line (equilibrium.solve_profile, ~line 119) so its verdict is
    the same one production would reach; reads only, no src change."""
    res_arr = np.asarray(res, dtype=float)
    m1i, m2i, fbi = E.T.reservoir_potentials(res_arr, chi)
    res_mu = (m1i, m2i)
    if warm is not None:
        x, y0 = warm
        if x.size > 2000:
            xs = np.linspace(0.0, L, n)
            y0 = np.vstack([np.interp(xs, x, y0[i]) for i in range(4)])
            x = xs
    else:
        x = np.linspace(0.0, L, n)
        y0 = E._guess(res_arr, seed, w, x)
    fun_jac = (lambda z, y: E._fun_jac(z, y, chi, KAPPA)) if E.USE_ANALYTIC_JAC else None
    bc_jac = (lambda ya, yb: E._bc_jac(ya, yb, surf, KAPPA)) if E.USE_ANALYTIC_JAC else None
    from scipy.integrate import solve_bvp
    sol = solve_bvp(
        lambda z, y: E._rhs(z, y, chi, res_arr, KAPPA, res_mu),
        lambda ya, yb: E._bc(ya, yb, chi, res_arr, surf, KAPPA),
        x, y0, fun_jac=fun_jac, bc_jac=bc_jac, tol=tol, max_nodes=max_nodes,
    )
    if not sol.success:
        return ("nonconv", f"status={sol.status} nodes={sol.x.size} msg={sol.message}")
    zf = np.linspace(0.0, L, 600)
    y = sol.sol(zf)
    reasons = []
    if y[0].min() < -1e-3:
        reasons.append(f"phi1min={y[0].min():.4f}")
    if y[1].min() < -1e-3:
        reasons.append(f"phi2min={y[1].min():.4f}")
    if (y[0] + y[1]).max() > 1.0:
        reasons.append(f"summax={(y[0]+y[1]).max():.4f}")
    if reasons:
        return ("reject", reasons)
    # film width: distance where phi1 has relaxed to within 1% of its wall excess
    wall_exc = abs(y[0][0] - res_arr[0])
    width = L
    if wall_exc > 1e-6:
        below = np.where(np.abs(y[0] - res_arr[0]) < 0.01 * wall_exc)[0]
        if below.size:
            width = zf[below[0]]
    return ("ok", (y[0][0], y[1][0]), width)


def _thick_autopsy(chi, phi2, surf, warm_thick, phi1, f_right_g):
    """At a failed (phi1, phi2), replay the thick branch's three attempts that
    _solve_branch_adaptive makes (warm@L=12, _guess@L=24, _guess@L=42) and print WHY each
    one is None (nonconv vs reject) or, if it converges, how wide the film is. This is the
    core question: is the thick branch un-converging, being rejected, or too wide for L."""
    res = (phi1, phi2)
    wall_seed = None
    if warm_thick is not None:
        _, wy = warm_thick
        wall_seed = (float(wy[0][0]), float(wy[1][0]))
    print(f"    [autopsy] thick branch at phi1={phi1:.4f} phi2={phi2:.4f} "
          f"wall_seed={wall_seed}", flush=True)
    # attempt 1: warm @ L=12 (what production tries first)
    r = _solve_profile_why(chi, res, surf, warm=warm_thick, L=12.0, n=300)
    print(f"    [autopsy] warm@L=12  -> {r[0]}: {r[1:]}", flush=True)
    # attempts 2,3: clean _guess seeded from wall, larger L (the fallback)
    for k in (2.0, 3.5):
        Lk, nk = 12.0 * k, int(round(300 * k))
        r = _solve_profile_why(chi, res, surf, seed=wall_seed, L=Lk, n=nk)
        print(f"    [autopsy] guess@L={Lk:.0f} (w=2.0) -> {r[0]}: {r[1:]}", flush=True)
    # extra probes: does a LONGER guess length-scale w help at L=42?
    for w in (6.0, 12.0):
        r = _solve_profile_why(chi, res, surf, seed=wall_seed, L=42.0, n=1050, w=w)
        print(f"    [autopsy] guess@L=42 w={w:.0f} -> {r[0]}: {r[1:]}", flush=True)
    # DECISIVE: warm@L=12 succeeded at the predictor guess but production judged thick=None,
    # so the secant must have probed a phi1 where the thick branch (via warm@L=12) fails.
    # Scan phi1 around the guess at THIS phi2, calling the production track_branches path,
    # and report thin/thick None-status + wall phi + why (nonconv/reject) for the thick.
    print(f"    [autopsy] phi1 scan at phi2={phi2:.4f} (thick via warm@L=12):", flush=True)
    fail_phi1 = None
    for phi1s in np.arange(phi1 - 0.004, phi1 + 0.0041, 0.0005):
        why = _solve_profile_why(chi, (float(phi1s), phi2), surf, warm=warm_thick,
                                 L=12.0, n=300)
        print(f"      phi1={phi1s:.4f}  thick warm@L=12 -> {why[0]}: {why[1:]}", flush=True)
        if why[0] == "nonconv" and phi1s > phi1 and fail_phi1 is None:
            fail_phi1 = float(phi1s)  # first failure to the RIGHT of the guess
    # CANDIDATE FIX test: at a phi1 that failed with tol=1e-8/max_nodes=30000, does a
    # looser tol or a bigger node budget converge? (failure was max_nodes exceeded.)
    if fail_phi1 is not None:
        print(f"    [autopsy] fix test at phi1={fail_phi1:.4f} phi2={phi2:.4f} "
              f"(thick warm@L=12):", flush=True)
        for tol, mn in [(1e-8, 100000), (1e-7, 30000), (1e-7, 100000), (1e-6, 100000)]:
            why = _solve_profile_why(chi, (fail_phi1, phi2), surf, warm=warm_thick,
                                     L=12.0, n=300, tol=tol, max_nodes=mn)
            print(f"      tol={tol:.0e} max_nodes={mn} -> {why[0]}: {why[1:]}", flush=True)
        # EXISTENCE test: warm keeps failing past phi1=0.0984 with runaway nodes, which
        # looks like the thick branch VANISHING (surface spinodal) rather than a solver
        # budget issue. Confirm from a COLD multi-start (find_states, no warm): does a
        # distinct thick state exist at all here? If find_states returns only ONE state,
        # the thick branch is physically gone -> the line truly terminates near phi2=0.026.
        dense = [(f_right_g(phi2), phi2), (0.97 * f_right_g(phi2), phi2)]
        for p1 in (float(phi1), fail_phi1, fail_phi1 + 0.002):
            st = E.find_states(chi, (p1, phi2), surf, KAPPA, dense_seeds=dense)
            walls = [(round(float(s.phi[0][0]), 4), round(float(s.phi[1][0]), 4))
                     for s in st]
            print(f"      cold find_states phi1={p1:.4f} -> n={len(st)} walls={walls}",
                  flush=True)


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
    # Stream every phi2 and every phi1 start point (no silent loop).
    band = np.arange(0.03, min(0.09, 0.85 * apex), 0.01)
    seeds = []
    for i, phi2 in enumerate(band):
        print(f"[lowend] seed {i+1}/{len(band)} 计算 phi2={phi2:.3f} ...", flush=True)
        r = V._seed_scan(chi, surf, float(phi2), f_left, f_right, progress="lowend")
        if r:
            seeds.append((float(phi2), r[0], r[1], r[2], r[3]))
        print(f"[lowend] seed {i+1}/{len(band)} phi2={phi2:.3f} DONE "
              f"{'phi1*=%.4f sep=%.3f' % (r[0], r[3]) if r else 'no pair'}", flush=True)
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
            # Autopsy the THICK branch at the failing point: is it non-converging, being
            # rejected, or just too wide for L=12 — and does a longer guess length-scale
            # rescue it? (warm_k is the last accepted thick profile = the warm start used.)
            _thick_autopsy(chi, phi2n, surf, warm_k, guess, f_right)
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
