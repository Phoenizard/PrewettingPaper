"""Equilibrium-DE pre-wetting solver (project_plan.md 3.3).

Solve the two Euler-Lagrange ODEs for the concentration profiles phi_1(z),
phi_2(z) against a wall at z=0, with a far-field reservoir at z=L:

    kappa_i phi_i''(z) = dW/dphi_i                    (z>0)
    kappa_i phi_i'(0)  = df_surf/dphi_i(0)            (wall BC)
    phi_i(L)           = phi_{i,inf}                  (far field)

State y = [phi1, phi2, phi1', phi2'] solved with scipy.integrate.solve_bvp.
Two initial guesses (flat = thin, wall-plateau = thick) give the thin/thick
profiles; the stable one has the smaller gamma. All thermodynamics reused from
thermo (dW, df_surf, W, f_surf).
"""
from __future__ import annotations

import json
import os
import sys
from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_bvp, simpson

import thermo as T

# Optimization toggles (see doc/plan). Both default on; flip off to bisect a
# regression. USE_ANALYTIC_JAC feeds solve_bvp the exact Jacobian (same equations,
# faster/steadier Newton). USE_WARM_START reuses a neighbour's converged profile
# as the initial guess, with a guard + cold-multistart fallback.
USE_ANALYTIC_JAC = True
USE_WARM_START = True


def _diag(**kw):
    """Emit one [PWDIAG] JSONL line to stderr when PW_DIAG is set (diagnosis only,
    no behaviour change, parallel-safe: stderr, no files)."""
    if os.environ.get("PW_DIAG"):
        print("[PWDIAG] " + json.dumps(kw), file=sys.stderr, flush=True)


@dataclass
class Profile:
    z: np.ndarray            # (N,) grid
    phi: np.ndarray          # (2, N) profiles [phi1(z), phi2(z)]
    gamma: float             # excess surface free energy
    ok: bool                 # solver converged
    kind: str                # 'thin' or 'thick'
    sol_x: np.ndarray = None  # raw solver mesh (for warm-starting the next point)
    sol_y: np.ndarray = None  # raw solver state (4, len(sol_x))


def _rhs(z, y, chi, res, kappa, res_mu=None):
    dW1, dW2 = T.dW(y[0], y[1], chi, res, res_mu=res_mu)
    return np.vstack([y[2], y[3], dW1 / kappa.k1, dW2 / kappa.k2])


def _fun_jac(z, y, chi, kappa):
    """d(rhs)/dy, shape (4, 4, m). d(dW_i)/dphi_j = hessian_fb (res-independent)."""
    f11, f12, f22 = T.hessian_fb(y[0], y[1], chi)
    m = y.shape[1]
    J = np.zeros((4, 4, m))
    J[0, 2] = 1.0
    J[1, 3] = 1.0
    J[2, 0] = f11 / kappa.k1
    J[2, 1] = f12 / kappa.k1
    J[3, 0] = f12 / kappa.k2
    J[3, 1] = f22 / kappa.k2
    return J


def _bc(ya, yb, chi, res, surf, kappa):
    g1, g2 = T.df_surf(ya[0], ya[1], surf)
    return np.array([
        kappa.k1 * ya[2] - g1,
        kappa.k2 * ya[3] - g2,
        yb[0] - res[0],
        yb[1] - res[1],
    ])


def _bc_jac(ya, yb, surf, kappa):
    """(d bc/d ya, d bc/d yb), each (4, 4). All entries constant (df_surf linear)."""
    Ja = np.zeros((4, 4))
    Jb = np.zeros((4, 4))
    Ja[0, 0] = -2.0 * surf.cbb1
    Ja[0, 1] = -surf.cbb12
    Ja[0, 2] = kappa.k1
    Ja[1, 0] = -surf.cbb12
    Ja[1, 1] = -2.0 * surf.cbb2
    Ja[1, 3] = kappa.k2
    Jb[2, 0] = 1.0
    Jb[3, 1] = 1.0
    return Ja, Jb


def _guess(res, seed, w, x):
    """Initial (4, m) guess: profile decays from wall value `seed` to `res`."""
    res = np.asarray(res, dtype=float)
    seed = np.asarray(seed, dtype=float)
    e = np.exp(-x / w)
    y = np.zeros((4, x.size))
    y[0] = res[0] + (seed[0] - res[0]) * e
    y[1] = res[1] + (seed[1] - res[1]) * e
    y[2] = -(seed[0] - res[0]) / w * e
    y[3] = -(seed[1] - res[1]) / w * e
    return y


def _gamma(sol, chi, res, surf, kappa, L, n_quad=801, res_mu=None, res_fb=None):
    z = np.linspace(0.0, L, n_quad)
    y = sol.sol(z)
    Wv = T.W(y[0], y[1], chi, res, res_mu=res_mu, res_fb=res_fb)
    grad = 0.5 * kappa.k1 * y[2] ** 2 + 0.5 * kappa.k2 * y[3] ** 2
    bulk = simpson(Wv + grad, x=z)
    return float(bulk + T.f_surf(y[0][0], y[1][0], surf))


def solve_profile(chi, res, surf, kappa, seed=None, w=2.0, L=12.0, n=300,
                  res_mu=None, res_fb=None, warm=None):
    """Solve one profile. Initial guess is either a wall->reservoir decay seeded
    near `seed` (target phi0), or, when `warm=(x, y)` is given, a neighbouring
    converged solution (resampled onto n nodes if its mesh grew large)."""
    if warm is not None:
        x, y0 = warm
        if x.size > 2000:
            xs = np.linspace(0.0, L, n)
            y0 = np.vstack([np.interp(xs, x, y0[i]) for i in range(4)])
            x = xs
    else:
        x = np.linspace(0.0, L, n)
        y0 = _guess(res, seed, w, x)
    fun_jac = (lambda z, y: _fun_jac(z, y, chi, kappa)) if USE_ANALYTIC_JAC else None
    bc_jac = (lambda ya, yb: _bc_jac(ya, yb, surf, kappa)) if USE_ANALYTIC_JAC else None
    sol = solve_bvp(
        lambda z, y: _rhs(z, y, chi, res, kappa, res_mu),
        lambda ya, yb: _bc(ya, yb, chi, res, surf, kappa),
        x, y0, fun_jac=fun_jac, bc_jac=bc_jac, tol=1e-8, max_nodes=30000,
    )
    if not sol.success:
        return None
    zf = np.linspace(0.0, L, 600)
    y = sol.sol(zf)
    if y[0].min() < -1e-3 or y[1].min() < -1e-3 or (y[0] + y[1]).max() > 1.0:
        return None
    gm = _gamma(sol, chi, res, surf, kappa, L, res_mu=res_mu, res_fb=res_fb)
    return Profile(z=zf, phi=y[:2], gamma=float(gm), ok=True, kind="?",
                   sol_x=sol.x, sol_y=sol.y)


def _distinct(profiles, tol=5e-3):
    """Keep profiles with distinct wall composition (phi1(0), phi2(0))."""
    out = []
    for p in profiles:
        phi0 = (p.phi[0][0], p.phi[1][0])
        if all(np.hypot(phi0[0] - q.phi[0][0], phi0[1] - q.phi[1][0]) > tol
               for q in out):
            out.append(p)
    return out


def _warm_ok(states, warm, band=0.05, sep=5e-3):
    """Accept a warm-started result only if it reproduces two distinct, correctly
    ordered branches that stayed close (wall phi1 within `band`) to last point's."""
    if len(states) != 2:
        return False
    states = sorted(states, key=lambda q: q.phi[0][0])
    if states[1].phi[0][0] - states[0].phi[0][0] <= sep:
        return False
    prev_thin_phi1 = warm[0][1][0][0]   # warm[0]=(x,y); y[0]=phi1 row; [0]=wall node
    prev_thick_phi1 = warm[1][1][0][0]
    return (abs(states[0].phi[0][0] - prev_thin_phi1) <= band
            and abs(states[1].phi[0][0] - prev_thick_phi1) <= band)


def find_states(chi, res, surf, kappa, dense_seeds=None, L=12.0, n=300, warm=None,
                distinct_tol=5e-3):
    """All distinct surface states for reservoir res, via targeted multi-start.

    Seeds: flat (thin basin) + wall-plateau targets in the dense basin. Pass
    `dense_seeds` (list of target phi0 near the dense coexisting phase, e.g. from
    the binodal) for a general topology; default targets suit T-a (phi1-rich).
    Returns states sorted by phi1(0) ascending: first = thinnest, last = thickest.

    `warm=(thin_state, thick_state)` (each an (x, y) from the previous point) tries
    two targeted warm-started solves first; the result is used only if it passes
    `_warm_ok`, otherwise we fall back to the full cold multi-start (== old path).

    `distinct_tol` is the wall-composition separation below which two profiles count
    as the same branch. The default 5e-3 is the production value; the pre-wetting
    endpoint fallback (verify.py) passes a smaller value so the near-merged thin/thick
    branches at the surface-critical end survive as two states.
    """
    res = np.asarray(res, dtype=float)
    m1i, m2i, fbi = T.reservoir_potentials(res, chi)
    res_mu = (m1i, m2i)

    if USE_WARM_START and warm is not None:
        cand = [solve_profile(chi, res, surf, kappa, L=L, n=n,
                              res_mu=res_mu, res_fb=fbi, warm=w) for w in warm]
        acc = _distinct([p for p in cand if p is not None], tol=distinct_tol)
        if _warm_ok(acc, warm, sep=distinct_tol):
            acc.sort(key=lambda q: q.phi[0][0])
            return acc
        # guard failed -> cold multi-start below

    if dense_seeds is None:
        dense_seeds = [(0.85, res[1]), (0.90, res[1]), (0.92, res[1] + 0.02)]
    seeds = [(res[0], res[1])] + list(dense_seeds)
    widths = (1.5, 2.5)
    found = []
    raw = []  # every converged candidate, before the distinctness filter (diag)
    for s in seeds:
        for w in widths:
            p = solve_profile(chi, res, surf, kappa, s, w=w, L=L, n=n,
                              res_mu=res_mu, res_fb=fbi)
            if p is None:
                continue
            phi0 = (p.phi[0][0], p.phi[1][0])
            raw.append((phi0, s, w))
            if all(np.hypot(phi0[0] - q.phi[0][0], phi0[1] - q.phi[1][0]) > distinct_tol
                   for q in found):
                found.append(p)
    found.sort(key=lambda q: q.phi[0][0])
    if os.environ.get("PW_DIAG"):
        # Split "thick never converged" (few raw candidates) from "branches merged
        # under the threshold" (>=2 raw but <2 distinct) at the empty phi2 lines.
        _diag(kind="find_states", res=[float(res[0]), float(res[1])],
              distinct_tol=distinct_tol, L=L, n_raw=len(raw), n_distinct=len(found),
              raw=[[float(p0[0]), float(p0[1]), [float(s[0]), float(s[1])], float(w)]
                   for (p0, s, w) in raw])
    return found


# ---------------------------------------------------------------------------
# Pre-wetting line by branch continuation (project note 3.3).
#
# Multi-start find_states rediscovers the thin/thick branches independently at each
# reservoir; near the ENDS of the pre-wetting line the two branches nearly merge and
# multi-start can no longer separate them, so the line drops points. Continuation
# instead TRACKS the two branches: from a converged (thin, thick) pair it warm-starts
# both at a neighbouring reservoir, so each branch is followed rather than re-found.
# This holds the pair apart down to the surface-critical terminus where they truly
# merge. Reuses solve_profile(warm=...) and Profile.sol_x/sol_y.
# ---------------------------------------------------------------------------

def _solve_branch_adaptive(chi, res, surf, kappa, warm, L, n, L_grow=(2.0, 3.5)):
    """Warm-solve one branch; if it fails, retry WITHOUT the warm start on a larger domain.

    The thick branch's interface widens as phi2 -> 0 and gets truncated by the fixed
    L=12, so solve_bvp fails. Simply resampling the warm profile onto the larger L does
    not help: the warm mesh only covers [0, 12], so beyond 12 np.interp pads a flat
    (wrong) tail — a bad guess. Instead we drop the warm start and rebuild a clean
    wall->reservoir decaying guess (_guess) on the larger domain, seeded from the warm
    profile's WALL composition (the one piece of the warm state that is still valid).
    n is scaled with L to keep the near-wall mesh density. Returns the first converged
    Profile, or None if every attempt fails."""
    p = solve_profile(chi, res, surf, kappa, L=L, n=n, warm=warm)
    if p is not None:
        return p
    wall_seed = None
    if warm is not None:
        _, wy = warm
        wall_seed = (float(wy[0][0]), float(wy[1][0]))  # wall (phi1, phi2) of the branch
    for k in L_grow:
        Lk, nk = L * k, int(round(n * k))
        p = solve_profile(chi, res, surf, kappa, seed=wall_seed, L=Lk, n=nk)
        if p is not None:
            return p
    return None


def track_branches(chi, res, surf, kappa, warm_thin, warm_thick, L=12.0, n=300):
    """Warm-solve the thin and thick branch at reservoir res from converged neighbours
    `warm_thin`/`warm_thick` (each an (x, y) mesh+state). Either return may be None if
    that branch failed to converge even after growing the domain."""
    thin = _solve_branch_adaptive(chi, res, surf, kappa, warm_thin, L, n)
    thick = _solve_branch_adaptive(chi, res, surf, kappa, warm_thick, L, n)
    return thin, thick


def pw_point(chi, phi2, surf, kappa, warm_thin, warm_thick, phi1_guess,
             L=12.0, n=300, merge_tol=3e-3, xtol=1e-6, max_it=40, h=3e-4):
    """Corrector at fixed phi2: find phi1_inf where gamma_thin == gamma_thick, tracking
    the two branches from warm starts via a secant root-find in phi1.

    Returns (phi1_star, thin_profile, thick_profile), or None when the branches merge or
    a branch collapses/fails to converge -- i.e. the pre-wetting line has terminated.
    `merge_tol` is the wall-composition gap below which the pair counts as merged.
    """
    def g(phi1):
        thin, thick = track_branches(chi, (phi1, phi2), surf, kappa,
                                     warm_thin, warm_thick, L=L, n=n)
        if thin is None or thick is None:
            return None, thin, thick
        # secant tracks each branch; if warm-thick collapsed onto thin (or they truly
        # merged) the wall compositions coincide -> treat as terminus.
        if abs(thin.phi[0][0] - thick.phi[0][0]) < merge_tol:
            return None, thin, thick
        return thin.gamma - thick.gamma, thin, thick

    x0 = phi1_guess
    g0, _, _ = g(x0)
    if g0 is None:
        return None
    x1 = phi1_guess + h
    g1, thin1, thick1 = g(x1)
    if g1 is None:
        return None
    for _ in range(max_it):
        if g1 == g0:
            break
        x2 = x1 - g1 * (x1 - x0) / (g1 - g0)
        step = x2 - x1
        if abs(step) > 0.02:           # damp runaway secant jumps
            x2 = x1 + (0.02 if step > 0 else -0.02)
        g2, thin2, thick2 = g(x2)
        if g2 is None:
            return None
        x0, g0 = x1, g1
        x1, g1, thin1, thick1 = x2, g2, thin2, thick2
        if abs(x1 - x0) < xtol:
            return x1, thin1, thick1
    return None
