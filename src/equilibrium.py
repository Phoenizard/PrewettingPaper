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

from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_bvp, simpson

import thermo as T


@dataclass
class Profile:
    z: np.ndarray            # (N,) grid
    phi: np.ndarray          # (2, N) profiles [phi1(z), phi2(z)]
    gamma: float             # excess surface free energy
    ok: bool                 # solver converged
    kind: str                # 'thin' or 'thick'


def _rhs(z, y, chi, res, kappa):
    dW1, dW2 = T.dW(y[0], y[1], chi, res)
    return np.vstack([y[2], y[3], dW1 / kappa.k1, dW2 / kappa.k2])


def _bc(ya, yb, chi, res, surf, kappa):
    g1, g2 = T.df_surf(ya[0], ya[1], surf)
    return np.array([
        kappa.k1 * ya[2] - g1,
        kappa.k2 * ya[3] - g2,
        yb[0] - res[0],
        yb[1] - res[1],
    ])


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


def _gamma(sol, chi, res, surf, kappa, L, n_quad=801):
    z = np.linspace(0.0, L, n_quad)
    y = sol.sol(z)
    Wv = T.W(y[0], y[1], chi, res)
    grad = 0.5 * kappa.k1 * y[2] ** 2 + 0.5 * kappa.k2 * y[3] ** 2
    bulk = simpson(Wv + grad, x=z)
    return float(bulk + T.f_surf(y[0][0], y[1][0], surf))


def solve_profile(chi, res, surf, kappa, seed, w=2.0, L=12.0, n=300):
    """Solve one profile whose wall value is seeded near `seed` (target phi0)."""
    x = np.linspace(0.0, L, n)
    y0 = _guess(res, seed, w, x)
    sol = solve_bvp(
        lambda z, y: _rhs(z, y, chi, res, kappa),
        lambda ya, yb: _bc(ya, yb, chi, res, surf, kappa),
        x, y0, tol=1e-8, max_nodes=30000,
    )
    if not sol.success:
        return None
    zf = np.linspace(0.0, L, 600)
    y = sol.sol(zf)
    if y[0].min() < -1e-3 or y[1].min() < -1e-3 or (y[0] + y[1]).max() > 1.0:
        return None
    gm = _gamma(sol, chi, res, surf, kappa, L)
    return Profile(z=zf, phi=y[:2], gamma=float(gm), ok=True, kind="?")


def find_states(chi, res, surf, kappa, dense_seeds=None, L=12.0, n=300):
    """All distinct surface states for reservoir res, via targeted multi-start.

    Seeds: flat (thin basin) + wall-plateau targets in the dense basin. Pass
    `dense_seeds` (list of target phi0 near the dense coexisting phase, e.g. from
    the binodal) for a general topology; default targets suit T-a (phi1-rich).
    Returns states sorted by phi1(0) ascending: first = thinnest, last = thickest.
    """
    res = np.asarray(res, dtype=float)
    if dense_seeds is None:
        dense_seeds = [(0.85, res[1]), (0.90, res[1]), (0.92, res[1] + 0.02)]
    seeds = [(res[0], res[1])] + list(dense_seeds)
    widths = (1.5, 2.5)
    found = []
    for s in seeds:
        for w in widths:
            p = solve_profile(chi, res, surf, kappa, s, w=w, L=L, n=n)
            if p is None:
                continue
            phi0 = (p.phi[0][0], p.phi[1][0])
            if all(np.hypot(phi0[0] - q.phi[0][0], phi0[1] - q.phi[1][0]) > 5e-3
                   for q in found):
                found.append(p)
    found.sort(key=lambda q: q.phi[0][0])
    return found
