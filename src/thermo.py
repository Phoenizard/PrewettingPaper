"""Bulk + surface thermodynamics for the regular ternary Flory-Huggins mixture.

All quantities use n1 = n2 = 1 and kT/nu = 1 (absorbed).  Every first
derivative is provided in closed form (derived in doc/note/ternary.md sec. 2.2)
so that downstream root-finding / quadrature is free of finite-difference noise.

Parameter bundles are plain namedtuples so functions stay vectorizable over
numpy arrays of phi.
"""
from __future__ import annotations

from collections import namedtuple

import numpy as np

# chi1s == chi13 (solute1-solvent), chi2s == chi23 (solute2-solvent)
Chi = namedtuple("Chi", ["chi12", "chi1s", "chi2s"])
Surf = namedtuple("Surf", ["w1", "w2", "cbb1", "cbb2", "cbb12"])
Kappa = namedtuple("Kappa", ["k1", "k2"])

# Below this volume fraction we treat a component as effectively absent; used
# only to keep the logs finite at the very edge of the simplex.
_EPS = 1e-300


def _phis(phi1, phi2):
    return 1.0 - phi1 - phi2


def _xlnx(x):
    """x*ln(x) with the x->0 limit (=0) handled; assumes x >= 0."""
    x = np.asarray(x, dtype=float)
    out = np.zeros_like(x)
    m = x > _EPS
    out[m] = x[m] * np.log(x[m])
    return out if out.ndim else float(out)


def f_b(phi1, phi2, chi: Chi):
    """Bulk free-energy density f_b(phi1, phi2).  Vectorized."""
    phi1 = np.asarray(phi1, dtype=float)
    phi2 = np.asarray(phi2, dtype=float)
    phis = _phis(phi1, phi2)
    entropy = _xlnx(phi1) + _xlnx(phi2) + _xlnx(phis)
    inter = chi.chi1s * phi1 * phis + chi.chi2s * phi2 * phis + chi.chi12 * phi1 * phi2
    return entropy + inter


def mu(phi1, phi2, chi: Chi):
    """Exchange chemical potentials (mu1, mu2) = (df_b/dphi1, df_b/dphi2).

    Closed form (n=1):
        mu1 = ln phi1 - ln phis + chi1s*(phis - phi1) + (chi12 - chi2s)*phi2
        mu2 = ln phi2 - ln phis + chi2s*(phis - phi2) + (chi12 - chi1s)*phi1
    """
    phi1 = np.asarray(phi1, dtype=float)
    phi2 = np.asarray(phi2, dtype=float)
    phis = _phis(phi1, phi2)
    lphi1 = np.log(np.clip(phi1, _EPS, None))
    lphi2 = np.log(np.clip(phi2, _EPS, None))
    lphis = np.log(np.clip(phis, _EPS, None))
    mu1 = lphi1 - lphis + chi.chi1s * (phis - phi1) + (chi.chi12 - chi.chi2s) * phi2
    mu2 = lphi2 - lphis + chi.chi2s * (phis - phi2) + (chi.chi12 - chi.chi1s) * phi1
    return mu1, mu2


def hessian_fb(phi1, phi2, chi: Chi):
    """Hessian [[f11, f12],[f12, f22]] of f_b, for spinodal (det=0) checks.

        f_ii = 1/phi_i + 1/phis - 2*chi_is
        f_12 = 1/phis + chi12 - chi1s - chi2s
    """
    phi1 = np.asarray(phi1, dtype=float)
    phi2 = np.asarray(phi2, dtype=float)
    phis = _phis(phi1, phi2)
    inv_s = 1.0 / phis
    f11 = 1.0 / phi1 + inv_s - 2.0 * chi.chi1s
    f22 = 1.0 / phi2 + inv_s - 2.0 * chi.chi2s
    f12 = inv_s + chi.chi12 - chi.chi1s - chi.chi2s
    return f11, f12, f22


def spinodal_det(phi1, phi2, chi: Chi):
    """det(Hessian) — zero locus is the spinodal."""
    f11, f12, f22 = hessian_fb(phi1, phi2, chi)
    return f11 * f22 - f12 * f12


def W(phi1, phi2, chi: Chi, res):
    """Grand-potential density relative to reservoir res=(phi1_inf, phi2_inf).

    W = f_b(phi) - f_b(res) - mu1_inf*(phi1-phi1_inf) - mu2_inf*(phi2-phi2_inf)
    i.e. the vertical distance of f_b above its tangent plane at the reservoir.
    W(res) = 0 and W is a local minimum there.
    """
    r1, r2 = res
    m1i, m2i = mu(r1, r2, chi)
    fbi = f_b(r1, r2, chi)
    return f_b(phi1, phi2, chi) - fbi - m1i * (phi1 - r1) - m2i * (phi2 - r2)


def dW(phi1, phi2, chi: Chi, res):
    """(dW/dphi1, dW/dphi2) = (mu1(phi) - mu1_inf, mu2(phi) - mu2_inf)."""
    r1, r2 = res
    m1i, m2i = mu(r1, r2, chi)
    m1, m2 = mu(phi1, phi2, chi)
    return m1 - m1i, m2 - m2i


def f_surf(p10, p20, surf: Surf):
    """Surface free energy evaluated at surface compositions (p10, p20)."""
    return (
        surf.w1 * p10
        + surf.w2 * p20
        + surf.cbb1 * p10 * p10
        + surf.cbb2 * p20 * p20
        + surf.cbb12 * p10 * p20
    )


def df_surf(p10, p20, surf: Surf):
    """(df_surf/dp10, df_surf/dp20)."""
    g1 = surf.w1 + 2.0 * surf.cbb1 * p10 + surf.cbb12 * p20
    g2 = surf.w2 + 2.0 * surf.cbb2 * p20 + surf.cbb12 * p10
    return g1, g2
