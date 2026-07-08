"""Bulk thermodynamics: free energy density, exchange chemical potentials,
Hessian. Pure numpy functions of (phi1, phi2) plus Physical params.

Formulas mirror both reference codes (identical there); see
doc/note/code_cross_comparison.md section 1. Dimensionless: kBT/nu = 1.
phi3 = 1 - phi1 - phi2 is the solvent.
"""

import numpy as np

EPS = 1e-12


def clip_composition(phi1, phi2, eps=EPS):
    p1 = np.clip(phi1, eps, 1.0 - eps)
    p2 = np.clip(phi2, eps, 1.0 - p1 - eps)
    p3 = np.clip(1.0 - p1 - p2, eps, 1.0)
    return p1, p2, p3


def free_energy(phi1, phi2, p):
    p1, p2, p3 = clip_composition(phi1, phi2)
    f = p1 * np.log(p1) / p.n1 + p2 * np.log(p2) / p.n2 + p3 * np.log(p3) / p.n3
    return f + p.chi_12 * p1 * p2 + p.chi_13 * p1 * p3 + p.chi_23 * p2 * p3


def chemical_potential(phi1, phi2, p):
    p1, p2, p3 = clip_composition(phi1, phi2)
    mu1 = (
        (np.log(p1) + 1.0) / p.n1
        - (np.log(p3) + 1.0) / p.n3
        + p.chi_12 * p2
        + p.chi_13 * (p3 - p1)
        - p.chi_23 * p2
    )
    mu2 = (
        (np.log(p2) + 1.0) / p.n2
        - (np.log(p3) + 1.0) / p.n3
        + p.chi_12 * p1
        - p.chi_13 * p1
        + p.chi_23 * (p3 - p2)
    )
    return mu1, mu2


def hessian(phi1, phi2, p):
    p1, p2, p3 = clip_composition(phi1, phi2)
    f11 = 1.0 / (p.n1 * p1) + 1.0 / (p.n3 * p3) - 2.0 * p.chi_13
    f12 = 1.0 / (p.n3 * p3) + p.chi_12 - p.chi_13 - p.chi_23
    f22 = 1.0 / (p.n2 * p2) + 1.0 / (p.n3 * p3) - 2.0 * p.chi_23
    return f11, f12, f22
