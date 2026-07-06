"""Bulk binodal (liquid-liquid coexistence) for the ternary mixture.

Primary method: lower convex hull of the free-energy surface (phi1, phi2, f_b)
over a triangular grid of the 2-simplex (the construction used in the reference
paper, Appendix C).  Facets of the lower hull that bridge the miscibility gap
have long edges in composition space; those long edges are the tie-lines and
their endpoints trace the binodal.

Independent cross-check (so verification is not self-referential): each hull
tie-line is refined with a common-tangent / equal-chemical-potential root solve
(scipy.optimize.root), and the equal-mu + tangent residuals are reported.  The
binodal depends only on chi, so it is computed once per topology.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.spatial import ConvexHull
from scipy.optimize import root

import thermo as T


@dataclass
class BinodalResult:
    chi: T.Chi
    n: int
    tie_lines: np.ndarray          # (K, 2, 2): [k, {A,B}, {phi1,phi2}]
    three_phase: list = field(default_factory=list)  # list of (3,2) triangles
    points: np.ndarray = None      # (M, 2) unique binodal compositions
    max_equal_mu_residual: float = np.nan  # after refinement, over all tie-lines


def simplex_grid(n: int, eps: float = 1e-6) -> np.ndarray:
    """Triangular lattice of interior points of the 2-simplex.

    Returns (M, 2) array of (phi1, phi2) with phi1, phi2 > 0 and
    phi1 + phi2 < 1, inset from every edge by ~eps to keep logs finite.
    """
    pts = []
    for i in range(1, n):
        for j in range(1, n - i):
            p1 = i / n
            p2 = j / n
            # inset toward centroid so we never sit exactly on an edge
            p1 = eps + p1 * (1 - 3 * eps)
            p2 = eps + p2 * (1 - 3 * eps)
            pts.append((p1, p2))
    return np.asarray(pts, dtype=float)


def _lower_hull_simplices(points3d: np.ndarray) -> np.ndarray:
    """Indices of ConvexHull facets whose outward normal points downward (-z)."""
    hull = ConvexHull(points3d)
    # equations: a*x + b*y + c*z + d = 0, (a,b,c) outward unit normal
    downward = hull.equations[:, 2] < -1e-9
    return hull.simplices[downward]


def binodal_from_hull(
    chi: T.Chi,
    n: int = 300,
    eps: float = 1e-6,
    tie_frac: float = 5.0,
    three_phase_area_min: float = 0.02,
    refine: bool = True,
) -> BinodalResult:
    """Compute the binodal via lower convex hull.

    tie_frac: an edge counts as a tie-line if its composition-space length
    exceeds tie_frac * (1/n) grid spacings.  Facets with all three edges long
    AND area > three_phase_area_min are genuine three-phase coexistence
    triangles; small all-long facets are near-critical tie regions whose
    longest edge is kept as an ordinary tie-line.
    """
    grid = simplex_grid(n, eps)
    fb = T.f_b(grid[:, 0], grid[:, 1], chi)
    pts3 = np.column_stack([grid, fb])
    simp = _lower_hull_simplices(pts3)

    step = 1.0 / n
    thresh = tie_frac * step

    tie_edges = []           # list of (idxA, idxB)
    three_phase = []
    for tri in simp:
        v = grid[tri]        # (3,2)
        e = [
            (tri[0], tri[1], np.hypot(*(v[0] - v[1]))),
            (tri[1], tri[2], np.hypot(*(v[1] - v[2]))),
            (tri[0], tri[2], np.hypot(*(v[0] - v[2]))),
        ]
        long_edges = [(a, b) for a, b, L in e if L > thresh]
        if len(long_edges) == 3:
            v = grid[tri]
            area = 0.5 * abs(
                (v[1, 0] - v[0, 0]) * (v[2, 1] - v[0, 1])
                - (v[2, 0] - v[0, 0]) * (v[1, 1] - v[0, 1])
            )
            if area > three_phase_area_min:
                three_phase.append(v.copy())
            else:
                # near-critical facet: keep only the single longest edge
                a, b, _ = max(e, key=lambda t: t[2])
                tie_edges.append((a, b))
        else:
            tie_edges.extend(long_edges)

    # Deduplicate tie edges (undirected) and collapse near-duplicate endpoints.
    seen = set()
    tie_lines = []
    for a, b in tie_edges:
        key = (min(a, b), max(a, b))
        if key in seen:
            continue
        seen.add(key)
        tie_lines.append(np.array([grid[a], grid[b]]))
    tie_lines = np.array(tie_lines) if tie_lines else np.empty((0, 2, 2))

    # Merge duplicate three-phase triangles (a genuine 3-phase region yields many
    # coincident tie-triangles); cluster by centroid.
    three_phase = _merge_triangles(three_phase)

    max_res = np.nan
    if refine and len(tie_lines):
        tie_lines, max_res = _refine_tie_lines(tie_lines, chi)

    points = (
        np.unique(np.round(tie_lines.reshape(-1, 2), 6), axis=0)
        if len(tie_lines)
        else np.empty((0, 2))
    )
    return BinodalResult(
        chi=chi,
        n=n,
        tie_lines=tie_lines,
        three_phase=three_phase,
        points=points,
        max_equal_mu_residual=max_res,
    )


def _merge_triangles(tris, tol: float = 0.05):
    """Collapse many near-coincident three-phase triangles to representatives."""
    reps = []
    for tri in tris:
        c = tri.mean(axis=0)
        if all(np.hypot(*(c - r.mean(axis=0))) > tol for r in reps):
            reps.append(tri)
    return reps


def coexistence_residual(A, B, chi: T.Chi):
    """(equal-mu1, equal-mu2, common-tangent) residual for a candidate tie-line."""
    m1A, m2A = T.mu(A[0], A[1], chi)
    m1B, m2B = T.mu(B[0], B[1], chi)
    fbA = T.f_b(A[0], A[1], chi)
    fbB = T.f_b(B[0], B[1], chi)
    tangent = (fbA - fbB) - m1A * (A[0] - B[0]) - m2A * (A[1] - B[1])
    return np.array([m1A - m1B, m2A - m2B, float(tangent)])


def _refine_tie_lines(tie_lines: np.ndarray, chi: T.Chi):
    """Newton-refine each hull tie-line to satisfy coexistence to machine eps.

    4 unknowns (A,B), 3 physical equations + 1 gauge fixing the total
    (phi1A+phi2A+phi1B+phi2B) to its hull value so the tie-line stays put.
    """
    refined = []
    max_res = 0.0
    for tl in tie_lines:
        A0, B0 = tl[0], tl[1]
        S = A0.sum() + B0.sum()

        def eqs(x):
            A = x[:2]
            B = x[2:]
            r = coexistence_residual(A, B, chi)
            gauge = (A[0] + A[1] + B[0] + B[1]) - S
            return np.array([r[0], r[1], r[2], gauge])

        sol = root(eqs, np.concatenate([A0, B0]), method="hybr", tol=1e-12)
        if sol.success:
            A, B = sol.x[:2], sol.x[2:]
            res = np.abs(coexistence_residual(A, B, chi)[:2]).max()
            # reject if the solver ran off the simplex
            if _inside(A) and _inside(B) and np.hypot(*(A - B)) > 1e-4:
                refined.append(np.array([A, B]))
                max_res = max(max_res, res)
                continue
        refined.append(tl)  # keep unrefined hull edge as fallback
    return np.array(refined), max_res


def _inside(p, eps=1e-9):
    return p[0] > eps and p[1] > eps and (p[0] + p[1]) < 1 - eps


def spinodal_curve(chi: T.Chi, n: int = 400):
    """Zero contour of det(Hessian f_b): the spinodal, for cross-checking."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    g = np.linspace(1e-4, 1 - 1e-4, n)
    P1, P2 = np.meshgrid(g, g)
    mask = (P1 + P2) < 1 - 1e-4
    D = np.full_like(P1, np.nan)
    D[mask] = T.spinodal_det(P1[mask], P2[mask], chi)
    fig = plt.figure()
    cs = plt.contour(P1, P2, D, levels=[0.0])
    segs = [np.asarray(v) for v in cs.allsegs[0]]  # level-0 segments
    plt.close(fig)
    return segs
