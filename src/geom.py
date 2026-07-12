"""Geometry tools for pre-wetting line measures.

All measures live in the raw (phi1_inf, phi2_inf) plane. Both axes are volume
fractions in [0, 1] with the same unit, so Euclidean distances and lengths are
physically meaningful without any rescaling. Do not normalise the axes.

Shared by every analysis angle (omega, chibb12, chibb11/22, chi): the caller
passes a case's pre-wetting points (fix_phi1 + fix_phi2 merged; the two are only
scan directions, not physical branches) and its binodal points, and gets back
length, distance-to-binodal, straight-line residual, and a branch split.
"""

import numpy as np
from scipy.spatial import cKDTree


def _as_points(points):
    """Coerce to an (n, 2) float array of (phi1, phi2). Empty input -> (0, 2)."""
    arr = np.asarray(points, dtype=float)
    if arr.size == 0:
        return arr.reshape(0, 2)
    if arr.ndim != 2 or arr.shape[1] != 2:
        raise ValueError(f"points must be (n, 2), got shape {arr.shape}")
    return arr


def fit_line_pca(points):
    """Best-fit line through the points by orthogonal (total) least squares.

    Uses PCA rather than y = k*x + b so a near-vertical pre-wetting line is fit
    just as stably as a horizontal one. Returns (centroid, direction) where
    direction is a unit vector along the line's principal axis. Needs >= 2
    distinct points; raises ValueError otherwise.
    """
    pts = _as_points(points)
    if len(pts) < 2:
        raise ValueError("need at least 2 points to fit a line")
    centroid = pts.mean(axis=0)
    centred = pts - centroid
    # Principal axis = eigenvector of the largest singular value.
    _, _, vh = np.linalg.svd(centred, full_matrices=False)
    direction = vh[0]
    norm = np.linalg.norm(direction)
    if norm == 0.0:
        raise ValueError("degenerate point set (all points coincide)")
    return centroid, direction / norm


def projection_span(points, line=None):
    """Length = extent of the points projected onto their best-fit line.

    Projects every point onto the fitted line's direction and returns the
    max-minus-min span of those projections. This measures how far the
    pre-wetting line reaches along its own direction and is insensitive to
    perpendicular scatter of individual points. Returns 0.0 for a single point.
    """
    pts = _as_points(points)
    if len(pts) < 2:
        return 0.0
    centroid, direction = line if line is not None else fit_line_pca(pts)
    t = (pts - centroid) @ direction
    return float(t.max() - t.min())


def residual_rms(points, line=None):
    """RMS of perpendicular distances from the points to their best-fit line.

    This is the branch/kink indicator: a clean single straight line has a small
    residual (near the data noise floor); a kinked or multi-branch point set
    fits one line poorly and gives a large residual. Returns 0.0 for < 3 points
    (two points always fit a line exactly).
    """
    pts = _as_points(points)
    if len(pts) < 3:
        return 0.0
    centroid, direction = line if line is not None else fit_line_pca(pts)
    centred = pts - centroid
    # Perpendicular component = full vector minus its projection on direction.
    proj = np.outer(centred @ direction, direction)
    perp = centred - proj
    d = np.linalg.norm(perp, axis=1)
    return float(np.sqrt(np.mean(d ** 2)))


def min_dist_to_set(points, target):
    """Per-point shortest Euclidean distance from each point to a target set.

    Used for the pre-wetting-line-to-binodal distance: measured on the raw
    pre-wetting points (not the fitted line), so the point count does not bias
    the result. The caller takes mean() and max() of the returned array.
    Returns an empty array if there are no points; raises if target is empty.
    """
    pts = _as_points(points)
    tgt = _as_points(target)
    if len(pts) == 0:
        return np.empty(0, dtype=float)
    if len(tgt) == 0:
        raise ValueError("target set is empty; cannot measure distance")
    tree = cKDTree(tgt)
    d, _ = tree.query(pts, k=1)
    return np.asarray(d, dtype=float)


def split_branches(points, link_tol=0.01):
    """Label pre-wetting points into branches by sweeping phi2 low -> high.

    Reimplements the reference algorithm (prewetting_project_clean/main.py
    _label_prewetting_branches): walk phi2 slices in increasing order; within a
    slice, link each active branch to its nearest unused point in phi1 if within
    link_tol, otherwise close it; leftover points in the slice open new branches.
    A single unbroken line yields one branch; a kink or a gap opens extra ones.

    Returns an int array of branch labels aligned with the input rows.
    Used only to flag/colour multi-branch cases; rigorous conclusions come from
    single-branch cases.
    """
    pts = _as_points(points)
    n = len(pts)
    if n == 0:
        return np.empty(0, dtype=int)

    # Group point indices by phi2 slice; sort slices ascending, points by phi1.
    by_slice = {}
    for idx, (phi1, phi2) in enumerate(pts):
        by_slice.setdefault(float(phi2), []).append((float(phi1), idx))
    phi2_sorted = sorted(by_slice)
    for p2 in phi2_sorted:
        by_slice[p2].sort()

    labels = np.full(n, -1, dtype=int)
    next_id = 0
    active = {}  # branch id -> last phi1 seen

    for p2 in phi2_sorted:
        slice_pts = by_slice[p2]
        used = set()
        for bid in list(active):
            last_phi1 = active[bid]
            best_j, best_d = None, None
            for j, (phi1, _idx) in enumerate(slice_pts):
                if j in used:
                    continue
                d = abs(phi1 - last_phi1)
                if best_d is None or d < best_d:
                    best_d, best_j = d, j
            if best_j is not None and best_d <= link_tol:
                phi1, idx = slice_pts[best_j]
                used.add(best_j)
                active[bid] = phi1
                labels[idx] = bid
            else:
                del active[bid]
        for j, (phi1, idx) in enumerate(slice_pts):
            if j in used:
                continue
            active[next_id] = phi1
            labels[idx] = next_id
            next_id += 1

    return labels
