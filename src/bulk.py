"""Bulk phase structure: demixed-region map, k-NN vote lookup, binodal points.

Mirrors ternary_ref_code src/thermodynamics.py (compute_bulk_demixed_map at
497-539, is_point_in_demixed_region at 542-568, regular_grid2 / hull_array /
compute_bulk_phase at 76-257 minus the tie-line/fsolve part, which only seeds
an initializer there and is not needed here).

The demixed map feeds the pre-wetting point filter (a crossing inside the
bulk two-phase region is spurious); binodal points are display-only for the
overlay figure.
"""

import numpy as np
from scipy.spatial import ConvexHull

import model
from logutil import log


def compute_demixed_map(p, bulk_cfg):
    eps = 1e-6
    n_grid = int(bulk_cfg.demixed_grid)
    span = np.linspace(eps, 1.0 - eps, n_grid)
    phi1_grid, phi2_grid = np.meshgrid(span, span, indexing="xy")
    mask = (phi1_grid + phi2_grid) < (1.0 - eps)

    x = phi1_grid[mask].ravel()
    y = phi2_grid[mask].ravel()
    f = model.free_energy(x, y, p)
    log(f"[demixed] grid {n_grid}x{n_grid} -> {len(x)} points, hull start")

    try:
        hull = ConvexHull(np.column_stack([x, y, f]))
    except Exception:
        log("[demixed] hull failed -> empty map")
        return {"x": x, "y": y, "demixed": np.zeros_like(x, dtype=bool)}

    eq = hull.equations
    lower = eq[:, 2] < -1e-12
    log(f"[demixed] hull done, {int(lower.sum())} lower facets, envelope eval start")
    if not np.any(lower):
        return {"x": x, "y": y, "demixed": np.zeros_like(x, dtype=bool)}

    nx_l, ny_l, nz_l, d_l = eq[lower, 0], eq[lower, 1], eq[lower, 2], eq[lower, 3]
    z_env = np.full_like(f, -np.inf, dtype=float)
    batch = 2000
    n_batch = (len(x) + batch - 1) // batch
    for bi, i in enumerate(range(0, len(x), batch)):
        j = min(i + batch, len(x))
        planes = -(nx_l[:, None] * x[i:j] + ny_l[:, None] * y[i:j] + d_l[:, None]) / nz_l[:, None]
        z_env[i:j] = np.max(planes, axis=0)
        if (bi + 1) % 10 == 0 or (bi + 1) == n_batch:
            log(f"[demixed] envelope batch {bi + 1}/{n_batch}")

    demixed = (f - z_env) > float(bulk_cfg.demix_tol)
    log(f"[demixed] done, {int(demixed.sum())}/{len(x)} points demixed")
    return {"x": x, "y": y, "demixed": demixed.astype(bool)}


def is_physical_state(phi1, phi2):
    return (
        phi1 > 1e-12
        and phi2 > 1e-12
        and (phi1 + phi2) < (1.0 - 1e-12)
        and phi1 < 1.0
        and phi2 < 1.0
    )


def is_point_in_demixed(phi1, phi2, demixed_map, bulk_cfg):
    if not is_physical_state(phi1, phi2):
        return False
    x, y, demixed = demixed_map["x"], demixed_map["y"], demixed_map["demixed"]
    if len(x) == 0:
        return False
    d2 = (x - float(phi1)) ** 2 + (y - float(phi2)) ** 2
    k = min(max(int(bulk_cfg.k_nearest), 1), len(d2))
    idx = np.argpartition(d2, k - 1)[:k]
    return float(np.mean(demixed[idx].astype(float))) >= float(bulk_cfg.vote_ratio)


def _regular_grid(n_grid):
    # Non-uniform grid: log-spaced near the low boundary, linear elsewhere.
    low, high = 1e-8, 1.0 - 1e-8
    refine = max(100, int(n_grid))
    mid = max(200, int(n_grid) * 2)
    span = np.unique(
        np.concatenate(
            (
                np.logspace(np.log10(max(low, 1e-12)), np.log10(1e-2), refine),
                np.linspace(1e-2, 0.9, mid),
                np.linspace(0.9, min(high, 1.0 - 1e-13), refine),
            )
        )
    )
    phi1_grid, phi2_grid = np.meshgrid(span, span, indexing="xy")
    phi1 = phi1_grid.ravel()
    phi2 = phi2_grid.ravel()
    mask = (phi1 + phi2) < (1.0 - 1e-13)
    return phi1[mask], phi2[mask]


def compute_binodal_points(p, bulk_cfg, vertex_jump_threshold=25):
    """Binodal point set via lower convex envelope + vertex-index jumps."""
    phi1, phi2 = _regular_grid(int(bulk_cfg.binodal_grid))
    f = model.free_energy(phi1, phi2, p)
    log(f"[binodal] grid -> {len(phi1)} points, hull start (may take minutes)")

    n_extra = 3
    max_f = float(np.max(f)) + 1000.0
    pts = np.zeros((len(phi1) + n_extra, 3))
    pts[0] = (0.0, 0.0, max_f)
    pts[1] = (1.0, 0.0, max_f)
    pts[2] = (0.0, 1.0, max_f)
    pts[n_extra:, 0] = phi1
    pts[n_extra:, 1] = phi2
    pts[n_extra:, 2] = f

    try:
        convex = ConvexHull(pts, incremental=True)
    except Exception:
        log("[binodal] hull failed -> no binodal points")
        return []
    log(f"[binodal] hull done, {len(convex.vertices)} vertices")

    hull_vertices = convex.vertices
    if len(hull_vertices) > n_extra:
        vertices = hull_vertices[n_extra:] - n_extra
    else:
        vertices = hull_vertices[hull_vertices >= n_extra] - n_extra
    vertices = vertices[(vertices >= 0) & (vertices < len(phi1))]
    if len(vertices) < 2:
        return []

    bc_x0, bc_y0, bc_x1, bc_y1 = [], [], [], []
    for i in range(len(vertices) - 1):
        idx_i = int(vertices[i])
        idx_j = int(vertices[i + 1])
        if (idx_j - idx_i) > vertex_jump_threshold:
            rho_1, rho_2 = float(phi1[idx_i]), float(phi2[idx_i])
            if abs(1.0 - (rho_1 + rho_2)) > 0.005 and rho_1 > 1e-2 and rho_2 > 1e-2:
                bc_x0.append(rho_1)
                bc_y0.append(rho_2)
            bc_x1.append(float(phi1[idx_j]))
            bc_y1.append(float(phi2[idx_j]))
    points = list(zip(bc_x0 + bc_x1, bc_y0 + bc_y1))
    log(f"[binodal] done, {len(points)} binodal points")
    return points
