"""Hysteresis line scan and crossing extraction.

One scan line: fix one far-field component, sweep the other over
linspace(phi_min, phi_max, n_scan) twice — thin branch forward (small
amplitude seed), thick branch backward (large amplitude seed) — with seed
continuation along the line (first point cold-starts from the branch guess,
later points reuse the last converged profile; on failure record NaN and
keep the old seed).

A pre-wetting candidate is a sign change of omega_diff = Omega_thin -
Omega_thick inside a contiguous run of valid points, where valid means
cs_gap = |cs_thick - cs_thin| > cs_threshold with both branches converged
(cs = cs1 + cs2 total). Linear interpolation locates the zero. All
crossings on the line are collected (filtering happens in pipeline.py).

Mechanisms deliberately not adopted (see doc/note/code_cross_comparison.md):
terminal boundary relaxation, auto expand of scan bounds, endpoint bisection
refinement, free-volume scaling of the initial guess.
"""

import numpy as np

_EPS = 1e-12


def _sign_from_wall(omega):
    return 1.0 if float(omega) < 0.0 else -1.0


def _project_feasible(phi1, phi2, eps=_EPS):
    phi1 = np.maximum(phi1, eps)
    phi2 = np.maximum(phi2, eps)
    s = phi1 + phi2
    over = s >= (1.0 - eps)
    if np.any(over):
        phi1[over] = phi1[over] * (1.0 - eps) / s[over]
        phi2[over] = phi2[over] * (1.0 - eps) / s[over]
    phi1 = np.minimum(phi1, 1.0 - 2 * eps)
    phi2 = np.minimum(phi2, 1.0 - 2 * eps - phi1)
    return phi1, phi2


def branch_guess(p, scan_cfg, mode):
    """Linear-decay profile phi_i(z) = phi_i_inf + sign_i * amp_i * (1 - z/L)."""
    z = np.linspace(0.0, float(p.L), int(p.N))
    layer = 1.0 - z / float(p.L)
    if mode == "thick":
        amp1, amp2 = scan_cfg.amp_phi1_thick, scan_cfg.amp_phi2_thick
    else:
        amp1, amp2 = scan_cfg.amp_phi1_thin, scan_cfg.amp_phi2_thin
    phi1 = p.phi1_inf + _sign_from_wall(p.omega_1) * amp1 * layer
    phi2 = p.phi2_inf + _sign_from_wall(p.omega_2) * amp2 * layer
    p1, p2 = _project_feasible(phi1, phi2)
    return np.concatenate([p1, p2])


def _sweep_branch(solver, scan_cfg, sweep_attr, values, mode):
    """One directional sweep. Mutates solver.p.<sweep_attr> point by point."""
    p = solver.p
    n = len(values)
    omega = np.full(n, np.nan)
    cs = np.full(n, np.nan)
    seed = None
    for i, val in enumerate(values):
        setattr(p, sweep_attr, float(val))
        U0 = seed if seed is not None else branch_guess(p, scan_cfg, mode)
        U, ok = solver.solve(U0)
        if ok:
            seed = U
            om, cs1, cs2 = solver.surface_metrics(U)
            omega[i] = om
            cs[i] = cs1 + cs2
    return omega, cs


def hysteresis_line(solver, scan_cfg, sweep_attr):
    """Bidirectional scan of solver.p.<sweep_attr> ('phi1_inf' or 'phi2_inf').

    Returns dict with scan values and per-branch Omega / total cs arrays.
    """
    values = np.linspace(scan_cfg.phi_min, scan_cfg.phi_max, int(scan_cfg.n_scan))

    omega_thin, cs_thin = _sweep_branch(solver, scan_cfg, sweep_attr, values, "thin")
    omega_thick_rev, cs_thick_rev = _sweep_branch(
        solver, scan_cfg, sweep_attr, values[::-1], "thick"
    )
    return {
        "values": values,
        "omega_thin": omega_thin,
        "cs_thin": cs_thin,
        "omega_thick": omega_thick_rev[::-1],
        "cs_thick": cs_thick_rev[::-1],
    }


def extract_crossings(line, crit):
    """All interpolated zero crossings of omega_diff within contiguous valid runs."""
    omega_diff = line["omega_thin"] - line["omega_thick"]
    cs_gap = np.abs(line["cs_thick"] - line["cs_thin"])
    values = line["values"]

    valid = np.isfinite(omega_diff) & np.isfinite(cs_gap) & (cs_gap > float(crit.cs_threshold))
    if int(valid.sum()) < int(crit.min_hysteresis_points):
        return []

    crossings = []
    idx = np.where(valid)[0]
    runs = np.split(idx, np.where(np.diff(idx) > 1)[0] + 1)
    for run in runs:
        if len(run) < 2:
            continue
        vals = omega_diff[run]
        for k in np.where(np.diff(np.sign(vals)) != 0)[0]:
            a, b = run[k], run[k + 1]
            va, vb = omega_diff[a], omega_diff[b]
            if vb == va:
                continue
            crossings.append(float(values[a] - va * (values[b] - values[a]) / (vb - va)))
    return crossings
