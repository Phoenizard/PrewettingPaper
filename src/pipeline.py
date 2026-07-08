"""Run one case end to end: bulk maps, both scan directions, pw_line.csv.

Compute only — no matplotlib here (plotting is a separate step, see
plotting.py / scripts/plot_case.py). Progress is printed per scan line and
flushed so parallel logs stream live.
"""

import csv
from pathlib import Path

import numpy as np

import bulk
import scan
from solver import NewtonSolver


def _scan_direction(newton, cfg, demixed_map, source, fixed_values):
    sweep_attr = "phi1_inf" if source == "fix_phi2" else "phi2_inf"
    fixed_attr = "phi2_inf" if source == "fix_phi2" else "phi1_inf"
    rows = []
    n = len(fixed_values)
    for i, fixed in enumerate(fixed_values):
        setattr(newton.p, fixed_attr, float(fixed))
        line = scan.hysteresis_line(newton, cfg.scan, sweep_attr)
        kept = []
        for eq in scan.extract_crossings(line, cfg.criterion):
            phi1, phi2 = (eq, float(fixed)) if source == "fix_phi2" else (float(fixed), eq)
            if not bulk.is_physical_state(phi1, phi2):
                continue
            if bulk.is_point_in_demixed(phi1, phi2, demixed_map, cfg.bulk):
                continue
            kept.append((phi1, phi2))
        rows += [(source, p1, p2) for p1, p2 in kept]
        print(f"[{source} {i + 1}/{n}] fixed={fixed:.6f} -> {len(kept)} point(s)", flush=True)
    return rows


def run_case(cfg, out_dir, max_lines=None):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"[bulk] demixed map, grid={cfg.bulk.demixed_grid}", flush=True)
    demixed_map = bulk.compute_demixed_map(cfg.physical, cfg.bulk)
    print(f"[bulk] binodal points, grid={cfg.bulk.binodal_grid}", flush=True)
    binodal = bulk.compute_binodal_points(cfg.physical, cfg.bulk)
    with open(out_dir / "binodal.csv", "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["phi1", "phi2"])
        writer.writerows(binodal)

    fixed_values = np.linspace(cfg.scan.phi_min, cfg.scan.phi_max, int(cfg.scan.n_lines))
    if max_lines is not None:
        fixed_values = fixed_values[: int(max_lines)]

    newton = NewtonSolver(cfg)
    rows = _scan_direction(newton, cfg, demixed_map, "fix_phi2", fixed_values)
    rows += _scan_direction(newton, cfg, demixed_map, "fix_phi1", fixed_values)

    with open(out_dir / "pw_line.csv", "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["source", "phi1_inf", "phi2_inf"])
        writer.writerows(rows)
    print(f"[case done] {len(rows)} pre-wetting points -> {out_dir / 'pw_line.csv'}", flush=True)
    return rows
