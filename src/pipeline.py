"""Run one case end to end: bulk maps, both scan directions, pw_line.csv.

Compute only — no matplotlib here (plotting is a separate step, see
plotting.py / scripts/plot_case.py). Progress is printed per scan line and
flushed so parallel logs stream live.
"""

import csv
import time
from pathlib import Path

import numpy as np

import bulk
import scan
from logutil import log
from solver import NewtonSolver


def _scan_direction(newton, cfg, demixed_map, source, fixed_values):
    sweep_attr = "phi1_inf" if source == "fix_phi2" else "phi2_inf"
    fixed_attr = "phi2_inf" if source == "fix_phi2" else "phi1_inf"
    rows = []
    n = len(fixed_values)
    log(f"[{source}] direction start, {n} lines x {cfg.scan.n_scan} points x 2 branches")
    t_dir = time.time()
    for i, fixed in enumerate(fixed_values):
        t_line = time.time()
        setattr(newton.p, fixed_attr, float(fixed))
        label = f"[{source} {i + 1}/{n} fixed={fixed:.6f}]"
        line = scan.hysteresis_line(newton, cfg.scan, sweep_attr, label=label)
        kept = []
        dropped = 0
        for eq in scan.extract_crossings(line, cfg.criterion):
            phi1, phi2 = (eq, float(fixed)) if source == "fix_phi2" else (float(fixed), eq)
            if not bulk.is_physical_state(phi1, phi2) or bulk.is_point_in_demixed(
                phi1, phi2, demixed_map, cfg.bulk
            ):
                dropped += 1
                continue
            kept.append((phi1, phi2))
        rows += [(source, p1, p2) for p1, p2 in kept]
        log(
            f"{label} line done in {time.time() - t_line:.1f}s -> "
            f"{len(kept)} point(s), {dropped} dropped, total {len(rows)}"
        )
    log(f"[{source}] direction done in {(time.time() - t_dir) / 60:.1f} min, {len(rows)} points")
    return rows


def run_case(cfg, out_dir, max_lines=None):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    t_case = time.time()
    demixed_map = bulk.compute_demixed_map(cfg.physical, cfg.bulk)
    binodal = bulk.compute_binodal_points(cfg.physical, cfg.bulk)
    log(f"[bulk] both maps done in {(time.time() - t_case) / 60:.1f} min")
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
    log(
        f"[case done] {len(rows)} pre-wetting points in "
        f"{(time.time() - t_case) / 60:.1f} min -> {out_dir / 'pw_line.csv'}"
    )
    return rows
