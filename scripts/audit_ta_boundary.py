"""Lightweight signed-side audit for displayed T-a re-entry curves.

This script reads existing binodal and prewetting CSV files only.  It does not
run the boundary-value solver.  A local tangent to the binodal defines the
normal; the normal is oriented toward the dilute origin, which is the one-phase
side for the displayed solvent-rich branch.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import cKDTree


CASES = [
    ("m0p5", -0.50, "strong"),
    ("m0p48", -0.48, "strong"),
    ("m0p46", -0.46, "strong"),
    ("m0p44", -0.44, "strong"),
    ("m0p1", -0.10, "weak"),
    ("m0p08", -0.08, "weak"),
    ("m0p06", -0.06, "weak"),
    ("m0p04", -0.04, "weak"),
]
NO_LINE_CONTROLS = [("m0p02", -0.02), ("m0p01", -0.01), ("m0p001", -0.001)]
CHI = "chi12_0__chi13_2p8__chi23_0"
CHIBB = "chibb11_0__chibb22_0__chibb12_0"


def read_xy(path: Path, x: str, y: str) -> np.ndarray:
    with path.open(newline="") as fh:
        rows = list(csv.DictReader(fh))
    if not rows:
        return np.empty((0, 2), dtype=float)
    return np.asarray([(float(r[x]), float(r[y])) for r in rows], dtype=float)


def case_dir(root: Path, code: str) -> Path:
    return root / CHI / f"om1_{code}__om2_{code}" / CHIBB


def physical_binodal_branch(points: np.ndarray) -> np.ndarray:
    """Keep the leading smooth solvent-rich branch, before hull edge vertices.

    The binodal CSV appends simplex-boundary convex-hull vertices after the
    physical branch.  The transition is a macroscopic jump, while neighbouring
    points on the physical branch are closely spaced.
    """
    if len(points) < 2:
        return points
    jumps = np.linalg.norm(np.diff(points, axis=0), axis=1)
    cut_candidates = np.where(jumps > 0.05)[0]
    if not len(cut_candidates):
        return points
    return points[: int(cut_candidates[0]) + 1]


def local_signed_distance(point: np.ndarray, binodal: np.ndarray, tree: cKDTree):
    # Nine nearby binodal points determine a local tangent by PCA.
    distances, indices = tree.query(point, k=min(9, len(binodal)))
    local = binodal[np.atleast_1d(indices)]
    centre = local.mean(axis=0)
    _, _, vh = np.linalg.svd(local - centre, full_matrices=False)
    tangent = vh[0]
    normal = np.asarray([-tangent[1], tangent[0]])
    nearest = binodal[int(np.atleast_1d(indices)[0])]
    # Orient toward the dilute origin: positive means the one-phase side.
    if float(normal @ (-nearest)) < 0.0:
        normal *= -1.0
    signed = float(normal @ (point - nearest))

    # Estimate local binodal sampling uncertainty from the same neighbourhood.
    local_tree = cKDTree(local)
    local_nn = local_tree.query(local, k=min(2, len(local)))[0]
    if local_nn.ndim == 2 and local_nn.shape[1] > 1:
        spacing = float(np.median(local_nn[:, 1]))
    else:
        spacing = 0.0
    tolerance = max(2.0 * spacing, 1.0e-4)
    return float(np.atleast_1d(distances)[0]), signed, tolerance, nearest


def longest_true_run(flags: list[bool]) -> int:
    best = current = 0
    for flag in flags:
        current = current + 1 if flag else 0
        best = max(best, current)
    return best


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-root", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    root = Path(args.data_root)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    binodal_all = read_xy(case_dir(root, "m0p3") / "binodal.csv", "phi1", "phi2")
    binodal = physical_binodal_branch(binodal_all)
    tree = cKDTree(binodal)
    point_rows: list[dict] = []
    summaries: list[dict] = []

    for code, value, side in CASES:
        points = read_xy(case_dir(root, code) / "pw_line.csv", "phi1_inf", "phi2_inf")
        records = []
        for point in points:
            distance, signed, tolerance, nearest = local_signed_distance(point, binodal, tree)
            if signed < -tolerance:
                classification = "two_phase_side"
            elif signed > tolerance:
                classification = "one_phase_side"
            else:
                classification = "boundary_tolerance"
            row = {
                "code": code,
                "wall_affinity": value,
                "side": side,
                "phi1_inf": float(point[0]),
                "phi2_inf": float(point[1]),
                "nearest_binodal_phi1": float(nearest[0]),
                "nearest_binodal_phi2": float(nearest[1]),
                "distance": distance,
                "signed_distance": signed,
                "tolerance": tolerance,
                "classification": classification,
            }
            records.append(row)
            point_rows.append(row)

        records.sort(key=lambda r: (r["phi2_inf"], r["phi1_inf"]))
        inside_flags = [r["classification"] == "two_phase_side" for r in records]
        n = len(records)
        n_inside = sum(inside_flags)
        n_boundary = sum(r["classification"] == "boundary_tolerance" for r in records)
        n_outside = sum(r["classification"] == "one_phase_side" for r in records)
        fraction = n_inside / n if n else 0.0
        longest = longest_true_run(inside_flags)
        if n < 5:
            status, reason = "invalid", "fewer_than_five_points"
        elif fraction > 0.05:
            status, reason = "invalid", "more_than_five_percent_on_two_phase_side"
        elif longest >= 2:
            status, reason = "invalid", "contiguous_two_phase_side_segment"
        else:
            status, reason = "valid", "one_phase_or_boundary_within_tolerance"
        summaries.append({
            "code": code,
            "wall_affinity": value,
            "side": side,
            "n_points": n,
            "n_one_phase": n_outside,
            "n_boundary": n_boundary,
            "n_two_phase": n_inside,
            "two_phase_fraction": fraction,
            "longest_two_phase_run": longest,
            "min_signed_distance": min((r["signed_distance"] for r in records), default=math.nan),
            "mean_signed_distance": float(np.mean([r["signed_distance"] for r in records])) if records else math.nan,
            "max_signed_distance": max((r["signed_distance"] for r in records), default=math.nan),
            "status": status,
            "reason": reason,
        })

    for code, value in NO_LINE_CONTROLS:
        points = read_xy(case_dir(root, code) / "pw_line.csv", "phi1_inf", "phi2_inf")
        summaries.append({
            "code": code, "wall_affinity": value, "side": "weak_control",
            "n_points": len(points), "n_one_phase": 0, "n_boundary": 0,
            "n_two_phase": 0, "two_phase_fraction": 0.0,
            "longest_two_phase_run": 0, "min_signed_distance": math.nan,
            "mean_signed_distance": math.nan, "max_signed_distance": math.nan,
            "status": "no_prewetting", "reason": "no_crossings_in_existing_result",
        })

    point_fields = list(point_rows[0]) if point_rows else []
    with (out / "point_classification.csv").open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=point_fields)
        writer.writeheader(); writer.writerows(point_rows)
    summary_fields = list(summaries[0])
    with (out / "curve_classification.csv").open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=summary_fields)
        writer.writeheader(); writer.writerows(summaries)
    (out / "classification.json").write_text(json.dumps(summaries, indent=2, allow_nan=True) + "\n")

    fig, axes = plt.subplots(1, 2, figsize=(8.0, 3.4))
    for ax, side, title in zip(axes, ("strong", "weak"), ("strong-attraction audit", "weak-attraction audit")):
        ax.scatter(binodal[:, 0], binodal[:, 1], s=1.0, c="0.65", label="binodal", zorder=1)
        for code, value, case_side in CASES:
            if case_side != side:
                continue
            rows = [r for r in point_rows if r["code"] == code]
            good = np.asarray([(r["phi1_inf"], r["phi2_inf"]) for r in rows if r["classification"] != "two_phase_side"])
            bad = np.asarray([(r["phi1_inf"], r["phi2_inf"]) for r in rows if r["classification"] == "two_phase_side"])
            if len(good):
                ax.scatter(good[:, 0], good[:, 1], s=11, label=f"{value:g} accepted")
            if len(bad):
                ax.scatter(bad[:, 0], bad[:, 1], s=15, marker="x", label=f"{value:g} rejected")
        ax.set_title(title); ax.set_xlabel(r"$\phi_{1,\infty}$"); ax.set_ylabel(r"$\phi_{2,\infty}$")
        ax.grid(alpha=0.15); ax.legend(frameon=False, fontsize=7, ncol=2)
    axes[0].set_xlim(0.072, 0.090); axes[0].set_ylim(-0.001, 0.035)
    axes[1].set_xlim(0.090, 0.190); axes[1].set_ylim(0.075, 0.203)
    fig.tight_layout()
    fig.savefig(out / "boundary_audit.png", dpi=220)
    fig.savefig(out / "boundary_audit.pdf")
    plt.close(fig)

    print(json.dumps(summaries, indent=2, allow_nan=True))


if __name__ == "__main__":
    main()
