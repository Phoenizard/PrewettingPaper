"""Compare EXP-TA-PANEL-01 results with the archived 150x150 T-a cases."""

import argparse
import csv
import sys
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from geom import min_dist_to_set, mst_length
from paperfig import configure, panel_label, save_pdf_png

CHI = "chi12_0__chi13_2p8__chi23_0"
CHIBB = "chibb11_0__chibb22_0__chibb12_0"
CASES = [
    ("om1_m0p18__om2_m0p3", r"$\omega_1=-0.18$, $\omega_2=-0.30$"),
    ("om1_m0p3__om2_m0p18", r"$\omega_1=-0.30$, $\omega_2=-0.18$"),
]


def read_points(path, xcol, ycol):
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh))
    points = np.array([[float(row[xcol]), float(row[ycol])] for row in rows], dtype=float)
    return points.reshape((-1, 2))


def validate(points, name):
    if not np.isfinite(points).all():
        raise ValueError(f"{name}: non-finite value")
    if ((points < 0.0).any() or (points > 1.0).any()
            or (points[:, 0] + points[:, 1] > 1.0 + 1e-12).any()):
        raise ValueError(f"{name}: composition outside simplex")
    if len(np.unique(points, axis=0)) != len(points):
        raise ValueError(f"{name}: duplicate points")


def measures(points, binodal):
    length, segments, gap, full = mst_length(points, gap_tol=0.015)
    distances = min_dist_to_set(points, binodal)
    return {
        "n_points": len(points),
        "pw_length": length,
        "n_segments": segments,
        "gap_total": gap,
        "full_length": full,
        "dist_mean": float(distances.mean()),
        "dist_min": float(distances.min()),
        "dist_max": float(distances.max()),
    }


def directed_nearest(source, target):
    distances, _ = cKDTree(target).query(source, k=1)
    return float(np.median(distances)), float(np.quantile(distances, 0.95)), float(distances.max())


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--archive-root", default="/root/autodl-fs/pw-space/data")
    parser.add_argument("--new-root", required=True)
    parser.add_argument("--out-dir", required=True)
    args = parser.parse_args()
    configure()
    archive_root = Path(args.archive_root)
    new_root = Path(args.new_root)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    records = []
    plots = []
    for case, title in CASES:
        old_dir = archive_root / CHI / case / CHIBB
        new_dir = new_root / CHI / case / CHIBB
        old = read_points(old_dir / "pw_line.csv", "phi1_inf", "phi2_inf")
        new = read_points(new_dir / "pw_line.csv", "phi1_inf", "phi2_inf")
        old_binodal = read_points(old_dir / "binodal.csv", "phi1", "phi2")
        new_binodal = read_points(new_dir / "binodal.csv", "phi1", "phi2")
        validate(old, f"{case} archive")
        validate(new, f"{case} recompute")
        old_m = measures(old, old_binodal)
        new_m = measures(new, new_binodal)
        new_to_old = directed_nearest(new, old)
        old_to_new = directed_nearest(old, new)
        for source, values in (("archive_150", old_m), ("recompute_400", new_m)):
            records.append({"case": case, "source": source, **values})
        records[-1].update({
            "new_to_old_median": new_to_old[0],
            "new_to_old_p95": new_to_old[1],
            "new_to_old_max": new_to_old[2],
            "old_to_new_median": old_to_new[0],
            "old_to_new_p95": old_to_new[1],
            "old_to_new_max": old_to_new[2],
        })
        plots.append((title, old, new, new_binodal))

    fields = sorted({key for record in records for key in record})
    fields = ["case", "source"] + [key for key in fields if key not in ("case", "source")]
    with open(out_dir / "comparison.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        writer.writerows(records)

    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.25), sharex=True, sharey=True)
    for index, (ax, (title, old, new, binodal)) in enumerate(zip(axes, plots)):
        ax.scatter(binodal[:, 0], binodal[:, 1], s=1.2, color="0.65", linewidth=0,
                   label="binodal", zorder=1)
        ax.scatter(old[:, 0], old[:, 1], s=19, facecolor="none", edgecolor="#2166ac",
                   linewidth=0.75, label="archive 150×150", zorder=2)
        ax.scatter(new[:, 0], new[:, 1], s=10, color="#b2182b", linewidth=0,
                   label="recomputed 400×400", zorder=3)
        ax.set_title(title)
        ax.set_xlabel(r"$\phi_{1,\infty}$")
        panel_label(ax, f"({chr(ord('a') + index)})")
        ax.legend(frameon=False, loc="best", handletextpad=0.35,
                  borderaxespad=0.2, labelspacing=0.3)
    axes[0].set_ylabel(r"$\phi_{2,\infty}$")
    all_new = np.vstack([item[2] for item in plots])
    xpad = 0.08 * np.ptp(all_new[:, 0])
    ypad = 0.08 * np.ptp(all_new[:, 1])
    axes[0].set_xlim(all_new[:, 0].min() - xpad, all_new[:, 0].max() + xpad)
    axes[0].set_ylim(all_new[:, 1].min() - ypad, all_new[:, 1].max() + ypad)
    fig.subplots_adjust(left=0.10, right=0.99, bottom=0.16, top=0.87, wspace=0.12)
    paths = save_pdf_png(fig, out_dir / "archive_vs_recompute.png")
    plt.close(fig)
    print("comparison:", out_dir / "comparison.csv")
    print("figure:", *paths)
    for record in records:
        print(record)


if __name__ == "__main__":
    main()
