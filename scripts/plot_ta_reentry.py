"""Publication figure and metrics table for EXP-TA-REENTRY-01."""

import argparse
import csv
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from paperfig import configure, panel_label, save_pdf_png


OUT_FIELDS = [
    "lambda", "prewetting_exists", "n_points", "n_segments", "pw_length",
    "dist_mean", "dist_min", "endpoint_1_phi1", "endpoint_1_phi2",
    "endpoint_2_phi1", "endpoint_2_phi2", "flag", "rel",
]


def read_rows(path):
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh))
    for row in rows:
        row["lambda"] = float(row["omega_1"])
        if row["lambda"] != float(row["omega_2"]):
            raise ValueError(f"non-diagonal case in input: {row['rel']}")
    return sorted(rows, key=lambda row: row["lambda"])


def read_points(path):
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)
        return np.asarray(
            [(float(row["phi1_inf"]), float(row["phi2_inf"])) for row in reader],
            dtype=float,
        ).reshape(-1, 2)


def farthest_pair(points):
    if len(points) == 0:
        return ("", "", "", "")
    if len(points) == 1:
        p = points[0]
        return (p[0], p[1], p[0], p[1])
    delta = points[:, None, :] - points[None, :, :]
    d2 = np.sum(delta * delta, axis=2)
    i, j = np.unravel_index(np.argmax(d2), d2.shape)
    return (points[i, 0], points[i, 1], points[j, 0], points[j, 1])


def build_metrics(rows, data_root):
    output = []
    for row in rows:
        points = read_points(data_root / row["rel"] / "pw_line.csv")
        endpoints = farthest_pair(points)
        exists = len(points) > 0
        output.append({
            "lambda": f"{row['lambda']:.6g}",
            "prewetting_exists": "1" if exists else "0",
            "n_points": row["n_points"],
            "n_segments": row["n_segments"],
            "pw_length": row["pw_length"],
            "dist_mean": row["dist_mean"],
            "dist_min": row["dist_min"],
            "endpoint_1_phi1": endpoints[0],
            "endpoint_1_phi2": endpoints[1],
            "endpoint_2_phi1": endpoints[2],
            "endpoint_2_phi2": endpoints[3],
            "flag": row["flag"],
            "rel": row["rel"],
        })
    return output


def write_metrics(rows, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=OUT_FIELDS)
        writer.writeheader()
        writer.writerows(rows)


def group_for(value):
    if value <= -0.44:
        return "strong"
    if value < -0.10:
        return "intermediate"
    return "weak"


def make_figure(rows, output):
    configure()
    x = np.asarray([float(row["lambda"]) for row in rows])
    exists = np.asarray([row["prewetting_exists"] == "1" for row in rows])
    length = np.asarray([float(row["pw_length"]) if row["pw_length"] else np.nan for row in rows])
    dmean = np.asarray([float(row["dist_mean"]) if row["dist_mean"] else np.nan for row in rows])
    dmin = np.asarray([float(row["dist_min"]) if row["dist_min"] else np.nan for row in rows])

    colors = {"strong": "#D55E00", "intermediate": "#0072B2", "weak": "#009E73"}
    markers = {"strong": "s", "intermediate": "o", "weak": "D"}
    groups = np.asarray([group_for(value) for value in x])

    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.15))
    ax = axes[0]
    ax.plot(x[exists], length[exists], color="0.45", lw=1.2, zorder=1)
    for group in ("strong", "intermediate", "weak"):
        take = exists & (groups == group)
        ax.scatter(x[take], length[take], s=31, marker=markers[group],
                   facecolor=colors[group], edgecolor="white", linewidth=0.55,
                   label=group, zorder=3)
    ax.scatter(x[~exists], np.zeros((~exists).sum()), marker="x", s=35,
               color="0.25", linewidth=1.2, label="no pre-wetting", zorder=4)
    ax.set_xlabel(r"diagonal wall affinity $\lambda$ ($\omega_1=\omega_2=\lambda$)")
    ax.set_ylabel(r"pre-wetting-line length $L$")
    ax.set_ylim(bottom=-0.004)
    panel_label(ax, "(a)")

    ax = axes[1]
    ax.plot(x[exists], dmean[exists], color="0.20", lw=1.35,
            marker="o", ms=4.2, label=r"mean distance $\bar d$", zorder=2)
    ax.plot(x[exists], dmin[exists], color="0.50", lw=1.25, ls="--",
            marker="^", ms=4.0, label=r"minimum distance $d_{\min}$", zorder=2)
    for group in ("strong", "intermediate", "weak"):
        take = exists & (groups == group)
        ax.scatter(x[take], dmean[take], s=26, marker=markers[group],
                   facecolor=colors[group], edgecolor="white", linewidth=0.45,
                   zorder=3)
    ax.scatter(x[~exists], np.zeros((~exists).sum()), marker="x", s=35,
               color="0.25", linewidth=1.2, zorder=4)
    ax.set_xlabel(r"diagonal wall affinity $\lambda$ ($\omega_1=\omega_2=\lambda$)")
    ax.set_ylabel("distance to binodal")
    ax.set_ylim(bottom=-0.0007)
    panel_label(ax, "(b)")
    ax.legend(frameon=False, loc="upper right")

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, frameon=False, ncol=4, loc="lower center",
               bbox_to_anchor=(0.5, -0.01))
    for ax in axes:
        ax.grid(alpha=0.18, linewidth=0.6)
        ax.tick_params(top=False, right=False)
    fig.subplots_adjust(left=0.10, right=0.985, bottom=0.26, top=0.92, wspace=0.31)
    return save_pdf_png(fig, output)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--measures", required=True)
    parser.add_argument("--data-root", required=True)
    parser.add_argument("--out-dir", required=True)
    args = parser.parse_args()
    out_dir = Path(args.out_dir)
    rows = read_rows(args.measures)
    metrics = build_metrics(rows, Path(args.data_root))
    write_metrics(metrics, out_dir / "reentry_metrics.csv")
    paths = make_figure(metrics, out_dir / "ta_reentry.png")
    print(f"cases: {len(metrics)}")
    print("figure:", *paths)


if __name__ == "__main__":
    main()
