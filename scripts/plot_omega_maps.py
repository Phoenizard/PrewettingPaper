"""Supplementary T-a extent heatmaps from the established measures table."""

import argparse
import csv
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from paperfig import configure, panel_label, save_pdf_png


def load_rows(path, stage):
    rows = []
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            if row["stage"] != stage:
                continue
            if any(float(row[key]) != 0.0 for key in
                   ("chi_bb_11", "chi_bb_22", "chi_bb_12")):
                continue
            rows.append(row)
    return rows


def make_grid(rows, field):
    omega1 = np.array(sorted({float(row["omega_1"]) for row in rows}))
    omega2 = np.array(sorted({float(row["omega_2"]) for row in rows}))
    i1 = {value: index for index, value in enumerate(omega1)}
    i2 = {value: index for index, value in enumerate(omega2)}
    grid = np.full((len(omega2), len(omega1)), np.nan)
    sampled = np.zeros_like(grid, dtype=bool)
    no_prewetting = np.zeros_like(grid, dtype=bool)
    for row in rows:
        y = i2[float(row["omega_2"])]
        x = i1[float(row["omega_1"])]
        sampled[y, x] = True
        if row["flag"].startswith("no_prewetting") or row.get(field, "") == "":
            no_prewetting[y, x] = True
        else:
            grid[y, x] = float(row[field])
    return omega1, omega2, grid, sampled, no_prewetting


def heatmap(ax, omega1, omega2, grid, sampled, no_prewetting,
            title, colorbar_label, cmap):
    ax.set_facecolor("0.92")
    image = ax.pcolormesh(omega1, omega2, np.ma.masked_invalid(grid),
                         shading="nearest", cmap=cmap)
    yy, xx = np.where(no_prewetting)
    if len(xx):
        ax.scatter(omega1[xx], omega2[yy], marker="s", s=92,
                   facecolor="white", edgecolor="0.35", linewidth=0.5, zorder=3)
        ax.scatter(omega1[xx], omega2[yy], marker="x", s=28,
                   color="0.20", linewidth=1.1, zorder=4)
    ax.set_xlabel(r"$\omega_1$")
    ax.set_ylabel(r"$\omega_2$")
    ax.set_title(title)
    ax.set_aspect("equal")
    colorbar = plt.colorbar(image, ax=ax, fraction=0.046, pad=0.04)
    colorbar.set_label(colorbar_label)
    return bool((~sampled).any()), bool(no_prewetting.any())


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--measures", required=True)
    parser.add_argument("--stage", default="chi12_0__chi13_2p8__chi23_0")
    parser.add_argument("--out-dir", required=True)
    args = parser.parse_args()
    configure()
    rows = load_rows(args.measures, args.stage)
    if not rows:
        raise SystemExit(f"no matching rows for stage {args.stage}")

    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.25))
    omega1, omega2, length, sampled_l, none_l = make_grid(rows, "pw_length")
    has_unsampled, has_none = heatmap(
        axes[0], omega1, omega2, length, sampled_l, none_l,
        r"pre-wetting-line length $L$", r"$L$", "viridis")
    omega1, omega2, distance, sampled_d, none_d = make_grid(rows, "dist_mean")
    unsampled_d, none_d_present = heatmap(
        axes[1], omega1, omega2, distance, sampled_d, none_d,
        r"mean distance to binodal $\bar d$", r"$\bar d$", "viridis_r")
    panel_label(axes[0], "(a)")
    panel_label(axes[1], "(b)")

    handles = []
    if has_none or none_d_present:
        handles.append(Line2D([], [], marker="x", ls="", color="0.20",
                              label="sampled; no pre-wetting"))
    if has_unsampled or unsampled_d:
        handles.append(Patch(facecolor="0.92", edgecolor="0.55",
                             label="not sampled"))
    if handles:
        fig.legend(handles=handles, loc="lower center", ncol=len(handles),
                   frameon=False, bbox_to_anchor=(0.5, -0.01))
        bottom = 0.22
    else:
        bottom = 0.16
    fig.subplots_adjust(left=0.09, right=0.98, bottom=bottom, top=0.87, wspace=0.34)
    paths = save_pdf_png(fig, Path(args.out_dir) / "omega_length_dist.png")
    plt.close(fig)
    print(f"T-a cases: {len(rows)}")
    print("figure:", *paths)


if __name__ == "__main__":
    main()
