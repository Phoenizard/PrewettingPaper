"""4x2 panel grid of T-a prewetting lines from the pw-space archive.

Top row: omega_2 = -0.30 fixed, one panel per omega_1 in {-0.46, -0.38, -0.30, -0.18}.
Bottom row: omega_1 = -0.30 fixed, one panel per omega_2 in the same four values.
Each panel: binodal (grey) + that case's prewetting line. All panels share the
same axis limits. Publication style: English labels, no embedded title.

Usage:
  python scripts/plot_pw_overlay.py [--data-root DIR] [--out PNG]
"""

import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.cm import viridis

CHI_DIR = "chi12_0__chi13_2p8__chi23_0"
CHIBB_DIR = "chibb11_0__chibb22_0__chibb12_0"
OM_VALUES = ["m0p46", "m0p38", "m0p3", "m0p18"]
OM_LABELS = ["-0.46", "-0.38", "-0.30", "-0.18"]
COLORS = [viridis(t) for t in (0.0, 0.25, 0.5, 0.7)]
MARKERS = ["o", "s", "^", "D"]


def read_xy(path, xcol, ycol):
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f))
    return ([float(r[xcol]) for r in rows], [float(r[ycol]) for r in rows])


def case_dir(root, om1, om2):
    return root / CHI_DIR / f"om1_{om1}__om2_{om2}" / CHIBB_DIR


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-root", default="/root/autodl-fs/pw-space/data")
    ap.add_argument("--out", default="out/analysis/omega/pw_overlay_Ta.png")
    args = ap.parse_args()
    root = Path(args.data_root)

    top = [(om, "m0p3", rf"$\omega_1 = {lab}$")
           for om, lab in zip(OM_VALUES, OM_LABELS)]
    bottom = [("m0p3", om, rf"$\omega_2 = {lab}$")
              for om, lab in zip(OM_VALUES, OM_LABELS)]
    row_notes = [r"$\omega_2 = -0.30$", r"$\omega_1 = -0.30$"]

    bx, by = read_xy(case_dir(root, "m0p3", "m0p3") / "binodal.csv",
                     "phi1", "phi2")

    fig, axes = plt.subplots(2, 4, figsize=(12.0, 6.2),
                             sharex=True, sharey=True)
    xs_all, ys_all = [], []
    for row, series, note in zip(axes, (top, bottom), row_notes):
        for ax, (om1, om2, lab), color, marker in zip(row, series,
                                                      COLORS, MARKERS):
            ax.scatter(bx, by, s=2, c="0.45", rasterized=True)
            px, py = read_xy(case_dir(root, om1, om2) / "pw_line.csv",
                             "phi1_inf", "phi2_inf")
            xs_all += px
            ys_all += py
            ax.scatter(px, py, s=8, marker=marker, color=color, linewidths=0)
            ax.set_title(lab, fontsize=10)
        row[0].text(0.05, 0.92, note, transform=row[0].transAxes,
                    va="top", fontsize=9, color="0.25")

    mx = 0.12 * (max(xs_all) - min(xs_all))
    my = 0.06 * (max(ys_all) - min(ys_all))
    axes[0, 0].set_xlim(min(xs_all) - mx, max(xs_all) + mx)
    axes[0, 0].set_ylim(min(ys_all) - my, max(ys_all) + my)
    axes[0, 0].xaxis.set_major_locator(plt.MaxNLocator(4))
    for ax in axes[1]:
        ax.set_xlabel(r"$\phi_{1,\infty}$")
    for ax in axes[:, 0]:
        ax.set_ylabel(r"$\phi_{2,\infty}$")
    fig.tight_layout()
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=300)
    print(f"[saved] {out}")


if __name__ == "__main__":
    main()
