"""Publication and supplementary T-a wall-affinity figures from archived data.

The main figure overlays four pre-wetting lines in each of two panels: varying
omega_1 at fixed omega_2, and varying omega_2 at fixed omega_1.  The original
case-by-case 2x4 layout is retained as a supplementary figure.
"""

import argparse
import csv
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.cm import viridis

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from paperfig import configure, panel_label, save_pdf_png

CHI_DIR = "chi12_0__chi13_2p8__chi23_0"
CHIBB_DIR = "chibb11_0__chibb22_0__chibb12_0"
OM_VALUES = ["m0p46", "m0p38", "m0p3", "m0p18"]
OM_FLOATS = [-0.46, -0.38, -0.30, -0.18]
COLORS = [viridis(t) for t in (0.08, 0.36, 0.64, 0.92)]
MARKERS = ["o", "s", "^", "D"]
BINODAL_COLOR = "0.55"


def read_xy(path, xcol, ycol):
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh))
    return ([float(r[xcol]) for r in rows], [float(r[ycol]) for r in rows])


def case_dir(root, om1, om2):
    return root / CHI_DIR / f"om1_{om1}__om2_{om2}" / CHIBB_DIR


def load_series(root):
    top = []
    bottom = []
    for encoded, value in zip(OM_VALUES, OM_FLOATS):
        px, py = read_xy(case_dir(root, encoded, "m0p3") / "pw_line.csv",
                         "phi1_inf", "phi2_inf")
        top.append((value, px, py))
        px, py = read_xy(case_dir(root, "m0p3", encoded) / "pw_line.csv",
                         "phi1_inf", "phi2_inf")
        bottom.append((value, px, py))
    bx, by = read_xy(case_dir(root, "m0p3", "m0p3") / "binodal.csv",
                     "phi1", "phi2")
    return bx, by, top, bottom


def common_limits(series):
    xs = [x for group in series for _, px, _ in group for x in px]
    ys = [y for group in series for _, _, py in group for y in py]
    dx = max(xs) - min(xs)
    dy = max(ys) - min(ys)
    return (min(xs) - 0.10 * dx, max(xs) + 0.10 * dx), \
           (min(ys) - 0.06 * dy, max(ys) + 0.06 * dy)


def draw_binodal(ax, bx, by):
    ax.scatter(bx, by, s=1.2, color=BINODAL_COLOR, linewidths=0,
               label="binodal", zorder=1)


def draw_overlay(ax, bx, by, series, varied):
    draw_binodal(ax, bx, by)
    for (value, px, py), color, marker in zip(series, COLORS, MARKERS):
        ax.scatter(px, py, s=14, marker=marker, facecolor=color,
                   edgecolor="white", linewidth=0.25,
                   label=rf"$\omega_{varied}={value:.2f}$", zorder=3)
    ax.set_xlabel(r"$\phi_{1,\infty}$")
    ax.legend(loc="best", frameon=False, handletextpad=0.35,
              borderaxespad=0.2, labelspacing=0.3)


def main_figure(bx, by, top, bottom, output):
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.25), sharex=True, sharey=True)
    draw_overlay(axes[0], bx, by, top, 1)
    draw_overlay(axes[1], bx, by, bottom, 2)
    axes[0].set_ylabel(r"$\phi_{2,\infty}$")
    axes[0].set_title(r"varying $\omega_1$; $\omega_2=-0.30$")
    axes[1].set_title(r"varying $\omega_2$; $\omega_1=-0.30$")
    panel_label(axes[0], "(a)")
    panel_label(axes[1], "(b)")
    xlim, ylim = common_limits((top, bottom))
    axes[0].set_xlim(*xlim)
    axes[0].set_ylim(*ylim)
    for ax in axes:
        ax.xaxis.set_major_locator(plt.MaxNLocator(4))
        ax.yaxis.set_major_locator(plt.MaxNLocator(5))
    fig.subplots_adjust(left=0.10, right=0.99, bottom=0.16, top=0.87, wspace=0.12)
    return save_pdf_png(fig, output)


def supplement_figure(bx, by, top, bottom, output):
    fig, axes = plt.subplots(2, 4, figsize=(10.8, 5.35), sharex=True, sharey=True)
    xlim, ylim = common_limits((top, bottom))
    for row, series, varied in zip(axes, (top, bottom), (1, 2)):
        for ax, (value, px, py), color, marker in zip(row, series, COLORS, MARKERS):
            draw_binodal(ax, bx, by)
            ax.scatter(px, py, s=12, marker=marker, color=color, linewidths=0,
                       zorder=3)
            ax.set_title(rf"$\omega_{varied}={value:.2f}$")
            ax.set_xlim(*xlim)
            ax.set_ylim(*ylim)
    for ax in axes[1]:
        ax.set_xlabel(r"$\phi_{1,\infty}$")
    for ax in axes[:, 0]:
        ax.set_ylabel(r"$\phi_{2,\infty}$")
    axes[0, 0].text(0.03, 0.96, r"$\omega_2=-0.30$", transform=axes[0, 0].transAxes,
                    va="top")
    axes[1, 0].text(0.03, 0.96, r"$\omega_1=-0.30$", transform=axes[1, 0].transAxes,
                    va="top")
    fig.tight_layout(pad=0.7, w_pad=0.5, h_pad=0.7)
    return save_pdf_png(fig, output)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-root", default="/root/autodl-fs/pw-space/data")
    ap.add_argument("--out", default="out/analysis/omega/pw_overlay_Ta.png")
    args = ap.parse_args()
    configure()
    bx, by, top, bottom = load_series(Path(args.data_root))
    output = Path(args.out)
    main_paths = main_figure(bx, by, top, bottom, output)
    supp = output.with_name(output.stem + "_supplement_2x4" + output.suffix)
    supp_paths = supplement_figure(bx, by, top, bottom, supp)
    plt.close("all")
    print("main:", *main_paths)
    print("supplement:", *supp_paths)


if __name__ == "__main__":
    main()
