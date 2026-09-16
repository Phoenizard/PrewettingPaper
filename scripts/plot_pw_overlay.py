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
import numpy as np
from scipy.interpolate import PchipInterpolator

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


def physical_binodal_branch(bx, by):
    points = np.column_stack([bx, by])
    jumps = np.linalg.norm(np.diff(points, axis=0), axis=1)
    cut = np.where(jumps > 0.05)[0]
    if len(cut):
        points = points[: int(cut[0]) + 1]
    return points[:, 0], points[:, 1]


def smooth_monotone_y(px, py, n=500):
    by_y = {}
    for x, y in zip(px, py):
        by_y.setdefault(float(y), []).append(float(x))
    y = np.asarray(sorted(by_y))
    x = np.asarray([np.mean(by_y[v]) for v in y])
    yd = np.linspace(y[0], y[-1], n)
    if len(y) == 2:
        return np.interp(yd, y, x), yd
    return PchipInterpolator(y, x)(yd), yd


def smooth_curve(px, py, n=240):
    """Smooth a nearly straight phase-boundary segment without scan-order kinks."""
    points = np.column_stack([px, py]).astype(float)
    centre = points.mean(axis=0)
    _, _, vh = np.linalg.svd(points - centre, full_matrices=False)
    direction = vh[0]
    if direction[1] < 0:
        direction *= -1
    normal = np.asarray([-direction[1], direction[0]])
    t = (points - centre) @ direction
    q = (points - centre) @ normal
    degree = min(3, len(points) - 1)
    coeff = np.polyfit(t, q, degree)
    td = np.linspace(t.min(), t.max(), n)
    qd = np.polyval(coeff, td)
    curve = centre + np.outer(td, direction) + np.outer(qd, normal)
    return curve[:, 0], curve[:, 1]


def case_dir(root, om1, om2):
    return root / CHI_DIR / f"om1_{om1}__om2_{om2}" / CHIBB_DIR


def pw_file(root, replacement_root, om1, om2):
    archive = case_dir(root, om1, om2) / "pw_line.csv"
    if replacement_root is None:
        return archive
    replacement = case_dir(replacement_root, om1, om2) / "pw_line.csv"
    return replacement if replacement.is_file() else archive


def load_series(root, replacement_root=None):
    top = []
    bottom = []
    for encoded, value in zip(OM_VALUES, OM_FLOATS):
        px, py = read_xy(pw_file(root, replacement_root, encoded, "m0p3"),
                         "phi1_inf", "phi2_inf")
        top.append((value, px, py))
        px, py = read_xy(pw_file(root, replacement_root, "m0p3", encoded),
                         "phi1_inf", "phi2_inf")
        bottom.append((value, px, py))
    bx, by = read_xy(case_dir(root, "m0p3", "m0p3") / "binodal.csv",
                     "phi1", "phi2")
    bx, by = physical_binodal_branch(bx, by)
    return bx, by, top, bottom


def common_limits(series):
    xs = [x for group in series for _, px, _ in group for x in px]
    ys = [y for group in series for _, _, py in group for y in py]
    dx = max(xs) - min(xs)
    dy = max(ys) - min(ys)
    return (min(xs) - 0.10 * dx, max(xs) + 0.10 * dx), \
           (min(ys) - 0.06 * dy, max(ys) + 0.06 * dy)


def draw_binodal(ax, bx, by):
    x, y = smooth_monotone_y(bx, by, n=500)
    ax.plot(x, y, color=BINODAL_COLOR, lw=1.3,
            label="binodal", zorder=1)


def draw_overlay(ax, bx, by, series, varied):
    draw_binodal(ax, bx, by)
    for (value, px, py), color, marker in zip(series, COLORS, MARKERS):
        x, y = smooth_curve(px, py)
        ax.plot(x, y, color=color, lw=1.8,
                label=rf"$\omega_{varied}={value:.2f}$", zorder=3)
        points = np.column_stack([px, py])
        centre = points.mean(axis=0)
        _, _, vh = np.linalg.svd(points - centre, full_matrices=False)
        ordered = points[np.argsort((points - centre) @ vh[0])]
        count = min(8, len(ordered))
        idx = np.unique(np.linspace(0, len(ordered) - 1, count).round().astype(int))
        ax.scatter(ordered[idx, 0], ordered[idx, 1], s=13, marker=marker,
                   facecolor=color, edgecolor="white", linewidth=0.25, zorder=4)
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
    ap.add_argument("--replacement-root",
                    help="optional case tree whose pw_line.csv files override the archive")
    ap.add_argument("--out", default="out/analysis/omega/pw_overlay_Ta.png")
    args = ap.parse_args()
    configure()
    replacement_root = Path(args.replacement_root) if args.replacement_root else None
    bx, by, top, bottom = load_series(Path(args.data_root), replacement_root)
    output = Path(args.out)
    main_paths = main_figure(bx, by, top, bottom, output)
    supp = output.with_name(output.stem + "_supplement_2x4" + output.suffix)
    supp_paths = supplement_figure(bx, by, top, bottom, supp)
    plt.close("all")
    print("main:", *main_paths)
    print("supplement:", *supp_paths)


if __name__ == "__main__":
    main()
