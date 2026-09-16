"""Publication phase-plane overlays for the T-a diagonal re-entry scan."""

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

CHI = "chi12_0__chi13_2p8__chi23_0"
CHIBB = "chibb11_0__chibb22_0__chibb12_0"
STRONG = [("m0p5", -0.50), ("m0p48", -0.48),
          ("m0p46", -0.46), ("m0p44", -0.44)]
WEAK = [("m0p1", -0.10), ("m0p08", -0.08),
        ("m0p06", -0.06), ("m0p04", -0.04)]


def read(path, x, y):
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh))
    return [float(r[x]) for r in rows], [float(r[y]) for r in rows]


def case_dir(root, code):
    return root / CHI / f"om1_{code}__om2_{code}" / CHIBB


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-root", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    root = Path(args.data_root)
    configure()
    bx, by = read(case_dir(root, "m0p3") / "binodal.csv", "phi1", "phi2")
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.25))
    groups = (STRONG, WEAK)
    titles = ("strong-attraction termination", "weak-attraction termination")
    legend_locations = ("upper left", "lower right")
    for ax, group, title, legend_loc in zip(axes, groups, titles, legend_locations):
        ax.scatter(bx, by, s=1.1, color="0.65", linewidths=0, label="binodal")
        colors = [viridis(v) for v in (0.08, 0.35, 0.62, 0.89)]
        for (code, value), color in zip(group, colors):
            px, py = read(case_dir(root, code) / "pw_line.csv", "phi1_inf", "phi2_inf")
            ax.scatter(px, py, s=14, color=color, edgecolor="white", linewidth=0.20,
                       label=rf"$\lambda={value:g}$", zorder=3)
        ax.set_title(title)
        ax.set_xlabel(r"$\phi_{1,\infty}$")
        ax.set_ylabel(r"$\phi_{2,\infty}$")
        ax.legend(frameon=False, markerscale=1.2, handletextpad=0.35,
                  labelspacing=0.25, borderaxespad=0.25, loc=legend_loc)
        ax.grid(alpha=0.15, linewidth=0.5)
    axes[0].set_xlim(0.072, 0.090)
    axes[0].set_ylim(-0.001, 0.035)
    axes[1].set_xlim(0.090, 0.190)
    axes[1].set_ylim(0.075, 0.203)
    axes[1].text(0.04, 0.08, r"no pre-wetting for $\lambda\geq-0.02$",
                 transform=axes[1].transAxes, ha="left", va="bottom", fontsize=8.5)
    for ax in axes:
        ax.xaxis.set_major_locator(plt.MaxNLocator(5))
        ax.yaxis.set_major_locator(plt.MaxNLocator(5))
    panel_label(axes[0], "(a)")
    panel_label(axes[1], "(b)")
    fig.subplots_adjust(left=0.09, right=0.99, bottom=0.16, top=0.87, wspace=0.28)
    print(*save_pdf_png(fig, Path(args.out)))


if __name__ == "__main__":
    main()
