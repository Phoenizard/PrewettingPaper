"""Diagnostic phase-plane overlays for the T-a diagonal re-entry scan."""

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
CASES = [
    ("m0p5", -0.50), ("m0p49", -0.49), ("m0p48", -0.48),
    ("m0p47", -0.47), ("m0p46", -0.46), ("m0p45", -0.45),
    ("m0p44", -0.44), ("m0p3", -0.30), ("m0p1", -0.10),
    ("m0p08", -0.08), ("m0p06", -0.06), ("m0p04", -0.04),
    ("m0p02", -0.02), ("m0p01", -0.01), ("m0p001", -0.001),
]


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
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.3), sharex=True, sharey=True)
    groups = (CASES[:7], CASES[8:])
    for ax, group, label in zip(axes, groups, ("strong side", "weak side")):
        ax.scatter(bx, by, s=1.1, color="0.65", linewidths=0, label="binodal")
        colors = [viridis(v) for v in [i / max(1, len(group) - 1) for i in range(len(group))]]
        for (code, value), color in zip(group, colors):
            px, py = read(case_dir(root, code) / "pw_line.csv", "phi1_inf", "phi2_inf")
            if px:
                ax.scatter(px, py, s=12, color=color, linewidths=0,
                           label=rf"$\lambda={value:g}$", zorder=3)
            else:
                ax.plot([], [], color=color, marker="x", ls="",
                        label=rf"$\lambda={value:g}$ (none)")
        ax.set_title(label)
        ax.set_xlabel(r"$\phi_{1,\infty}$")
        ax.legend(frameon=False, markerscale=1.25, handletextpad=0.35,
                  labelspacing=0.25, borderaxespad=0.2, loc="best")
        ax.set_xlim(0.065, 0.205)
        ax.set_ylim(-0.004, 0.205)
        ax.grid(alpha=0.15, linewidth=0.5)
    axes[0].set_ylabel(r"$\phi_{2,\infty}$")
    panel_label(axes[0], "(a)")
    panel_label(axes[1], "(b)")
    fig.subplots_adjust(left=0.09, right=0.99, bottom=0.16, top=0.89, wspace=0.10)
    print(*save_pdf_png(fig, Path(args.out)))


if __name__ == "__main__":
    main()
