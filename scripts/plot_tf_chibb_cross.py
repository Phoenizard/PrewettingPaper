"""Publication and supplementary T-f surface-interaction figures from archived data."""

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.cm import cividis

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from paperfig import configure, panel_label, save_pdf_png

CHI_DIR = "chi12_m8p5__chi13_0__chi23_0"
OM_DIR = "om1_0p25__om2_m0p375"
VALUES = [-0.1, -0.05, 0.0, 0.05, 0.1]
COLORS = [cividis(t) for t in (0.05, 0.27, 0.50, 0.73, 0.95)]
MARKERS = ["o", "s", "^", "D", "P"]


def encode(value):
    return f"{value:g}".replace("-", "m").replace(".", "p")


def chibb_dir(chibb11, chibb22, chibb12=0.0):
    return (f"chibb11_{encode(chibb11)}__chibb22_{encode(chibb22)}"
            f"__chibb12_{encode(chibb12)}")


def read_cols(path, cols):
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh))
    return [[float(r[c]) for r in rows] for c in cols]


def load_case(root, chibb11, chibb22):
    directory = root / CHI_DIR / OM_DIR / chibb_dir(chibb11, chibb22)
    px, py = read_cols(directory / "pw_line.csv", ["phi1_inf", "phi2_inf"])
    return {"chibb11": chibb11, "chibb22": chibb22, "px": px, "py": py}


def load_gapless_binodal(path):
    branches = defaultdict(lambda: {"a1": [], "a2": [], "b1": [], "b2": []})
    columns = {"a1": "phi1_a", "a2": "phi2_a",
               "b1": "phi1_b", "b2": "phi2_b"}
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            branch = int(row["branch"])
            for key, column in columns.items():
                branches[branch][key].append(float(row[column]))
    return [branches[k] for k in sorted(branches)]


def draw_binodal(ax, branches):
    first = True
    for branch in branches:
        for xkey, ykey in (("a1", "a2"), ("b1", "b2")):
            ax.plot(branch[xkey], branch[ykey], color="0.42", lw=1.35,
                    label="binodal" if first else None, zorder=1)
            first = False


def limits(branches, cases):
    vals = []
    for branch in branches:
        vals.extend(branch["a1"] + branch["a2"] + branch["b1"] + branch["b2"])
    for case in cases:
        vals.extend(case["px"] + case["py"])
    lo, hi = min(vals), max(vals)
    pad = 0.045 * (hi - lo)
    return lo - pad, hi + pad


def style_axes(ax, lim):
    ax.plot(lim, lim, "--", lw=0.8, color="0.80", zorder=0)
    ax.set_xlim(*lim)
    ax.set_ylim(*lim)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel(r"$\phi_{1,\infty}$")


def main_figure(branches, row1, row2, lim, output):
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.35), sharex=True, sharey=True)
    for ax, cases, key, subscript in ((axes[0], row1, "chibb11", "11"),
                                      (axes[1], row2, "chibb22", "22")):
        draw_binodal(ax, branches)
        for case, color, marker in zip(cases, COLORS, MARKERS):
            ax.scatter(case["px"], case["py"], s=15, marker=marker,
                       color=color, edgecolor="white", linewidth=0.25,
                       label=rf"$\chi_{{bb,{subscript}}}={case[key]:g}$", zorder=3)
        style_axes(ax, lim)
        ax.legend(loc="best", frameon=False, handletextpad=0.3,
                  borderaxespad=0.2, labelspacing=0.25)
    axes[0].set_ylabel(r"$\phi_{2,\infty}$")
    axes[0].set_title(r"varying $\chi_{bb,11}$")
    axes[1].set_title(r"varying $\chi_{bb,22}$")
    panel_label(axes[0], "(a)")
    panel_label(axes[1], "(b)")
    fig.subplots_adjust(left=0.10, right=0.99, bottom=0.15, top=0.88, wspace=0.12)
    return save_pdf_png(fig, output)


def supplementary_figure(branches, row1, row2, lim, output):
    fig, axes = plt.subplots(2, 5, figsize=(12.2, 5.2), sharex=True, sharey=True)
    for row, cases, key, subscript in ((axes[0], row1, "chibb11", "11"),
                                       (axes[1], row2, "chibb22", "22")):
        for ax, case, color, marker in zip(row, cases, COLORS, MARKERS):
            draw_binodal(ax, branches)
            ax.scatter(case["px"], case["py"], s=12, marker=marker,
                       color=color, linewidth=0, zorder=3)
            style_axes(ax, lim)
            ax.set_title(rf"$\chi_{{bb,{subscript}}}={case[key]:g}$")
    for ax in axes[0]:
        ax.set_xlabel("")
    axes[0, 0].set_ylabel(r"$\phi_{2,\infty}$")
    axes[1, 0].set_ylabel(r"$\phi_{2,\infty}$")
    fig.tight_layout(pad=0.6, w_pad=0.35, h_pad=0.6)
    return save_pdf_png(fig, output)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-root", default="/root/autodl-fs/pw-space/data")
    ap.add_argument("--out-dir", default="out/analysis/Tf_omega_cross")
    ap.add_argument("--binodal-tie",
                    default="out/analysis/Tf_omega_cross/critical/binodal_tie.csv")
    args = ap.parse_args()
    configure()
    root = Path(args.data_root)
    out = Path(args.out_dir)
    row1 = [load_case(root, value, 0.0) for value in VALUES]
    row2 = [load_case(root, 0.0, value) for value in VALUES]
    branches = load_gapless_binodal(Path(args.binodal_tie))
    lim = limits(branches, row1 + row2)
    main_paths = main_figure(branches, row1, row2, lim, out / "chibb_cross_main.png")
    supp_paths = supplementary_figure(branches, row1, row2, lim,
                                      out / "chibb_cross_supplement_2x5.png")
    plt.close("all")
    print("main:", *main_paths)
    print("supplement:", *supp_paths)


if __name__ == "__main__":
    main()
