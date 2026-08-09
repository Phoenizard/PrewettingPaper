"""T-f (chi12=-8.5, chi13=0, chi23=0) prewetting lines along the chibb cross.

Companion to plot_tf_omega_cross.py. There the wall affinities omega1, omega2
were swept at zero surface-enhanced interaction; here omega is pinned at the
centre of that cross, (omega1, omega2) = (0.25, -0.375), and the surface-enhanced
interactions are swept instead. The layout mirrors the omega figure exactly: the
top row walks the parameter belonging to solute 1, the bottom row walks the one
belonging to solute 2, and column 3 of both rows is the shared all-zero centre.

Controlled variables. Every panel holds chi = (-8.5, 0, 0), omega = (0.25,
-0.375) and chibb12 = 0; the top row additionally holds chibb22 = 0 and the
bottom row chibb11 = 0. chibb12 = 0 is not a choice: every archived case has at
most one non-zero chibb component, so no data exists in which chibb12 sits at a
non-zero value while chibb11 or chibb22 is swept.

All prewetting points of a case are drawn as one set: the source column of
pw_line.csv records only which scan direction found a point, which is a
numerical detail, not a physical distinction.

Both axes share one square window and equal aspect, so the phi1 = phi2 diagonal
(drawn faint in every panel) sits at 45 degrees: this topology is symmetric
under swapping the two solutes together with omega1 and omega2, and the figure
has to show that symmetry undistorted.

Usage:
  python scripts/plot_tf_chibb_cross.py [--data-root DIR] [--out-dir DIR]
"""

import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CHI_DIR = "chi12_m8p5__chi13_0__chi23_0"
OM_DIR = "om1_0p25__om2_m0p375"

CENTRE_OM1 = 0.25
CENTRE_OM2 = -0.375

ROW1_CHIBB11 = [-0.1, -0.05, 0.0, 0.05, 0.1]
ROW2_CHIBB22 = [-0.1, -0.05, 0.0, 0.05, 0.1]
FIXED_CHIBB12 = 0.0

BINODAL_COLOR = "0.55"
PW_COLOR = "#d62728"
DIAG_COLOR = "0.80"
PAD_FRAC = 0.05
AXIS_EPS = 1e-6


def encode(value):
    """Encode a parameter value the way the archive directory names do."""
    text = f"{value:g}"
    return text.replace("-", "m").replace(".", "p")


def chibb_dir(chibb11, chibb22, chibb12):
    return (f"chibb11_{encode(chibb11)}"
            f"__chibb22_{encode(chibb22)}"
            f"__chibb12_{encode(chibb12)}")


def read_cols(path, cols):
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f))
    return [[float(r[c]) for r in rows] for c in cols]


def load_case(root, chibb11, chibb22, chibb12=FIXED_CHIBB12):
    case_dir = root / CHI_DIR / OM_DIR / chibb_dir(chibb11, chibb22, chibb12)
    bx, by = read_cols(case_dir / "binodal.csv", ["phi1", "phi2"])
    px, py = read_cols(case_dir / "pw_line.csv", ["phi1_inf", "phi2_inf"])
    return {"chibb11": chibb11, "chibb22": chibb22, "chibb12": chibb12,
            "bx": bx, "by": by, "px": px, "py": py}


def draw(ax, case, label):
    ax.plot(case["bx"], case["by"], ".", ms=1.0, color=BINODAL_COLOR, zorder=1)
    if case["px"]:
        ax.plot(case["px"], case["py"], "o", ms=2.8, color=PW_COLOR, lw=0,
                zorder=3)
    ax.set_title(label, fontsize=10)


def square_window(cases):
    """One window used for both axes, so the phi1 = phi2 diagonal is at 45 deg.

    Framed on the closed binodal loop and the prewetting points. binodal.csv
    also carries a dense row of points at phi2 < 1e-6 running the whole phi1
    range, with no mirror row on the phi1 = 0 side: on that edge the mixture is
    binary with chi13 = 0, which has no coexistence, so those are degenerate
    solutions of the binodal solve. Framing on them would shrink the loop to a
    corner. They are still drawn -- they just fall outside the window.
    """
    vals = []
    for c in cases:
        loop = [(x, y) for x, y in zip(c["bx"], c["by"]) if y > AXIS_EPS]
        vals += [v for xy in loop for v in xy] + c["px"] + c["py"]
    lo, hi = min(vals), max(vals)
    pad = (hi - lo) * PAD_FRAC
    return lo - pad, hi + pad


def style_axes(ax, lim):
    ax.plot(lim, lim, "--", lw=0.7, color=DIAG_COLOR, zorder=0)
    ax.set_xlim(*lim)
    ax.set_ylim(*lim)
    ax.set_aspect("equal", adjustable="box")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-root", default="/root/autodl-fs/pw-space/data")
    ap.add_argument("--out-dir", default="out/analysis/Tf_omega_cross")
    args = ap.parse_args()

    root = Path(args.data_root)
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    row1 = [load_case(root, b11, 0.0) for b11 in ROW1_CHIBB11]
    row2 = [load_case(root, 0.0, b22) for b22 in ROW2_CHIBB22]
    row1_labels = [rf"$\chi_{{bb,11}} = {b11:g}$" for b11 in ROW1_CHIBB11]
    row2_labels = [rf"$\chi_{{bb,22}} = {b22:g}$" for b22 in ROW2_CHIBB22]
    lim = square_window(row1 + row2)

    fig, axes = plt.subplots(2, 5, figsize=(15.0, 7.0),
                             sharex=True, sharey=True)
    for ax_row, cases, labels in zip(axes, (row1, row2),
                                     (row1_labels, row2_labels)):
        for ax, case, label in zip(ax_row, cases, labels):
            draw(ax, case, label)
            style_axes(ax, lim)
    for ax in axes[1]:
        ax.set_xlabel(r"$\phi_{1,\infty}$")
    axes[0][0].set_ylabel(r"$\chi_{bb,22} = 0$ fixed"
                          "\n" r"$\phi_{2,\infty}$")
    axes[1][0].set_ylabel(r"$\chi_{bb,11} = 0$ fixed"
                          "\n" r"$\phi_{2,\infty}$")
    fig.suptitle(rf"$\chi = (-8.5, 0, 0)$,  "
                 rf"$\omega_1 = {CENTRE_OM1:g}$,  "
                 rf"$\omega_2 = {CENTRE_OM2:g}$,  "
                 rf"$\chi_{{bb,12}} = {FIXED_CHIBB12:g}$ throughout",
                 fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(out / "chibb_cross_2x5.png", dpi=160)
    plt.close(fig)

    singles = ([(f"row1_chibb11_{encode(c['chibb11'])}", c,
                 rf"$\chi_{{bb,11}} = {c['chibb11']:g}$, "
                 rf"$\chi_{{bb,22}} = {c['chibb22']:g}$, "
                 rf"$\chi_{{bb,12}} = {c['chibb12']:g}$")
                for c in row1]
               + [(f"row2_chibb22_{encode(c['chibb22'])}", c,
                   rf"$\chi_{{bb,11}} = {c['chibb11']:g}$, "
                   rf"$\chi_{{bb,22}} = {c['chibb22']:g}$, "
                   rf"$\chi_{{bb,12}} = {c['chibb12']:g}$")
                  for c in row2])
    for tag, case, label in singles:
        fig, ax = plt.subplots(figsize=(5.0, 5.0))
        draw(ax, case, label)
        style_axes(ax, lim)
        ax.set_xlabel(r"$\phi_{1,\infty}$")
        ax.set_ylabel(r"$\phi_{2,\infty}$")
        fig.tight_layout()
        fig.savefig(out / f"{tag}.png", dpi=170)
        plt.close(fig)
        print(f"{tag}: n_pw = {len(case['px'])}", flush=True)

    print(f"grid written to {out/'chibb_cross_2x5.png'}", flush=True)


if __name__ == "__main__":
    main()
