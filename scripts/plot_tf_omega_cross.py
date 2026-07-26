"""T-f (chi12=-8.5, chi13=0, chi23=0) prewetting lines along the omega cross.

The archive samples omega for this topology as a cross through the centre point
(omega1, omega2) = (0.25, -0.375), all at chibb11 = chibb22 = chibb12 = 0 --
the only surface-coupling setting with broad omega coverage. This draws a 2x5
grid: the top row walks omega1 at fixed omega2, the bottom row walks omega2 at
fixed omega1, and column 3 of both rows is the shared centre case.

All prewetting points of a case are drawn as one set: the source column of
pw_line.csv records only which scan direction found a point, which is a
numerical detail, not a physical distinction.

Usage:
  python scripts/plot_tf_omega_cross.py [--data-root DIR] [--out-dir DIR]
"""

import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CHI_DIR = "chi12_m8p5__chi13_0__chi23_0"
CHIBB_DIR = "chibb11_0__chibb22_0__chibb12_0"

CENTRE_OM1 = 0.25
CENTRE_OM2 = -0.375

ROW1_OM1 = [0.22, 0.24, 0.25, 0.26, 0.28]
ROW2_OM2 = [-0.405, -0.385, -0.375, -0.365, -0.345]

BINODAL_COLOR = "0.55"
PW_COLOR = "#d62728"
PAD_FRAC = 0.05


def encode(value):
    """Encode a parameter value the way the archive directory names do."""
    text = f"{value:g}"
    return text.replace("-", "m").replace(".", "p")


def om_dir(om1, om2):
    return f"om1_{encode(om1)}__om2_{encode(om2)}"


def read_cols(path, cols):
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f))
    return [[float(r[c]) for r in rows] for c in cols]


def load_case(root, om1, om2):
    case_dir = root / CHI_DIR / om_dir(om1, om2) / CHIBB_DIR
    bx, by = read_cols(case_dir / "binodal.csv", ["phi1", "phi2"])
    px, py = read_cols(case_dir / "pw_line.csv", ["phi1_inf", "phi2_inf"])
    return {"om1": om1, "om2": om2, "bx": bx, "by": by, "px": px, "py": py}


def draw(ax, case, label):
    ax.plot(case["bx"], case["by"], ".", ms=1.0, color=BINODAL_COLOR, zorder=1)
    if case["px"]:
        ax.plot(case["px"], case["py"], "o", ms=2.8, color=PW_COLOR, lw=0,
                zorder=3)
    ax.set_title(label, fontsize=10)


def shared_limits(cases):
    xs, ys = [], []
    for c in cases:
        xs += c["bx"] + c["px"]
        ys += c["by"] + c["py"]
    x0, x1, y0, y1 = min(xs), max(xs), min(ys), max(ys)
    px, py = (x1 - x0) * PAD_FRAC, (y1 - y0) * PAD_FRAC
    return (x0 - px, x1 + px), (y0 - py, y1 + py)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-root", default="/root/autodl-fs/pw-space/data")
    ap.add_argument("--out-dir", default="out/analysis/Tf_omega_cross")
    args = ap.parse_args()

    root = Path(args.data_root)
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    row1 = [load_case(root, om1, CENTRE_OM2) for om1 in ROW1_OM1]
    row2 = [load_case(root, CENTRE_OM1, om2) for om2 in ROW2_OM2]
    row1_labels = [rf"$\omega_1 = {om1:g}$" for om1 in ROW1_OM1]
    row2_labels = [rf"$\omega_2 = {om2:g}$" for om2 in ROW2_OM2]
    xlim, ylim = shared_limits(row1 + row2)

    fig, axes = plt.subplots(2, 5, figsize=(16.0, 7.0),
                             sharex=True, sharey=True)
    for ax_row, cases, labels in zip(axes, (row1, row2),
                                     (row1_labels, row2_labels)):
        for ax, case, label in zip(ax_row, cases, labels):
            draw(ax, case, label)
            ax.set_xlim(*xlim)
            ax.set_ylim(*ylim)
    for ax in axes[1]:
        ax.set_xlabel(r"$\phi_{1,\infty}$")
    axes[0][0].set_ylabel(rf"$\omega_2 = {CENTRE_OM2:g}$ fixed"
                          "\n" r"$\phi_{2,\infty}$")
    axes[1][0].set_ylabel(rf"$\omega_1 = {CENTRE_OM1:g}$ fixed"
                          "\n" r"$\phi_{2,\infty}$")
    fig.tight_layout()
    fig.savefig(out / "omega_cross_2x5.png", dpi=160)
    plt.close(fig)

    singles = ([(f"row1_om1_{encode(c['om1'])}", c,
                 rf"$\omega_1 = {c['om1']:g}$, $\omega_2 = {c['om2']:g}$")
                for c in row1]
               + [(f"row2_om2_{encode(c['om2'])}", c,
                   rf"$\omega_1 = {c['om1']:g}$, $\omega_2 = {c['om2']:g}$")
                  for c in row2])
    for tag, case, label in singles:
        fig, ax = plt.subplots(figsize=(5.4, 4.8))
        draw(ax, case, label)
        ax.set_xlim(*xlim)
        ax.set_ylim(*ylim)
        ax.set_xlabel(r"$\phi_{1,\infty}$")
        ax.set_ylabel(r"$\phi_{2,\infty}$")
        fig.tight_layout()
        fig.savefig(out / f"{tag}.png", dpi=170)
        plt.close(fig)
        print(f"{tag}: n_pw = {len(case['px'])}", flush=True)

    print(f"grid written to {out/'omega_cross_2x5.png'}", flush=True)


if __name__ == "__main__":
    main()
