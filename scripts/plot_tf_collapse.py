"""Do the two arms of the T-f omega cross collapse onto omega1 + omega2?

Falsification test for the claim that in the topology chi12 = -8.5, chi13 = chi23 = 0
the wall acts on the prewetting film only through its net affinity for the two
solutes. If that holds, the arm that varies omega1 at fixed omega2 and the arm that
varies omega2 at fixed omega1 must fall on a single curve when plotted against
omega1 + omega2. Points with the same net affinity but a different split between the
solutes test how much the split still matters.

Reads the archive-wide measure table; computes nothing new.

Usage:
  python scripts/plot_tf_collapse.py [--measures out/analysis/measures.csv]
                                     [--out-dir out/analysis/Tf]
"""

import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

STAGE = "chi12_m8p5__chi13_0__chi23_0"
CENTRE_OM1 = 0.25
CENTRE_OM2 = -0.375
TOL = 1e-9

# label -> (colour, marker, filled)
STYLE = {
    "arm A: omega2 fixed, omega1 varied": ("#1f77b4", "o", True),
    "arm B: omega1 fixed, omega2 varied": ("#d62728", "s", True),
    "mirror of arm A": ("#1f77b4", "o", False),
    "mirror of arm B": ("#d62728", "s", False),
    "off-cross points": ("#2ca02c", "D", True),
}


def load(path):
    rows = []
    with open(path, newline="") as f:
        for r in csv.DictReader(f):
            if r["stage"] != STAGE:
                continue
            if any(float(r[k]) != 0.0
                   for k in ("chi_bb_11", "chi_bb_22", "chi_bb_12")):
                continue
            if r["flag"]:
                continue
            rows.append({
                "om1": float(r["omega_1"]),
                "om2": float(r["omega_2"]),
                "L": float(r["pw_length"]),
                "d": float(r["dist_mean"]),
            })
    return rows


def classify(r):
    om1, om2 = r["om1"], r["om2"]
    if abs(om2 - CENTRE_OM2) < TOL and om1 > 0:
        return "arm A: omega2 fixed, omega1 varied"
    if abs(om1 - CENTRE_OM1) < TOL:
        return "arm B: omega1 fixed, omega2 varied"
    if abs(om1 - CENTRE_OM2) < TOL and om2 > 0:
        return "mirror of arm A"
    if abs(om2 - CENTRE_OM1) < TOL:
        return "mirror of arm B"
    return "off-cross points"


def panel(ax, groups, field, ylabel):
    for label, (color, marker, filled) in STYLE.items():
        pts = groups.get(label, [])
        if not pts:
            continue
        xs = [p["om1"] + p["om2"] for p in pts]
        ys = [p[field] for p in pts]
        ax.plot(xs, ys, marker, ms=7, lw=0, label=label,
                color=color if filled else "none",
                markerfacecolor=color if filled else "none",
                markeredgecolor=color, markeredgewidth=1.4)
    ax.set_xlabel(r"$\omega_1 + \omega_2$")
    ax.set_ylabel(ylabel)
    ax.grid(alpha=0.25)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--measures", default="out/analysis/measures.csv")
    ap.add_argument("--out-dir", default="out/analysis/Tf")
    args = ap.parse_args()

    rows = load(args.measures)
    groups = {}
    for r in rows:
        groups.setdefault(classify(r), []).append(r)

    print(f"{'group':<38}{'om1':>8}{'om2':>9}{'sum':>9}"
          f"{'L':>10}{'d_mean':>10}")
    for label in STYLE:
        for p in sorted(groups.get(label, []), key=lambda q: q["om1"] + q["om2"]):
            print(f"{label:<38}{p['om1']:8.3f}{p['om2']:9.3f}"
                  f"{p['om1'] + p['om2']:9.3f}{p['L']:10.5f}{p['d']:10.5f}")

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(12.4, 5.0))
    panel(axes[0], groups, "L", r"prewetting line length $L$")
    panel(axes[1], groups, "d", r"mean distance to binodal $\bar{d}$")
    axes[0].legend(fontsize=8, loc="best")
    fig.tight_layout()
    fig.savefig(out / "tf_collapse.png", dpi=170)
    plt.close(fig)
    print(f"\nwritten {out / 'tf_collapse.png'}", flush=True)


if __name__ == "__main__":
    main()
