"""Survey figures for the T-f topology (chi12=-8.5, chi13=0, chi23=0).

Picks 9 cases from the pw-space archive that span every control variable this
stage scans (omega1, omega2, chibb12, chibb22, chibb11) around the anchor
(omega1, omega2) = (0.25, -0.375), chibb = 0. Writes one PNG per case plus a
3x3 contact sheet. All prewetting points of a case are drawn as one set: the
source column of pw_line.csv records only which scan direction found a point,
which is a numerical detail, not a physical distinction.

Usage:
  python scripts/plot_tf_survey.py [--data-root DIR] [--out-dir DIR]
"""

import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CHI_DIR = "chi12_m8p5__chi13_0__chi23_0"
BASE_CHIBB = "chibb11_0__chibb22_0__chibb12_0"

# (tag, om_dir, chibb_dir, panel label, control-variable group)
CASES = [
    ("01_om1_0p22", "om1_0p22__om2_m0p375", BASE_CHIBB,
     r"$\omega_1 = 0.22$, $\omega_2 = -0.375$", "omega1 scan"),
    ("02_om1_0p25_anchor", "om1_0p25__om2_m0p375", BASE_CHIBB,
     r"$\omega_1 = 0.25$, $\omega_2 = -0.375$ (anchor)", "omega1 scan"),
    ("03_om1_0p28", "om1_0p28__om2_m0p375", BASE_CHIBB,
     r"$\omega_1 = 0.28$, $\omega_2 = -0.375$", "omega1 scan"),
    ("04_om2_m0p345", "om1_0p25__om2_m0p345", BASE_CHIBB,
     r"$\omega_1 = 0.25$, $\omega_2 = -0.345$", "omega2 scan"),
    ("05_om2_m0p405", "om1_0p25__om2_m0p405", BASE_CHIBB,
     r"$\omega_1 = 0.25$, $\omega_2 = -0.405$", "omega2 scan"),
    ("06_chibb12_p0p2", "om1_0p25__om2_m0p375",
     "chibb11_0__chibb22_0__chibb12_0p2",
     r"anchor, $\chi_{bb,12} = +0.2$", "chibb12 sweep"),
    ("07_chibb12_m0p1", "om1_0p25__om2_m0p375",
     "chibb11_0__chibb22_0__chibb12_m0p1",
     r"anchor, $\chi_{bb,12} = -0.1$", "chibb12 sweep"),
    ("08_chibb22_p0p1", "om1_0p25__om2_m0p375",
     "chibb11_0__chibb22_0p1__chibb12_0",
     r"anchor, $\chi_{bb,22} = +0.1$", "chibb22 sweep"),
    ("09_chibb11_m0p1", "om1_0p25__om2_m0p375",
     "chibb11_m0p1__chibb22_0__chibb12_0",
     r"anchor, $\chi_{bb,11} = -0.1$", "chibb11 sweep"),
]

PW_COLOR = "#d62728"


def read_rows(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def draw(ax, case_dir, label):
    b = read_rows(case_dir / "binodal.csv")
    ax.plot([float(r["phi1"]) for r in b], [float(r["phi2"]) for r in b],
            ".", ms=1.0, color="0.55", zorder=1)

    pw = read_rows(case_dir / "pw_line.csv")
    xs = [float(r["phi1_inf"]) for r in pw]
    ys = [float(r["phi2_inf"]) for r in pw]
    n = len(xs)
    if xs:
        ax.plot(xs, ys, "o", ms=2.6, color=PW_COLOR, lw=0, zorder=3)
    ax.set_title(f"{label}   (n = {n})", fontsize=9)
    ax.set_xlabel(r"$\phi_{1,\infty}$")
    ax.set_ylabel(r"$\phi_{2,\infty}$")
    return n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-root", default="/root/autodl-fs/pw-space/data")
    ap.add_argument("--out-dir", default="out/analysis/Tf_survey")
    args = ap.parse_args()

    root = Path(args.data_root) / CHI_DIR
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    for tag, om, chibb, label, group in CASES:
        fig, ax = plt.subplots(figsize=(5.4, 4.6))
        n = draw(ax, root / om / chibb, label)
        fig.tight_layout()
        fig.savefig(out / f"{tag}.png", dpi=170)
        plt.close(fig)
        print(f"{tag}: {group}, n={n}", flush=True)

    fig, axes = plt.subplots(3, 3, figsize=(14.0, 11.5))
    for ax, (tag, om, chibb, label, group) in zip(axes.ravel(), CASES):
        draw(ax, root / om / chibb, label)
    fig.tight_layout()
    fig.savefig(out / "00_contact_sheet.png", dpi=150)
    plt.close(fig)
    print("contact sheet written", flush=True)


if __name__ == "__main__":
    main()
