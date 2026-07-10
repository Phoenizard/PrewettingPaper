"""Angle 1 -- wall affinity omega1/omega2 vs the prewetting line.

Sweeps the T-a omega grid (chibb=0), extracts PW-line metrics from the overlay
PNGs, and produces a (omega1, omega2) -> phi2_max phase map with the existence
boundary marked. Answers: does making the wall more attractive (omega more
negative) grow the prewetting window, and where does prewetting vanish.

Run: /opt/miniconda3/envs/numenv/bin/python scripts/angle1_omega.py
"""

import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
from pwpix import metrics, find_overlay  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
CHI_TA = "chi12_0p0__chi13_2p8__chi23_0p0"
CHIBB0 = "chibb11_0p00__chibb22_0p00__chibb12_0p00"
OUT = ROOT / "out" / "analysis"
OUT.mkdir(parents=True, exist_ok=True)


def decode(tok):
    return float(tok.replace("m", "-").replace("p", "."))


def collect():
    base = ROOT / "data" / CHI_TA
    rows = []
    for om_dir in sorted(base.iterdir()):
        if not om_dir.name.startswith("om1"):
            continue
        case = om_dir / CHIBB0
        ov = find_overlay(case)
        if ov is None:
            continue
        o1, o2 = om_dir.name.split("__")
        om1 = decode(o1.replace("om1_", ""))
        om2 = decode(o2.replace("om2_", ""))
        m = metrics(ov)
        rows.append(dict(om1=om1, om2=om2, **m))
    return rows


def main():
    rows = collect()
    om1s = sorted({r["om1"] for r in rows})
    om2s = sorted({r["om2"] for r in rows})
    grid = np.full((len(om2s), len(om1s)), np.nan)
    exist = np.zeros_like(grid, dtype=bool)
    idx1 = {v: i for i, v in enumerate(om1s)}
    idx2 = {v: i for i, v in enumerate(om2s)}
    for r in rows:
        i, j = idx2[r["om2"]], idx1[r["om1"]]
        exist[i, j] = r["exists"]
        if r["exists"]:
            grid[i, j] = r["phi2_max"]

    fig, ax = plt.subplots(figsize=(7.2, 6.0))
    # phi2_max heatmap where PW exists
    im = ax.imshow(
        grid, origin="lower", aspect="auto",
        extent=[min(om1s), max(om1s), min(om2s), max(om2s)],
        cmap="viridis", interpolation="nearest",
    )
    cb = fig.colorbar(im, ax=ax)
    cb.set_label(r"$\phi_{2,\infty}$ extent of prewetting line  (top of PW line)")

    # hatch the no-prewetting region
    no_pw = ~exist
    ax.contourf(
        om1s, om2s, no_pw.astype(float), levels=[0.5, 1.5],
        colors="none", hatches=["xxx"], alpha=0,
    )
    # existence boundary
    ax.contour(om1s, om2s, exist.astype(float), levels=[0.5],
               colors="crimson", linewidths=2.0)

    ax.set_xlabel(r"$\omega_1$  (wall affinity for solute 1)")
    ax.set_ylabel(r"$\omega_2$  (wall affinity for solute 2)")
    ax.set_title("Angle 1: wall affinity controls the prewetting window\n"
                 r"T-a  ($\chi_{12}{=}0,\ \chi_{13}{=}2.8,\ \chi_{23}{=}0$),  "
                 r"$\chi_{bb}{=}0$")
    ax.text(0.02, 0.02,
            "hatched = no prewetting;  red = existence boundary\n"
            "more negative omega = wall more attractive = larger PW window",
            transform=ax.transAxes, fontsize=8, va="bottom",
            bbox=dict(boxstyle="round", fc="white", ec="0.7", alpha=0.85))

    fig.tight_layout()
    out = OUT / "angle1_omega_phasemap.png"
    fig.savefig(out, dpi=150)
    print("wrote", out)

    # also dump the numeric table
    import csv
    tbl = OUT / "angle1_omega_metrics.csv"
    with open(tbl, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["om1", "om2", "exists", "n_branch", "phi2_max",
                    "phi1_min", "phi1_max", "pw_coverage"])
        for r in sorted(rows, key=lambda r: (r["om1"], r["om2"])):
            w.writerow([r["om1"], r["om2"], int(r["exists"]), r["n_branch"],
                        f"{r['phi2_max']:.4f}" if r["exists"] else "",
                        f"{r['phi1_min']:.4f}" if r["exists"] else "",
                        f"{r['phi1_max']:.4f}" if r["exists"] else "",
                        f"{r['pw_coverage']:.5f}"])
    print("wrote", tbl)

    n_exist = sum(r["exists"] for r in rows)
    print(f"cases: {len(rows)}, prewetting present: {n_exist}, absent: {len(rows)-n_exist}")


if __name__ == "__main__":
    main()
