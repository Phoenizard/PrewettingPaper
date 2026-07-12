"""Topic-1 (omega) figures: pw_length and dist_mean over the (omega1, omega2) grid.

Reads the master measures table produced by extract_measures.py, filters to one
chi topology with chi_bb == 0 (the omega topic), and draws two heatmaps on the
(omega_1, omega_2) grid: the headline pre-wetting-line length pw_length and the
headline distance-to-binodal dist_mean. Cases with no pre-wetting are drawn as
empty (white) cells so the extinction region is visible.

These are analysis figures, not per-case maps: they go under out/ (never
doc/analysis/). All labels are English.

Usage:
  python scripts/plot_omega_maps.py --measures <measures.csv> \
      --stage chi12_0__chi13_2p8__chi23_0 --out-dir <out/analysis/omega>
"""

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def _load(measures_path, stage):
    rows = []
    with open(measures_path, newline="") as fh:
        for r in csv.DictReader(fh):
            if r["stage"] != stage:
                continue
            # omega topic: surface cross/self terms all zero
            if any(float(r[k]) != 0.0 for k in
                   ("chi_bb_11", "chi_bb_22", "chi_bb_12")):
                continue
            rows.append(r)
    return rows


def _grid(rows, field):
    """Build a (omega2 x omega1) grid of `field`; NaN where no pre-wetting."""
    om1 = sorted({float(r["omega_1"]) for r in rows})
    om2 = sorted({float(r["omega_2"]) for r in rows})
    i1 = {v: i for i, v in enumerate(om1)}
    i2 = {v: i for i, v in enumerate(om2)}
    grid = np.full((len(om2), len(om1)), np.nan)
    for r in rows:
        val = r.get(field, "")
        if val == "" or r["flag"] == "no_prewetting":
            continue
        grid[i2[float(r["omega_2"])], i1[float(r["omega_1"])]] = float(val)
    return np.array(om1), np.array(om2), grid


def _heatmap(ax, om1, om2, grid, title, cbar_label):
    ax.set_facecolor("white")
    im = ax.pcolormesh(om1, om2, grid, shading="nearest", cmap="viridis")
    ax.set_xlabel(r"$\omega_1$")
    ax.set_ylabel(r"$\omega_2$")
    ax.set_title(title, fontsize=10)
    ax.set_aspect("equal")
    cb = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label(cbar_label)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--measures", required=True)
    ap.add_argument("--stage", default="chi12_0__chi13_2p8__chi23_0")
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()

    rows = _load(args.measures, args.stage)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))
    om1, om2, gl = _grid(rows, "pw_length")
    _heatmap(axes[0], om1, om2, gl,
             r"pre-wetting line length $L$", r"$L$")
    om1, om2, gd = _grid(rows, "dist_mean")
    _heatmap(axes[1], om1, om2, gd,
             r"distance to binodal $\bar d$", r"$\bar d$")
    fig.suptitle(f"omega topic ({args.stage}, chibb=0)", fontsize=11)
    fig.tight_layout()
    out = out_dir / "omega_length_dist.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)

    n_pw = sum(1 for r in rows if r["flag"] != "no_prewetting")
    print(f"omega topic cases: {len(rows)}  with pre-wetting: {n_pw}")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
