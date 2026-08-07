"""Mark the T-f wetting transition point on the binodal, with the pre-wetting line.

Draws what scripts/tf_wetting_point.py computed (wetting_points.csv) on the same
frame as scripts/plot_tf_critical.py uses: the archived binodal points, the tie
lines, the gapless tie-line binodal loop, and the archived pre-wetting points for
this omega setting. This script draws data, it computes nothing.

The wetting transition point is where the wall stops holding a thin layer against
the coexisting phase and lets it spread; it is the end of the pre-wetting line
that sits on the binodal, so the square should land where the red points meet the
loop.

Usage:
  python scripts/plot_tf_wetting.py [--in-dir DIR] [--tie-dir DIR]
      [--field-dir DIR] [--out-dir DIR] [--om1 0.28] [--om2 -0.375]
"""

import argparse
import csv
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BINODAL_COLOR = "0.55"      # archived binodal points
LOOP_COLOR = "#1f77b4"      # recomputed loop
TIE_COLOR = "0.75"
PW_COLOR = "#d62728"
WET_COLOR = "#9467bd"
DIAG_COLOR = "0.80"

PAD_FRAC = 0.06
AXIS_EPS = 1e-6
TIE_EVERY = 6


def read_rows(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def load_archive(field_dir):
    binodal = [(float(r["phi1"]), float(r["phi2"]))
               for r in read_rows(field_dir / "binodal.csv")]
    # binodal.csv also carries a dense degenerate row at phi2 < 1e-6 spanning
    # the whole phi1 range; it is not part of the loop.
    loop = np.array([q for q in binodal if q[1] > AXIS_EPS])
    pw = np.array([(float(r["phi1_inf"]), float(r["phi2_inf"]))
                   for r in read_rows(field_dir / "pw_line.csv")])
    return loop, pw


def load_tie_lines(path):
    branches = {}
    for r in read_rows(path):
        branches.setdefault(int(r["branch"]), []).append(
            [float(r["phi1_a"]), float(r["phi2_a"]),
             float(r["phi1_b"]), float(r["phi2_b"]), float(r["length"])])
    return {k: np.array(v) for k, v in branches.items()}


def closed_loop(branches):
    b0, b1 = branches[0], branches[1]
    parts = [b0[:, 0:2], b0[::-1, 2:4], b1[:, 2:4], b1[::-1, 0:2]]
    return np.vstack(parts)


def square_window(*point_sets):
    vals = np.concatenate([np.asarray(s).ravel() for s in point_sets if len(s)])
    lo, hi = float(vals.min()), float(vals.max())
    pad = (hi - lo) * PAD_FRAC
    return lo - pad, hi + pad


def style_axes(ax, xlim, ylim):
    lo = min(xlim[0], ylim[0])
    hi = max(xlim[1], ylim[1])
    ax.plot([lo, hi], [lo, hi], "--", lw=0.7, color=DIAG_COLOR, zorder=0)
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel(r"$\phi_{1,\infty}$")
    ax.set_ylabel(r"$\phi_{2,\infty}$")


def classify(scan_rows):
    """Split the loop points into the ones where the wall is coated and the rest.

    A far field counts as completely wet when either its own uncoated state has
    ceased to exist, or that state sits above the coated one (Delta > 0). Where
    Delta is not available because the coexisting phase is itself coated, the
    wall at this far field is the bare one.
    """
    coated, bare = [], []
    for r in scan_rows:
        q = (float(r["phi1"]), float(r["phi2"]))
        delta = float(r["delta"])
        wet = (int(r["thin_ok"]) == 0 and r["reason"] == "coated") or \
              (np.isfinite(delta) and delta > 0.0)
        (coated if wet else bare).append(q)
    return np.array(coated), np.array(bare)


def draw(ax, archive_loop, loop, branches, pw, wet, coated_pts, bare_pts,
         tie_every):
    ax.plot(archive_loop[:, 0], archive_loop[:, 1], ".", ms=2.0,
            color=BINODAL_COLOR, zorder=1, label="binodal (archived points)")
    for b in branches.values():
        for row in b[::tie_every]:
            ax.plot([row[0], row[2]], [row[1], row[3]], "-", lw=0.5,
                    color=TIE_COLOR, zorder=2)
    ax.plot(loop[:, 0], loop[:, 1], "-", lw=1.0, color=LOOP_COLOR, zorder=3,
            label="binodal (tie-line continuation)")
    if len(bare_pts):
        ax.plot(bare_pts[:, 0], bare_pts[:, 1], "o", ms=3.0, color=LOOP_COLOR,
                lw=0, zorder=5, label="wall stays bare")
    if len(coated_pts):
        ax.plot(coated_pts[:, 0], coated_pts[:, 1], "o", ms=3.0,
                color=WET_COLOR, lw=0, zorder=5,
                label="wall completely wet")
    if len(pw):
        ax.plot(pw[:, 0], pw[:, 1], "o", ms=2.6, color=PW_COLOR, lw=0, zorder=4,
                label="pre-wetting points")
    if len(wet):
        ax.plot(wet[:, 0], wet[:, 1], "s", ms=11, color="#111111",
                markeredgecolor="k", markeredgewidth=0.6, lw=0, zorder=6,
                label="wetting transition point")


def annotate_wet(ax, wet):
    for phi1, phi2 in wet:
        ax.annotate(f"({phi1:.4f}, {phi2:.4f})", (phi1, phi2),
                    textcoords="offset points", xytext=(9, 7), fontsize=8,
                    color=WET_COLOR)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in-dir", default="out/analysis/Tf_omega_cross/wetting")
    ap.add_argument("--tie-dir", default="out/analysis/Tf_omega_cross/critical")
    ap.add_argument("--field-dir", default="out/analysis/Tf_omega_cross/field")
    ap.add_argument("--out-dir")
    ap.add_argument("--om1", type=float, default=0.28)
    ap.add_argument("--om2", type=float, default=-0.375)
    args = ap.parse_args()

    in_dir = Path(args.in_dir)
    tie_dir = Path(args.tie_dir)
    field_dir = Path(args.field_dir)
    out_dir = Path(args.out_dir) if args.out_dir else in_dir
    for d in (in_dir, tie_dir, field_dir, out_dir):
        if not d.is_dir():
            raise SystemExit(f"directory does not exist: {d}")

    scan_rows = read_rows(in_dir / "wetting_scan.csv")
    coated_pts, bare_pts = classify(scan_rows)
    wet_rows = read_rows(in_dir / "wetting_points.csv")
    wet = np.array([(float(r["phi1_w"]), float(r["phi2_w"])) for r in wet_rows]) \
        if wet_rows else np.empty((0, 2))
    branches = load_tie_lines(tie_dir / "binodal_tie.csv")
    loop = closed_loop(branches)
    archive_loop, pw = load_archive(field_dir)

    print(f"loop points: completely wet {len(coated_pts)}, bare {len(bare_pts)}",
          flush=True)
    print(f"wetting transition points: {len(wet)}", flush=True)
    for r in wet_rows:
        print(f"  ({float(r['phi1_w']):.6f}, {float(r['phi2_w']):.6f})  "
              f"nearest pre-wetting point at {float(r['nearest_pw_dist']):.6f}, "
              f"tangent angle {float(r['tangent_angle_deg']):.2f} deg", flush=True)

    title = (rf"T-f: $\chi_{{12}} = -8.5$, $\chi_{{13}} = \chi_{{23}} = 0$, "
             rf"$\omega_1 = {args.om1:g}$, $\omega_2 = {args.om2:g}$")

    lim = square_window(loop, pw)
    fig, ax = plt.subplots(figsize=(6.4, 6.4))
    draw(ax, archive_loop, loop, branches, pw, wet, coated_pts, bare_pts,
         TIE_EVERY)
    annotate_wet(ax, wet)
    style_axes(ax, lim, lim)
    ax.set_title(title, fontsize=10)
    ax.legend(fontsize=7.5, loc="upper right", framealpha=0.9)
    fig.tight_layout()
    out_png = out_dir / "tf_wetting_overlay.png"
    fig.savefig(out_png, dpi=180)
    plt.close(fig)
    print(f"wrote {out_png}", flush=True)


if __name__ == "__main__":
    main()
