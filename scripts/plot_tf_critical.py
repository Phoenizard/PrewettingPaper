"""Mark the T-f binodal critical points on the loop, with the pre-wetting line.

Reads what scripts/tf_critical_point.py wrote (critical_points.csv,
binodal_tie.csv) and draws the recomputed loop against the archived
binodal.csv, so where the archived detection dropped points is visible directly.
The pre-wetting points are the archived ones for a single omega setting -- this
draws data, it does not compute any.

The tie lines are drawn every few steps: they run nearly along the phi1 = phi2
direction and shrink to nothing at the critical points, which is what makes the
critical points readable off the figure rather than asserted.

Also prints the pre-wetting points ranked by distance to the nearer critical
point, which is the input the profile step needs.

Usage:
  python scripts/plot_tf_critical.py [--data-root DIR] [--in-dir DIR]
      [--out-dir DIR] [--om1 0.28] [--om2 -0.375]
"""

import argparse
import csv
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CHI_DIR = "chi12_m8p5__chi13_0__chi23_0"
CHIBB_DIR = "chibb11_0__chibb22_0__chibb12_0"

BINODAL_COLOR = "0.55"      # archived binodal points
LOOP_COLOR = "#1f77b4"      # recomputed loop
TIE_COLOR = "0.75"
PW_COLOR = "#d62728"
CRIT_COLOR = "#2ca02c"
DIAG_COLOR = "0.80"

PAD_FRAC = 0.06
AXIS_EPS = 1e-6
TIE_EVERY = 6
N_RANKED = 10


def encode(value):
    return f"{value:g}".replace("-", "m").replace(".", "p")


def read_rows(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def load_archive(root, om1, om2):
    case = root / CHI_DIR / f"om1_{encode(om1)}__om2_{encode(om2)}" / CHIBB_DIR
    binodal = [(float(r["phi1"]), float(r["phi2"]))
               for r in read_rows(case / "binodal.csv")]
    # binodal.csv also carries a dense degenerate row at phi2 < 1e-6 spanning
    # the whole phi1 range; it is not part of the loop.
    loop = np.array([q for q in binodal if q[1] > AXIS_EPS])
    pw = np.array([(float(r["phi1_inf"]), float(r["phi2_inf"]))
                   for r in read_rows(case / "pw_line.csv")])
    return loop, pw


def load_tie_lines(path):
    rows = read_rows(path)
    branches = {}
    for r in rows:
        branches.setdefault(int(r["branch"]), []).append(
            [float(r["phi1_a"]), float(r["phi2_a"]),
             float(r["phi1_b"]), float(r["phi2_b"]), float(r["length"])])
    return {k: np.array(v) for k, v in branches.items()}


def closed_loop(branches):
    """Stitch the two traced branches into one closed curve.

    Each branch runs from the symmetric tie line out to a critical point, so its
    a-ends trace one arc and its b-ends trace the other. Walking a-ends out,
    b-ends back, then the same on the mirror branch closes the loop.
    """
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


def draw(ax, archive_loop, loop, branches, pw, crit, tie_every):
    ax.plot(archive_loop[:, 0], archive_loop[:, 1], ".", ms=2.0,
            color=BINODAL_COLOR, zorder=1, label="binodal (archived points)")
    for b in branches.values():
        for row in b[::tie_every]:
            ax.plot([row[0], row[2]], [row[1], row[3]], "-", lw=0.5,
                    color=TIE_COLOR, zorder=2)
    ax.plot(loop[:, 0], loop[:, 1], "-", lw=1.4, color=LOOP_COLOR, zorder=3,
            label="binodal (tie-line continuation)")
    if len(pw):
        ax.plot(pw[:, 0], pw[:, 1], "o", ms=2.6, color=PW_COLOR, lw=0, zorder=4,
                label="pre-wetting points")
    ax.plot(crit[:, 0], crit[:, 1], "*", ms=14, color=CRIT_COLOR,
            markeredgecolor="k", markeredgewidth=0.6, lw=0, zorder=5,
            label="critical point")


def annotate_crit(ax, crit):
    for phi1, phi2 in crit:
        ax.annotate(f"({phi1:.4f}, {phi2:.4f})", (phi1, phi2),
                    textcoords="offset points", xytext=(9, 7), fontsize=8,
                    color=CRIT_COLOR)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-root", default="/root/autodl-fs/pw-space/data")
    ap.add_argument("--in-dir", default="out/analysis/Tf_omega_cross/critical")
    ap.add_argument("--out-dir")
    ap.add_argument("--om1", type=float, default=0.28)
    ap.add_argument("--om2", type=float, default=-0.375)
    args = ap.parse_args()

    in_dir = Path(args.in_dir)
    out_dir = Path(args.out_dir) if args.out_dir else in_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    crit = np.array([(float(r["phi1_c"]), float(r["phi2_c"]))
                     for r in read_rows(in_dir / "critical_points.csv")])
    branches = load_tie_lines(in_dir / "binodal_tie.csv")
    loop = closed_loop(branches)
    archive_loop, pw = load_archive(Path(args.data_root), args.om1, args.om2)

    title = (rf"T-f: $\chi_{{12}} = -8.5$, $\chi_{{13}} = \chi_{{23}} = 0$, "
             rf"$\omega_1 = {args.om1:g}$, $\omega_2 = {args.om2:g}$")

    # how far the archived points sit from the recomputed loop
    d_arch = np.min(np.hypot(archive_loop[:, None, 0] - loop[None, :, 0],
                             archive_loop[:, None, 1] - loop[None, :, 1]),
                    axis=1)
    print(f"archived loop points: {len(archive_loop)}, distance to the "
          f"recomputed loop  max {d_arch.max():.2e}  mean {d_arch.mean():.2e}",
          flush=True)

    # which critical point the pre-wetting line is near, and the ranking
    d_pw = np.hypot(pw[:, None, 0] - crit[None, :, 0],
                    pw[:, None, 1] - crit[None, :, 1])
    which = int(np.argmin(d_pw.min(axis=0)))
    order = np.argsort(d_pw[:, which])
    target = crit[which]
    print(f"\npre-wetting points: {len(pw)}; nearer critical point "
          f"({target[0]:.6f}, {target[1]:.6f})", flush=True)
    print(f"{'rank':>4} {'phi1_inf':>10} {'phi2_inf':>10} {'distance':>10}",
          flush=True)
    for rank, i in enumerate(order[:N_RANKED]):
        print(f"{rank:>4} {pw[i, 0]:>10.6f} {pw[i, 1]:>10.6f} "
              f"{d_pw[i, which]:>10.6f}", flush=True)

    # full view
    lim = square_window(loop, pw)
    fig, ax = plt.subplots(figsize=(6.4, 6.4))
    draw(ax, archive_loop, loop, branches, pw, crit, TIE_EVERY)
    annotate_crit(ax, crit)
    style_axes(ax, lim, lim)
    ax.set_title(title, fontsize=10)
    ax.legend(fontsize=7.5, loc="upper right", framealpha=0.9)
    fig.tight_layout()
    fig.savefig(out_dir / "tf_critical_overlay.png", dpi=180)
    plt.close(fig)
    print(f"\nwrote {out_dir / 'tf_critical_overlay.png'}", flush=True)

    # zoom on the critical point the pre-wetting line runs toward
    near_pw = pw[order[:30]]
    box = np.vstack([near_pw, target[None, :]])
    half = max(box[:, 0].ptp(), box[:, 1].ptp()) * (0.5 + PAD_FRAC) + 0.005
    cx, cy = 0.5 * (box[:, 0].min() + box[:, 0].max()), \
             0.5 * (box[:, 1].min() + box[:, 1].max())
    fig, ax = plt.subplots(figsize=(6.4, 6.4))
    draw(ax, archive_loop, loop, branches, pw, crit, max(TIE_EVERY // 3, 1))
    annotate_crit(ax, crit)
    style_axes(ax, (cx - half, cx + half), (cy - half, cy + half))
    ax.set_title(title + "  (zoom)", fontsize=10)
    fig.tight_layout()
    fig.savefig(out_dir / "tf_critical_zoom.png", dpi=180)
    plt.close(fig)
    print(f"wrote {out_dir / 'tf_critical_zoom.png'}", flush=True)


if __name__ == "__main__":
    main()
