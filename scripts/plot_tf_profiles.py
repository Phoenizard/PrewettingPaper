"""Plot the thin / thick surface profiles and the scan lines behind them.

Reads what scripts/tf_profiles.py wrote: profiles.csv (long table of phi_1(z),
phi_2(z) per point and branch), profile_summary.csv, scanline_<id>.csv.

Two figures:
  tf_profiles.png   one panel per far field; phi_1 and phi_2 against z, thin
                    state solid, thick state dashed, far-field levels dotted.
                    x is cut where the profile has relaxed back to the far field;
                    tf_profiles_fullz.png keeps the whole box as a check.
  tf_scanlines.png  the hysteresis line through each far field: gamma of both
                    branches (top row) and total adsorption cs (bottom row)
                    against the swept phi_2_inf, with the pre-wetting point marked.

Usage:
  python scripts/plot_tf_profiles.py --in-dir DIR [--out-dir DIR] [--title STR]
"""

import argparse
import csv
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

C1 = "#1f77b4"   # solute 1
C2 = "#d62728"   # solute 2
RELAX_TOL = 1e-3
PAD_FRAC = 0.25


def read_rows(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def load_profiles(path):
    grouped = defaultdict(lambda: {"z": [], "phi1": [], "phi2": []})
    for r in read_rows(path):
        d = grouped[(r["point_id"], r["branch"])]
        d["z"].append(float(r["z"]))
        d["phi1"].append(float(r["phi1"]))
        d["phi2"].append(float(r["phi2"]))
    return grouped


def relaxation_z(curves, phi1_inf, phi2_inf):
    """Largest z where either profile still deviates from the far field."""
    z_max = 0.0
    for d in curves:
        for k in range(len(d["z"])):
            dev = max(abs(d["phi1"][k] - phi1_inf), abs(d["phi2"][k] - phi2_inf))
            if dev > RELAX_TOL:
                z_max = max(z_max, d["z"][k])
    return z_max


def draw_profiles(ax, prof, summary, point_id, zlim):
    phi1_inf = float(summary[point_id]["thin"]["phi1_inf"])
    phi2_inf = float(summary[point_id]["thin"]["phi2_inf"])
    for branch, style in (("thin", "-"), ("thick", "--")):
        d = prof.get((point_id, branch))
        if d is None:
            continue
        ax.plot(d["z"], d["phi1"], style, color=C1, lw=1.6)
        ax.plot(d["z"], d["phi2"], style, color=C2, lw=1.6)
    ax.axhline(phi1_inf, color=C1, ls=":", lw=0.8)
    ax.axhline(phi2_inf, color=C2, ls=":", lw=0.8)
    ax.set_xlim(0.0, zlim)
    ax.set_xlabel(r"$z$")
    ax.set_title(rf"{point_id}: $\phi_{{1,\infty}} = {phi1_inf:.4f}$, "
                 rf"$\phi_{{2,\infty}} = {phi2_inf:.4f}$", fontsize=10)


def profile_figure(prof, summary, point_ids, out_path, zlim_by_point, title):
    n = len(point_ids)
    fig, axes = plt.subplots(1, n, figsize=(4.6 * n, 4.2), sharey=True)
    axes = [axes] if n == 1 else list(axes)
    for ax, pid in zip(axes, point_ids):
        draw_profiles(ax, prof, summary, pid, zlim_by_point[pid])
    axes[0].set_ylabel(r"$\phi_i(z)$")
    handles = [
        plt.Line2D([], [], color=C1, ls="-", label=r"$\phi_1$ thin"),
        plt.Line2D([], [], color=C2, ls="-", label=r"$\phi_2$ thin"),
        plt.Line2D([], [], color=C1, ls="--", label=r"$\phi_1$ thick"),
        plt.Line2D([], [], color=C2, ls="--", label=r"$\phi_2$ thick"),
    ]
    axes[-1].legend(handles=handles, fontsize=8, loc="upper right")
    if title:
        fig.suptitle(title, fontsize=11)
        fig.tight_layout(rect=(0, 0, 1, 0.95))
    else:
        fig.tight_layout()
    fig.savefig(out_path, dpi=170)
    plt.close(fig)
    print(f"wrote {out_path}", flush=True)


def scanline_figure(in_dir, summary, point_ids, out_path, title):
    n = len(point_ids)
    fig, axes = plt.subplots(2, n, figsize=(4.6 * n, 6.4), sharex="col")
    axes = axes.reshape(2, n)
    for col, pid in enumerate(point_ids):
        rows = read_rows(in_dir / f"scanline_{pid}.csv")
        x = [float(r["phi2_inf"]) for r in rows]
        pw = float(summary[pid]["thin"]["pw_phi2_archive"])
        for row, (key_thin, key_thick, ylabel) in enumerate((
            ("gamma_thin", "gamma_thick", r"$\gamma$"),
            ("cs_thin", "cs_thick", r"$c_s = c_{s,1} + c_{s,2}$"),
        )):
            ax = axes[row][col]
            ax.plot(x, [float(r[key_thin]) for r in rows], "-", color=C1,
                    lw=1.4, label="thin branch")
            ax.plot(x, [float(r[key_thick]) for r in rows], "--", color=C2,
                    lw=1.4, label="thick branch")
            ax.axvline(pw, color="0.4", ls=":", lw=1.0)
            if col == 0:
                ax.set_ylabel(ylabel)
            if row == 0:
                ax.set_title(rf"{pid}: $\phi_{{1,\infty}} = "
                             rf"{float(summary[pid]['thin']['phi1_inf']):.4f}$",
                             fontsize=10)
            else:
                ax.set_xlabel(r"$\phi_{2,\infty}$")
    axes[0][-1].legend(fontsize=8)
    if title:
        fig.suptitle(title, fontsize=11)
        fig.tight_layout(rect=(0, 0, 1, 0.95))
    else:
        fig.tight_layout()
    fig.savefig(out_path, dpi=170)
    plt.close(fig)
    print(f"wrote {out_path}", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in-dir", required=True)
    ap.add_argument("--out-dir")
    ap.add_argument("--title", default="")
    args = ap.parse_args()

    in_dir = Path(args.in_dir)
    out_dir = Path(args.out_dir) if args.out_dir else in_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    prof = load_profiles(in_dir / "profiles.csv")
    summary = defaultdict(dict)
    for r in read_rows(in_dir / "profile_summary.csv"):
        summary[r["point_id"]][r["branch"]] = r
    point_ids = sorted(summary)

    zlim_cut, zlim_full = {}, {}
    for pid in point_ids:
        curves = [prof[(pid, b)] for b in ("thin", "thick") if (pid, b) in prof]
        z_full = max(max(d["z"]) for d in curves)
        phi1_inf = float(summary[pid]["thin"]["phi1_inf"])
        phi2_inf = float(summary[pid]["thin"]["phi2_inf"])
        z_rel = relaxation_z(curves, phi1_inf, phi2_inf)
        zlim_cut[pid] = min(z_full, max(z_rel * (1.0 + PAD_FRAC), 0.5))
        zlim_full[pid] = z_full

    profile_figure(prof, summary, point_ids, out_dir / "tf_profiles.png",
                   zlim_cut, args.title)
    profile_figure(prof, summary, point_ids, out_dir / "tf_profiles_fullz.png",
                   zlim_full, args.title)
    scanline_figure(in_dir, summary, point_ids, out_dir / "tf_scanlines.png",
                    args.title)


if __name__ == "__main__":
    main()
