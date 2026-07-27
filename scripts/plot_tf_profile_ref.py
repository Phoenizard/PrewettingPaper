"""Plot the two surface states found by scripts/tf_profile_ref.py.

Reads profiles.csv and states.csv from that script's output directory and draws
one panel per far field: phi_1(z) and phi_2(z) for the low-adsorption state
(thin film, solid) and the high-adsorption state (thick film, dashed), with the
far-field levels as dotted lines. The x axis is cut where both profiles have
relaxed back to the far field; tf_profiles_fullz.png keeps the whole box.

Also prints the numbers the T-f note's two claims turn on: the wall contact
values, the width over which each state deviates from the far field, and how
close phi_1 and phi_2 run inside the thick film.

Usage:
  python scripts/plot_tf_profile_ref.py --in-dir DIR [--out-dir DIR]
"""

import argparse
import csv
from collections import defaultdict
from pathlib import Path

import numpy as np
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
        d = grouped[(r["point_id"], r["state"])]
        d["z"].append(float(r["z"]))
        d["phi1"].append(float(r["phi1"]))
        d["phi2"].append(float(r["phi2"]))
    return {k: {a: np.array(v) for a, v in d.items()} for k, d in grouped.items()}


def deviation_width(d, phi1_inf, phi2_inf, tol=RELAX_TOL):
    """Largest z where either profile still deviates from the far field."""
    dev = np.maximum(np.abs(d["phi1"] - phi1_inf), np.abs(d["phi2"] - phi2_inf))
    idx = np.where(dev > tol)[0]
    return float(d["z"][idx[-1]]) if len(idx) else 0.0


def report(prof, states, point_ids):
    print(f"{'pt':>3} {'state':>6} {'gamma':>13} {'cs':>10} {'phi1(0)':>9} "
          f"{'phi2(0)':>9} {'width':>7} {'max|p1-p2|':>11} {'in film':>9}", flush=True)
    for pid in point_ids:
        for state in sorted(states[pid]):
            s = states[pid][state]
            d = prof[(pid, state)]
            p1i, p2i = float(s["phi1_inf"]), float(s["phi2_inf"])
            w = deviation_width(d, p1i, p2i)
            inside = d["z"] <= w if w > 0 else np.zeros_like(d["z"], dtype=bool)
            diff = np.abs(d["phi1"] - d["phi2"])
            # away from the wall layer: outer half of the deviating region
            outer = inside & (d["z"] > 0.5 * w)
            print(f"{pid:>3} {state:>6} {float(s['gamma']):>13.8f} {float(s['cs']):>10.6f} "
                  f"{float(s['phi1_0']):>9.6f} {float(s['phi2_0']):>9.6f} {w:>7.3f} "
                  f"{diff[inside].max() if inside.any() else float('nan'):>11.6f} "
                  f"{diff[outer].max() if outer.any() else float('nan'):>9.6f}", flush=True)
        g = [float(states[pid][s]["gamma"]) for s in sorted(states[pid])]
        c = [float(states[pid][s]["cs"]) for s in sorted(states[pid])]
        if len(g) == 2:
            print(f"{pid:>3} -> gamma gap {abs(g[0] - g[1]):.3e} "
                  f"(relative {abs(g[0] - g[1]) / abs(g[0]):.3e}), "
                  f"cs jump {abs(c[0] - c[1]):.6f}", flush=True)


def draw(ax, prof, states, pid, zlim):
    s0 = states[pid][sorted(states[pid])[0]]
    p1i, p2i = float(s0["phi1_inf"]), float(s0["phi2_inf"])
    for state, style in zip(sorted(states[pid]), ("-", "--")):
        d = prof[(pid, state)]
        ax.plot(d["z"], d["phi1"], style, color=C1, lw=1.7)
        ax.plot(d["z"], d["phi2"], style, color=C2, lw=1.7)
    ax.axhline(p1i, color=C1, ls=":", lw=0.9)
    ax.axhline(p2i, color=C2, ls=":", lw=0.9)
    ax.set_xlim(0.0, zlim)
    ax.set_xlabel(r"$z$")
    ax.set_title(rf"{pid}: $\phi_{{1,\infty}} = {p1i:.4f}$, "
                 rf"$\phi_{{2,\infty}} = {p2i:.4f}$", fontsize=10)


def figure(prof, states, point_ids, zlims, out_path, title):
    n = len(point_ids)
    fig, axes = plt.subplots(1, n, figsize=(4.7 * n, 4.3), sharey=True)
    axes = [axes] if n == 1 else list(axes)
    for ax, pid in zip(axes, point_ids):
        draw(ax, prof, states, pid, zlims[pid])
    axes[0].set_ylabel(r"$\phi_i(z)$")
    handles = [
        plt.Line2D([], [], color=C1, ls="-", label=r"$\phi_1$ thin state"),
        plt.Line2D([], [], color=C2, ls="-", label=r"$\phi_2$ thin state"),
        plt.Line2D([], [], color=C1, ls="--", label=r"$\phi_1$ thick state"),
        plt.Line2D([], [], color=C2, ls="--", label=r"$\phi_2$ thick state"),
    ]
    axes[-1].legend(handles=handles, fontsize=8, loc="upper right")
    if title:
        fig.suptitle(title, fontsize=11)
        fig.tight_layout(rect=(0, 0, 1, 0.94))
    else:
        fig.tight_layout()
    fig.savefig(out_path, dpi=170)
    plt.close(fig)
    print(f"wrote {out_path}", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in-dir", required=True)
    ap.add_argument("--out-dir")
    ap.add_argument("--title",
                    default=r"T-f: $\chi_{12} = -8.5$, $\chi_{13} = \chi_{23} = 0$, "
                            r"$\omega_1 = 0.25$, $\omega_2 = -0.375$")
    args = ap.parse_args()

    in_dir = Path(args.in_dir)
    out_dir = Path(args.out_dir) if args.out_dir else in_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    prof = load_profiles(in_dir / "profiles.csv")
    states = defaultdict(dict)
    for r in read_rows(in_dir / "states.csv"):
        states[r["point_id"]][r["state"]] = r
    point_ids = sorted(states)

    report(prof, states, point_ids)

    zlim_cut, zlim_full = {}, {}
    for pid in point_ids:
        s0 = states[pid][sorted(states[pid])[0]]
        p1i, p2i = float(s0["phi1_inf"]), float(s0["phi2_inf"])
        widths = [deviation_width(prof[(pid, s)], p1i, p2i) for s in states[pid]]
        z_full = max(prof[(pid, s)]["z"].max() for s in states[pid])
        zlim_cut[pid] = min(z_full, max(max(widths) * (1.0 + PAD_FRAC), 0.5))
        zlim_full[pid] = z_full

    figure(prof, states, point_ids, zlim_cut, out_dir / "tf_profiles.png", args.title)
    figure(prof, states, point_ids, zlim_full, out_dir / "tf_profiles_fullz.png", args.title)


if __name__ == "__main__":
    main()
