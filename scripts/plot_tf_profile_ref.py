"""Publication T-f thin/thick profiles from the existing archived solutions."""

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from paperfig import configure, panel_label, save_pdf_png

C1 = "#2166ac"
C2 = "#b2182b"
RELAX_TOL = 1e-3
PAD_FRAC = 0.25


def read_rows(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def load_profiles(path):
    grouped = defaultdict(lambda: {"z": [], "phi1": [], "phi2": []})
    for row in read_rows(path):
        values = grouped[(row["point_id"], row["state"])]
        for key in ("z", "phi1", "phi2"):
            values[key].append(float(row[key]))
    return {key: {name: np.asarray(vals) for name, vals in data.items()}
            for key, data in grouped.items()}


def deviation_width(data, phi1_inf, phi2_inf, tol=RELAX_TOL):
    dev = np.maximum(np.abs(data["phi1"] - phi1_inf),
                     np.abs(data["phi2"] - phi2_inf))
    indices = np.where(dev > tol)[0]
    return float(data["z"][indices[-1]]) if len(indices) else 0.0


def report(profiles, states, point_ids):
    print("point,state,gamma,cs,phi1_0,phi2_0,width", flush=True)
    for point in point_ids:
        for state in sorted(states[point]):
            meta = states[point][state]
            data = profiles[(point, state)]
            p1 = float(meta["phi1_inf"])
            p2 = float(meta["phi2_inf"])
            width = deviation_width(data, p1, p2)
            print(f"{point},{state},{float(meta['gamma']):.10g},"
                  f"{float(meta['cs']):.10g},{float(meta['phi1_0']):.10g},"
                  f"{float(meta['phi2_0']):.10g},{width:.10g}", flush=True)


def draw(ax, profiles, states, point, zlim):
    ordered = sorted(states[point])
    meta = states[point][ordered[0]]
    p1 = float(meta["phi1_inf"])
    p2 = float(meta["phi2_inf"])
    for state, style in zip(ordered, ("-", "--")):
        data = profiles[(point, state)]
        ax.plot(data["z"], data["phi1"], style, color=C1, lw=1.6)
        ax.plot(data["z"], data["phi2"], style, color=C2, lw=1.6)
    ax.axhline(p1, color=C1, ls=":", lw=0.9, alpha=0.8)
    ax.axhline(p2, color=C2, ls=":", lw=0.9, alpha=0.8)
    ax.set_xlim(0.0, zlim)
    ax.set_xlabel(r"$z$")
    ax.text(0.97, 0.96,
            rf"$\phi_{{1,\infty}}={p1:.3f}$" "\n" rf"$\phi_{{2,\infty}}={p2:.3f}$",
            transform=ax.transAxes, ha="right", va="top", fontsize=8.5)


def figure(profiles, states, point_ids, zlims, output):
    fig, axes = plt.subplots(1, len(point_ids), figsize=(7.2, 2.75), sharey=True)
    axes = [axes] if len(point_ids) == 1 else list(axes)
    for index, (ax, point) in enumerate(zip(axes, point_ids)):
        draw(ax, profiles, states, point, zlims[point])
        panel_label(ax, f"({chr(ord('a') + index)}) {point}")
    axes[0].set_ylabel(r"$\phi_i(z)$")
    handles = [
        plt.Line2D([], [], color=C1, ls="-", label=r"$\phi_1$, thin"),
        plt.Line2D([], [], color=C2, ls="-", label=r"$\phi_2$, thin"),
        plt.Line2D([], [], color=C1, ls="--", label=r"$\phi_1$, thick"),
        plt.Line2D([], [], color=C2, ls="--", label=r"$\phi_2$, thick"),
    ]
    fig.legend(handles=handles, loc="upper center", ncol=4, frameon=False,
               bbox_to_anchor=(0.55, 1.01), columnspacing=1.0, handlelength=2.2)
    fig.subplots_adjust(left=0.09, right=0.995, bottom=0.19, top=0.78, wspace=0.12)
    return save_pdf_png(fig, output)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in-dir", required=True)
    ap.add_argument("--out-dir")
    args = ap.parse_args()
    configure()
    in_dir = Path(args.in_dir)
    out_dir = Path(args.out_dir) if args.out_dir else in_dir
    profiles = load_profiles(in_dir / "profiles.csv")
    states = defaultdict(dict)
    for row in read_rows(in_dir / "states.csv"):
        states[row["point_id"]][row["state"]] = row
    point_ids = sorted(states)
    report(profiles, states, point_ids)
    zlims = {}
    for point in point_ids:
        meta = states[point][sorted(states[point])[0]]
        p1 = float(meta["phi1_inf"])
        p2 = float(meta["phi2_inf"])
        widths = [deviation_width(profiles[(point, state)], p1, p2)
                  for state in states[point]]
        full = max(profiles[(point, state)]["z"].max() for state in states[point])
        zlims[point] = min(full, max(max(widths) * (1.0 + PAD_FRAC), 0.5))
    paths = figure(profiles, states, point_ids, zlims, out_dir / "tf_profiles.png")
    plt.close("all")
    print("figure:", *paths)


if __name__ == "__main__":
    main()
