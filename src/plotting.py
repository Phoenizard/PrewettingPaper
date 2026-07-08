"""Render overlay.png from pw_line.csv + binodal.csv (plot step only)."""

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def _read_rows(path):
    with open(path, newline="") as fh:
        reader = csv.reader(fh)
        next(reader)
        return list(reader)


def render_overlay(case_dir, title=None):
    case_dir = Path(case_dir)
    binodal = [(float(a), float(b)) for a, b in _read_rows(case_dir / "binodal.csv")]
    pw = [(s, float(a), float(b)) for s, a, b in _read_rows(case_dir / "pw_line.csv")]

    fig, ax = plt.subplots(figsize=(6, 6))
    if binodal:
        bx, by = zip(*binodal)
        ax.scatter(bx, by, s=4, c="0.65", label="binodal")
    for source, marker, color in (
        ("fix_phi2", "o", "tab:blue"),
        ("fix_phi1", "^", "tab:red"),
    ):
        pts = [(a, b) for s, a, b in pw if s == source]
        if pts:
            xs, ys = zip(*pts)
            ax.scatter(xs, ys, s=10, marker=marker, c=color, label=f"pre-wetting ({source})")
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.0)
    ax.set_aspect("equal")
    ax.set_xlabel(r"$\phi_{1,\infty}$")
    ax.set_ylabel(r"$\phi_{2,\infty}$")
    if title:
        ax.set_title(title, fontsize=8)
    ax.legend(loc="upper right", fontsize=8)
    fig.tight_layout()
    out_path = case_dir / "overlay.png"
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    return out_path
