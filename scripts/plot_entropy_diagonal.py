"""Mixing-entropy part of f_b along the equal-solute line phi1 = phi2.

Plots x ln x + y ln y + (1-x-y) ln(1-x-y) with x = y = phi, so the curve is
2 phi ln phi + (1 - 2 phi) ln(1 - 2 phi) for phi in (0, 0.5). No interaction
term: this is the entropy alone, the part that opposes phase separation.

Usage:
  python scripts/plot_entropy_diagonal.py [--out-dir DIR]
"""

import argparse
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", default="out/analysis/Tf_omega_cross")
    args = ap.parse_args()
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    phi = np.linspace(1e-6, 0.5 - 1e-6, 4000)
    s = 2.0 * phi * np.log(phi) + (1.0 - 2.0 * phi) * np.log(1.0 - 2.0 * phi)

    i = int(np.argmin(s))
    fig, ax = plt.subplots(figsize=(5.6, 4.4))
    ax.plot(phi, s, "-", color="#1f77b4", lw=1.6)
    ax.plot([phi[i]], [s[i]], "o", ms=4.5, color="#d62728", zorder=3)
    ax.annotate(rf"minimum at $\phi = {phi[i]:.3f}$, value ${s[i]:.3f}$",
                xy=(phi[i], s[i]), xytext=(0.30, s[i] + 0.10),
                fontsize=9, color="#d62728")
    ax.axhline(0.0, color="0.85", lw=0.8, zorder=0)
    ax.set_xlabel(r"$\phi$  (with $\phi_1 = \phi_2 = \phi$)")
    ax.set_ylabel(r"$2\phi\ln\phi + (1-2\phi)\ln(1-2\phi)$")
    ax.set_title("Mixing entropy along the equal-solute line", fontsize=10)
    ax.set_xlim(0.0, 0.5)
    fig.tight_layout()
    path = out / "entropy_diagonal.png"
    fig.savefig(path, dpi=170)
    plt.close(fig)
    print(f"min at phi = {phi[i]:.6f}, value = {s[i]:.6f}", flush=True)
    print(f"written to {path}", flush=True)


if __name__ == "__main__":
    main()
