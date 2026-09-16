"""Publication figures from the lightweight T-a boundary audit."""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.cm import viridis
from scipy.interpolate import PchipInterpolator

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from paperfig import configure, save_pdf_png

CHI = "chi12_0__chi13_2p8__chi23_0"
CHIBB = "chibb11_0__chibb22_0__chibb12_0"
DISPLAY = [("m0p5", -0.50), ("m0p48", -0.48), ("m0p46", -0.46), ("m0p44", -0.44),
           ("m0p1", -0.10), ("m0p08", -0.08), ("m0p06", -0.06), ("m0p04", -0.04)]


def case_dir(root: Path, code: str) -> Path:
    return root / CHI / f"om1_{code}__om2_{code}" / CHIBB


def read_xy(path: Path, x: str, y: str) -> np.ndarray:
    with path.open(newline="") as fh:
        rows = list(csv.DictReader(fh))
    return np.asarray([(float(r[x]), float(r[y])) for r in rows], dtype=float)


def physical_binodal(points: np.ndarray) -> np.ndarray:
    jumps = np.linalg.norm(np.diff(points, axis=0), axis=1)
    cut = np.where(jumps > 0.05)[0]
    return points[: int(cut[0]) + 1] if len(cut) else points


def pchip_xy(points: np.ndarray, n: int = 240) -> tuple[np.ndarray, np.ndarray]:
    # The solvent-rich branch and the accepted prewetting curves are monotone in
    # phi2 over the displayed ranges.  Average duplicate slices before PCHIP.
    by_y: dict[float, list[float]] = {}
    for x, y in points:
        by_y.setdefault(float(y), []).append(float(x))
    y = np.asarray(sorted(by_y))
    x = np.asarray([np.mean(by_y[v]) for v in y])
    if len(y) == 2:
        yd = np.linspace(y[0], y[-1], n)
        return np.interp(yd, y, x), yd
    yd = np.linspace(y[0], y[-1], n)
    return PchipInterpolator(y, x)(yd), yd


def smooth_binodal(points: np.ndarray, ylim: tuple[float, float], n: int = 500) -> tuple[np.ndarray, np.ndarray]:
    """Return a smooth local fit to the solvent-rich bulk coexistence branch."""
    margin = 0.006
    local = points[(points[:, 1] >= ylim[0] - margin) &
                   (points[:, 1] <= ylim[1] + margin)]
    by_y: dict[float, list[float]] = {}
    for x, y in local:
        by_y.setdefault(float(y), []).append(float(x))
    y = np.asarray(sorted(by_y))
    x = np.asarray([np.mean(by_y[v]) for v in y])
    degree = min(3, len(y) - 1)
    fit = np.polynomial.Polynomial.fit(y, x, degree)
    yd = np.linspace(max(ylim[0], y.min()), min(ylim[1], y.max()), n)
    return fit(yd), yd


def load_status(path: Path) -> dict[str, str]:
    with path.open(newline="") as fh:
        return {r["code"]: r["status"] for r in csv.DictReader(fh)}


def draw(data_root: Path, audit_csv: Path, out: Path, markers: bool) -> None:
    configure()
    status = load_status(audit_csv)
    raw_binodal = read_xy(case_dir(data_root, "m0p3") / "binodal.csv", "phi1", "phi2")
    binodal = physical_binodal(raw_binodal)
    xlim = (0.072, 0.108)
    ylim = (-0.001, 0.041)
    bx, by = smooth_binodal(binodal, ylim)

    fig, ax = plt.subplots(figsize=(4.65, 3.45))
    ax.fill_betweenx(by, bx, xlim[1], color="0.90", zorder=0)
    ax.plot(bx, by, color="0.08", lw=2.2,
            label="bulk coexistence boundary", zorder=5)
    ax.set_xlim(*xlim); ax.set_ylim(*ylim)
    ax.set_title("strengthening both wall attractions")
    ax.set_xlabel(r"$\phi_{1,\infty}$"); ax.set_ylabel(r"$\phi_{2,\infty}$")
    ax.grid(alpha=0.13, linewidth=0.5)
    ax.text(0.56, 0.66, "uniform bulk", transform=ax.transAxes,
            ha="center", color="0.20", fontsize=8)
    ax.text(0.83, 0.84, "bulk phase\nseparation", transform=ax.transAxes,
            ha="center", va="center", color="0.25", fontsize=8)

    valid_strong = [(c, v) for c, v in DISPLAY[:4] if status.get(c) == "valid"]
    colors = [viridis(v) for v in np.linspace(0.18, 0.82, max(len(valid_strong), 1))]
    for (code, value), color in zip(valid_strong, colors):
        points = read_xy(case_dir(data_root, code) / "pw_line.csv", "phi1_inf", "phi2_inf")
        x, y = pchip_xy(points)
        ax.plot(x, y, color=color, lw=2.2,
                label=rf"$\omega_1=\omega_2={value:g}$", zorder=3)
        if markers:
            count = min(9, len(points))
            idx = np.unique(np.linspace(0, len(points) - 1, count).round().astype(int))
            ordered = points[np.argsort(points[:, 1])]
            ax.scatter(ordered[idx, 0], ordered[idx, 1], s=17, color=color,
                       edgecolor="white", linewidth=0.25, zorder=4)

    ax.legend(frameon=False, loc="upper left", handlelength=1.7,
              handletextpad=0.45, labelspacing=0.28, borderaxespad=0.35,
              fontsize=8)
    fig.subplots_adjust(left=0.16, right=0.98, bottom=0.17, top=0.86)
    save_pdf_png(fig, out)
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-root", required=True)
    ap.add_argument("--audit", required=True)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()
    out = Path(args.out_dir); out.mkdir(parents=True, exist_ok=True)
    draw(Path(args.data_root), Path(args.audit), out / "ta_boundary_curves_markers.pdf", True)
    draw(Path(args.data_root), Path(args.audit), out / "ta_boundary_curves_only.pdf", False)


if __name__ == "__main__":
    main()
