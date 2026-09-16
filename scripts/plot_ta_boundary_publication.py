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
from paperfig import configure, panel_label, save_pdf_png

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


def load_status(path: Path) -> dict[str, str]:
    with path.open(newline="") as fh:
        return {r["code"]: r["status"] for r in csv.DictReader(fh)}


def draw(data_root: Path, audit_csv: Path, out: Path, markers: bool) -> None:
    configure()
    status = load_status(audit_csv)
    raw_binodal = read_xy(case_dir(data_root, "m0p3") / "binodal.csv", "phi1", "phi2")
    binodal = physical_binodal(raw_binodal)
    bx, by = pchip_xy(binodal, 500)

    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.25))
    xlims = ((0.072, 0.108), (0.090, 0.190))
    ylims = ((-0.001, 0.041), (0.075, 0.203))
    titles = ("strong wall attraction", "weak wall attraction")
    for ax, xlim, ylim, title in zip(axes, xlims, ylims, titles):
        visible = (by >= max(ylim[0], by.min())) & (by <= min(ylim[1], by.max()))
        ax.fill_betweenx(by[visible], bx[visible], xlim[1], color="0.94", zorder=0)
        ax.plot(bx[visible], by[visible], color="0.55", lw=1.5, label="binodal", zorder=2)
        ax.set_xlim(*xlim); ax.set_ylim(*ylim); ax.set_title(title)
        ax.set_xlabel(r"$\phi_{1,\infty}$"); ax.set_ylabel(r"$\phi_{2,\infty}$")
        ax.grid(alpha=0.14, linewidth=0.5)
    axes[0].text(0.69, 0.89, "two-phase", transform=axes[0].transAxes,
                 color="0.45", fontsize=8)
    axes[1].text(0.73, 0.88, "two-phase", transform=axes[1].transAxes,
                 color="0.45", fontsize=8)

    valid_strong = [(c, v) for c, v in DISPLAY[:4] if status.get(c) == "valid"]
    colors = [viridis(v) for v in np.linspace(0.18, 0.82, max(len(valid_strong), 1))]
    for (code, value), color in zip(valid_strong, colors):
        points = read_xy(case_dir(data_root, code) / "pw_line.csv", "phi1_inf", "phi2_inf")
        x, y = pchip_xy(points)
        axes[0].plot(x, y, color=color, lw=2.0,
                     label=rf"$\omega_1=\omega_2={value:g}$", zorder=3)
        if markers:
            count = min(9, len(points))
            idx = np.unique(np.linspace(0, len(points) - 1, count).round().astype(int))
            ordered = points[np.argsort(points[:, 1])]
            axes[0].scatter(ordered[idx, 0], ordered[idx, 1], s=15, color=color,
                            edgecolor="white", linewidth=0.25, zorder=4)

    axes[1].text(0.50, 0.42,
                 "no pre-wetting in the\none-phase region",
                 transform=axes[1].transAxes, ha="center", va="center", fontsize=9)
    axes[0].legend(frameon=False, loc="upper left", handlelength=1.4,
                   handletextpad=0.45, labelspacing=0.25, borderaxespad=0.25)
    axes[1].legend(frameon=False, loc="upper left", handlelength=1.4,
                   handletextpad=0.45, borderaxespad=0.25)
    panel_label(axes[0], "(a)"); panel_label(axes[1], "(b)")
    fig.subplots_adjust(left=0.09, right=0.99, bottom=0.16, top=0.87, wspace=0.28)
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
