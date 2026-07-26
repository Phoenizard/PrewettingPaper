"""Render a research-grade phase map per case: binodal + pre-wetting line.

Style (agreed): clean paper style, white background; binodal as thin grey
scatter; pre-wetting points as small solid dots (the two scan directions
fix_phi1/fix_phi2 are NOT distinguished). The minimum spanning tree over the
points is overlaid: solid blue edges are the segment arclength that makes up
pw_length; dashed red edges are the gaps (MST edges longer than gap_tol, cut and
excluded from pw_length). A corner note shows L (pw_length), d-bar (dist_mean),
and seg (n_segments). No branch splitting. Parameters go in the title (LaTeX).
Axes are the raw (phi1_inf, phi2_inf) plane, [0,1] x [0,1], aspect equal.

The rendered PNG is written back into each case directory (measure_map.png),
so it lives next to the case's CSVs as a paper asset and a visual-inspection aid.

Usage:
  python scripts/plot_all_cases.py --case-dir <one case dir>            # single
  python scripts/plot_all_cases.py --data-root <pw-space/data>          # all
      [--gap-tol 0.01] [--out-name measure_map.png]
"""

import argparse
import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

import numpy as np
from scipy.spatial.distance import pdist, squareform
from scipy.sparse.csgraph import minimum_spanning_tree

import geom
from cases import iter_cases

SINGLE_COLOR = "#1f3b73"  # deep blue for pre-wetting points
BINODAL_COLOR = "0.6"
MST_COLOR = "#1f3b73"     # segment edges (arclength contributing to pw_length)
GAP_COLOR = "#c02a2a"     # cut bridge edges (gaps, excluded from pw_length)


def _read_csv_rows(path):
    with open(path, newline="") as fh:
        reader = csv.reader(fh)
        next(reader, None)
        return [row for row in reader if row]


def _title(params):
    p = params
    return (
        rf"$\chi_{{12}}={p['chi_12']:g}\ \chi_{{13}}={p['chi_13']:g}\ "
        rf"\chi_{{23}}={p['chi_23']:g}\ \ "
        rf"\omega_1={p['omega_1']:g}\ \omega_2={p['omega_2']:g}\ \ "
        rf"\chi^{{bb}}_{{11}}={p['chi_bb_11']:g}\ "
        rf"\chi^{{bb}}_{{22}}={p['chi_bb_22']:g}\ "
        rf"\chi^{{bb}}_{{12}}={p['chi_bb_12']:g}$"
    )


def render_case(case_dir, params, gap_tol, out_name):
    case_dir = Path(case_dir)
    binodal = [(float(a), float(b)) for a, b in _read_csv_rows(case_dir / "binodal.csv")]
    pts = [(float(r[1]), float(r[2])) for r in _read_csv_rows(case_dir / "pw_line.csv")]

    fig, ax = plt.subplots(figsize=(5, 5))
    if binodal:
        bx, by = zip(*binodal)
        # scatter, not plot: binodal points are not path-ordered and may form
        # several disjoint loops; connecting them with a line draws false edges.
        ax.scatter(bx, by, s=2, color=BINODAL_COLOR, label="binodal", linewidths=0)

    if pts:
        parr = np.asarray(pts, dtype=float)
        xs, ys = parr[:, 0], parr[:, 1]
        ax.scatter(xs, ys, s=10, color=SINGLE_COLOR, zorder=3,
                   label="pre-wetting")

        # Overlay the MST: kept edges (<= gap_tol) are the arclength that makes
        # up pw_length; cut edges (> gap_tol) are gaps between segments.
        length = dist_mean = None
        n_segments = 0
        if len(parr) >= 2:
            length, n_segments, _gap, _full = geom.mst_length(parr, gap_tol)
            mst = minimum_spanning_tree(squareform(pdist(parr))).tocoo()
            for i, j, w in zip(mst.row, mst.col, mst.data):
                seg = w <= gap_tol
                ax.plot([xs[i], xs[j]], [ys[i], ys[j]],
                        color=MST_COLOR if seg else GAP_COLOR,
                        lw=1.0 if seg else 0.8,
                        ls="-" if seg else "--", zorder=2)
        if binodal:
            d = geom.min_dist_to_set(parr, binodal)
            dist_mean = float(d.mean())
        note = []
        if length is not None:
            note.append(rf"$L={length:.3f}$")
        if dist_mean is not None:
            note.append(rf"$\bar d={dist_mean:.3f}$")
        note.append(rf"seg$={n_segments}$")
        ax.text(0.03, 0.97, "  ".join(note), transform=ax.transAxes,
                fontsize=8, va="top", ha="left")

    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.0)
    ax.set_aspect("equal")
    ax.set_xlabel(r"$\phi_{1,\infty}$")
    ax.set_ylabel(r"$\phi_{2,\infty}$")
    ax.set_title(_title(params), fontsize=7)
    ax.legend(loc="upper right", fontsize=7, frameon=False)
    fig.tight_layout()
    out_path = case_dir / out_name
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    return out_path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--case-dir", help="render a single case directory")
    ap.add_argument("--data-root", help="render every case under this data root")
    ap.add_argument("--gap-tol", type=float, default=0.015,
                    help="MST edges longer than this are drawn as gaps (dashed); "
                         "0.015 calibrated on the omega topic")
    ap.add_argument("--out-name", default="measure_map.png")
    args = ap.parse_args()

    if args.case_dir:
        from cases import parse_case_rel
        case_dir = Path(args.case_dir)
        # rel is the last three path components (chi/om/chibb).
        rel = "/".join(case_dir.parts[-3:])
        params = parse_case_rel(rel)
        out = render_case(case_dir, params, args.gap_tol, args.out_name)
        print(f"wrote {out}")
        return

    if not args.data_root:
        ap.error("provide either --case-dir or --data-root")

    n = 0
    for case in iter_cases(args.data_root):
        render_case(case.pw_path.parent, case.params,
                    args.gap_tol, args.out_name)
        n += 1
        if n % 100 == 0:
            print(f"  rendered {n} cases")
    print(f"rendered {n} cases")


if __name__ == "__main__":
    main()
