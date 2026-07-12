"""Render a research-grade phase map per case: binodal + pre-wetting line.

Style (agreed): clean paper style, white background; binodal as one thin grey
line; pre-wetting points as small solid dots; the two scan directions
(fix_phi1, fix_phi2) are NOT distinguished (same line, same marker); points are
coloured by branch ONLY when the case is multi-branch, otherwise a single
colour. Parameters go in the title (LaTeX). Axes are the raw
(phi1_inf, phi2_inf) plane, [0,1] x [0,1], aspect equal.

The rendered PNG is written back into each case directory (measure_map.png),
so it lives next to the case's CSVs as a paper asset and a visual-inspection aid.

Usage:
  python scripts/plot_all_cases.py --case-dir <one case dir>            # single
  python scripts/plot_all_cases.py --data-root <pw-space/data>          # all
      [--residual-threshold 0.003] [--out-name measure_map.png]
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

import geom
from cases import iter_cases

SINGLE_COLOR = "#1f3b73"  # deep blue for single-branch pre-wetting points
BINODAL_COLOR = "0.6"


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


def render_case(case_dir, params, residual_threshold, out_name):
    case_dir = Path(case_dir)
    binodal = [(float(a), float(b)) for a, b in _read_csv_rows(case_dir / "binodal.csv")]
    pts = [(float(r[1]), float(r[2])) for r in _read_csv_rows(case_dir / "pw_line.csv")]

    fig, ax = plt.subplots(figsize=(5, 5))
    if binodal:
        bx, by = zip(*binodal)
        ax.plot(bx, by, lw=1.0, color=BINODAL_COLOR, label="binodal")

    if pts:
        labels = geom.split_branches(pts)
        n_branch = int(labels.max()) + 1 if len(labels) else 0
        rms = geom.residual_rms(pts) if len(pts) >= 3 else 0.0
        multi = n_branch > 1 or rms > residual_threshold
        xs, ys = zip(*pts)
        if multi:
            cmap = plt.get_cmap("viridis")
            ncol = max(n_branch, 1)
            for bid in range(ncol):
                sel = [i for i, b in enumerate(labels) if b == bid]
                if sel:
                    ax.scatter([xs[i] for i in sel], [ys[i] for i in sel],
                               s=12, color=cmap(bid / max(ncol - 1, 1)),
                               label=f"PW branch {bid}")
        else:
            ax.scatter(xs, ys, s=12, color=SINGLE_COLOR, label="pre-wetting")

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
    ap.add_argument("--residual-threshold", type=float, default=0.003)
    ap.add_argument("--out-name", default="measure_map.png")
    args = ap.parse_args()

    if args.case_dir:
        from params import parse_case_rel
        case_dir = Path(args.case_dir)
        # rel is the last three path components (chi/om/chibb).
        rel = "/".join(case_dir.parts[-3:])
        params = parse_case_rel(rel)
        out = render_case(case_dir, params, args.residual_threshold, args.out_name)
        print(f"wrote {out}")
        return

    if not args.data_root:
        ap.error("provide either --case-dir or --data-root")

    n = 0
    for case in iter_cases(args.data_root):
        render_case(case.pw_path.parent, case.params,
                    args.residual_threshold, args.out_name)
        n += 1
        if n % 100 == 0:
            print(f"  rendered {n} cases")
    print(f"rendered {n} cases")


if __name__ == "__main__":
    main()
