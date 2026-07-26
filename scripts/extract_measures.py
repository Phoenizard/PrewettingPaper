"""Extract pre-wetting occurrence measures for every case in a data/ tree.

Generic (not omega-specific): produces one master table covering all cases, one
row each. Each analysis angle (omega, chibb12, chibb11/22, chi) filters this
table by its own parameters. Rigorous conclusions are drawn on single-branch
cases only; multi-branch cases are flagged (is_single = 0) for separate study.

Headline measures (paper), raw (phi1_inf, phi2_inf) plane, no rescaling:
  pw_length   total pre-wetting-line length = MST arclength within clusters,
              gaps excluded (bend-robust; replaces the old PCA single-line span)
  dist_mean   mean over points of shortest distance to the binodal (typical depth
              into the one-phase region)

Diagnostics:
  n_points    number of pre-wetting points (fix_phi1 + fix_phi2 merged)
  n_segments  connected components after cutting MST edges > gap_tol (clusters)
  gap_total   summed length of the cut bridge edges (how far segments pull apart)
  full_length MST arclength WITHOUT the cut (extent including gaps); vs pw_length
              this says whether the gap carries physical meaning
  dist_max    max shortest distance to the binodal (deepest reach)
  dist_min    min shortest distance to the binodal (closest hug)
  residual_rms  RMS perpendicular distance to the best-fit line (straightness)
  pca_span    old PCA single-line span, kept to show why MST wins on bent lines
  n_branch    split_branches count (unreliable; informational only)
  is_single   informational only; conclusions use the merged whole, no filtering

Usage:
  python scripts/extract_measures.py --data-root <pw-space/data> \
      --out <analysis/prewetting_measures/measures.csv> \
      [--gap-tol 0.01] [--residual-threshold 0.003]

gap_tol cuts MST edges longer than it (segment separators); it must exceed the
natural inter-point spacing of a single scan line or a straight line fragments.
Calibrate it from the actual point spacing before drawing conclusions.
"""

import argparse
import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import geom
from cases import iter_cases

FIELDS = [
    "stage", "rel",
    "chi_12", "chi_13", "chi_23",
    "omega_1", "omega_2",
    "chi_bb_11", "chi_bb_22", "chi_bb_12",
    # headline
    "pw_length", "dist_mean",
    # diagnostics
    "n_points", "n_segments", "gap_total", "full_length",
    "dist_max", "dist_min", "residual_rms", "pca_span",
    "n_branch", "is_single",
    "flag",
]


def _read_csv_rows(path):
    """Read a CSV skipping the header; tolerant of CRLF line endings."""
    with open(path, newline="") as fh:
        reader = csv.reader(fh)
        next(reader, None)
        return [row for row in reader if row]


def _load_pw_points(pw_path):
    """Return list of (phi1, phi2); the two scan directions are merged."""
    return [(float(r[1]), float(r[2])) for r in _read_csv_rows(pw_path)]


def _load_binodal(binodal_path):
    return [(float(r[0]), float(r[1])) for r in _read_csv_rows(binodal_path)]


def measure_case(case, residual_threshold, gap_tol):
    """Compute the measure row for one case. Never raises on data content;
    degenerate cases get a flag and blank measures."""
    p = case.params
    row = {
        "stage": case.rel.split("/")[0],
        "rel": case.rel,
        "chi_12": p["chi_12"], "chi_13": p["chi_13"], "chi_23": p["chi_23"],
        "omega_1": p["omega_1"], "omega_2": p["omega_2"],
        "chi_bb_11": p["chi_bb_11"], "chi_bb_22": p["chi_bb_22"],
        "chi_bb_12": p["chi_bb_12"],
        "pw_length": "", "dist_mean": "",
        "n_points": 0, "n_segments": 0, "gap_total": "", "full_length": "",
        "dist_max": "", "dist_min": "", "residual_rms": "", "pca_span": "",
        "n_branch": 0, "is_single": 0,
        "flag": "",
    }

    pts = _load_pw_points(case.pw_path)
    row["n_points"] = len(pts)
    if not pts:
        row["flag"] = "no_prewetting"
        return row

    labels = geom.split_branches(pts)
    n_branch = int(labels.max()) + 1 if len(labels) else 0
    row["n_branch"] = n_branch

    # Headline length + gap diagnostics from the bend-robust MST.
    length, n_segments, gap_total, full_length = geom.mst_length(pts, gap_tol)
    row["pw_length"] = round(length, 8)
    row["n_segments"] = n_segments
    row["gap_total"] = round(gap_total, 8)
    row["full_length"] = round(full_length, 8)

    if len(pts) < 2:
        # A single point: no line to fit; still measure distance below.
        row["residual_rms"] = 0.0
        row["pca_span"] = 0.0
        row["flag"] = "too_few_points"
    else:
        line = geom.fit_line_pca(pts)
        row["pca_span"] = round(geom.projection_span(pts, line), 8)
        rms = geom.residual_rms(pts, line)
        row["residual_rms"] = round(rms, 8)
        row["is_single"] = int(n_branch == 1 and rms <= residual_threshold)

    binodal = _load_binodal(case.binodal_path)
    if binodal:
        d = geom.min_dist_to_set(pts, binodal)
        row["dist_mean"] = round(float(d.mean()), 8)
        row["dist_max"] = round(float(d.max()), 8)
        row["dist_min"] = round(float(d.min()), 8)
    else:
        row["flag"] = (row["flag"] + ";" if row["flag"] else "") + "no_binodal"

    return row


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-root", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--gap-tol", type=float, default=0.015,
                    help="cut MST edges longer than this (segment separators); "
                         "0.015 calibrated on the omega topic: ~10x the natural "
                         "in-line point spacing (~0.0014), below the 0.02 grid step")
    ap.add_argument("--residual-threshold", type=float, default=0.003)
    args = ap.parse_args()

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    rows = []
    for case in iter_cases(args.data_root):
        rows.append(measure_case(case, args.residual_threshold, args.gap_tol))

    with open(out_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)

    n_multiseg = sum(1 for r in rows if isinstance(r["n_segments"], int)
                     and r["n_segments"] > 1)
    n_empty = sum(1 for r in rows if r["flag"] and "no_prewetting" in r["flag"])
    print(f"cases: {len(rows)}  multi-segment (n_segments>1): {n_multiseg}  "
          f"no-prewetting: {n_empty}  gap_tol={args.gap_tol}")
    print(f"wrote {out_path}")


if __name__ == "__main__":
    main()
