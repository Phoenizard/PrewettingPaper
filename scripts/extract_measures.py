"""Extract pre-wetting occurrence measures for every case in a data/ tree.

Generic (not omega-specific): produces one master table covering all cases, one
row each. Each analysis angle (omega, chibb12, chibb11/22, chi) filters this
table by its own parameters. Rigorous conclusions are drawn on single-branch
cases only; multi-branch cases are flagged (is_single = 0) for separate study.

Measures per case, all in the raw (phi1_inf, phi2_inf) plane (no rescaling):
  length      projection span of the pre-wetting points on their best-fit line
  dist_mean   mean over points of shortest distance to the binodal
  dist_max    max  over points of shortest distance to the binodal
  residual_rms  RMS perpendicular distance to the best-fit line (kink indicator)
  n_points    number of pre-wetting points (fix_phi1 + fix_phi2 merged)
  n_branch    branch count from split_branches
  is_single   1 if single-branch (n_branch == 1 and residual within threshold)

Usage:
  python scripts/extract_measures.py --data-root <pw-space/data> \
      --out <analysis/prewetting_measures/measures.csv> \
      [--residual-threshold 0.003]

The residual threshold defaults to a placeholder; calibrate it from the actual
residual distribution (see the calibration step) before drawing conclusions.
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
    "n_points", "n_branch", "is_single",
    "residual_rms", "length", "dist_mean", "dist_max",
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


def measure_case(case, residual_threshold):
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
        "n_points": 0, "n_branch": 0, "is_single": 0,
        "residual_rms": "", "length": "", "dist_mean": "", "dist_max": "",
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

    if len(pts) < 2:
        # A single point: length 0, no line to fit; still measure distance.
        row["length"] = 0.0
        row["residual_rms"] = 0.0
        row["flag"] = "too_few_points"
    else:
        line = geom.fit_line_pca(pts)
        row["length"] = round(geom.projection_span(pts, line), 8)
        rms = geom.residual_rms(pts, line)
        row["residual_rms"] = round(rms, 8)
        row["is_single"] = int(n_branch == 1 and rms <= residual_threshold)

    binodal = _load_binodal(case.binodal_path)
    if binodal:
        d = geom.min_dist_to_set(pts, binodal)
        row["dist_mean"] = round(float(d.mean()), 8)
        row["dist_max"] = round(float(d.max()), 8)
    else:
        row["flag"] = (row["flag"] + ";" if row["flag"] else "") + "no_binodal"

    return row


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-root", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--residual-threshold", type=float, default=0.003)
    args = ap.parse_args()

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    rows = []
    for case in iter_cases(args.data_root):
        rows.append(measure_case(case, args.residual_threshold))

    with open(out_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)

    n_single = sum(r["is_single"] for r in rows)
    n_multi = sum(1 for r in rows if r["n_branch"] > 1)
    n_empty = sum(1 for r in rows if r["flag"] and "no_prewetting" in r["flag"])
    print(f"cases: {len(rows)}  single-branch: {n_single}  "
          f"multi-branch: {n_multi}  no-prewetting: {n_empty}")
    print(f"wrote {out_path}")


if __name__ == "__main__":
    main()
