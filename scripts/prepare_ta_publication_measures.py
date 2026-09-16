"""Merge the T-a boundary audit into a publication-only measures table."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def key(row):
    return (float(row["omega_1"]), float(row["omega_2"]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", required=True)
    ap.add_argument("--reentry", required=True)
    ap.add_argument("--audit", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    with open(args.audit, newline="") as fh:
        audit = {float(r["wall_affinity"]): r for r in csv.DictReader(fh)
                 if r["side"] in {"strong", "weak"}}
    with open(args.reentry, newline="") as fh:
        rerun = {key(r): r for r in csv.DictReader(fh)}
    with open(args.base, newline="") as fh:
        reader = csv.DictReader(fh); fields = reader.fieldnames; rows = list(reader)

    stage = "chi12_0__chi13_2p8__chi23_0"
    for row in rows:
        if row["stage"] != stage or any(float(row[k]) != 0.0 for k in
                                        ("chi_bb_11", "chi_bb_22", "chi_bb_12")):
            continue
        w1, w2 = key(row)
        if w1 != w2 or w1 not in audit:
            continue
        verdict = audit[w1]["status"]
        if verdict == "valid" and (w1, w2) in rerun:
            source = rerun[(w1, w2)]
            for field in ("pw_length", "dist_mean", "n_points", "n_segments",
                          "gap_total", "full_length", "dist_max", "dist_min",
                          "residual_rms", "pca_span", "n_branch", "is_single"):
                if field in row and field in source:
                    row[field] = source[field]
            row["flag"] = ""
        else:
            row["pw_length"] = "0.0"
            row["dist_mean"] = ""
            for field in ("dist_max", "dist_min"):
                if field in row:
                    row[field] = ""
            row["flag"] = "no_prewetting_boundary_audit"

    out = Path(args.out); out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader(); writer.writerows(rows)


if __name__ == "__main__":
    main()
