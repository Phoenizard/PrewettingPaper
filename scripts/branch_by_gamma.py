"""Recompute the thermodynamics behind a case's pre-wetting points.

Branch identity is a physical question, not a geometric one: two pre-wetting
points belong to the same branch when the same pair of surface states (one thin,
one thick) carries the equal-gamma transition continuously between them. gamma
equality holds at every pre-wetting point (that is the transition), so it cannot
tell branches apart; the wall adsorption cs = cs1 + cs2 (from surface_metrics)
does — it tracks the actual surface state and jumps when the state pair changes.

For each pre-wetting point (phi1_inf, phi2_inf) this sets the far field to that
point, solves the thin branch and the thick branch (branch_guess seeds), and
records (gamma_thin, gamma_thick, cs_thin, cs_thick). Downstream, walk the line
and split a branch wherever cs_thin or cs_thick jumps.

Usage:
  python scripts/branch_by_gamma.py --case-dir <case dir> \
      [--config config/base.yaml] [--out <csv>]
"""

import argparse
import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np

import params
import scan
from cases import parse_case_rel
from solver import NewtonSolver


def _read_pw(pw_path):
    with open(pw_path, newline="") as fh:
        reader = csv.reader(fh)
        next(reader, None)
        return [(r[0], float(r[1]), float(r[2])) for r in reader if r]


def _solve_branch(cfg, phi1_inf, phi2_inf, mode):
    """Solve one surface branch at a fixed far field. Returns (gamma, cs, ok)."""
    newton = NewtonSolver(cfg)
    newton.p.phi1_inf = float(phi1_inf)
    newton.p.phi2_inf = float(phi2_inf)
    U0 = scan.branch_guess(newton.p, cfg.scan, mode)
    U, ok = newton.solve(U0)
    if not ok:
        return float("nan"), float("nan"), False
    gamma, cs1, cs2 = newton.surface_metrics(U)
    return float(gamma), float(cs1 + cs2), True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--case-dir", required=True)
    ap.add_argument("--config", default=str(ROOT / "config" / "base.yaml"))
    ap.add_argument("--out")
    args = ap.parse_args()

    case_dir = Path(args.case_dir)
    rel = "/".join(case_dir.parts[-3:])
    overrides = parse_case_rel(rel)
    cfg = params.apply_case(params.load_config(args.config), overrides)

    pw = _read_pw(case_dir / "pw_line.csv")
    out_path = Path(args.out) if args.out else case_dir / "branch_gamma.csv"

    rows = []
    for source, phi1_inf, phi2_inf in pw:
        g_thin, cs_thin, ok_thin = _solve_branch(cfg, phi1_inf, phi2_inf, "thin")
        g_thick, cs_thick, ok_thick = _solve_branch(cfg, phi1_inf, phi2_inf, "thick")
        rows.append({
            "source": source,
            "phi1_inf": phi1_inf, "phi2_inf": phi2_inf,
            "gamma_thin": g_thin, "gamma_thick": g_thick,
            "cs_thin": cs_thin, "cs_thick": cs_thick,
            "ok_thin": int(ok_thin), "ok_thick": int(ok_thick),
        })

    fields = ["source", "phi1_inf", "phi2_inf", "gamma_thin", "gamma_thick",
              "cs_thin", "cs_thick", "ok_thin", "ok_thick"]
    with open(out_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

    ok = sum(r["ok_thin"] and r["ok_thick"] for r in rows)
    print(f"{len(rows)} points, both branches converged on {ok}")
    print(f"wrote {out_path}")


if __name__ == "__main__":
    main()
