"""Extract thin / thick surface profiles phi_1(z), phi_2(z) at pre-wetting points.

A pre-wetting point is a far field (phi1_inf, phi2_inf) where the thin-film and
thick-film surface states have equal gamma. Both states therefore exist there,
and their profiles are what shows the shape of each state.

The profiles must come from the same scan path that produced pw_line.csv, seed
continuation included: cold-starting the thick branch from branch_guess is known
to die on T-f (thick guess amplitude 0.3/0.3, see PROGRESS 2026-07-08). So this
runs the full hysteresis line through the point (fixed phi1_inf, sweeping
phi2_inf) with the profile-collection hook, takes the converged profiles at the
scan node nearest the crossing as seeds, and re-solves at the exact far field.

Compute only, no matplotlib (plotting is scripts/plot_tf_profiles.py).

Usage:
  python scripts/tf_profiles.py --case-rel <chi/om/chibb> --pw-line <pw_line.csv> \
      --out-dir <dir> [--config config/base.yaml] [--far-field p1,p2 ...]
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
from logutil import log
from solver import NewtonSolver

POINT_IDS = ["A", "B", "C"]


def _read_pw(pw_path):
    with open(pw_path, newline="") as fh:
        reader = csv.DictReader(fh)
        return [(float(r["phi1_inf"]), float(r["phi2_inf"])) for r in reader if r]


def _pick_far_fields(pw_points, n=3):
    """Start / middle / far end of the pre-wetting line, ordered along phi1_inf."""
    ordered = sorted(pw_points)
    idx = np.linspace(0, len(ordered) - 1, n).round().astype(int)
    return [ordered[i] for i in idx]


def _solve_at(cfg, phi1_inf, phi2_inf, seed):
    """Re-solve one branch at the exact far field, warm-started from seed."""
    newton = NewtonSolver(cfg)
    newton.p.phi1_inf = float(phi1_inf)
    newton.p.phi2_inf = float(phi2_inf)
    U, ok = newton.solve(seed)
    gamma, cs1, cs2 = newton.surface_metrics(U)
    return newton.z, U, ok, float(gamma), float(cs1), float(cs2)


def _nearest_index(values, target, store):
    """Index of the scan node nearest target that has a stored profile."""
    have = sorted(store)
    if not have:
        return None
    return min(have, key=lambda i: abs(values[i] - target))


def run_point(cfg, point_id, phi1_inf, phi2_inf, out_dir):
    label = f"[{point_id} phi1_inf={phi1_inf:.6f}]"
    newton = NewtonSolver(cfg)
    newton.p.phi1_inf = float(phi1_inf)

    store_thin, store_thick = {}, {}
    line = scan.hysteresis_line(
        newton, cfg.scan, "phi2_inf", label=label,
        store_thin=store_thin, store_thick=store_thick,
    )
    values = line["values"]

    with open(out_dir / f"scanline_{point_id}.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["phi2_inf", "gamma_thin", "gamma_thick", "cs_thin", "cs_thick"])
        for i, v in enumerate(values):
            w.writerow([v, line["omega_thin"][i], line["omega_thick"][i],
                        line["cs_thin"][i], line["cs_thick"][i]])

    crossings = scan.extract_crossings(line, cfg.criterion)
    recomputed = min(crossings, key=lambda c: abs(c - phi2_inf)) if crossings else float("nan")
    log(f"{label} {len(crossings)} crossing(s); nearest={recomputed:.6f} "
        f"archive={phi2_inf:.6f}")

    rows, profiles = [], []
    for branch, store in (("thin", store_thin), ("thick", store_thick)):
        i = _nearest_index(values, phi2_inf, store)
        if i is None:
            log(f"{label} {branch}: no converged profile on the whole line")
            rows.append({
                "point_id": point_id, "phi1_inf": phi1_inf, "phi2_inf": phi2_inf,
                "branch": branch, "converged": 0, "gamma": float("nan"),
                "cs1": float("nan"), "cs2": float("nan"),
                "phi1_0": float("nan"), "phi2_0": float("nan"),
                "seed_phi2_inf": float("nan"),
                "pw_phi2_recomputed": recomputed, "pw_phi2_archive": phi2_inf,
            })
            continue
        z, U, ok, gamma, cs1, cs2 = _solve_at(cfg, phi1_inf, phi2_inf, store[i])
        N = len(z)
        phi1, phi2 = U[:N], U[N:]
        rows.append({
            "point_id": point_id, "phi1_inf": phi1_inf, "phi2_inf": phi2_inf,
            "branch": branch, "converged": int(ok), "gamma": gamma,
            "cs1": cs1, "cs2": cs2,
            "phi1_0": float(phi1[0]), "phi2_0": float(phi2[0]),
            "seed_phi2_inf": float(values[i]),
            "pw_phi2_recomputed": recomputed, "pw_phi2_archive": phi2_inf,
        })
        profiles += [(point_id, phi1_inf, phi2_inf, branch, z[k], phi1[k], phi2[k])
                     for k in range(N)]
        log(f"{label} {branch}: ok={ok} gamma={gamma:.8f} cs={cs1 + cs2:.6f} "
            f"phi1(0)={phi1[0]:.6f} phi2(0)={phi2[0]:.6f}")

    return rows, profiles


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--case-rel", required=True)
    ap.add_argument("--pw-line")
    ap.add_argument("--far-field", action="append", default=[],
                    help="'phi1,phi2'; repeatable. Overrides --pw-line picking.")
    ap.add_argument("--config", default=str(ROOT / "config" / "base.yaml"))
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()

    cfg = params.apply_case(params.load_config(args.config),
                            params.parse_case_rel(args.case_rel))
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if args.far_field:
        far_fields = [tuple(float(x) for x in s.split(",")) for s in args.far_field]
    elif args.pw_line:
        far_fields = _pick_far_fields(_read_pw(args.pw_line))
    else:
        ap.error("give --pw-line or --far-field")

    p = cfg.physical
    log(f"case {args.case_rel}: chi=({p.chi_12}, {p.chi_13}, {p.chi_23}) "
        f"omega=({p.omega_1}, {p.omega_2}) "
        f"chibb=({p.chi_bb_11}, {p.chi_bb_22}, {p.chi_bb_12})")

    all_rows, all_profiles = [], []
    for point_id, (phi1_inf, phi2_inf) in zip(POINT_IDS, far_fields):
        rows, profiles = run_point(cfg, point_id, phi1_inf, phi2_inf, out_dir)
        all_rows += rows
        all_profiles += profiles

    fields = ["point_id", "phi1_inf", "phi2_inf", "branch", "converged", "gamma",
              "cs1", "cs2", "phi1_0", "phi2_0", "seed_phi2_inf",
              "pw_phi2_recomputed", "pw_phi2_archive"]
    with open(out_dir / "profile_summary.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(all_rows)

    with open(out_dir / "profiles.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["point_id", "phi1_inf", "phi2_inf", "branch", "z", "phi1", "phi2"])
        w.writerows(all_profiles)

    log(f"wrote {out_dir}/profiles.csv ({len(all_profiles)} rows), "
        f"profile_summary.csv ({len(all_rows)} rows)")


if __name__ == "__main__":
    main()
