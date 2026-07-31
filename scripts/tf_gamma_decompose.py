"""Decompose gamma into its four pieces on the saved T-f profiles.

Purpose: test the first TODO in
out/analysis/Tf_omega_cross/omega_existance_analysis.md, which claims the wall
term f_surf is numerically far smaller than the two integral terms, so that the
thin/thick balance can be discussed from the non-wall terms alone.

Reads the already-solved profiles (no re-solving) and splits

    gamma = int W dz + int f_grad dz + f_surf

term by term, mirroring solver.Solver.surface_metrics exactly. Here chibb = 0,
so f_surf = omega_1 phi_1(0) + omega_2 phi_2(0).

The W integral is reported both as a whole and split into the bulk-mixing part
and the chi_12 solute-solute part, since the note discusses those separately.
"""

import csv
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import model  # noqa: E402


@dataclass
class Physical:
    """The subset of params.Physical that model.py reads, same defaults.

    Declared here rather than imported so this script does not pull in the yaml
    config loader; the values match config/base.yaml.
    """

    kappa_11: float = 1.0
    kappa_22: float = 1.0
    kappa_12: float = 0.0
    n1: float = 1.0
    n2: float = 1.0
    n3: float = 1.0
    chi_12: float = 0.0
    chi_13: float = 0.0
    chi_23: float = 0.0
    omega_1: float = 0.0
    omega_2: float = 0.0
    phi1_inf: float = 0.0
    phi2_inf: float = 0.0

PROFILES = ROOT / "out/analysis/Tf_omega_cross/profiles/profiles.csv"
OUT_CSV = ROOT / "out/analysis/Tf_omega_cross/profiles/gamma_terms.csv"

# T-f cross centre: the case the profile figure was made from.
CHI_12, CHI_13, CHI_23 = -8.5, 0.0, 0.0
OMEGA_1, OMEGA_2 = 0.25, -0.375


def load_profiles(path):
    """Group the saved profile rows by (point_id, state), keeping z order."""
    groups = {}
    with open(path) as fh:
        for row in csv.DictReader(fh):
            key = (row["point_id"], row["state"])
            rec = groups.setdefault(
                key,
                {
                    "phi1_inf": float(row["phi1_inf"]),
                    "phi2_inf": float(row["phi2_inf"]),
                    "z": [],
                    "phi1": [],
                    "phi2": [],
                },
            )
            rec["z"].append(float(row["z"]))
            rec["phi1"].append(float(row["phi1"]))
            rec["phi2"].append(float(row["phi2"]))
    for rec in groups.values():
        for k in ("z", "phi1", "phi2"):
            rec[k] = np.asarray(rec[k])
    return groups


def decompose(rec):
    """Split gamma for one profile. Mirrors solver.surface_metrics."""
    p = Physical(
        chi_12=CHI_12,
        chi_13=CHI_13,
        chi_23=CHI_23,
        omega_1=OMEGA_1,
        omega_2=OMEGA_2,
        phi1_inf=rec["phi1_inf"],
        phi2_inf=rec["phi2_inf"],
    )
    z, phi1, phi2 = rec["z"], rec["phi1"], rec["phi2"]

    # Full W, exactly as the solver builds it.
    f = model.free_energy(phi1, phi2, p)
    f_inf = model.free_energy(p.phi1_inf, p.phi2_inf, p)
    mu1_inf, mu2_inf = model.chemical_potential(p.phi1_inf, p.phi2_inf, p)
    W = f - f_inf - mu1_inf * (phi1 - p.phi1_inf) - mu2_inf * (phi2 - p.phi2_inf)

    # Same construction with chi_12 switched off, so the difference isolates
    # the solute-solute contribution to W (tangent-plane subtraction included).
    p0 = Physical(
        chi_12=0.0,
        chi_13=CHI_13,
        chi_23=CHI_23,
        omega_1=OMEGA_1,
        omega_2=OMEGA_2,
        phi1_inf=rec["phi1_inf"],
        phi2_inf=rec["phi2_inf"],
    )
    f0 = model.free_energy(phi1, phi2, p0)
    f0_inf = model.free_energy(p0.phi1_inf, p0.phi2_inf, p0)
    mu1_0, mu2_0 = model.chemical_potential(p0.phi1_inf, p0.phi2_inf, p0)
    W_mix = f0 - f0_inf - mu1_0 * (phi1 - p0.phi1_inf) - mu2_0 * (phi2 - p0.phi2_inf)
    W_chi12 = W - W_mix

    h = z[1] - z[0]
    dphi1 = np.gradient(phi1, h)
    dphi2 = np.gradient(phi2, h)
    f_grad = 0.5 * p.kappa_11 * dphi1**2 + 0.5 * p.kappa_22 * dphi2**2

    # Diagnostic, not a piece of the gamma sum: the chi_12 coupling written on
    # the deviations from the far field rather than on phi itself.
    int_chi12_dev = np.trapezoid(
        p.chi_12 * (phi1 - p.phi1_inf) * (phi2 - p.phi2_inf), x=z
    )

    int_W = np.trapezoid(W, x=z)
    int_W_mix = np.trapezoid(W_mix, x=z)
    int_W_chi12 = np.trapezoid(W_chi12, x=z)
    int_grad = np.trapezoid(f_grad, x=z)
    f_surf = p.omega_1 * phi1[0] + p.omega_2 * phi2[0]

    return {
        "int_W": int_W,
        "int_W_mix": int_W_mix,
        "int_W_chi12": int_W_chi12,
        "int_chi12_dev": int_chi12_dev,
        "int_grad": int_grad,
        "f_surf": f_surf,
        "gamma": int_W + int_grad + f_surf,
    }


def main():
    groups = load_profiles(PROFILES)
    fields = [
        "point_id",
        "state",
        "int_W",
        "int_W_mix",
        "int_W_chi12",
        "int_chi12_dev",
        "int_grad",
        "f_surf",
        "gamma",
    ]
    rows = []
    for (pid, state), rec in sorted(groups.items()):
        terms = decompose(rec)
        rows.append({"point_id": pid, "state": state, **terms})

    with open(OUT_CSV, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        for r in rows:
            writer.writerow({k: (r[k] if isinstance(r[k], str) else f"{r[k]:.6f}") for k in fields})

    head = (
        f"{'pt':<4}{'state':<9}{'int_W':>11}{'W_mix':>11}{'W_chi12':>11}"
        f"{'chi12_dev':>11}{'grad':>11}{'f_surf':>11}{'gamma':>11}"
    )
    print(head)
    print("-" * len(head))
    for r in rows:
        print(
            f"{r['point_id']:<4}{r['state']:<9}{r['int_W']:11.4f}{r['int_W_mix']:11.4f}"
            f"{r['int_W_chi12']:11.4f}{r['int_chi12_dev']:11.4f}{r['int_grad']:11.4f}"
            f"{r['f_surf']:11.4f}{r['gamma']:11.4f}"
        )

    print("\nthin -> thick changes (thick minus thin):")
    print(
        f"{'pt':<4}{'d int_W':>11}{'d chi12_dev':>13}{'d grad':>11}"
        f"{'d f_surf':>11}{'d gamma':>11}"
    )
    by_pt = {}
    for r in rows:
        by_pt.setdefault(r["point_id"], {})[r["state"]] = r
    for pid, states in sorted(by_pt.items()):
        if "state0" in states and "state1" in states:
            a, b = states["state0"], states["state1"]
            print(
                f"{pid:<4}{b['int_W'] - a['int_W']:11.4f}"
                f"{b['int_chi12_dev'] - a['int_chi12_dev']:13.4f}"
                f"{b['int_grad'] - a['int_grad']:11.4f}"
                f"{b['f_surf'] - a['f_surf']:11.4f}{b['gamma'] - a['gamma']:11.4f}"
            )

    print(f"\nwrote {OUT_CSV.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
