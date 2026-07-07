#!/usr/bin/env python
"""At a chosen reservoir phi_inf on the om2=-0.40 line, list EVERY surface state
(exact solve_bvp multi-start) with its wall composition, gamma, and film thickness L,
to test the pre-wetting picture directly:

  - a genuine pre-wetting point should have exactly TWO states (thin, thick) with EQUAL
    gamma (that is the definition of the transition);
  - the disputed low-phi2 region shows THREE states — print all three gammas so we can
    see whether any two are equal (a real transition) or none are (then a line drawn
    there, as in the reference, would NOT be a pre-wetting transition).

Film thickness L = Gibbs adsorption of solute 1 = INT_0^inf (phi1(z) - phi1_inf) dz,
the standard thermodynamic thickness conjugate to gamma. Exact profiles, no ansatz.

  python scripts/diagnostics/pw_check_transition.py [<chi> <om> <chibb>] [phi2_inf ...]
Default: om2=-0.40 at a mid (agreed) point 0.05 and the disputed low end 0.026, 0.023.
Run on the server.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.dirname(HERE)
ROOT = os.path.dirname(SCRIPTS)
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, os.path.join(ROOT, "src"))

import numpy as np  # noqa: E402
import cases  # noqa: E402
import binodal as B  # noqa: E402
import equilibrium as E  # noqa: E402
import verify as V  # noqa: E402

KAPPA = E.T.Kappa(1.0, 1.0)
_trapz = getattr(np, "trapezoid", getattr(np, "trapz", None))


def film_L(prof, phi1_inf):
    """Gibbs adsorption thickness INT (phi1(z) - phi1_inf) dz."""
    return float(_trapz(prof.phi[0] - phi1_inf, prof.z))


def main(argv):
    args = [a for a in argv if not _isfloat(a)]
    floats = [float(a) for a in argv if _isfloat(a)]
    rel = tuple(args[:3]) if len(args) >= 3 else (
        "chi12_0p0__chi13_2p8__chi23_0p0", "om1_m0p30__om2_m0p40",
        "chibb11_0p00__chibb22_0p00__chibb12_0p00")
    chi, surf = cases.parse_case(*rel)
    bd = B.binodal_from_hull(chi)
    f_left, f_right, apex = V._branches(bd)

    # Reservoirs to test: the continuation's OWN transition points (pw_line.csv), NOT the
    # binodal left flank — the pre-wetting point sits inside the flank, and using f_left
    # picks a phi1_inf where no surface states exist. Read (phi1_inf, phi2_inf) rows and
    # test the ones nearest the requested phi2 values.
    dbroot = os.environ.get("RESULTS_DB", os.path.join(ROOT, "database"))
    pw_csv = os.path.join(dbroot, *rel, "pw_line.csv")
    if not os.path.exists(pw_csv):
        pw_csv = os.path.join(ROOT, "out", *rel, "pw_line.csv")
    pw = np.loadtxt(pw_csv, delimiter=",", skiprows=1, ndmin=2)  # (n,2) phi1_inf,phi2_inf
    want = floats if floats else [0.06, 0.05, 0.04, 0.03, 0.026]
    # pick the pw-line row whose phi2_inf is closest to each requested value
    reservoirs = []
    for w in want:
        k = int(np.argmin(np.abs(pw[:, 1] - w)))
        reservoirs.append((float(pw[k, 0]), float(pw[k, 1])))

    print(f"{'/'.join(rel[1:2])}  (exact multi-start states at pw-line transition points)",
          flush=True)
    for phi1_inf, phi2_inf in reservoirs:
        fr = f_right(phi2_inf)
        dense = [(fr, phi2_inf), (0.97 * fr, phi2_inf), (0.9, phi2_inf),
                 (0.7, phi2_inf), (0.5, phi2_inf)]
        st = E.find_states(chi, (phi1_inf, phi2_inf), surf, KAPPA, dense_seeds=dense)
        print(f"\nphi_inf=({phi1_inf:.4f},{phi2_inf:.4f})  ->  {len(st)} state(s):",
              flush=True)
        rows = []
        for s in st:
            L = film_L(s, phi1_inf)
            rows.append((s.phi[0][0], s.phi[1][0], s.gamma, L))
            print(f"   wall=({s.phi[0][0]:.4f},{s.phi[1][0]:.4f})  "
                  f"gamma={s.gamma:+.6f}  L={L:.3f}", flush=True)
        # pairwise gamma differences: which (if any) two states are equal-energy?
        if len(rows) >= 2:
            print("   pairwise |dgamma| (equal-energy pair => pre-wetting):", flush=True)
            for i in range(len(rows)):
                for j in range(i + 1, len(rows)):
                    dg = abs(rows[i][2] - rows[j][2])
                    print(f"     states {i}&{j}: |dgamma|={dg:.6f}", flush=True)
    print("\nDONE", flush=True)


def _isfloat(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


if __name__ == "__main__":
    main(sys.argv[1:])
