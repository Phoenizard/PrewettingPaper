#!/usr/bin/env python
"""On-line pre-wetting test for a single reservoir (phi1_inf, phi2_inf).

Given a reservoir, decide whether it sits ON the pre-wetting line: i.e. whether the
thinnest and thickest surface states have equal excess free energy gamma (the defining
condition of the first-order thin<->thick transition).

Definition (on-line test):
    is_prewetting = (>= 2 distinct surface states) AND
                    |gamma(thinnest) - gamma(thickest)| < tol

We compare ONLY the thinnest (st[0]) and thickest (st[-1]) states, never any middle
state. At low phi2 a THIRD (middle) surface state appears (wall phi1 ≈ 0.34 / 0.68 /
0.92); it is a higher-gamma barrier between the two equal-depth wells and does not take
part in the transition. Comparing thin vs thick — never the middle — is exactly the
correctness point diagnosed for the continuation solver.

The surface states come from equilibrium.find_states (exact solve_bvp multi-start),
sorted by wall phi1 ascending, so st[0] is thinnest and st[-1] is thickest. The dense
multi-start seeds fan across the phi1-rich flank so all three states are found.

Film thickness reported per state is the Gibbs adsorption L = INT (phi1(z)-phi1_inf) dz,
the standard thermodynamic thickness conjugate to gamma.

  python scripts/is_prewetting.py <chi_dir> <om_dir> <chibb_dir> <phi1_inf> <phi2_inf>

Exit code 0 if the reservoir is on the pre-wetting line, 1 otherwise (so it composes in
shell pipelines). Run on the server per project rule (exact BVP solves).
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = HERE
ROOT = os.path.dirname(SCRIPTS)
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, os.path.join(ROOT, "src"))

import numpy as np  # noqa: E402
import thermo as T  # noqa: E402
import binodal as B  # noqa: E402
import equilibrium as E  # noqa: E402
import cases  # noqa: E402
import verify as V  # noqa: E402

KAPPA = T.Kappa(1.0, 1.0)
_trapz = getattr(np, "trapezoid", getattr(np, "trapz", None))  # numpy>=2 renamed trapz

# gamma tolerance for the equal-depth condition. 3e-3 sits well above the _gamma simpson
# noise (tol 1e-8 solves) and far below the ~0.05 gamma gap up to the middle-state
# barrier; measured transition points have |dgamma| in 0..3.4e-4. Stored pw_line.csv
# points converge on phi1 (xtol=1e-6), not gamma, so their |dgamma| may reach a few e-3;
# loosen via the tol argument if a known on-line point reports just over the bar.
DEFAULT_TOL = 3e-3


def film_L(prof, phi1_inf):
    """Gibbs adsorption thickness INT (phi1(z) - phi1_inf) dz."""
    return float(_trapz(prof.phi[0] - phi1_inf, prof.z))


def _dense_seeds(f_right, phi2_inf):
    """Multi-start targets across the phi1-rich flank (thin, middle, thick basins)."""
    fr = f_right(phi2_inf)
    return [(fr, phi2_inf), (0.97 * fr, phi2_inf), (0.9, phi2_inf),
            (0.7, phi2_inf), (0.5, phi2_inf)]


def is_prewetting(chi, surf, phi1_inf, phi2_inf, f_right=None, tol=DEFAULT_TOL,
                  distinct_tol=5e-3, L=12.0, n=300, dense_seeds=None):
    """Return (on_line: bool, info: dict) for reservoir (phi1_inf, phi2_inf).

    on_line is True iff there are >= 2 distinct surface states and the thinnest and
    thickest have equal gamma within `tol`. `info` reports every state so a caller can
    also inspect the (weaker) "structure supports pre-wetting" question (n_states >= 2).

    `f_right` is the binodal dense-flank interpolator phi1(phi2); if None it is computed
    from the binodal. Pass it (or `dense_seeds`) to avoid recomputing the binodal per call.
    """
    if dense_seeds is None:
        if f_right is None:
            _, f_right, _ = V._branches(B.binodal_from_hull(chi))
        dense_seeds = _dense_seeds(f_right, phi2_inf)
    st = E.find_states(chi, (phi1_inf, phi2_inf), surf, KAPPA,
                       dense_seeds=dense_seeds, L=L, n=n, distinct_tol=distinct_tol)
    states = [(float(s.phi[0][0]), float(s.phi[1][0]), float(s.gamma),
               film_L(s, phi1_inf)) for s in st]
    dgamma = abs(states[0][2] - states[-1][2]) if len(states) >= 2 else None
    on_line = len(states) >= 2 and dgamma < tol
    info = {
        "phi1_inf": float(phi1_inf),
        "phi2_inf": float(phi2_inf),
        "n_states": len(states),
        "states": states,                       # (wall_phi1, wall_phi2, gamma, L) ascending
        "dgamma_thin_thick": dgamma,            # |gamma(st[0]) - gamma(st[-1])|
        "has_middle": len(states) >= 3,
        "tol": tol,
        "is_online": on_line,
    }
    return on_line, info


def _print(info):
    print(f"reservoir=({info['phi1_inf']:.4f},{info['phi2_inf']:.4f})  "
          f"n_states={info['n_states']}"
          + ("  [middle state present]" if info["has_middle"] else ""))
    for i, (p1, p2, g, L) in enumerate(info["states"]):
        role = ""
        if info["n_states"] >= 2:
            role = " (thin)" if i == 0 else " (thick)" if i == info["n_states"] - 1 \
                else " (middle)"
        print(f"  state{i}{role}: wall=({p1:.4f},{p2:.4f})  gamma={g:+.6f}  L={L:.3f}")
    dg = info["dgamma_thin_thick"]
    if dg is not None:
        print(f"  |dgamma(thin,thick)|={dg:.6f}  tol={info['tol']:.1e}")
    print(f"  PRE-WETTING (on line): {info['is_online']}")


def main(argv):
    if len(argv) < 5:
        print(__doc__)
        return 2
    chi_dir, om_dir, chibb_dir = argv[:3]
    phi1_inf, phi2_inf = float(argv[3]), float(argv[4])
    chi, surf = cases.parse_case(chi_dir, om_dir, chibb_dir)
    on_line, info = is_prewetting(chi, surf, phi1_inf, phi2_inf)
    _print(info)
    return 0 if on_line else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
