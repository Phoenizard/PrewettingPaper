#!/usr/bin/env python
"""Free-energy-vs-film-thickness curve gamma(L) at a fixed reservoir, the standard
pre-wetting picture: a double (here possibly triple) well whose local minima ARE the
thin / (middle) / thick surface states.

Method (project note ternary.md sections 7-8, no ODE solve):
  First integral (note 7):        sum_i (kappa_i/2)(phi_i')^2 = W(phi)
  => along the interfacial path from the reservoir phi_inf out to a wall composition
     phi0, the excess free energy and the physical film thickness are pure path integrals
       gamma(phi0) = f_surf(phi0) + INT_{phi_inf}^{phi0} sqrt(2 kappa_eff W) ds     (note 8.1)
       L(phi0)     =               INT_{phi_inf}^{phi0}   ds / sqrt(2 W / kappa_eff)
  ds is arclength along the composition path; we use the linear-path ansatz (note 8.2)
  phi(t) = phi_inf + t (phi0 - phi_inf), t in [0,1], so kappa_eff is kappa projected on
  the fixed path direction. Sweeping the wall endpoint phi0 along the dense (thick) flank
  traces gamma vs L; the minima of gamma(L) are the surface states, and pre-wetting is the
  first-order transition between two equally deep minima.

We sweep phi0 along a ray from the reservoir toward the phi1-rich (wall-favoured) corner,
parameterised by the wall phi1 value p10 in [phi1_inf, p10_max]; phi2 at the wall follows
the same ray. For each p10 we integrate along the linear path and record (L, gamma, p10).

  python scripts/diagnostics/gamma_vs_L.py [<chi_dir> <om_dir> <chibb_dir>] [phi2_inf]
Default: om2=-0.40, phi2_inf swept at a few values across the 3-state region.
Run on the server (cheap, but keep it on the compute box per project rule).
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
import thermo as T  # noqa: E402
import verify as V  # noqa: E402

KAPPA = T.Kappa(1.0, 1.0)


def gamma_and_L(chi, surf, res, phi0, npath=4000):
    """Integrate gamma and L along the LINEAR path phi_inf -> phi0. Returns (L, gamma) or
    (None, None) if the path dips into W<0 (path leaves the metastable well; not a valid
    single-interface profile)."""
    res = np.asarray(res, float)
    phi0 = np.asarray(phi0, float)
    res_mu = T.mu(res[0], res[1], chi)
    res_fb = T.f_b(res[0], res[1], chi)
    t = np.linspace(0.0, 1.0, npath)
    path = res[None, :] + t[:, None] * (phi0 - res)[None, :]         # (npath, 2)
    d = phi0 - res                                                    # path direction
    kappa_eff = KAPPA.k1 * d[0] ** 2 + KAPPA.k2 * d[1] ** 2           # kappa along the ray
    seglen = np.hypot(d[0], d[1])                                     # |d|; ds = seglen*dt
    Wv = np.array([T.W(p[0], p[1], chi, res, res_mu=res_mu, res_fb=res_fb) for p in path])
    if np.any(Wv[1:] < -1e-9):        # allow tiny negative at the reservoir end only
        return None, None
    Wv = np.clip(Wv, 0.0, None)
    # gamma = f_surf + INT sqrt(2 kappa_eff W) ds ; but kappa_eff above already has |d|^2,
    # and ds = seglen dt, so use the per-component form consistent with note 8.1:
    integ_g = np.sqrt(2.0 * kappa_eff * Wv)                           # d(gamma)/dt
    # L = INT ds / sqrt(2 W / kappa_eff_unit); with the linear ansatz the arclength speed
    # is seglen, and the local inverse-gradient is sqrt(kappa_eff)/sqrt(2 W) per unit t:
    with np.errstate(divide="ignore"):
        integ_L = np.where(Wv > 1e-12, np.sqrt(kappa_eff) / np.sqrt(2.0 * Wv), 0.0)
    gamma = T.f_surf(phi0[0], phi0[1], surf) + np.trapz(integ_g, t)
    L = np.trapz(integ_L, t)
    return float(L), float(gamma)


def curve(chi, surf, res, p10_max=0.97, n=140):
    """Sweep the wall endpoint p10 from just above phi1_inf to p10_max along a ray whose
    phi2 tracks the dense flank; return arrays (p10, L, gamma)."""
    r1, r2 = float(res[0]), float(res[1])
    # ray in phi2: from reservoir phi2 up toward a modest wall phi2 (dense corner). Keep it
    # simple: hold phi2 wall = phi2_inf (phi1-rich wall) — matches T-a dense_seeds style.
    out = []
    for p10 in np.linspace(r1 + 0.02, p10_max, n):
        p20 = r2  # phi1-rich wall ray; phi2 stays near reservoir
        L, g = gamma_and_L(chi, surf, res, (p10, p20))
        if L is not None:
            out.append((p10, L, g))
    return np.array(out) if out else np.empty((0, 3))


def main(argv):
    args = [a for a in argv if not _isfloat(a)]
    floats = [float(a) for a in argv if _isfloat(a)]
    rel = tuple(args[:3]) if len(args) >= 3 else (
        "chi12_0p0__chi13_2p8__chi23_0p0", "om1_m0p30__om2_m0p40",
        "chibb11_0p00__chibb22_0p00__chibb12_0p00")
    chi, surf = cases.parse_case(*rel)
    bd = B.binodal_from_hull(chi)
    f_left, f_right, apex = V._branches(bd)
    # phi2_inf values to probe: default a few spanning the 3-state low-end region.
    phi2s = floats if floats else [0.05, 0.035, 0.026, 0.02]

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7, 6))
    for phi2_inf in phi2s:
        phi1_inf = float(f_left(phi2_inf))
        arr = curve(chi, surf, (phi1_inf, phi2_inf))
        if not len(arr):
            print(f"phi2_inf={phi2_inf:.3f}: empty", flush=True)
            continue
        p10, L, g = arr[:, 0], arr[:, 1], arr[:, 2]
        # count local minima of gamma(L)
        mins = [i for i in range(1, len(g) - 1) if g[i] < g[i-1] and g[i] < g[i+1]]
        print(f"phi2_inf={phi2_inf:.3f}  phi1_inf={phi1_inf:.4f}  n_pts={len(g)}  "
              f"local_minima={len(mins)} at L={[round(float(L[i]),2) for i in mins]}",
              flush=True)
        ax.plot(L, g, "-", label=f"phi2_inf={phi2_inf:.3f} ({len(mins)} min)")
        ax.plot(L[mins], g[mins], "o", ms=6)
    ax.set_xlabel("film thickness L")
    ax.set_ylabel("gamma")
    ax.set_title(f"{rel[1]}  gamma vs L (path-integral)")
    ax.legend(fontsize=8)
    outdir = os.path.join(ROOT, "tmp", "gamma_vs_L")
    os.makedirs(outdir, exist_ok=True)
    png = os.path.join(outdir, f"{rel[1]}.png")
    fig.savefig(png, dpi=130, bbox_inches="tight")
    print(f"wrote {png}", flush=True)
    print("DONE", flush=True)


def _isfloat(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


if __name__ == "__main__":
    main(sys.argv[1:])
