#!/usr/bin/env python
"""Trace one case's pre-wetting line by continuation at several phi2 step sizes.
Tells a real merge terminus (finer steps reach no further) from a march that stepped
too coarsely and quit early (finer steps reach further -> real points were missed).
Reuses the production solver (verify._seed_scan/_continue). Run on the server.
  python scripts/diagnostics/pw_probe.py [<chi_dir> <om_dir> <chibb_dir>]   # default: om2=-0.40,-0.42
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.dirname(HERE)
ROOT = os.path.dirname(SCRIPTS)
sys.path.insert(0, SCRIPTS)                       # import verify
sys.path.insert(0, os.path.join(ROOT, "src"))

import numpy as np  # noqa: E402
import cases  # noqa: E402
import binodal as B  # noqa: E402
import verify as V  # noqa: E402


def probe(rel, dphi2_list=(0.005, 0.002, 0.001)):
    chi, surf = cases.parse_case(*rel)
    bd = B.binodal_from_hull(chi)
    f_left, f_right, apex = V._branches(bd)
    seeds = []
    for phi2 in np.arange(0.03, min(0.09, 0.85 * apex), 0.01):
        r = V._seed_scan(chi, surf, float(phi2), f_left, f_right)
        if r:
            seeds.append((float(phi2), r[0], r[1], r[2], r[3]))
    if not seeds:
        print(f"{'/'.join(rel)}: no seed on band")
        return None
    best = max(seeds, key=lambda s: s[4])
    seed = (best[0], best[1], best[2], best[3])
    print(f"{'/'.join(rel[1:2])}  apex_phi2={apex:.4f}  "
          f"seed phi2={best[0]:.3f} phi1*={best[1]:.4f} sep={best[4]:.3f}  "
          f"seeds@phi2={[round(s[0],3) for s in seeds]}")
    print("  dphi2    down_reach  up_reach   n_pts")
    finest = None
    for dp in dphi2_list:
        down = V._continue(chi, surf, seed, apex, -1, dphi2=dp)
        up = V._continue(chi, surf, seed, apex, +1, dphi2=dp)
        line = sorted([(best[1], best[0])] + down + up, key=lambda p: p[1])
        print(f"  {dp:.4f}   {line[0][1]:.4f}      {line[-1][1]:.4f}     {len(line)}")
        finest = line
    return bd, finest


def main(argv):
    if len(argv) >= 3:
        rels = [tuple(argv[:3])]
    else:
        chib = "chibb11_0p00__chibb22_0p00__chibb12_0p00"
        rels = [("chi12_0p0__chi13_2p8__chi23_0p0", "om1_m0p30__om2_m0p40", chib),
                ("chi12_0p0__chi13_2p8__chi23_0p0", "om1_m0p30__om2_m0p42", chib)]
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    for rel in rels:
        res = probe(rel)
        if not res:
            continue
        bd, line = res
        la = np.array(line)
        fig, ax = plt.subplots(figsize=(6, 6))
        ax.scatter(bd.points[:, 0], bd.points[:, 1], s=6, c="#9db4c8", label="binodal")
        ax.plot(la[:, 0], la[:, 1], "-o", ms=3, color="#5a3d99",
                label="PW continuation (finest)")
        ax.set_xlim(0, 0.3)
        ax.set_ylim(0, 0.15)
        ax.set_xlabel("phi1")
        ax.set_ylabel("phi2")
        ax.set_title(rel[1])
        ax.legend(fontsize=8)
        outdir = os.path.join(ROOT, "out", "diagnostics", rel[1])
        os.makedirs(outdir, exist_ok=True)
        png = os.path.join(outdir, "pw_probe.png")
        fig.savefig(png, dpi=120, bbox_inches="tight")
        plt.close(fig)
        print(f"  wrote {png}\n")


if __name__ == "__main__":
    main(sys.argv[1:])
