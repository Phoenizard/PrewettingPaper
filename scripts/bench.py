#!/usr/bin/env python
"""Benchmark one prewetting case: wall time + solve_bvp call count, per layer.

Toggles equilibrium.USE_ANALYTIC_JAC / USE_WARM_START so each optimization layer
can be measured on its own. Prints only; writes nothing.

Usage:
  conda run -n numenv python scripts/bench.py                         # T-a template, both layers on
  conda run -n numenv python scripts/bench.py <chi_dir> <om_dir> <chibb_dir>
  --no-jac    disable analytic Jacobian
  --no-warm   disable warm-start
"""
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, HERE)

import equilibrium as E  # noqa: E402
import binodal as B  # noqa: E402
import cases  # noqa: E402
import verify  # noqa: E402  (reuse prewetting_line; import does not run main)

DEFAULT = ("chi12_0p0__chi13_2p8__chi23_0p0",
           "om1_m0p30__om2_m0p30",
           "chibb11_0p00__chibb22_0p00__chibb12_0p00")


def main(argv):
    E.USE_ANALYTIC_JAC = "--no-jac" not in argv
    E.USE_WARM_START = "--no-warm" not in argv
    rel = [a for a in argv if not a.startswith("--")]
    rel = tuple(rel) if len(rel) == 3 else DEFAULT

    # count solve_bvp calls by wrapping the name in equilibrium's namespace
    orig = E.solve_bvp
    n_calls = [0]

    def counting(*a, **k):
        n_calls[0] += 1
        return orig(*a, **k)

    E.solve_bvp = counting

    chi, surf = cases.parse_case(*rel)
    binodal = B.binodal_from_hull(chi)
    t0 = time.perf_counter()
    pw = verify.prewetting_line(chi, surf, binodal)
    dt = time.perf_counter() - t0

    print(f"case             {'/'.join(rel)}")
    print(f"USE_ANALYTIC_JAC {E.USE_ANALYTIC_JAC}")
    print(f"USE_WARM_START   {E.USE_WARM_START}")
    print(f"wall_seconds     {dt:.1f}")
    print(f"solve_bvp_calls  {n_calls[0]}")
    print(f"pw_points        {len(pw)}")


if __name__ == "__main__":
    main(sys.argv[1:])
