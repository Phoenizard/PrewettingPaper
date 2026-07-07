#!/usr/bin/env python
"""Compute the bulk binodal for the 6 project topologies and save phase maps.

Run directly (see CLAUDE.md — everything runs in the numenv conda env):

    conda run -n numenv python scripts/binodal_check.py

This is a plain driver, not a package: it puts src/ on sys.path and imports the
modules directly (no prewet. prefix), then reuses binodal.binodal_from_hull and
plotting.plot_phase_map. Topologies T-a..T-f follow doc/note/project_plan.md
section 2.1; Chi is (chi12, chi1s, chi2s) with chi1s=chi_{1s}, chi2s=chi_{2s}.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "src"))

import binodal  # noqa: E402
import thermo  # noqa: E402
import plotting  # noqa: E402

# (chi12, chi1s, chi2s) per project_plan.md section 2.1
TOPOLOGIES = {
    "T-a": thermo.Chi(chi12=0.0, chi1s=2.8, chi2s=0.0),
    "T-b": thermo.Chi(chi12=0.0, chi1s=2.8, chi2s=2.6),
    "T-c": thermo.Chi(chi12=2.3, chi1s=2.8, chi2s=2.6),
    "T-d": thermo.Chi(chi12=2.6, chi1s=2.6, chi2s=2.6),
    "T-e": thermo.Chi(chi12=2.8, chi1s=2.8, chi2s=2.8),
    "T-f": thermo.Chi(chi12=-8.5, chi1s=0.0, chi2s=0.0),
}

OUT_DIR = os.path.join(ROOT, "tmp", "binodal_check")


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    for name, chi in TOPOLOGIES.items():
        res = binodal.binodal_from_hull(chi)
        path = os.path.join(OUT_DIR, f"binodal_{name}.png")
        plotting.plot_phase_map(
            path,
            res,
            title=f"{name}: chi12={chi.chi12}, chi1s={chi.chi1s}, chi2s={chi.chi2s}",
        )
        print(
            f"{name}: {len(res.tie_lines)} tie-lines, "
            f"{len(res.three_phase)} three-phase region(s), "
            f"max equal-mu residual={res.max_equal_mu_residual:.2e} -> {path}"
        )


if __name__ == "__main__":
    main()
