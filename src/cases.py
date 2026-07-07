"""Case naming <-> parameters, and result/ <-> out/ path mapping.

We adopt result/'s directory naming verbatim (it is filesystem-safe, self-
describing, sortable, and 1:1 with the parameters), so every verification output
maps 1:1 to the result/ case it checks. Layout mirrored under out/:

    result/ <chi_dir>/<om_dir>/<chibb_dir>/overlay_omega1_omega2.png  (theirs)
    out/    <chi_dir>/<om_dir>/<chibb_dir>/overlay.png , pw_line.csv  (ours)

Encoding (as in result/): value -> digits with '.'->'p', negative prefixed 'm'.
  chi_dir   : chi12_0p0__chi13_2p8__chi23_0p0     (chi13=chi1s, chi23=chi2s)
  om_dir    : om1_m0p30__om2_m0p30                (om_i = omega_{b,i})
  chibb_dir : chibb11_0p00__chibb22_0p00__chibb12_0p00
"""
from __future__ import annotations

import os

import thermo as T

RESULT_ROOT = "result"
VERIFY_ROOT = "out"  # our results live directly under out/<chi>/<om>/<chibb>/ (no stage layer)


def decode(tok: str) -> float:
    """'0p0'->0.0, '2p8'->2.8, 'm0p30'->-0.30, 'm8p5'->-8.5."""
    neg = tok.startswith("m")
    if neg:
        tok = tok[1:]
    v = float(tok.replace("p", "."))
    return -v if neg else v


def _fields(dirname: str) -> dict:
    """'chi12_0p0__chi13_2p8' -> {'chi12':0.0,'chi13':2.8} (last '_' splits key/val)."""
    out = {}
    for part in dirname.split("__"):
        key, _, val = part.rpartition("_")
        out[key] = decode(val)
    return out


def parse_case(chi_dir: str, om_dir: str, chibb_dir: str):
    """Return (Chi, Surf) for a result/ leaf. Kappa is fixed (1,1) per result/."""
    c = _fields(chi_dir)
    o = _fields(om_dir)
    b = _fields(chibb_dir)
    chi = T.Chi(chi12=c["chi12"], chi1s=c["chi13"], chi2s=c["chi23"])
    surf = T.Surf(w1=o["om1"], w2=o["om2"],
                  cbb1=b["chibb11"], cbb2=b["chibb22"], cbb12=b["chibb12"])
    return chi, surf


def iter_cases(root: str = RESULT_ROOT):
    """Yield ((chi_dir, om_dir, chibb_dir), Chi, Surf) for every result/ leaf, sorted."""
    for chi_dir in sorted(os.listdir(root)):
        cp = os.path.join(root, chi_dir)
        if not os.path.isdir(cp):
            continue
        for om_dir in sorted(os.listdir(cp)):
            op = os.path.join(cp, om_dir)
            if not os.path.isdir(op):
                continue
            for chibb_dir in sorted(os.listdir(op)):
                if not os.path.isdir(os.path.join(op, chibb_dir)):
                    continue
                chi, surf = parse_case(chi_dir, om_dir, chibb_dir)
                yield (chi_dir, om_dir, chibb_dir), chi, surf


def verify_dir(rel: tuple, root: str = ".") -> str:
    """out/<chi_dir>/<om_dir>/<chibb_dir> for a case key `rel`."""
    return os.path.join(root, VERIFY_ROOT, *rel)


def result_overlay(rel: tuple, root: str = ".") -> str:
    """Path to the matching result/ overlay PNG (theirs) for visual comparison."""
    return os.path.join(root, RESULT_ROOT, *rel, "overlay_omega1_omega2.png")
