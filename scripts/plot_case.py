"""Render overlay.png for one case directory (plot step only).

Usage: python scripts/plot_case.py --case-dir out/verify/<chi>/<om>/<chibb>
"""

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import plotting


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--case-dir", required=True)
    ap.add_argument("--title", default=None)
    args = ap.parse_args()

    case_dir = Path(args.case_dir)
    title = args.title or "/".join(case_dir.parts[-3:])
    out = plotting.render_overlay(case_dir, title=title)
    print(f"[plot done] {out}", flush=True)


if __name__ == "__main__":
    main()
