"""Run one verification case (compute only).

Usage:
  python scripts/run_case.py --case-rel <chi_dir>/<om_dir>/<chibb_dir> \
      [--config config/base.yaml] [--out-root out/verify] [--max-lines N]

--max-lines truncates the number of fixed-value scan lines per direction
(dry-run use).
"""

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import params
import pipeline


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--case-rel", required=True)
    ap.add_argument("--config", default=str(ROOT / "config" / "base.yaml"))
    ap.add_argument("--out-root", default=str(ROOT / "out" / "verify"))
    ap.add_argument("--max-lines", type=int, default=None)
    args = ap.parse_args()

    cfg = params.load_config(args.config)
    cfg = params.apply_case(cfg, params.parse_case_rel(args.case_rel))
    print(f"[case] {args.case_rel}", flush=True)
    print(f"[params] {cfg.physical}", flush=True)
    pipeline.run_case(cfg, Path(args.out_root) / args.case_rel, max_lines=args.max_lines)


if __name__ == "__main__":
    main()
