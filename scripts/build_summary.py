"""Aggregate out/<chi>/<om>/<chibb>/pw_line.csv into out/SUMMARY.csv.

One row per case: point counts per direction and covered ranges.
"""

import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main():
    verify_root = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "out"
    rows = []
    for pw_path in sorted(verify_root.glob("*/*/*/pw_line.csv")):
        rel = pw_path.parent.relative_to(verify_root)
        with open(pw_path, newline="") as fh:
            reader = csv.reader(fh)
            next(reader)
            pts = [(s, float(a), float(b)) for s, a, b in reader]
        row = {"case": str(rel)}
        for source in ("fix_phi2", "fix_phi1"):
            sub = [(a, b) for s, a, b in pts if s == source]
            row[f"n_{source}"] = len(sub)
            if sub:
                xs, ys = zip(*sub)
                row[f"{source}_phi1_range"] = f"{min(xs):.4f}..{max(xs):.4f}"
                row[f"{source}_phi2_range"] = f"{min(ys):.4f}..{max(ys):.4f}"
            else:
                row[f"{source}_phi1_range"] = ""
                row[f"{source}_phi2_range"] = ""
        rows.append(row)

    if not rows:
        print(f"[summary] no pw_line.csv under {verify_root}", flush=True)
        return
    out_path = verify_root / "SUMMARY.csv"
    with open(out_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"[summary] {len(rows)} case(s) -> {out_path}", flush=True)


if __name__ == "__main__":
    main()
