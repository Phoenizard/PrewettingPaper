"""Enumerate leaf cases in a pw-space data/ tree.

Layout (four levels): data_root / chi.. / om.. / chibb.. / {pw_line.csv, binodal.csv}.
A leaf case is any directory holding pw_line.csv. Enumeration does no parameter
filtering; downstream callers filter the resulting measure table by whatever
parameters their analysis angle needs (omega analysis keeps chibb == 0, etc.).
"""

from collections import namedtuple
from pathlib import Path

from params import parse_case_rel

Case = namedtuple("Case", ["rel", "params", "pw_path", "binodal_path"])


def iter_cases(data_root):
    """Yield one Case per leaf directory, sorted for reproducible order.

    rel is the three-level relative path 'chi../om../chibb..'; params is the
    decoded parameter dict from parse_case_rel; pw_path and binodal_path point
    at the two CSVs. Directories without pw_line.csv are skipped.
    """
    data_root = Path(data_root)
    for pw_path in sorted(data_root.rglob("pw_line.csv")):
        case_dir = pw_path.parent
        rel = case_dir.relative_to(data_root).as_posix()
        params = parse_case_rel(rel)
        yield Case(rel, params, pw_path, case_dir / "binodal.csv")
