"""Enumerate leaf cases in a pw-space data/ tree.

Layout (four levels): data_root / chi.. / om.. / chibb.. / {pw_line.csv, binodal.csv}.
A leaf case is any directory holding pw_line.csv. Enumeration does no parameter
filtering; downstream callers filter the resulting measure table by whatever
parameters their analysis angle needs (omega analysis keeps chibb == 0, etc.).
"""

from collections import namedtuple
from pathlib import Path

Case = namedtuple("Case", ["rel", "params", "pw_path", "binodal_path"])

# Directory-name token -> parameter key. Encoding: m = minus, p = decimal point.
# Kept local so this measure layer needs no config loader (no yaml dependency);
# mirrors params.parse_case_rel.
_CASE_KEYS = {
    "chi12": "chi_12", "chi13": "chi_13", "chi23": "chi_23",
    "om1": "omega_1", "om2": "omega_2",
    "chibb11": "chi_bb_11", "chibb22": "chi_bb_22", "chibb12": "chi_bb_12",
}


def parse_case_rel(rel):
    """'chi../om../chibb..' three-level path -> dict of decoded parameters."""
    parts = Path(rel).parts
    if len(parts) != 3:
        raise ValueError(f"case rel must have 3 levels (chi/om/chibb): {rel}")
    params = {}
    for part in parts:
        for token in part.split("__"):
            key, _, value = token.partition("_")
            if key not in _CASE_KEYS:
                raise ValueError(f"unknown token {token!r} in {rel}")
            params[_CASE_KEYS[key]] = float(value.replace("m", "-").replace("p", "."))
    return params


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
