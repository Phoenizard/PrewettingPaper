"""Parameter dataclasses and the single-yaml config loader.

Every physical / numerical / criterion parameter lives in config/base.yaml.
A case (chi topology x omega x chibb) overrides only physical fields, parsed
from an unvalidate_data-style relative path (encoding: m = minus, p = point).
"""

from dataclasses import dataclass, replace
from pathlib import Path

import yaml


@dataclass
class Physical:
    L: float = 10.0
    N: int = 1000
    kappa_11: float = 1.0
    kappa_22: float = 1.0
    kappa_12: float = 0.0
    n1: float = 1.0
    n2: float = 1.0
    n3: float = 1.0
    chi_12: float = 0.0
    chi_13: float = 0.0
    chi_23: float = 0.0
    omega_1: float = 0.0
    omega_2: float = 0.0
    chi_bb_11: float = 0.0
    chi_bb_22: float = 0.0
    chi_bb_12: float = 0.0
    phi1_inf: float = 0.0
    phi2_inf: float = 0.0


@dataclass
class Solver:
    tol: float = 1e-8
    max_iter: int = 500
    min_alpha: float = 1e-3


@dataclass
class Scan:
    phi_min: float = 1e-4
    phi_max: float = 0.2
    n_scan: int = 400
    n_lines: int = 400
    amp_phi1_thin: float = 0.015
    amp_phi2_thin: float = 0.004
    amp_phi1_thick: float = 0.3
    amp_phi2_thick: float = 0.3


@dataclass
class Criterion:
    cs_threshold: float = 0.1
    min_hysteresis_points: int = 2


@dataclass
class Bulk:
    demixed_grid: int = 451
    demix_tol: float = 1e-5
    binodal_grid: int = 801
    k_nearest: int = 8
    vote_ratio: float = 0.5


@dataclass
class Config:
    physical: Physical
    solver: Solver
    scan: Scan
    criterion: Criterion
    bulk: Bulk


_SECTIONS = {
    "physical": Physical,
    "solver": Solver,
    "scan": Scan,
    "criterion": Criterion,
    "bulk": Bulk,
}


def load_config(path):
    raw = yaml.safe_load(Path(path).read_text()) or {}
    kwargs = {}
    for name, cls in _SECTIONS.items():
        section = raw.get(name) or {}
        unknown = set(section) - set(cls.__dataclass_fields__)
        if unknown:
            raise KeyError(f"unknown keys in [{name}]: {sorted(unknown)}")
        kwargs[name] = cls(**section)
    return Config(**kwargs)


_CASE_KEYS = {
    "chi12": "chi_12",
    "chi13": "chi_13",
    "chi23": "chi_23",
    "om1": "omega_1",
    "om2": "omega_2",
    "chibb11": "chi_bb_11",
    "chibb22": "chi_bb_22",
    "chibb12": "chi_bb_12",
}


def _decode_value(text):
    return float(text.replace("m", "-").replace("p", "."))


def parse_case_rel(rel):
    """'chi12_0p0__chi13_2p8__chi23_0p0/om1_m0p42__om2_m0p3/chibb11_0p0__...'
    -> dict of Physical overrides."""
    parts = Path(rel).parts
    if len(parts) != 3:
        raise ValueError(f"case rel must have 3 levels (chi/om/chibb): {rel}")
    overrides = {}
    for part in parts:
        for token in part.split("__"):
            key, _, value = token.partition("_")
            if key not in _CASE_KEYS:
                raise ValueError(f"unknown token {token!r} in {rel}")
            overrides[_CASE_KEYS[key]] = _decode_value(value)
    return overrides


def apply_case(cfg, overrides):
    return Config(
        physical=replace(cfg.physical, **overrides),
        solver=cfg.solver,
        scan=cfg.scan,
        criterion=cfg.criterion,
        bulk=cfg.bulk,
    )
