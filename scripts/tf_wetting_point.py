"""Locate the wetting transition point on the T-f binodal for one omega setting.

The wetting transition point is a property of the binodal and the wall alone --
it does not use the pre-wetting line, which is kept aside as an independent check.

On a tie line the two coexisting phases are R (the far field / reservoir) and P
(the conjugate phase that could coat the wall). Three quantities per tie line:

  gamma(R)  wall surface free energy with far field R, the solution that stays
            near R (no coating)
  gamma(P)  the same with far field P
  sigma     free planar interfacial free energy between R and P, no wall

The wall is coated once gamma(R) = gamma(P) + sigma; walking the closed binodal
loop, the wetting transition point is the sign change of

  Delta = gamma(R) - gamma(P) - sigma.

gamma comes from the production solver (src/solver.py, NewtonSolver): its
surface_metrics returns exactly this quantity. Every loop point is a reservoir
once and a conjugate once, so one wall solve per loop point suffices.

sigma uses the planar-interface first integral. With kappa_11 = kappa_22 = 1 and
kappa_12 = 0 the profile obeys W = (phi1'^2 + phi2'^2) / 2, hence

  sigma = integral of sqrt(2 W) along arc length in the (phi1, phi2) plane,

minimised over paths joining the two coexisting compositions. Writing it this way
removes the position of the interface along z from the problem, so there is no
near-singular translation mode to fight. The path is parametrised as a graph over
the tie line (offset v transverse to it, endpoints pinned at v = 0), which also
removes the reparametrisation freedom; the objective and its gradient are both
analytic.

Reads binodal_tie.csv (the gapless tie-line binodal written by
scripts/tf_critical_point.py) and pw_line.csv (archived pre-wetting points, used
only for the checks at the end). Writes wetting_scan.csv and wetting_points.csv.

Usage:
  python scripts/tf_wetting_point.py [--om1 0.28] [--om2 -0.375]
      [--tie-dir DIR] [--field-dir DIR] [--out-dir DIR]
      [--max-lines N] [--path-nodes 401]
"""

import argparse
import csv
import sys
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import model                      # noqa: E402
from scan import branch_guess     # noqa: E402
from solver import NewtonSolver   # noqa: E402

CHI_12, CHI_13, CHI_23 = -8.5, 0.0, 0.0


# Parameters are declared here rather than read from config/base.yaml, because
# src/params.py needs pyyaml and numenv does not have it. Values are the ones in
# config/base.yaml; chi and omega are the T-f case. Same arrangement as
# scripts/tf_critical_point.py.
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
    chi_12: float = CHI_12
    chi_13: float = CHI_13
    chi_23: float = CHI_23
    omega_1: float = 0.28
    omega_2: float = -0.375
    chi_bb_11: float = 0.0
    chi_bb_22: float = 0.0
    chi_bb_12: float = 0.0
    phi1_inf: float = 0.0
    phi2_inf: float = 0.0


@dataclass
class SolverCfg:
    tol: float = 1e-8
    max_iter: int = 500
    min_alpha: float = 1e-3


@dataclass
class ScanCfg:
    amp_phi1_thin: float = 0.015
    amp_phi2_thin: float = 0.004
    amp_phi1_thick: float = 0.3
    amp_phi2_thick: float = 0.3


@dataclass
class Config:
    physical: Physical
    solver: SolverCfg
    scan: ScanCfg

COATED_FRAC = 0.3     # profile counted as coated if it reaches within this
                      # fraction of the tie-line length of the conjugate phase
V_BOUND_FRAC = 0.4    # transverse offset allowed, as a fraction of tie length
FINE_NODES = 20001    # nodes for the straight-path cross-check quadrature
N_TANGENT_PW = 15     # pre-wetting points used for the local direction fit

_LOG_FH = None


def log(msg):
    line = str(msg)
    print(line, flush=True)
    if _LOG_FH is not None:
        _LOG_FH.write(line + "\n")
        _LOG_FH.flush()


# --------------------------------------------------------------------------
# input


def read_rows(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def load_tie_lines(path):
    branches = {}
    for r in read_rows(path):
        branches.setdefault(int(r["branch"]), []).append(
            [float(r["phi1_a"]), float(r["phi2_a"]),
             float(r["phi1_b"]), float(r["phi2_b"]), float(r["length"])])
    return {k: np.array(v) for k, v in sorted(branches.items())}


def build_loop(branches, max_lines=None):
    """Ordered loop points, each carrying its coexisting partner.

    Same stitching as plot_tf_critical.closed_loop: branch 0 a-ends out, branch 0
    b-ends back, branch 1 b-ends out, branch 1 a-ends back. Every tie line
    therefore contributes both of its endpoints exactly once.
    """
    b0, b1 = branches[0], branches[1]
    if max_lines is not None:
        b0, b1 = b0[:max_lines], b1[:max_lines]

    entries = []

    def add(row, branch, row_idx, flip):
        a = (row[0], row[1])
        b = (row[2], row[3])
        pt, conj = (b, a) if flip else (a, b)
        entries.append({"phi1": pt[0], "phi2": pt[1],
                        "phi1_conj": conj[0], "phi2_conj": conj[1],
                        "branch": branch, "row": row_idx, "length": row[4]})

    for i, row in enumerate(b0):
        add(row, 0, i, flip=False)
    for i in range(len(b0) - 1, -1, -1):
        add(b0[i], 0, i, flip=True)
    for i, row in enumerate(b1):
        add(row, 1, i, flip=True)
    for i in range(len(b1) - 1, -1, -1):
        add(b1[i], 1, i, flip=False)
    return entries


# --------------------------------------------------------------------------
# wall surface free energy


def build_cfg(om1, om2):
    return Config(physical=Physical(omega_1=float(om1), omega_2=float(om2)),
                  solver=SolverCfg(), scan=ScanCfg())


def _sign(omega):
    return 1.0 if float(omega) < 0.0 else -1.0


def _feasible(phi1, phi2, eps=1e-12):
    phi1 = np.clip(phi1, eps, 1.0 - 2 * eps)
    phi2 = np.clip(phi2, eps, 1.0 - 2 * eps - phi1)
    return np.concatenate([phi1, phi2])


def guess_bank(p, z, conj):
    """Initial profiles spanning thin, intermediate and coated surface states.

    This stage is known to carry three surface states at low phi2 (thin,
    intermediate, thick), so a single small-amplitude seed can land on the
    intermediate one and report it as the thin branch. The bank mirrors the one
    in scripts/tf_profile_ref.py (linear and exponential seeds, both the
    wall-sign and the both-up variants) and adds seeds that start from the
    coexisting phase itself, which is what a coated wall looks like.
    """
    layer = 1.0 - z / float(p.L)
    s1, s2 = _sign(p.omega_1), _sign(p.omega_2)
    out = {}

    for tag, a1, a2 in [("lin_thin", 0.015, 0.004),
                        ("lin_0p09", 0.09, 0.02),
                        ("lin_0p3", 0.3, 0.3)]:
        out[tag] = _feasible(p.phi1_inf + s1 * a1 * layer,
                             p.phi2_inf + s2 * a2 * layer)
    for tag, a1, a2 in [("lin_bothup_0p09", 0.09, 0.02),
                        ("lin_bothup_0p3", 0.3, 0.3)]:
        out[tag] = _feasible(p.phi1_inf + a1 * layer, p.phi2_inf + a2 * layer)

    decay = np.exp(-z / 1.0)
    for tag, a1, a2, sa, sb in [("exp_thin", 0.015, 0.004, s1, s2),
                                ("exp_0p09", 0.09, 0.02, s1, s2),
                                ("exp_bothup", 0.09, 0.02, 1.0, 1.0)]:
        out[tag] = _feasible(p.phi1_inf + sa * a1 * decay,
                             p.phi2_inf + sb * a2 * decay)

    # a slab of the coexisting phase against the wall, three widths
    for width in (1.0, 2.0, 3.0):
        ramp = np.clip((z - width) / 1.0, 0.0, 1.0)
        out[f"coat_{width:g}"] = _feasible(
            conj[0] + (p.phi1_inf - conj[0]) * ramp,
            conj[1] + (p.phi2_inf - conj[1]) * ramp)
    return out


def wall_states(solver, entry, dedupe_tol=1e-6, only=None):
    """Every distinct wall solution at this far field, classified.

    A state counts as coated when its profile reaches the coexisting
    composition; at coexistence such a state has gamma = gamma(conjugate) +
    sigma and carries no information the other two quantities do not already
    give. The thin branch is the lowest-gamma state that is not coated.
    """
    p = solver.p
    p.phi1_inf = float(entry["phi1"])
    p.phi2_inf = float(entry["phi2"])
    conj = (entry["phi1_conj"], entry["phi2_conj"])
    n = solver.N

    bank = guess_bank(p, solver.z, conj)
    if only is not None:
        bank = {k: v for k, v in bank.items() if k in only}
    found = []
    for tag, U0 in bank.items():
        U, ok = solver.solve(U0)
        if not ok:
            continue
        gamma, cs1, cs2 = solver.surface_metrics(U)
        cs = cs1 + cs2
        if any(abs(cs - q["cs"]) < dedupe_tol for q in found):
            continue
        phi1, phi2 = U[:n], U[n:]
        found.append({
            "from_guess": tag,
            "gamma": float(gamma),
            "cs": float(cs),
            "max_dev": float(np.hypot(phi1 - p.phi1_inf,
                                      phi2 - p.phi2_inf).max()),
            "to_conj": float(np.hypot(phi1 - conj[0], phi2 - conj[1]).min()),
        })
    for q in found:
        q["coated"] = q["to_conj"] < COATED_FRAC * entry["length"]
    return sorted(found, key=lambda q: q["cs"])


def wall_gamma(solver, scan_cfg, entry):
    states = wall_states(solver, entry)
    out = {"ok": False, "reason": "no_convergence", "gamma": float("nan"),
           "max_dev": float("nan"), "cs": float("nan"), "n_states": len(states),
           "n_uncoated": 0, "from_guess": "", "gamma_coated": float("nan")}
    if not states:
        return out

    coated = [q for q in states if q["coated"]]
    if coated:
        out["gamma_coated"] = min(q["gamma"] for q in coated)

    bare = [q for q in states if not q["coated"]]
    out["n_uncoated"] = len(bare)
    if not bare:
        out["reason"] = "coated"
        return out

    best = min(bare, key=lambda q: q["gamma"])
    out.update({"ok": True, "reason": "", "gamma": best["gamma"],
                "max_dev": best["max_dev"], "cs": best["cs"],
                "from_guess": best["from_guess"]})
    return out


# --------------------------------------------------------------------------
# planar interfacial free energy


class TiePath:
    """sqrt(2 W) arc-length functional along a path drawn over one tie line."""

    def __init__(self, a, b, p, n_nodes):
        self.a = np.asarray(a, dtype=float)
        self.b = np.asarray(b, dtype=float)
        self.p = p
        d = self.b - self.a
        self.length = float(np.hypot(d[0], d[1]))
        t_hat = d / self.length
        self.n_hat = np.array([-t_hat[1], t_hat[0]])
        self.u = np.linspace(0.0, 1.0, int(n_nodes))
        self.base = self.a[None, :] + self.u[:, None] * d[None, :]
        self.f_a = float(model.free_energy(self.a[0], self.a[1], p))
        mu1_a, mu2_a = model.chemical_potential(self.a[0], self.a[1], p)
        self.mu_a = np.array([float(mu1_a), float(mu2_a)])

    def points(self, v_free):
        v = np.concatenate([[0.0], np.asarray(v_free, dtype=float), [0.0]])
        return self.base + v[:, None] * self.n_hat[None, :]

    def _W(self, pts):
        f = model.free_energy(pts[:, 0], pts[:, 1], self.p)
        w = (f - self.f_a
             - self.mu_a[0] * (pts[:, 0] - self.a[0])
             - self.mu_a[1] * (pts[:, 1] - self.a[1]))
        return np.maximum(w, 0.0)

    def value(self, v_free):
        pts = self.points(v_free)
        seg = pts[1:] - pts[:-1]
        seg_len = np.hypot(seg[:, 0], seg[:, 1])
        mid = 0.5 * (pts[1:] + pts[:-1])
        return float(np.sum(np.sqrt(2.0 * self._W(mid)) * seg_len))

    def value_grad(self, v_free):
        pts = self.points(v_free)
        seg = pts[1:] - pts[:-1]
        seg_len = np.hypot(seg[:, 0], seg[:, 1])
        mid = 0.5 * (pts[1:] + pts[:-1])
        s = np.sqrt(2.0 * self._W(mid))
        total = float(np.sum(s * seg_len))

        mu1, mu2 = model.chemical_potential(mid[:, 0], mid[:, 1], self.p)
        grad_w = np.stack([mu1 - self.mu_a[0], mu2 - self.mu_a[1]], axis=1)
        safe = np.maximum(s, 1e-30)
        grad_s = np.where(s[:, None] > 0.0, grad_w / safe[:, None], 0.0)
        t_seg = seg / seg_len[:, None]

        d_pts = np.zeros_like(pts)
        d_pts[:-1] += 0.5 * grad_s * seg_len[:, None] - s[:, None] * t_seg
        d_pts[1:] += 0.5 * grad_s * seg_len[:, None] + s[:, None] * t_seg
        return total, d_pts[1:-1] @ self.n_hat


def sigma_tie(a, b, p, n_nodes):
    path = TiePath(a, b, p, n_nodes)
    v0 = np.zeros(len(path.u) - 2)
    bound = V_BOUND_FRAC * path.length
    res = minimize(path.value_grad, v0, jac=True, method="L-BFGS-B",
                   bounds=[(-bound, bound)] * len(v0),
                   options={"maxiter": 500, "ftol": 1e-14, "gtol": 1e-10})
    v = res.x
    return {"sigma": float(res.fun), "straight": path.value(v0),
            "max_offset": float(np.abs(v).max()) if len(v) else 0.0,
            "length": path.length, "nit": int(res.nit)}


# --------------------------------------------------------------------------


def pca_direction(points):
    q = np.asarray(points, dtype=float)
    q = q - q.mean(axis=0)
    _, _, vt = np.linalg.svd(q, full_matrices=False)
    return vt[0]


def angle_between(u, v):
    u = np.asarray(u) / np.linalg.norm(u)
    v = np.asarray(v) / np.linalg.norm(v)
    c = abs(float(np.dot(u, v)))
    return float(np.degrees(np.arccos(min(1.0, max(-1.0, c)))))


BARE_TAGS = ("lin_thin", "exp_thin", "lin_0p09")


def omega_scan(branches, p_ref, omegas, out_csv, n_tie=40, path_nodes=401):
    """Does a second zero of S appear inside the arc as the wall weakens?

    S = gamma(R) - gamma(P) - sigma is evaluated on a subsample of tie lines for
    several wall strengths. sigma is a bulk quantity and does not depend on
    omega, so it is computed once and reused. Only the uncoated wall state is
    needed at each end; where it does not exist the wall is coated outright and
    S counts as positive.
    """
    rows_all = []
    for branch, rows in branches.items():
        for i, row in enumerate(rows):
            rows_all.append((branch, i, row))
    pick = [rows_all[k] for k in
            np.unique(np.linspace(0, len(rows_all) - 1, n_tie).astype(int))]

    log(f"--- omega scan: {len(pick)} tie lines x {len(omegas)} wall settings ---")
    t0 = time.time()
    sigma = {}
    for branch, i, row in pick:
        sigma[(branch, i)] = sigma_tie((row[0], row[1]), (row[2], row[3]),
                                       p_ref, path_nodes)["sigma"]
    log(f"sigma on the sampled tie lines: {time.time() - t0:.1f}s")

    out = []
    for om1, om2 in omegas:
        cfg = build_cfg(om1, om2)
        solver = NewtonSolver(cfg)
        t1 = time.time()
        recs = []
        for branch, i, row in pick:
            sg = sigma[(branch, i)]
            length = row[4]
            ends = []
            for a, b in (((row[0], row[1]), (row[2], row[3])),
                         ((row[2], row[3]), (row[0], row[1]))):
                entry = {"phi1": a[0], "phi2": a[1], "phi1_conj": b[0],
                         "phi2_conj": b[1], "length": length}
                st = wall_states(solver, entry, only=BARE_TAGS)
                bare = [q for q in st if not q["coated"]]
                ends.append(min(q["gamma"] for q in bare) if bare else None)
            g_a, g_b = ends
            if g_a is None or g_b is None:
                # no uncoated state at one end: that end is coated outright
                s_max = float("inf")
                s_a = s_b = float("nan")
            else:
                s_a = g_a - g_b - sg
                s_b = g_b - g_a - sg
                s_max = max(s_a, s_b)
            recs.append((branch, i, length, sg, s_a, s_b, s_max))
            out.append([om1, om2, branch, i, length, sg, s_a, s_b, s_max])

        # The two traced branches are mirror images of each other, and with
        # omega_1 different from omega_2 the mirror is not a symmetry of the
        # wall, so S differs between them at the same tie-line length. They must
        # be walked separately or the interleaving invents sign changes.
        n_wet = sum(1 for r in recs if r[6] > 0)
        parts = []
        for branch in sorted({r[0] for r in recs}):
            finite = sorted((r[2], r[6]) for r in recs
                            if r[0] == branch and np.isfinite(r[6]))
            cross = []
            for (l0, s0), (l1, s1) in zip(finite, finite[1:]):
                if np.sign(s0) != np.sign(s1) and s0 != s1:
                    cross.append(l0 + (l1 - l0) * s0 / (s0 - s1))
            parts.append(f"branch {branch}: {len(cross)} zero(s)"
                         + (f" at tie length {['%.4f' % c for c in cross]}"
                            if cross else ""))
        log(f"omega = ({om1:+.6f}, {om2:+.6f}): wet on {n_wet}/{len(recs)} tie "
            f"lines; " + "; ".join(parts) + f"  [{time.time() - t1:.1f}s]")

    with open(out_csv, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["om1", "om2", "branch", "row", "tie_length", "sigma",
                    "S_a", "S_b", "S_max"])
        for r in out:
            w.writerow([f"{r[0]:.6f}", f"{r[1]:.6f}", r[2], r[3]] +
                       [f"{v:.8f}" for v in r[4:]])
    log(f"wrote {out_csv}")


def verify_pw(pw_csv, om1, om2, n_sample=9):
    """At an archived pre-wetting point the two surface states must have equal
    gamma. This checks the state enumeration and gamma against the archived
    line; the line is used to verify, never to locate anything."""
    cfg = build_cfg(om1, om2)
    solver = NewtonSolver(cfg)
    pw = [(float(r["phi1_inf"]), float(r["phi2_inf"])) for r in read_rows(pw_csv)]
    pick = [pw[i] for i in np.linspace(0, len(pw) - 1, n_sample).astype(int)]
    log(f"verifying at {len(pick)} of {len(pw)} archived pre-wetting points, "
        f"omega = ({om1}, {om2})")
    worst = 0.0
    for p1, p2 in pick:
        entry = {"phi1": p1, "phi2": p2, "phi1_conj": p1, "phi2_conj": p2,
                 "length": 0.0}
        got = wall_states(solver, entry)
        if len(got) < 2:
            log(f"  ({p1:.6f}, {p2:.6f}): only {len(got)} state found")
            continue
        thin, thick = got[0], got[-1]
        gap = thin["gamma"] - thick["gamma"]
        worst = max(worst, abs(gap))
        log(f"  ({p1:.6f}, {p2:.6f}): {len(got)} states, "
            f"gamma thin {thin['gamma']:+.8f} (cs {thin['cs']:.4f}), "
            f"thick {thick['gamma']:+.8f} (cs {thick['cs']:.4f}), "
            f"gap {gap:+.2e}")
    log(f"worst gamma gap on the archived pre-wetting line {worst:.2e}")
    return worst


def verify_states(states_csv, om1, om2):
    """Re-solve at far fields whose surface states are already on record.

    profiles/states.csv holds the two states found at three pre-wetting points
    of the omega = (0.25, -0.375) case, produced by scripts/tf_profile_ref.py
    and already used in the gamma decomposition. Reproducing those gamma values
    with the bank here is a direct check that the state enumeration is sound.
    """
    cfg = build_cfg(om1, om2)
    solver = NewtonSolver(cfg)
    rows = read_rows(states_csv)
    log(f"verifying against {states_csv} at omega = ({om1}, {om2})")
    worst = 0.0
    by_point = {}
    for r in rows:
        by_point.setdefault((float(r["phi1_inf"]), float(r["phi2_inf"])), []).append(r)
    for (p1, p2), recs in by_point.items():
        entry = {"phi1": p1, "phi2": p2, "phi1_conj": p1, "phi2_conj": p2,
                 "length": 0.0}
        got = wall_states(solver, entry)
        log(f"  far field ({p1:.6f}, {p2:.6f}): on record {len(recs)} states, "
            f"found {len(got)}")
        for rec in recs:
            g_ref, cs_ref = float(rec["gamma"]), float(rec["cs"])
            near = min(got, key=lambda q: abs(q["cs"] - cs_ref))
            d_gamma = abs(near["gamma"] - g_ref)
            worst = max(worst, d_gamma)
            log(f"    {rec['state']}: gamma on record {g_ref:+.10f}, "
                f"found {near['gamma']:+.10f}, difference {d_gamma:.2e}; "
                f"cs {cs_ref:.6f} vs {near['cs']:.6f}")
    log(f"worst gamma difference {worst:.2e}")
    return worst


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--om1", type=float, default=0.28)
    ap.add_argument("--om2", type=float, default=-0.375)
    ap.add_argument("--tie-dir", default="out/analysis/Tf_omega_cross/critical")
    ap.add_argument("--field-dir", default="out/analysis/Tf_omega_cross/field")
    ap.add_argument("--out-dir", default="out/analysis/Tf_omega_cross/wetting")
    ap.add_argument("--max-lines", type=int)
    ap.add_argument("--path-nodes", type=int, default=401)
    ap.add_argument("--verify-states", action="store_true",
                    help="only check the state enumeration against "
                         "profiles/states.csv, then stop")
    ap.add_argument("--omega-scan", action="store_true",
                    help="only scan the wall strength for a second zero of S "
                         "inside the arc, then stop")
    ap.add_argument("--verify-pw", action="store_true",
                    help="only check that thin and thick gamma agree on the "
                         "archived pre-wetting line, then stop")
    args = ap.parse_args()

    tie_dir = ROOT / args.tie_dir
    field_dir = ROOT / args.field_dir
    out_dir = ROOT / args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    global _LOG_FH
    _LOG_FH = open(out_dir / "tf_wetting_point.log", "w")
    t_start = time.time()

    log(f"T-f wetting transition point: chi = ({CHI_12}, {CHI_13}, {CHI_23}), "
        f"omega = ({args.om1}, {args.om2}), chibb = 0")
    log(f"tie lines from {tie_dir / 'binodal_tie.csv'}")

    if args.omega_scan:
        branches = load_tie_lines(tie_dir / "binodal_tie.csv")
        base = build_cfg(args.om1, args.om2).physical
        fracs = [1.0, 0.5, 0.25, 0.125, 0.05]
        omegas = [(args.om1 * f, args.om2 * f) for f in fracs]
        omega_scan(branches, base, omegas, out_dir / "omega_scan.csv")
        log("EXIT=0")
        return

    if args.verify_pw:
        verify_pw(field_dir / "pw_line.csv", args.om1, args.om2)
        log("EXIT=0")
        return

    if args.verify_states:
        verify_states(ROOT / "out/analysis/Tf_omega_cross/profiles/states.csv",
                      0.25, -0.375)
        log("EXIT=0")
        return

    cfg = build_cfg(args.om1, args.om2)
    p = cfg.physical
    log(f"grid L = {p.L}, N = {p.N}; kappa = ({p.kappa_11}, {p.kappa_22}, "
        f"{p.kappa_12}); n = ({p.n1}, {p.n2}, {p.n3})")
    log(f"path nodes {args.path_nodes}, coated fraction {COATED_FRAC}")

    branches = load_tie_lines(tie_dir / "binodal_tie.csv")
    entries = build_loop(branches, args.max_lines)
    n_rows = sum(len(b) for b in branches.values())
    if args.max_lines:
        n_rows = min(args.max_lines, len(branches[0])) + \
                 min(args.max_lines, len(branches[1]))
    log(f"tie lines used {n_rows}, loop points {len(entries)}")
    if args.max_lines:
        log("WARNING: --max-lines truncates the loop, so the two arcs are cut "
            "and stitched at a seam. Sign changes reported across that seam are "
            "an artefact of the truncation, not wetting transitions.")

    # ---- sigma, one per tie line -----------------------------------------
    log("")
    log("--- sigma along each tie line ---")
    sigma_of = {}
    t0 = time.time()
    for branch, rows in branches.items():
        rows_used = rows[:args.max_lines] if args.max_lines else rows
        for i, row in enumerate(rows_used):
            out = sigma_tie((row[0], row[1]), (row[2], row[3]), p,
                            args.path_nodes)
            sigma_of[(branch, i)] = out
            if i % 10 == 0 or i == len(rows_used) - 1:
                log(f"tie b{branch} {i + 1}/{len(rows_used)}  "
                    f"len={out['length']:.6f}  sigma={out['sigma']:.8f}  "
                    f"bow={out['max_offset']:.2e}  nit={out['nit']}")
    log(f"sigma done in {time.time() - t0:.1f}s")

    # cross-check on the symmetric tie line: by the solute exchange symmetry of
    # this stage its path lies on phi1 = phi2, so the straight-path quadrature
    # is the answer and can be refined independently.
    sym = branches[0][0]
    fine = TiePath((sym[0], sym[1]), (sym[2], sym[3]), p, FINE_NODES)
    sigma_fine = fine.value(np.zeros(FINE_NODES - 2))
    sigma_sym = sigma_of[(0, 0)]["sigma"]
    rel = abs(sigma_sym - sigma_fine) / sigma_fine
    log("")
    log(f"symmetric tie line ({sym[0]:.6f}, {sym[1]:.6f})-({sym[2]:.6f}, "
        f"{sym[3]:.6f})")
    log(f"  optimised path, {args.path_nodes} nodes : {sigma_sym:.10f}")
    log(f"  straight path,  {FINE_NODES} nodes      : {sigma_fine:.10f}")
    log(f"  relative difference                     : {rel:.2e}")
    if rel > 1e-3:
        raise SystemExit("sigma cross-check failed: relative difference > 1e-3")

    # ---- gamma, one wall solve per loop point ----------------------------
    log("")
    log("--- wall surface free energy at each loop point ---")
    solver = NewtonSolver(cfg)
    t0 = time.time()
    n_ok = n_coated = n_fail = 0
    for i, e in enumerate(entries):
        res = wall_gamma(solver, cfg.scan, e)
        e.update(res)
        if res["ok"]:
            n_ok += 1
        elif res["reason"] == "coated":
            n_coated += 1
        else:
            n_fail += 1
        if i % 20 == 0 or i == len(entries) - 1:
            log(f"loop {i + 1}/{len(entries)}  "
                f"phi=({e['phi1']:.6f}, {e['phi2']:.6f})  "
                f"gamma={res['gamma']:.8f}  ok={int(res['ok'])}  "
                f"states={res['n_states']} bare={res['n_uncoated']}  "
                f"max_dev={res['max_dev']:.4f}  from={res['from_guess']}")
    log(f"wall solves done in {time.time() - t0:.1f}s: "
        f"uncoated {n_ok}, coated {n_coated}, not converged {n_fail}")

    # ---- Delta around the loop -------------------------------------------
    pts = np.array([[e["phi1"], e["phi2"]] for e in entries])
    conj = np.array([[e["phi1_conj"], e["phi2_conj"]] for e in entries])
    d = np.hypot(pts[None, :, 0] - conj[:, None, 0],
                 pts[None, :, 1] - conj[:, None, 1])
    partner = np.argmin(d, axis=1)
    log(f"partner lookup: worst mismatch {d[np.arange(len(entries)), partner].max():.2e}")

    gamma = np.array([e["gamma"] if e["ok"] else np.nan for e in entries])
    sigma = np.array([sigma_of[(e["branch"], e["row"])]["sigma"] for e in entries])
    gamma_conj = gamma[partner]
    delta = gamma - gamma_conj - sigma

    # ---- sign change -----------------------------------------------------
    n = len(entries)
    found = []
    for i in range(n):
        j = (i + 1) % n
        di, dj = delta[i], delta[j]
        if not (np.isfinite(di) and np.isfinite(dj)) or di == dj:
            continue
        if np.sign(di) == np.sign(dj):
            continue
        t = di / (di - dj)
        q = pts[i] + t * (pts[j] - pts[i])
        found.append({"phi1_w": float(q[0]), "phi2_w": float(q[1]),
                      "i": i, "j": j, "delta_left": float(di),
                      "delta_right": float(dj),
                      "sigma_at_point": float(sigma[i] + t * (sigma[j] - sigma[i]))})
    for w in found:
        w["kind"] = "interior"

    # If Delta never reaches zero at a generic point of the arc, the wall is
    # completely wet all the way along it and the only zeros of Delta are the
    # ends of that arc, where the two phases merge and gamma(R) - gamma(P) and
    # sigma vanish together. Those ends are the bulk critical points, so the
    # wetting transition sits on them.
    if not found:
        pos = delta[np.isfinite(delta)] > 0
        if pos.any() and bool(np.all(delta[np.isfinite(delta) & (delta > 0)] > 0)):
            crit = [(float(r["phi1_c"]), float(r["phi2_c"]))
                    for r in read_rows(tie_dir / "critical_points.csv")]
            for c in crit:
                k = int(np.argmin(np.hypot(pts[:, 0] - c[0], pts[:, 1] - c[1])))
                found.append({"phi1_w": c[0], "phi2_w": c[1], "i": k,
                              "j": (k + 1) % n, "delta_left": float("nan"),
                              "delta_right": float("nan"),
                              "sigma_at_point": float(sigma[k]),
                              "kind": "critical"})

    log("")
    log(f"--- wetting transition points: {len(found)} ---")
    for w in found:
        log(f"  ({w['phi1_w']:.6f}, {w['phi2_w']:.6f})  kind {w['kind']}  "
            f"at loop index {w['i']}  Delta {w['delta_left']:+.6f} -> "
            f"{w['delta_right']:+.6f}")

    # how Delta approaches zero: the ratio to sigma should grow toward the arc
    # end, since the wall preference is linear in the composition difference
    # between the phases while sigma is cubic in it.
    wet_idx = np.where(np.isfinite(delta) & (delta > 0))[0]
    if len(wet_idx):
        log(f"completely wet arc: loop index {wet_idx.min()}..{wet_idx.max()}, "
            f"Delta {delta[wet_idx].min():+.5f}..{delta[wet_idx].max():+.5f}, "
            f"Delta/sigma {(delta[wet_idx] / sigma[wet_idx]).min():.3f}.."
            f"{(delta[wet_idx] / sigma[wet_idx]).max():.3f}")

    # ---- checks ----------------------------------------------------------
    log("")
    log("--- checks ---")
    log(f"sigma range [{np.nanmin(sigma):.6f}, {np.nanmax(sigma):.6f}], "
        f"all positive: {bool(np.all(sigma > 0))}")

    # The Young bound |gamma(R) - gamma(P)| <= sigma is a statement about the
    # partially wet side only. Where Delta > 0 the wall is already coated and the
    # uncoated profile is metastable, so exceeding the bound there is expected
    # and is in fact the signature of complete wetting.
    both = np.isfinite(gamma) & np.isfinite(gamma_conj)
    excess = np.abs(gamma - gamma_conj) - sigma
    partial = both & (delta < 0)
    wet = both & (delta > 0)
    n_viol = int(np.sum(excess[partial] > 0))
    log(f"Young bound on the partially wet side (Delta < 0): {n_viol} of "
        f"{int(partial.sum())} exceed it, worst excess "
        f"{(excess[partial].max() if partial.any() else float('nan')):+.2e}")
    log(f"completely wet side (Delta > 0): {int(wet.sum())} points, "
        f"as expected all beyond the bound: "
        f"{bool(np.all(excess[wet] > 0)) if wet.any() else 'n/a'}")

    # The solver's own coated state must satisfy gamma_coated = gamma(P) + sigma:
    # at coexistence such a profile is a slab of the conjugate phase against the
    # wall plus one free interface. This checks sigma and the whole bookkeeping
    # at once, without going through the sign of Delta.
    g_coat = np.array([e["gamma_coated"] for e in entries])
    predicted = gamma_conj + sigma
    ok_pair = np.isfinite(g_coat) & np.isfinite(predicted)
    if ok_pair.any():
        resid = np.abs(g_coat - predicted)[ok_pair]
        log(f"coated state vs gamma(P) + sigma: {int(ok_pair.sum())} points, "
            f"max difference {resid.max():.2e}, mean {resid.mean():.2e}")
    else:
        log("coated state vs gamma(P) + sigma: no point has both")

    pw = np.array([(float(r["phi1_inf"]), float(r["phi2_inf"]))
                   for r in read_rows(field_dir / "pw_line.csv")])
    log(f"pre-wetting points loaded: {len(pw)}")

    rows_out = []
    for w in found:
        q = np.array([w["phi1_w"], w["phi2_w"]])
        dist = np.hypot(pw[:, 0] - q[0], pw[:, 1] - q[1])
        k = int(np.argmin(dist))
        order = np.argsort(dist)[:N_TANGENT_PW]
        pw_dir = pca_direction(pw[order])
        tangent = pts[w["j"]] - pts[w["i"]]
        ang = angle_between(pw_dir, tangent)
        log(f"  point ({q[0]:.6f}, {q[1]:.6f}): nearest pre-wetting point "
            f"({pw[k, 0]:.6f}, {pw[k, 1]:.6f}) at {dist[k]:.6f}; "
            f"angle between the pre-wetting line end and the binodal tangent "
            f"{ang:.2f} deg")
        rows_out.append({
            "phi1_w": f"{q[0]:.8f}", "phi2_w": f"{q[1]:.8f}",
            "delta_left": f"{w['delta_left']:.8f}",
            "delta_right": f"{w['delta_right']:.8f}",
            "loop_index_left": w["i"], "loop_index_right": w["j"],
            "sigma_at_point": f"{w['sigma_at_point']:.8f}",
            "nearest_pw_phi1": f"{pw[k, 0]:.8f}",
            "nearest_pw_phi2": f"{pw[k, 1]:.8f}",
            "nearest_pw_dist": f"{dist[k]:.8f}",
            "tangent_angle_deg": f"{ang:.4f}",
            "kind": w["kind"],
        })

    # ---- write -----------------------------------------------------------
    scan_path = out_dir / "wetting_scan.csv"
    with open(scan_path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["loop_index", "branch", "row", "phi1", "phi2",
                    "phi1_conj", "phi2_conj", "gamma", "gamma_conj", "sigma",
                    "delta", "thin_ok", "reason", "n_states", "n_uncoated",
                    "from_guess", "max_dev", "cs_total", "gamma_coated"])
        for i, e in enumerate(entries):
            w.writerow([i, e["branch"], e["row"],
                        f"{e['phi1']:.8f}", f"{e['phi2']:.8f}",
                        f"{e['phi1_conj']:.8f}", f"{e['phi2_conj']:.8f}",
                        f"{gamma[i]:.8f}", f"{gamma_conj[i]:.8f}",
                        f"{sigma[i]:.8f}", f"{delta[i]:.8f}",
                        int(bool(e["ok"])), e["reason"],
                        e["n_states"], e["n_uncoated"], e["from_guess"],
                        f"{e['max_dev']:.8f}", f"{e['cs']:.8f}",
                        f"{e['gamma_coated']:.8f}"])
    log("")
    log(f"wrote {scan_path}")

    pts_path = out_dir / "wetting_points.csv"
    fields = ["phi1_w", "phi2_w", "delta_left", "delta_right",
              "loop_index_left", "loop_index_right", "sigma_at_point",
              "nearest_pw_phi1", "nearest_pw_phi2", "nearest_pw_dist",
              "tangent_angle_deg", "kind"]
    with open(pts_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows_out)
    log(f"wrote {pts_path}")
    log(f"total {time.time() - t_start:.1f}s")
    log("EXIT=0")


if __name__ == "__main__":
    main()
