"""Surface profiles phi_1(z), phi_2(z) at a KNOWN pre-wetting far field.

Not a search: the far field (phi1_inf, phi2_inf) is already known to be a
pre-wetting point, so there is nothing to scan for. At that far field the
Euler-Lagrange system

    kappa_i phi_i'' = partial W / partial phi_i        (z > 0)
    kappa_i phi_i'(0) = partial f_surf / partial phi_i(0)
    phi_i(L) = phi_i_inf

has two solutions, the thin-film and the thick-film state, and they have equal
gamma. Each initial guess just selects which solution Newton lands on; the
converged profile is a genuine solution either way. This solves at the fixed far
field from several initial guesses, keeps the distinct solutions, and reports
gamma / cs / the wall contact values for each.

The numerics below are copied from reference/prewetting_project_clean
(src/criterion.py, src/solver.py, src/initial_guess.py) rather than imported:
that project is a read-only reference, and copying keeps this script
self-contained and pinned to the code the group vouches for. kappa_12 and the
chi_bb terms are zero for this case, so the reference's simpler forms apply
verbatim.

Usage:
  python scripts/tf_profile_ref.py --out-dir DIR [--far-field p1,p2 ...]
"""

import argparse
import csv
from pathlib import Path

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

# ---------------------------------------------------------------- case setup

CASE = dict(
    L=10.0, N=1000,
    kappa_11=1.0, kappa_22=1.0, kappa_12=0.0,
    n1=1.0, n2=1.0, n3=1.0,
    chi_12=-8.5, chi_13=0.0, chi_23=0.0,
    omega_1=0.25, omega_2=-0.375,
    chi_bb_11=0.0, chi_bb_22=0.0, chi_bb_12=0.0,
    phi1_inf=0.0, phi2_inf=0.0,
)

# Far fields on the archived pre-wetting line of this case
# (pw-space data/chi12_m8p5__chi13_0__chi23_0/om1_0p25__om2_m0p375/chibb11_0__chibb22_0__chibb12_0).
DEFAULT_FAR_FIELDS = [
    ("A", 0.12985964912280704, 0.14731505515757057),
    ("B", 0.17995989974937346, 0.11878945672644850),
    ("C", 0.20000000000000000, 0.11413811988344350),
]

MAX_ITER = 500
TOL = 1e-8


class Params:
    def __init__(self, **kw):
        self.__dict__.update(kw)


# ------------------------------- copied from reference/src/criterion.py -----

def get_free_energy(p1, p2, params, eps=1e-12):
    p3 = 1.0 - p1 - p2
    p1_c = np.clip(p1, eps, 1.0 - eps)
    p2_c = np.clip(p2, eps, 1.0 - p1_c - eps)
    p3_c = np.clip(p3, eps, 1.0)

    f_ideal = (p1_c * np.log(p1_c)) / params.n1 + \
              (p2_c * np.log(p2_c)) / params.n2 + \
              (p3_c * np.log(p3_c)) / params.n3
    f_inter = params.chi_12 * p1_c * p2_c + \
              params.chi_13 * p1_c * p3_c + \
              params.chi_23 * p2_c * p3_c
    return f_ideal + f_inter


def get_chemical_potential(p1, p2, params, eps=1e-12):
    p3 = 1.0 - p1 - p2
    p1_c = np.clip(p1, eps, 1.0 - eps)
    p2_c = np.clip(p2, eps, 1.0 - p1_c - eps)
    p3_c = np.clip(p3, eps, 1.0)

    mu1 = (np.log(p1_c) + 1) / params.n1 - (np.log(p3_c) + 1) / params.n3 + \
          params.chi_12 * p2_c + params.chi_13 * (p3_c - p1_c) - params.chi_23 * p2_c
    mu2 = (np.log(p2_c) + 1) / params.n2 - (np.log(p3_c) + 1) / params.n3 + \
          params.chi_12 * p1_c - params.chi_13 * p1_c + params.chi_23 * (p3_c - p2_c)
    return mu1, mu2


# --------------------------------- copied from reference/src/solver.py ------

class TernaryNewtonSolver:
    def __init__(self, p):
        self.p = p
        self.N = int(p.N)
        self.L = float(p.L)
        self.h = self.L / (self.N - 1)
        self.z = np.linspace(0.0, self.L, self.N)
        self.A = self._build_matrix_A()

    def _build_matrix_A(self):
        N, h = self.N, self.h

        def make_block(kappa):
            k = float(kappa) / (h**2)
            main_diag = 2.0 * k * np.ones(N)
            off_diag = -k * np.ones(N - 1)
            M = sp.diags([main_diag, off_diag, off_diag], [0, 1, -1], shape=(N, N), format="lil")
            M[0, 0] = 2.0 * k
            M[0, 1] = -2.0 * k
            return M

        A11 = make_block(self.p.kappa_11)
        A12 = make_block(self.p.kappa_12)
        A22 = make_block(self.p.kappa_22)
        return sp.bmat([[A11, A12], [A12, A22]], format="csr")

    def get_bulk_derivatives(self, phi1, phi2):
        p = self.p
        p1 = np.clip(phi1, 1e-12, 1.0 - 1e-12)
        p2 = np.clip(phi2, 1e-12, 1.0 - p1 - 1e-12)
        p3 = np.clip(1.0 - p1 - p2, 1e-12, 1.0)

        f11 = 1.0 / (p.n1 * p1) + 1.0 / (p.n3 * p3) - 2.0 * p.chi_13
        f12 = 1.0 / (p.n3 * p3) + p.chi_12 - p.chi_13 - p.chi_23
        f22 = 1.0 / (p.n2 * p2) + 1.0 / (p.n3 * p3) - 2.0 * p.chi_23
        return f11, f12, f12, f22

    def _build_rhs(self, U):
        N = self.N
        p1, p2 = U[:N], U[N:]

        mu1_inf, mu2_inf = get_chemical_potential(self.p.phi1_inf, self.p.phi2_inf, self.p)
        m1, m2 = get_chemical_potential(p1, p2, self.p)

        b1 = mu1_inf - m1
        b2 = mu2_inf - m2

        hs1 = self.p.omega_1 + 2.0 * self.p.chi_bb_11 * p1[0] + self.p.chi_bb_12 * p2[0]
        hs2 = self.p.omega_2 + 2.0 * self.p.chi_bb_22 * p2[0] + self.p.chi_bb_12 * p1[0]
        b1[0] -= (2.0 / self.h) * hs1
        b2[0] -= (2.0 / self.h) * hs2

        k_f = 1.0 / (self.h**2)
        b1[N - 1] += k_f * (self.p.kappa_11 * self.p.phi1_inf)
        b2[N - 1] += k_f * (self.p.kappa_22 * self.p.phi2_inf)
        return np.concatenate([b1, b2])

    def residual(self, U):
        return self.A.dot(U) - self._build_rhs(U)

    def solve(self, U_init, max_iter=MAX_ITER, tol=TOL):
        N = self.N
        U = U_init.copy()
        for _ in range(max_iter):
            p1, p2 = U[:N], U[N:]
            res = self.residual(U)
            if float(np.max(np.abs(res))) < tol:
                return U, True, "converged"

            f11, f12, f21, f22 = self.get_bulk_derivatives(p1, p2)
            j_bulk = sp.diags(
                [np.concatenate([f11, f22]), f12, f21],
                [0, N, -N], shape=(2 * N, 2 * N), format="csr",
            )
            try:
                dU = spla.spsolve(self.A + j_bulk, -res)
            except Exception:
                return U, False, "linear_solve_failed"

            alpha = 1.0
            accepted = False
            while alpha > 1e-3:
                U_new = U + alpha * dU
                p1n, p2n = U_new[:N], U_new[N:]
                if np.all(p1n > 0.0) and np.all(p2n > 0.0) and np.all(1.0 - p1n - p2n > 0.0):
                    U = U_new
                    accepted = True
                    break
                alpha *= 0.5
            if not accepted:
                return U, False, "line_search_failed"
        return U, False, "max_iter_reached"

    def calculate_omega(self, U):
        N = self.N
        phi1, phi2 = U[:N], U[N:]

        f = get_free_energy(phi1, phi2, self.p)
        f_inf = get_free_energy(self.p.phi1_inf, self.p.phi2_inf, self.p)
        mu1_inf, mu2_inf = get_chemical_potential(self.p.phi1_inf, self.p.phi2_inf, self.p)

        W = f - f_inf - mu1_inf * (phi1 - self.p.phi1_inf) - mu2_inf * (phi2 - self.p.phi2_inf)
        dphi1 = np.gradient(phi1, self.h)
        dphi2 = np.gradient(phi2, self.h)
        f_grad = 0.5 * self.p.kappa_11 * dphi1**2 + 0.5 * self.p.kappa_22 * dphi2**2
        omega_bulk = np.trapezoid(W + f_grad, x=self.z)

        p1_0, p2_0 = phi1[0], phi2[0]
        j_0 = (self.p.omega_1 * p1_0 + self.p.omega_2 * p2_0
               + self.p.chi_bb_11 * p1_0**2 + self.p.chi_bb_22 * p2_0**2
               + self.p.chi_bb_12 * p1_0 * p2_0)
        return omega_bulk + j_0

    def surface_excess(self, U):
        N = self.N
        cs1 = np.trapezoid(U[:N] - self.p.phi1_inf, x=self.z)
        cs2 = np.trapezoid(U[N:] - self.p.phi2_inf, x=self.z)
        return float(cs1), float(cs2)


# -------------------------- copied from reference/src/initial_guess.py ------

def _project_strict_feasible(phi1, phi2, eps=1e-12):
    phi1 = np.maximum(phi1, eps)
    phi2 = np.maximum(phi2, eps)
    s = phi1 + phi2
    over = s >= (1.0 - eps)
    if np.any(over):
        phi1[over] = phi1[over] * (1.0 - eps) / s[over]
        phi2[over] = phi2[over] * (1.0 - eps) / s[over]
    phi1 = np.minimum(phi1, 1.0 - 2 * eps)
    phi2 = np.minimum(phi2, 1.0 - 2 * eps - phi1)
    return phi1, phi2


def _sign_from_wall_pref(omega):
    return +1.0 if float(omega) < 0 else -1.0


def guess_linear(p, amp1, amp2, sign1=None, sign2=None):
    z = np.linspace(0.0, p.L, int(p.N))
    layer = 1.0 - z / max(float(p.L), 1e-12)
    s1 = _sign_from_wall_pref(p.omega_1) if sign1 is None else sign1
    s2 = _sign_from_wall_pref(p.omega_2) if sign2 is None else sign2
    phi1 = p.phi1_inf + s1 * amp1 * layer
    phi2 = p.phi2_inf + s2 * amp2 * layer
    return np.concatenate(_project_strict_feasible(phi1, phi2))


def guess_exponential(p, amp1, amp2, ell_ratio, sign1=None, sign2=None):
    z = np.linspace(0.0, p.L, int(p.N))
    ell = max(ell_ratio * p.L, 1e-12)
    layer = (np.exp(-z / ell) - np.exp(-p.L / ell)) / (1.0 - np.exp(-p.L / ell))
    s1 = _sign_from_wall_pref(p.omega_1) if sign1 is None else sign1
    s2 = _sign_from_wall_pref(p.omega_2) if sign2 is None else sign2
    phi1 = p.phi1_inf + s1 * amp1 * layer
    phi2 = p.phi2_inf + s2 * amp2 * layer
    return np.concatenate(_project_strict_feasible(phi1, phi2))


# ------------------------------------------------------------------ driver

def guess_bank(p):
    """Initial guesses, all from the reference's own two guess families.

    The wall sign rule puts phi_1 below the far field (omega_1 > 0) and phi_2
    above it (omega_2 < 0). The last entries lift both solutes together, which
    is the shape a chi_12 = -8.5 thick film is expected to take.
    """
    return [
        ("lin_thin", guess_linear(p, 0.015, 0.004)),
        ("lin_thick_0p3", guess_linear(p, 0.3, 0.3)),
        ("lin_thick_0p09", guess_linear(p, 0.09, 0.02)),
        ("exp_thin", guess_exponential(p, 0.02, 0.005, 0.05)),
        ("exp_thick", guess_exponential(p, 0.10, 0.02, 0.30)),
        ("lin_bothup_0p09", guess_linear(p, 0.09, 0.02, sign1=+1.0, sign2=+1.0)),
        ("lin_bothup_0p3", guess_linear(p, 0.3, 0.3, sign1=+1.0, sign2=+1.0)),
        ("exp_bothup", guess_exponential(p, 0.10, 0.02, 0.30, sign1=+1.0, sign2=+1.0)),
    ]


def run_far_field(tag, phi1_inf, phi2_inf):
    p = Params(**CASE)
    p.phi1_inf = float(phi1_inf)
    p.phi2_inf = float(phi2_inf)
    slv = TernaryNewtonSolver(p)

    found = []
    for name, U0 in guess_bank(p):
        U, ok, status = slv.solve(U0)
        if not ok:
            print(f"  {tag} {name:16s} FAILED ({status})", flush=True)
            continue
        gamma = float(slv.calculate_omega(U))
        cs1, cs2 = slv.surface_excess(U)
        rec = {"guess": name, "gamma": gamma, "cs1": cs1, "cs2": cs2,
               "cs": cs1 + cs2, "phi1_0": float(U[:slv.N][0]),
               "phi2_0": float(U[slv.N:][0]), "U": U}
        print(f"  {tag} {name:16s} ok  gamma={gamma:+.8f} cs={rec['cs']:+.6f} "
              f"phi1(0)={rec['phi1_0']:.6f} phi2(0)={rec['phi2_0']:.6f}", flush=True)
        found.append(rec)

    # Distinct solutions of the same boundary-value problem, by total adsorption.
    distinct = []
    for rec in sorted(found, key=lambda r: r["cs"]):
        if not distinct or abs(rec["cs"] - distinct[-1]["cs"]) > 1e-6:
            distinct.append(rec)
    return slv.z, found, distinct


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--far-field", action="append", default=[],
                    help="'tag,phi1,phi2'; repeatable. Default: the three archived points.")
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if args.far_field:
        far_fields = []
        for s in args.far_field:
            tag, p1, p2 = s.split(",")
            far_fields.append((tag, float(p1), float(p2)))
    else:
        far_fields = DEFAULT_FAR_FIELDS

    print(f"case chi=({CASE['chi_12']}, {CASE['chi_13']}, {CASE['chi_23']}) "
          f"omega=({CASE['omega_1']}, {CASE['omega_2']}) chibb=0", flush=True)

    attempts, states, profiles = [], [], []
    for tag, phi1_inf, phi2_inf in far_fields:
        print(f"far field {tag}: phi1_inf={phi1_inf:.8f} phi2_inf={phi2_inf:.8f}", flush=True)
        z, found, distinct = run_far_field(tag, phi1_inf, phi2_inf)
        for rec in found:
            attempts.append({"point_id": tag, "phi1_inf": phi1_inf, "phi2_inf": phi2_inf,
                             **{k: rec[k] for k in
                                ("guess", "gamma", "cs1", "cs2", "cs", "phi1_0", "phi2_0")}})
        print(f"  -> {len(distinct)} distinct solution(s): "
              f"{[round(r['cs'], 6) for r in distinct]}", flush=True)
        for k, rec in enumerate(distinct):
            state = f"state{k}"
            states.append({"point_id": tag, "phi1_inf": phi1_inf, "phi2_inf": phi2_inf,
                           "state": state, "from_guess": rec["guess"],
                           "gamma": rec["gamma"], "cs1": rec["cs1"], "cs2": rec["cs2"],
                           "cs": rec["cs"], "phi1_0": rec["phi1_0"], "phi2_0": rec["phi2_0"]})
            N = len(z)
            phi1, phi2 = rec["U"][:N], rec["U"][N:]
            profiles += [(tag, phi1_inf, phi2_inf, state, z[i], phi1[i], phi2[i])
                         for i in range(N)]

    def dump(path, rows, fields):
        with open(path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=fields)
            w.writeheader()
            w.writerows(rows)
        print(f"wrote {path} ({len(rows)} rows)", flush=True)

    dump(out_dir / "solve_attempts.csv", attempts,
         ["point_id", "phi1_inf", "phi2_inf", "guess", "gamma", "cs1", "cs2",
          "cs", "phi1_0", "phi2_0"])
    dump(out_dir / "states.csv", states,
         ["point_id", "phi1_inf", "phi2_inf", "state", "from_guess", "gamma",
          "cs1", "cs2", "cs", "phi1_0", "phi2_0"])
    with open(out_dir / "profiles.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["point_id", "phi1_inf", "phi2_inf", "state", "z", "phi1", "phi2"])
        w.writerows(profiles)
    print(f"wrote {out_dir / 'profiles.csv'} ({len(profiles)} rows)", flush=True)


if __name__ == "__main__":
    main()
