"""Single-profile solver: fixed uniform grid, Newton with analytic Jacobian.

Discretization and iteration mirror the reference codes (identical backbone
there); Jacobian includes the wall chi_bb terms and all kappa_12 terms
(the more complete ternary_ref_code form). See doc/note/reference_method.md
sections 1-2 and doc/note/code_cross_comparison.md section 2.

Unknown vector U = [phi1(z_0..z_{N-1}), phi2(...)], length 2N.
EL equation: kappa_ii phi_i'' = mu_i(phi(z)) - mu_{i,inf}, written as
F(U) = A U - b(U) = 0.
Boundary z=0: Neumann row (2k, -2k) with wall field
hs_i = omega_i + 2 chi_bb_ii phi_i(0) + chi_bb_12 phi_j(0), b_i[0] -= (2/h) hs_i.
Boundary z=L: Dirichlet phi_i(L) = phi_i_inf via the RHS ghost term.
"""

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

import model

_trapz = np.trapezoid if hasattr(np, "trapezoid") else np.trapz


class NewtonSolver:
    def __init__(self, cfg):
        self.p = cfg.physical
        self.s = cfg.solver
        self.N = int(self.p.N)
        self.h = float(self.p.L) / (self.N - 1)
        self.z = np.linspace(0.0, float(self.p.L), self.N)
        self.A = self._build_matrix_A()

    def _build_matrix_A(self):
        N, h = self.N, self.h

        def block(kappa):
            k = float(kappa) / h**2
            M = sp.diags(
                [2.0 * k * np.ones(N), -k * np.ones(N - 1), -k * np.ones(N - 1)],
                [0, 1, -1],
                format="lil",
            )
            M[0, 0] = 2.0 * k
            M[0, 1] = -2.0 * k
            return M

        A11 = block(self.p.kappa_11)
        A12 = block(self.p.kappa_12)
        A22 = block(self.p.kappa_22)
        return sp.bmat([[A11, A12], [A12, A22]], format="csr")

    def _build_rhs(self, U):
        p, N = self.p, self.N
        p1, p2 = U[:N], U[N:]

        mu1_inf, mu2_inf = model.chemical_potential(p.phi1_inf, p.phi2_inf, p)
        m1, m2 = model.chemical_potential(p1, p2, p)
        b1 = mu1_inf - m1
        b2 = mu2_inf - m2

        hs1 = p.omega_1 + 2.0 * p.chi_bb_11 * p1[0] + p.chi_bb_12 * p2[0]
        hs2 = p.omega_2 + 2.0 * p.chi_bb_22 * p2[0] + p.chi_bb_12 * p1[0]
        b1[0] -= (2.0 / self.h) * hs1
        b2[0] -= (2.0 / self.h) * hs2

        k_f = 1.0 / self.h**2
        b1[N - 1] += k_f * (p.kappa_11 * p.phi1_inf + p.kappa_12 * p.phi2_inf)
        b2[N - 1] += k_f * (p.kappa_12 * p.phi1_inf + p.kappa_22 * p.phi2_inf)
        return np.concatenate([b1, b2])

    def residual(self, U):
        return self.A.dot(U) - self._build_rhs(U)

    def _jacobian(self, U):
        p, N = self.p, self.N
        f11, f12, f22 = model.hessian(U[:N], U[N:], p)
        j_bulk = sp.diags(
            [np.concatenate([f11, f22]), f12, f12],
            [0, N, -N],
            shape=(2 * N, 2 * N),
            format="csr",
        )
        # d(-b_i[0])/d(phi(0)) from the wall field hs_i.
        fac = 2.0 / self.h
        j_wall = sp.coo_matrix(
            (
                [fac * 2.0 * p.chi_bb_11, fac * p.chi_bb_12,
                 fac * 2.0 * p.chi_bb_22, fac * p.chi_bb_12],
                ([0, 0, N, N], [0, N, N, 0]),
            ),
            shape=(2 * N, 2 * N),
        ).tocsr()
        return self.A + j_bulk + j_wall

    def solve(self, U_init):
        N = self.N
        U = U_init.copy()
        for _ in range(int(self.s.max_iter)):
            res = self.residual(U)
            if np.max(np.abs(res)) < float(self.s.tol):
                return U, True

            try:
                dU = spla.spsolve(self._jacobian(U), -res)
            except Exception:
                return U, False

            alpha = 1.0
            while alpha > float(self.s.min_alpha):
                U_new = U + alpha * dU
                p1, p2 = U_new[:N], U_new[N:]
                if np.all(p1 > 0.0) and np.all(p2 > 0.0) and np.all(1.0 - p1 - p2 > 0.0):
                    U = U_new
                    break
                alpha *= 0.5
            else:
                return U, False
        return U, False

    def surface_metrics(self, U):
        """Return (Omega, cs1, cs2) for a converged profile."""
        p, N = self.p, self.N
        phi1, phi2 = U[:N], U[N:]

        f = model.free_energy(phi1, phi2, p)
        f_inf = model.free_energy(p.phi1_inf, p.phi2_inf, p)
        mu1_inf, mu2_inf = model.chemical_potential(p.phi1_inf, p.phi2_inf, p)
        W = f - f_inf - mu1_inf * (phi1 - p.phi1_inf) - mu2_inf * (phi2 - p.phi2_inf)

        dphi1 = np.gradient(phi1, self.h)
        dphi2 = np.gradient(phi2, self.h)
        f_grad = (
            0.5 * p.kappa_11 * dphi1**2
            + 0.5 * p.kappa_22 * dphi2**2
            + p.kappa_12 * dphi1 * dphi2
        )

        j_0 = (
            p.omega_1 * phi1[0]
            + p.omega_2 * phi2[0]
            + p.chi_bb_11 * phi1[0] ** 2
            + p.chi_bb_22 * phi2[0] ** 2
            + p.chi_bb_12 * phi1[0] * phi2[0]
        )

        omega = _trapz(W + f_grad, x=self.z) + j_0
        cs1 = _trapz(phi1 - p.phi1_inf, x=self.z)
        cs2 = _trapz(phi2 - p.phi2_inf, x=self.z)
        return float(omega), float(cs1), float(cs2)
