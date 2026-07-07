"""Fixed-grid damped-Newton solver for the equilibrium pre-wetting profile.

An independent alternative to the scipy solve_bvp path in equilibrium.py. The adaptive
mesh in solve_bvp blows past max_nodes on wide (thick-film) profiles at low phi2, so the
thick surface branch is lost exactly where the pre-wetting transition needs it. A FIXED
grid cannot exhaust a node budget: a thick film of width ~O(1..5) sits fine on a box of
length L with N points. We solve the same Euler-Lagrange system as equilibrium.py, on a
fixed grid, by a damped Newton iteration with a feasibility line search — the discretisation
is the standard second-order finite-difference one used by the reference implementation,
but written here from our own thermodynamics (thermo.dW / hessian_fb / df_surf / W /
f_surf); no reference code is imported or copied.

Model (kappa_12 = 0, matching our verification cases):
    kappa_i phi_i''(z) = dW/dphi_i(phi)             interior 0 < z < L
    kappa_i phi_i'(0)  = df_surf/dphi_i(phi(0))      wall (Robin), z = 0
    phi_i(L)           = phi_{i,inf}                 far field (Dirichlet), z = L

Grid z_j = j h, j = 0..N-1, h = L/(N-1). Unknown U = [phi1_0..phi1_{N-1},
phi2_0..phi2_{N-1}] (length 2N). Residual F(U) = A U - b(U):
  A = block-diag of per-species (-kappa_i * d2/dz2) with a Neumann ghost-node row at the
      wall (main 2k, first off -2k, k = kappa_i/h^2) and the diffusion row kept at the far
      boundary;
  b(U): interior/far rows carry -(dW/dphi_i) = mu_{i,inf} - mu_i(phi); the wall row carries
      the -(2/h) df_surf source; the far row carries the +k * kappa_i phi_{i,inf} Dirichlet
      term. This is algebraically the same F used by equilibrium._rhs/_bc, discretised.

Jacobian J = A + bulk-Hessian block (d2W/dphi = hessian_fb, res-independent) + wall-BC
block ((2/h) [[2 cbb1, cbb12],[cbb12, 2 cbb2]] at node 0). Newton step dU solves J dU = -F;
a backtracking line search keeps phi1, phi2, phi_s = 1 - phi1 - phi2 strictly positive.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

import thermo as T

_EPS = 1e-12


@dataclass
class FDProfile:
    z: np.ndarray            # (N,) grid
    phi: np.ndarray          # (2, N) profiles [phi1(z), phi2(z)]
    gamma: float             # excess surface free energy
    ok: bool                 # Newton converged
    kind: str = "?"          # 'thin' / 'thick' label (set by caller)
    U: np.ndarray = None     # raw (2N,) solution, for warm-starting the next point


def _laplacian_block(kappa, N, h):
    """Per-species (-kappa d2/dz2) with a Neumann ghost-node wall row (row 0) and the plain
    diffusion row kept at the far boundary (row N-1). Matches the reference stencil so the
    wall Robin BC enters as a source term, not a matrix row replacement."""
    k = kappa / (h * h)
    main = 2.0 * k * np.ones(N)
    off = -k * np.ones(N - 1)
    M = sp.diags([main, off, off], [0, 1, -1], shape=(N, N), format="lil")
    M[0, 0] = 2.0 * k
    M[0, 1] = -2.0 * k          # ghost-node Neumann: phi_{-1} eliminated via phi'(0)
    return M


class FDSolver:
    """Fixed grid (L, N); build A and the constant wall-Jacobian once, reuse per solve."""

    def __init__(self, chi, surf, kappa, L=40.0, N=2000):
        self.chi = chi
        self.surf = surf
        self.kappa = kappa
        self.L = float(L)
        self.N = int(N)
        self.h = self.L / (self.N - 1)
        self.z = np.linspace(0.0, self.L, self.N)
        self._fac = 2.0 / self.h                       # wall source prefactor (2/h)
        self._kf1 = kappa.k1 / (self.h * self.h)       # far Dirichlet term, species 1
        self._kf2 = kappa.k2 / (self.h * self.h)       # far Dirichlet term, species 2
        A1 = _laplacian_block(kappa.k1, self.N, self.h)
        A2 = _laplacian_block(kappa.k2, self.N, self.h)
        self.A = sp.block_diag([A1, A2], format="csr")
        self._wall_jac = self._build_wall_jac()

    def _build_wall_jac(self):
        """Constant wall-BC Jacobian block: (2/h) d(df_surf)/dphi at node 0.
        d(df_surf)/dphi = [[2 cbb1, cbb12],[cbb12, 2 cbb2]] (df_surf is linear in phi)."""
        N = self.N
        rows = [0, 0, N, N]
        cols = [0, N, 0, N]
        s = self.surf
        data = [self._fac * 2.0 * s.cbb1, self._fac * s.cbb12,
                self._fac * s.cbb12, self._fac * 2.0 * s.cbb2]
        return sp.coo_matrix((data, (rows, cols)), shape=(2 * N, 2 * N)).tocsr()

    def _rhs(self, U, res, res_mu):
        """b(U): -(dW/dphi) on all rows, minus the wall source at node 0, plus the far
        Dirichlet term at node N-1. F = A U - b."""
        N = self.N
        p1 = U[:N]
        p2 = U[N:]
        dW1, dW2 = T.dW(p1, p2, self.chi, res, res_mu=res_mu)
        b1 = -dW1
        b2 = -dW2
        g1, g2 = T.df_surf(p1[0], p2[0], self.surf)     # wall force at node 0
        b1[0] -= self._fac * g1
        b2[0] -= self._fac * g2
        b1[N - 1] += self._kf1 * res[0]                 # phi1(L) = phi1_inf
        b2[N - 1] += self._kf2 * res[1]                 # phi2(L) = phi2_inf
        return np.concatenate([b1, b2])

    def _bulk_jac(self, U):
        """Bulk Hessian block d2W/dphi = hessian_fb(phi) (res-independent); assembled as
        diag(f11,f22) with f12 on the +-N off-diagonals (species coupling)."""
        N = self.N
        f11, f12, f22 = T.hessian_fb(U[:N], U[N:], self.chi)
        main = np.concatenate([f11, f22])
        return sp.diags([main, f12 * np.ones(N), f12 * np.ones(N)], [0, N, -N],
                        shape=(2 * N, 2 * N), format="csr")

    def solve(self, U_init, res, max_iter=200, tol=1e-9):
        """Damped Newton from initial guess U_init (2N,). Returns (U, converged)."""
        N = self.N
        res = (float(res[0]), float(res[1]))
        res_mu = T.mu(res[0], res[1], self.chi)
        U = np.clip(U_init.astype(float).copy(), _EPS, 1.0 - 2.0 * _EPS)
        for _ in range(max_iter):
            F = self.A.dot(U) - self._rhs(U, res, res_mu)
            if float(np.max(np.abs(F))) < tol:
                return U, True
            J = self.A + self._bulk_jac(U) + self._wall_jac
            try:
                dU = spla.spsolve(J.tocsr(), -F)
            except Exception:
                return U, False
            alpha = 1.0
            stepped = False
            while alpha > 1e-4:
                Un = U + alpha * dU
                p1, p2 = Un[:N], Un[N:]
                if np.all(p1 > _EPS) and np.all(p2 > _EPS) and np.all(1.0 - p1 - p2 > _EPS):
                    U = Un
                    stepped = True
                    break
                alpha *= 0.5
            if not stepped:
                return U, False
        return U, False

    def gamma(self, U, res):
        """Excess surface free energy: INT_0^L [W + (1/2) sum kappa_i (dphi_i/dz)^2] dz
        + f_surf(phi(0)). Trapezoid on the fixed grid (our thermo.W / f_surf)."""
        N = self.N
        p1 = U[:N]
        p2 = U[N:]
        res_mu = T.mu(res[0], res[1], self.chi)
        res_fb = T.f_b(res[0], res[1], self.chi)
        Wv = T.W(p1, p2, self.chi, res, res_mu=res_mu, res_fb=res_fb)
        d1 = np.gradient(p1, self.z)
        d2 = np.gradient(p2, self.z)
        grad = 0.5 * self.kappa.k1 * d1 * d1 + 0.5 * self.kappa.k2 * d2 * d2
        bulk = np.trapezoid(Wv + grad, self.z)
        return float(bulk + T.f_surf(p1[0], p2[0], self.surf))

    def make_profile(self, U, res, kind="?"):
        p1 = np.clip(U[:self.N], 0.0, 1.0)
        p2 = np.clip(U[self.N:], 0.0, 1.0)
        return FDProfile(z=self.z.copy(), phi=np.vstack([p1, p2]),
                         gamma=self.gamma(U, res), ok=True, kind=kind,
                         U=U.copy())


def guess_plateau(res, wall, N, film, L):
    """Initial (2N,) guess: a wall PLATEAU of composition `wall`=(w1,w2) sustained over a
    film of thickness `film`, then a smooth drop to the reservoir `res` — the shape a thick
    surface film takes. A large `film` seeds the thick basin; a small `film` the thin basin.
    Uses a tanh shoulder so the guess is smooth (helps Newton)."""
    z = np.linspace(0.0, L, N)
    res = np.asarray(res, float)
    wall = np.asarray(wall, float)
    edge = 0.5 * (1.0 - np.tanh((z - film) / max(0.5, 0.15 * L)))   # 1 at wall -> 0 far
    p1 = res[0] + (wall[0] - res[0]) * edge
    p2 = res[1] + (wall[1] - res[1]) * edge
    p1 = np.clip(p1, _EPS, 1.0 - 2.0 * _EPS)
    p2 = np.clip(p2, _EPS, 1.0 - p1 - _EPS)
    return np.concatenate([p1, p2])
