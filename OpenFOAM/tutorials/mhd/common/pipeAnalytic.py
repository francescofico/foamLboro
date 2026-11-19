"""Reference solutions for fully developed, inductionless MHD flow in a circular
pipe in a uniform transverse magnetic field.

Non-dimensional form. Lengths are scaled by the radius a, the velocity by
G a^2/(rho nu), and the induced field b so that the equations read

    lap(u) + Ha db/dy + 1 = 0
    lap(b) + Ha du/dy     = 0

with B along y. Polar coordinates in the cross-section: y = r sin(theta),
z = r cos(theta). The current is j = curl(b e_x), so j_y = db/dz, j_z = -db/dy.

Wall conditions at r = 1 (u = 0 always):
    insulating wall                b = 0                                   (Gold 1962)
    thin wall, c = s_w t_w/(s a)   db/dr + b/c = 0
    thick wall, a < r < a + t_w    db_m/dr + kappa_m b_m = 0 for each Fourier mode,
                                   kappa_m = (sigma/sigma_w) |m| coth(|m| ln(1+t)),
                                   kappa_0 = (sigma/sigma_w) / ln(1+t)

The thick-wall condition comes from solving Laplace's equation for b exactly
in the annular wall, with b and (1/sigma) db/dn continuous at r = 1 and b = 0
outside. As t -> 0 at fixed c it tends to the thin-wall condition.

PipeSolution uses Chebyshev collocation in r (with the parity trick: no
point at r = 0) and a Fourier expansion in theta. The y-derivative couples
neighbouring Fourier modes, which gives a block-tridiagonal system.
gold_insulating() is an independent, closed-form Bessel-series solution for
the insulating pipe, accurate up to Ha ~ 20 in double precision, which is
used to check PipeSolution.
"""
import numpy as np
from scipy.special import ive


def _cheb(N):
    """Chebyshev points x_k = cos(pi k/N) and differentiation matrix (Trefethen)."""
    k = np.arange(N + 1)
    x = np.cos(np.pi*k/N)
    c = np.ones(N + 1)
    c[0] = c[-1] = 2.0
    c *= (-1.0)**k
    X = x[:, None] - x[None, :]
    D = (c[:, None]/c[None, :])/(X + np.eye(N + 1))
    D -= np.diag(D.sum(axis=1))
    return x, D


def _bary_matrix(x, xi):
    """Barycentric interpolation matrix from Chebyshev points x to points xi."""
    N = len(x) - 1
    w = (-1.0)**np.arange(N + 1)
    w[0] *= 0.5
    w[-1] *= 0.5
    d = xi[:, None] - x[None, :]
    exact = np.isclose(d, 0.0, atol=1e-15)
    d[exact] = 1.0
    M = w[None, :]/d
    M /= M.sum(axis=1, keepdims=True)
    rows = exact.any(axis=1)
    M[rows] = exact[rows].astype(float)
    return M


class PipeSolution:
    def __init__(self, Ha, wall="insulating", c=0.0, sigma_ratio=1.0, thickness=0.0,
                 n_radial=None, n_modes=None):
        self.Ha, self.wall = float(Ha), wall
        self.c, self.s, self.t = float(c), float(sigma_ratio), float(thickness)
        N = n_radial if n_radial else 2*int(max(40, 4*np.sqrt(Ha) + 20)) + 1   # odd
        M = n_modes if n_modes else int(max(48, 3*Ha))
        self.N, self.M = N, M
        self.x, D = _cheb(N)
        self.P = (N + 1)//2                       # points with r > 0, r[0] = 1
        self.r = self.x[:self.P]
        self._D, self._D2 = D, D @ D
        self._solve()

    # ------------------------------------------------------------------ setup
    def _fold(self, A, m):
        """Restrict an operator on [-1, 1] to r > 0 for a function of parity (-1)^m."""
        P, N = self.P, self.N
        p = 1.0 if m % 2 == 0 else -1.0
        return A[:P, :P] + p*A[:P, N - np.arange(P)]

    def _kappa(self, m):
        m = abs(m)
        if self.wall == "thin":
            return 1.0/self.c
        if self.wall == "thick":
            L = np.log1p(self.t)
            return (1.0/self.s)*(1.0/L if m == 0 else m/np.tanh(m*L))
        raise ValueError(self.wall)

    def _solve(self):
        Ha, P, M = self.Ha, self.P, self.M
        r = self.r
        R = np.diag(1.0/r)
        modes = np.arange(-M, M + 1)
        K = len(modes)
        n = 2*P
        A = np.zeros((K, n, n), complex)
        Lo = np.zeros((K, n, n), complex)       # coupling to mode m-1
        Up = np.zeros((K, n, n), complex)       # coupling to mode m+1
        F = np.zeros((K, n), complex)
        for i, m in enumerate(modes):
            D1, D2 = self._fold(self._D, m), self._fold(self._D2, m)
            lap = D2 + R @ D1 - m**2*R @ R
            # d/dy couplings: (1/2i)[f'_{m-1} - (m-1) f_{m-1}/r] - (1/2i)[f'_{m+1} + (m+1) f_{m+1}/r]
            Dm = self._fold(self._D, m - 1)
            Dp = self._fold(self._D, m + 1)
            cm = (Dm - (m - 1)*R)/(2j)
            cp = -(Dp + (m + 1)*R)/(2j)
            # equation rows: u-equation rows 0..P-1, b-equation rows P..2P-1
            A[i, :P, :P] = lap
            A[i, P:, P:] = lap
            Lo[i, :P, P:] = Ha*cm               # Ha db/dy in u-equation
            Up[i, :P, P:] = Ha*cp
            Lo[i, P:, :P] = Ha*cm               # Ha du/dy in b-equation
            Up[i, P:, :P] = Ha*cp
            if m == 0:
                F[i, :P] = -1.0
            # wall row r = 1 (index 0): u = 0
            for X in (A, Lo, Up):
                X[i, 0, :] = 0.0
            A[i, 0, 0] = 1.0
            F[i, 0] = 0.0
            # wall row for b (index P)
            for X in (A, Lo, Up):
                X[i, P, :] = 0.0
            if self.wall == "insulating":
                A[i, P, P] = 1.0
            else:
                A[i, P, P:] = D1[0, :]
                A[i, P, P] += self._kappa(m)
            F[i, P] = 0.0
        # block Thomas algorithm
        Cp = np.zeros_like(A)
        Fp = np.zeros_like(F)
        inv = np.linalg.inv(A[0])
        Cp[0] = inv @ Up[0]
        Fp[0] = inv @ F[0]
        for i in range(1, K):
            den = A[i] - Lo[i] @ Cp[i - 1]
            inv = np.linalg.inv(den)
            Cp[i] = inv @ Up[i]
            Fp[i] = inv @ (F[i] - Lo[i] @ Fp[i - 1])
        X = np.zeros_like(F)
        X[-1] = Fp[-1]
        for i in range(K - 2, -1, -1):
            X[i] = Fp[i] - Cp[i] @ X[i + 1]
        self.modes = modes
        self.u_m = X[:, :P]                    # (K, P), values at r[0..P-1]
        self.b_m = X[:, P:]
        # mean velocity: ubar = 2 int_0^1 u_0(r) r dr
        g, w = np.polynomial.legendre.leggauss(200)
        rq = 0.5*(g + 1.0)
        u0 = self._radial(self.u_m[M], 0, rq).real
        self.ubar = 2.0*np.sum(0.5*w*u0*rq)

    # ------------------------------------------------------------- evaluation
    def _radial(self, vals, m, rq):
        """Interpolate mode values (at r > 0 points) to radii rq, using parity."""
        N, P = self.N, self.P
        p = 1.0 if m % 2 == 0 else -1.0
        full = np.empty(N + 1, complex)
        full[:P] = vals
        full[N - np.arange(P)] = p*vals
        return _bary_matrix(self.x, np.asarray(rq, float)) @ full

    def _field(self, coef, y, z):
        y, z = np.broadcast_arrays(np.asarray(y, float), np.asarray(z, float))
        r = np.hypot(y, z).ravel()
        th = np.arctan2(y, z).ravel()
        N, P = self.N, self.P
        B = _bary_matrix(self.x, np.minimum(r, 1.0))      # (npts, N+1)
        out = np.zeros(r.shape, complex)
        for i, m in enumerate(self.modes):
            p = 1.0 if m % 2 == 0 else -1.0
            full = np.empty(N + 1, complex)
            full[:P] = coef[i]
            full[N - np.arange(P)] = p*coef[i]
            out += (B @ full)*np.exp(1j*m*th)
        return out.real.reshape(y.shape)

    def velocity(self, y, z):
        """u/ubar at points (y, z) with r <= 1."""
        return self._field(self.u_m, y, z)/self.ubar

    def induced_field(self, y, z):
        """b/ubar (same scaling as u)."""
        return self._field(self.b_m, y, z)/self.ubar

    @property
    def gradp(self):
        """Non-dimensional pressure gradient G a^2/(rho nu ubar)."""
        return 1.0/self.ubar


def gold_insulating(Ha, y, z, n_fft=4096):
    """Exact solution for the insulating pipe (the Gold 1962 problem), u in units
    of G a^2/(rho nu). W = u + b solves lap W + Ha dW/dy = -1 with W = 0 at r = 1.
    W = -y/Ha + exp(-Ha y/2) V, where lap V = (Ha/2)^2 V, so V is a series in
    I_m(Ha r/2). By symmetry u(y, z) = (W(y, z) + W(-y, z))/2.
    Round-off grows like exp(Ha)*1e-16, so keep Ha below about 20."""
    k = Ha/2.0
    th = 2*np.pi*np.arange(n_fft)/n_fft
    # boundary data for V: exp(k sin th) sin th / Ha
    Vb = np.exp(k*np.sin(th))*np.sin(th)/Ha
    cm = np.fft.fft(Vb)/n_fft
    mm = np.fft.fftfreq(n_fft, 1.0/n_fft).astype(int)

    def W(yy, zz):
        r = np.hypot(yy, zz)
        t = np.arctan2(yy, zz)
        V = np.zeros(r.shape, complex)
        for c_, m in zip(cm, mm):
            den = ive(abs(m), k)
            if abs(c_) < 1e-17*np.abs(cm).max() or den == 0.0:
                continue
            ratio = ive(abs(m), k*r)/den*np.exp(k*(r - 1.0))
            V += c_*ratio*np.exp(1j*m*t)
        return -yy/Ha + np.exp(-k*yy)*V.real

    y, z = np.broadcast_arrays(np.asarray(y, float), np.asarray(z, float))
    return 0.5*(W(y, z) + W(-y, z))
