"""Analytical solutions for fully developed, inductionless MHD flows.

Non-dimensional form. Lengths are scaled by the half-width a along B, the
velocity by G a^2/(rho nu) (G = -dp/dx), and the potential by B a G a^2/(rho nu).

Coordinates: x is streamwise, y is along B (Hartmann direction), z is normal
to B (side-wall direction). The duct occupies -1 < y < 1 and -b < z < b.

Governing equations (u = u(y, z), phi = phi(y, z)):

    lap(u) - Ha^2 (u - dphi/dz) + 1 = 0
    lap(phi) = du/dz

Side walls (z = +-b): no slip, insulating (dphi/dz = 0).
Hartmann walls (y = +-1): no slip. The electrical condition is written, mode by
mode, as dphi/dy + kappa_k phi = 0, with

    insulating        kappa = 0                              (Shercliff 1953)
    thin wall         kappa = c alpha^2                      (Hunt 1965)
    thick wall        kappa = s alpha tanh(alpha t)          (exact, see below)
    perfect conductor phi = 0

For the thick wall, s = sigma_w/sigma and t = t_w/a. The wall covers -b < z < b,
its outer face and its ends are insulating, and Laplace's equation is solved
inside it exactly. That gives the Robin condition above. As t -> 0 at fixed
c = s t, it reduces to Hunt's thin-wall condition.
"""
import numpy as np


# ----------------------------------------------------------------------------
# 1D Hartmann channel
# ----------------------------------------------------------------------------
def hartmann_profile(eta, Ha):
    """u/ubar for the Hartmann channel. The shape is the same for insulating
    and short-circuited walls."""
    eta = np.asarray(eta, dtype=float)
    # (cosh Ha - cosh Ha eta) / (cosh Ha - sinh Ha / Ha), written to avoid overflow
    num = 1.0 - np.exp(-Ha*(1.0 - np.abs(eta)))*(1.0 + np.exp(-2.0*Ha*np.abs(eta)))/(1.0 + np.exp(-2.0*Ha))
    den = 1.0 - np.tanh(Ha)/Ha
    return num/den


def hartmann_gradp(Ha, electrodes="insulating"):
    """Non-dimensional pressure gradient G a^2/(rho nu ubar)."""
    if electrodes == "insulating":
        return Ha**2/(Ha/np.tanh(Ha) - 1.0)
    if electrodes == "shorted":
        return Ha**2/(1.0 - np.tanh(Ha)/Ha)
    raise ValueError(electrodes)


# ----------------------------------------------------------------------------
# Rectangular duct (Shercliff / Hunt / thick Hartmann walls)
# ----------------------------------------------------------------------------
class DuctSolution:
    """Fourier-series solution in z, exact exponentials in y."""

    def __init__(self, Ha, wall="insulating", c=0.0, sigma_ratio=0.0,
                 thickness=0.0, b=1.0, n_modes=600):
        self.Ha, self.b, self.wall = float(Ha), float(b), wall
        self.c, self.s, self.t = float(c), float(sigma_ratio), float(thickness)
        k = np.arange(n_modes)
        self.alpha = (k + 0.5)*np.pi/self.b
        self.f = 2.0*(-1.0)**k/(self.alpha*self.b)   # Fourier coefficients of 1
        self._solve()

    def _kappa(self):
        a = self.alpha
        if self.wall == "insulating":
            return np.zeros_like(a)
        if self.wall == "thin":
            return self.c*a**2
        if self.wall == "thick":
            return self.s*a*np.tanh(a*self.t)
        if self.wall == "conducting":
            return np.full_like(a, np.inf)
        raise ValueError(self.wall)

    def _solve(self):
        Ha, a, f = self.Ha, self.alpha, self.f
        r = np.sqrt(Ha**2 + 4.0*a**2)
        lam = np.stack([(r + Ha)/2.0, 2.0*a**2/(r + Ha)])      # stable roots
        s = lam**2 - a**2                                       # s1 > 0 > s2
        T = np.tanh(lam)
        kap = self._kappa()
        # Unknowns B_i = A_i cosh(lam_i):  u_k(y) = f/a^2 + sum B_i E_i(y)
        # Row 1: u(1) = 0
        # Row 2: phi'(1) + kappa phi(1) = 0, with
        #        phi(y) = f/a^3 - a sum B_i/s_i E_i(y)
        M = np.zeros((len(a), 2, 2))
        rhs = np.zeros((len(a), 2))
        M[:, 0, :] = 1.0
        rhs[:, 0] = -f/a**2
        finite = np.isfinite(kap)
        kf = np.where(finite, kap, 0.0)
        for i in range(2):
            M[:, 1, i] = np.where(finite, -a/s[i]*(lam[i]*T[i] + kf), -a/s[i])
        rhs[:, 1] = np.where(finite, -kf*f/a**3, -f/a**3)
        B = np.linalg.solve(M, rhs[..., None])[..., 0].T           # (2, n)
        self.lam, self.sroot, self.B = lam, s, B
        Uk = 2.0*f/a**2 + (2.0*B*T/lam).sum(axis=0)                # int u_k dy
        self.ubar = (Uk*2.0*(-1.0)**np.arange(len(a))/a).sum()/(4.0*self.b)

    @staticmethod
    def _E(lam, y):
        """cosh(lam y)/cosh(lam) without overflow."""
        y = np.abs(y)
        return (np.exp(lam*(y - 1.0)) + np.exp(-lam*(y + 1.0)))/(1.0 + np.exp(-2.0*lam))

    def _modes(self, y):
        y = np.asarray(y, dtype=float)[..., None]
        a = self.alpha
        E = [self._E(self.lam[i], y) for i in range(2)]
        uk = self.f/a**2 + self.B[0]*E[0] + self.B[1]*E[1]
        pk = self.f/a**3 - a*(self.B[0]/self.sroot[0]*E[0] + self.B[1]/self.sroot[1]*E[1])
        return uk, pk

    def velocity(self, y, z):
        """u/ubar at points (y, z), for |y| <= 1."""
        y, z = np.broadcast_arrays(np.asarray(y, float), np.asarray(z, float))
        uk, _ = self._modes(y)
        return (uk*np.cos(self.alpha*z[..., None])).sum(-1)/self.ubar

    def potential(self, y, z):
        """phi/(ubar B a) in the fluid. It is odd in z, so its zero-mean gauge is fixed."""
        y, z = np.broadcast_arrays(np.asarray(y, float), np.asarray(z, float))
        _, pk = self._modes(y)
        return (pk*np.sin(self.alpha*z[..., None])).sum(-1)/self.ubar

    def wall_potential(self, y, z):
        """phi/(ubar B a) inside a thick wall, 1 <= |y| <= 1 + t."""
        if self.wall != "thick":
            raise ValueError("wall_potential needs wall='thick'")
        y, z = np.broadcast_arrays(np.asarray(y, float), np.asarray(z, float))
        _, pk = self._modes(np.ones_like(y))
        a, t = self.alpha, self.t
        d = (1.0 + t - np.abs(y))[..., None]                      # distance to outer face
        # cosh(a d)/cosh(a t), written to avoid overflow
        ratio = np.exp(a*(d - t))*(1.0 + np.exp(-2.0*a*d))/(1.0 + np.exp(-2.0*a*t))
        return (pk*ratio*np.sin(a*z[..., None])).sum(-1)/self.ubar

    @property
    def gradp(self):
        """Non-dimensional pressure gradient G a^2/(rho nu ubar)."""
        return 1.0/self.ubar


# ----------------------------------------------------------------------------
# Buoyancy-driven fully developed duct flow (vertical duct, zero net flow)
# ----------------------------------------------------------------------------
def buoyant_channel_profile(eta, Ha):
    """Hartmann channel driven by a forcing f = eta (linear across the channel,
    along B), zero net flow: u'' - Ha^2 u + eta = 0, u(+-1) = 0.
    Returns u per unit forcing amplitude: (eta - sinh(Ha eta)/sinh(Ha))/Ha^2."""
    eta = np.asarray(eta, dtype=float)
    s = np.sign(eta)*np.exp(-Ha*(1.0 - np.abs(eta)))*(1.0 - np.exp(-2.0*Ha*np.abs(eta)))/(1.0 - np.exp(-2.0*Ha))
    return (eta - s)/Ha**2


class BuoyantDuctSolution:
    """Same duct and wall models as DuctSolution, driven by a body force that is
    linear across the duct instead of a uniform pressure gradient:

        lap(u) - Ha^2 (u - dphi/dz) + f = 0,     f = y  (forcing="y", grad T along B)
                                                 f = z  (forcing="z", grad T normal to B)

    The net flow is zero. u is returned per unit forcing amplitude F (u = F u_hat),
    and phi per unit B F a (phi = B F a phi_hat).

    forcing="y": u odd in y, even in z; u ~ cos(alpha z), alpha = (k+1/2) pi/b
    forcing="z": u even in y, odd in z; u ~ sin(alpha z), alpha = k pi/b
    The thick-wall Robin coefficient (sigma_w/sigma) alpha tanh(alpha t) holds for
    both z bases, because both satisfy dphi/dz = 0 at the slab ends.
    """

    def __init__(self, Ha, forcing, wall="insulating", c=0.0, sigma_ratio=0.0,
                 thickness=0.0, b=1.0, n_modes=600):
        if forcing not in ("y", "z"):
            raise ValueError(forcing)
        self.Ha, self.b, self.wall, self.forcing = float(Ha), float(b), wall, forcing
        self.c, self.s, self.t = float(c), float(sigma_ratio), float(thickness)
        if forcing == "y":
            k = np.arange(n_modes)
            self.alpha = (k + 0.5)*np.pi/self.b
            self.f = 2.0*(-1.0)**k/(self.alpha*self.b)           # cos coefficients of 1
        else:
            k = np.arange(1, n_modes + 1)
            self.alpha = k*np.pi/self.b
            self.f = 2.0*(-1.0)**(k + 1)/self.alpha              # sin coefficients of z
        self.odd_y = forcing == "y"
        self._solve()

    _kappa = DuctSolution._kappa

    def _solve(self):
        Ha, a, f = self.Ha, self.alpha, self.f
        r = np.sqrt(Ha**2 + 4.0*a**2)
        lam = np.stack([(r + Ha)/2.0, 2.0*a**2/(r + Ha)])
        s = lam**2 - a**2
        Hp = lam/np.tanh(lam) if self.odd_y else lam*np.tanh(lam)   # H'(1)
        dPy = 1.0 if self.odd_y else 0.0                              # d/dy of y-profile of f at 1
        kap = self._kappa()
        finite = np.isfinite(kap)
        kf = np.where(finite, kap, 0.0)
        M = np.zeros((len(a), 2, 2))
        rhs = np.zeros((len(a), 2))
        M[:, 0, :] = 1.0
        rhs[:, 0] = -f/a**2
        for i in range(2):
            M[:, 1, i] = np.where(finite, -a/s[i]*(Hp[i] + kf), -a/s[i])
        rhs[:, 1] = np.where(finite, -f*(dPy + kf)/a**3, -f/a**3)
        self.B = np.linalg.solve(M, rhs[..., None])[..., 0].T
        self.lam, self.sroot = lam, s

    def _H(self, lam, y):
        """cosh(lam y)/cosh(lam) (even) or sinh(lam y)/sinh(lam) (odd), stable."""
        ay = np.abs(y)
        if not self.odd_y:
            return (np.exp(lam*(ay - 1.0)) + np.exp(-lam*(ay + 1.0)))/(1.0 + np.exp(-2.0*lam))
        small = lam < 1e-3
        big = np.exp(lam*(ay - 1.0))*(-np.expm1(-2.0*lam*ay))/(-np.expm1(-2.0*lam))
        out = np.where(small, ay*(1.0 + lam**2*(ay**2 - 1.0)/6.0), big)
        return np.sign(y)*out

    def _modes(self, y):
        y = np.asarray(y, dtype=float)[..., None]
        a = self.alpha
        Py = y if self.odd_y else np.ones_like(y)
        H = [self._H(self.lam[i], y) for i in range(2)]
        uk = self.f*Py/a**2 + self.B[0]*H[0] + self.B[1]*H[1]
        pk = self.f*Py/a**3 - a*(self.B[0]/self.sroot[0]*H[0] + self.B[1]/self.sroot[1]*H[1])
        return uk, pk

    def _zbasis(self, z):
        az = self.alpha*np.asarray(z, float)[..., None]
        # (u basis, phi basis); for forcing z the potential carries a minus sign
        return (np.cos(az), np.sin(az)) if self.odd_y else (np.sin(az), -np.cos(az))

    def velocity(self, y, z):
        y, z = np.broadcast_arrays(np.asarray(y, float), np.asarray(z, float))
        uk, _ = self._modes(y)
        return (uk*self._zbasis(z)[0]).sum(-1)

    def potential(self, y, z):
        y, z = np.broadcast_arrays(np.asarray(y, float), np.asarray(z, float))
        _, pk = self._modes(y)
        return (pk*self._zbasis(z)[1]).sum(-1)

    def wall_potential(self, y, z):
        """Potential inside a thick Hartmann wall, 1 <= |y| <= 1 + t."""
        if self.wall != "thick":
            raise ValueError("wall_potential needs wall='thick'")
        y, z = np.broadcast_arrays(np.asarray(y, float), np.asarray(z, float))
        sgn = np.sign(y)[..., None] if self.odd_y else 1.0
        _, pk = self._modes(np.ones_like(y))
        a, t = self.alpha, self.t
        d = (1.0 + t - np.abs(y))[..., None]
        ratio = np.exp(a*(d - t))*(1.0 + np.exp(-2.0*a*d))/(1.0 + np.exp(-2.0*a*t))
        return (sgn*pk*ratio*self._zbasis(z)[1]).sum(-1)
