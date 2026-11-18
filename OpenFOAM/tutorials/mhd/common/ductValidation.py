"""Comparison of a fully developed duct solution with mhdAnalytic.DuctSolution.

Shared by the shercliff_duct, hunt_duct and thick_wall_duct tutorials.
"""
import numpy as np

from foamCase import (latest_time, pressure_gradient_history, read_internal_field,
                      report)


def load(case, region=None):
    """Return the latest time directory, cell centres, volumes, U and PotE."""
    t = latest_time(case)
    d = t/region if region else t
    C = read_internal_field(d/"C")
    n = len(C)
    V = read_internal_field(d/"V", n) if (d/"V").exists() else np.ones(n)
    U = read_internal_field(d/"U", n)
    phi = read_internal_field(d/"PotE", n)
    return t, C, V, U, phi


def compare(case, sol, log, title, tol, region=None, B=None, png="validation.png",
            extra_rows=(), extra_plot=None):
    """Print the comparison table, write a figure, return True on pass."""
    t, C, V, U, phi = load(case, region)
    y, z, u = C[:, 1], C[:, 2], U[:, 0]
    ubar = (u*V).sum()/V.sum()
    ua = sol.velocity(y, z)
    un = u/ubar
    e_inf = np.abs(un - ua).max()/ua.max()
    e_l2 = np.sqrt(((un - ua)**2*V).sum()/(ua**2*V).sum())

    # Potential, scaled by ubar B a. The numerical gauge is arbitrary, so
    # compare after removing the volume mean (the exact potential is odd in z).
    B = B if B is not None else sol.Ha
    pn = phi/(ubar*B)
    pn = pn - (pn*V).sum()/V.sum()
    pa = sol.potential(y, z)
    e_phi = np.abs(pn - pa).max()/np.abs(pa).max()

    hist = pressure_gradient_history(log)
    gp = hist[-1]
    tail = hist[int(0.9*len(hist)):]
    drift = (tail.max() - tail.min())/abs(gp)

    # Rows closest to the duct centre lines
    y0 = np.abs(y).min()
    z0 = np.abs(z).min()
    row = np.isclose(np.abs(y), y0) & (y > 0)
    col = np.isclose(np.abs(z), z0) & (z > 0)

    print(f"\n{title}, time {t.name}, {len(u)} fluid cells\n")
    rows = [("pressure gradient G a^2/(rho nu ubar)", gp, sol.gradp),
            ("u_max / ubar", un.max(), ua.max())]
    if np.any(row):
        rows.append((f"u_max along z at y = {y0:.3g}", un[row].max(), ua[row].max()))
    rows += list(extra_rows)
    ok = report(rows, tol)
    for name, val in [("velocity error (max, rel. to u_max)", e_inf),
                      ("velocity error (L2, relative)", e_l2),
                      ("potential error (max, rel. to |phi|max)", e_phi),
                      ("steady-state drift of G, last 10%", drift)]:
        flag = "  <-- above tolerance" if val > tol else ""
        print(f"{name:38s} {val:14.3e}{flag}")
    ok = ok and e_inf <= tol and drift <= 0.1*tol

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return ok

    s = np.linspace(-1, 1, 1601)
    fig, ax = plt.subplots(2, 2, figsize=(11, 8.5))
    i = np.argsort(z[row])
    ax[0, 0].plot(s, sol.velocity(np.full_like(s, y0), s), "k-", label="analytical")
    ax[0, 0].plot(z[row][i], un[row][i], "o", ms=3, mfc="none", label="OpenFOAM")
    ax[0, 0].set(xlabel="z / a", ylabel="u / ubar", title=f"across the side layers, y = {y0:.3g}")
    ax[0, 0].legend()
    j = np.argsort(y[col])
    ax[0, 1].plot(s, sol.velocity(s, np.full_like(s, z0)), "k-")
    ax[0, 1].plot(y[col][j], un[col][j], "o", ms=3, mfc="none")
    ax[0, 1].set(xlabel="y / a", ylabel="u / ubar", title=f"across the Hartmann layers, z = {z0:.3g}")
    ax[1, 0].plot(s, sol.potential(np.full_like(s, y0), s), "k-")
    ax[1, 0].plot(z[row][i], pn[row][i], "o", ms=3, mfc="none")
    ax[1, 0].set(xlabel="z / a", ylabel="phi / (ubar B a)", title=f"potential, y = {y0:.3g}")
    if extra_plot:
        extra_plot(ax[1, 1])
    else:
        err = (un - ua)/ua.max()
        m = np.abs(err).max()
        sc = ax[1, 1].tricontourf(z, y, err, np.linspace(-m, m, 31), cmap="RdBu_r")
        fig.colorbar(sc, ax=ax[1, 1], label="(u - u_exact)/u_max")
        ax[1, 1].set(xlabel="z / a", ylabel="y / a", title="velocity error", aspect="equal")
    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(case/png, dpi=130)
    print(f"\nWrote {png}")
    return ok
