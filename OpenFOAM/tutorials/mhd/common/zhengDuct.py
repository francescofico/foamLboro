"""MHD mixed convection in a vertical square duct (Zheng, Chen & Ni 2025, JFM 1018 A48),
cases I (insulating walls) and II (thick conducting Hartmann walls, insulating
side walls). Shared by epotFoamBuo/mixed_convection_insulating and
cEpotBuoyantFoam/mixed_convection_conducting.

Mapping to the paper (l = side length = 2a):
    xbar = (z + 1)/2   (across the heated side walls, normal to B)
    ybar = (y + 1)/2   (along B)
    Ha_paper = 2 Ha    (Ha here is half-width based)
    Gr/Re = g beta dT rho / G  (buoyancy over pressure gradient; no length scale)
    wbar = u/ubar      (their flow rate Q = 1 over the unit square equals ubar = 1)
    T = 1/2 - xbar = -z/2

The paper assumes inertialess, conduction-dominated flow. In a fully developed
periodic duct u.grad(u) = u.grad(T) = 0, so the OpenFOAM case satisfies the
same equations exactly. These are linear, so the reference is the
superposition of the pressure-driven duct solution (DuctSolution) and the
buoyancy-driven one (BuoyantDuctSolution, forcing in z). The buoyant part
carries no net flow, so with ubar = 1 the pressure gradient G equals the
pressure-driven value for every Gr/Re.
"""
import sys

import numpy as np

import buoyantDuct
from foamCase import (pressure_gradient_history, read_internal_field,
                      read_parameters, report, time_dirs)
from mhdAnalytic import BuoyantDuctSolution, DuctSolution


def references(Ha, conducting, sigma_ratio=1.0, thickness=0.0):
    if conducting:
        return (DuctSolution(Ha, "thick", sigma_ratio=sigma_ratio, thickness=thickness),
                BuoyantDuctSolution(Ha, "z", "thick", sigma_ratio=sigma_ratio, thickness=thickness))
    return DuctSolution(Ha), BuoyantDuctSolution(Ha, "z")


def parse(case, description, defaults, conducting):
    """Reuse the buoyant-duct options; --gr is replaced by --gr-re (Gr/Re of the paper)."""
    old = read_parameters(case/"system"/"caseParameters")
    argv = list(sys.argv[1:])
    gr_re = old.get("grRe", defaults["grRe"])
    if "--gr-re" in argv:
        i = argv.index("--gr-re")
        gr_re = float(argv[i + 1])
        del argv[i:i + 2]
    sys.argv[1:] = argv + ["--heating", "normal"]
    d = dict(defaults, Gr=0.0, heating="normal")
    a = buoyantDuct.parse(case, description + "\n\n    --gr-re R      Gr/Re of the paper (default %g)" % gr_re,
                          d, conducting)
    a.gr_re = gr_re
    return a


def write(case, a, title):
    u1, _ = references(a.ha, a.conducting,
                       getattr(a, "sigma_ratio", 1.0), getattr(a, "thickness", 0.0))
    G = u1.gradp                                  # G a^2/(rho nu ubar), independent of Gr/Re
    a.gr = a.gr_re*G                              # g beta dT = (Gr/Re) G
    extra = [
        ("grRe", a.gr_re, "Gr/Re of Zheng et al. (2025): buoyancy over pressure gradient"),
        ("haPaper", 2.0*a.ha, "Hartmann number based on the side length l = 2a (paper convention)"),
        ("gradpRef", G, "reference G a^2/(rho nu ubar); setup uses it to set Gr = (Gr/Re) G"),
    ]
    buoyantDuct.write(case, a, title, mean_velocity=True, extra_entries=extra)
    print(f"Gr/Re = {a.gr_re:g}, Ha_paper = {2*a.ha:g}; reference G = {G:.6g}, buoyancy coefficient {a.gr:.6g}")


def validate(case, title, tol, log, conducting):
    par = read_parameters(case/"system"/"caseParameters")
    Ha, R = par["Ha"], par["grRe"]
    s, tw = par.get("sigmaRatio", 1.0), par.get("wallThickness", 0.0)
    u1, uz = references(Ha, conducting, s, tw)
    G = u1.gradp
    F = -0.5*R*G                                  # body force F*z (T = -z/2)

    times = [t for t in time_dirs(case) if t.name != "0"]
    t, tp = times[-1], (times[-2] if len(times) > 1 else None)
    sub = (lambda d: d/"fluid") if conducting else (lambda d: d)
    C = read_internal_field(sub(t)/"C")
    n = len(C)
    V = read_internal_field(sub(t)/"V", n)
    U = read_internal_field(sub(t)/"U", n)
    P = read_internal_field(sub(t)/"PotE", n)
    T = read_internal_field(sub(t)/"T", n)
    y, z, u = C[:, 1], C[:, 2], U[:, 0]
    ubar = (u*V).sum()/V.sum()
    ua = u1.velocity(y, z) + F*uz.velocity(y, z)
    pa = Ha*(u1.potential(y, z) + F*uz.potential(y, z))
    pn = P - (P*V).sum()/V.sum()
    e_u = np.abs(u/ubar - ua).max()/np.abs(ua).max()
    e_l2 = np.sqrt(((u/ubar - ua)**2*V).sum()/(ua**2*V).sum())
    e_p = np.abs(pn - pa).max()/np.abs(pa).max()
    e_T = np.abs(T + 0.5*z).max()/0.5
    hist = pressure_gradient_history(log)
    gp = hist[-1]
    drift = np.abs(u - read_internal_field(sub(tp)/"U", n)[:, 0]).max()/np.abs(ua).max() if tp else 0.0

    # quantities in the paper's normalisation, along ybar = 0.5 (y = 0)
    y0 = np.abs(y).min()
    row = np.isclose(y, y0)
    o = np.argsort(z[row])
    xb = 0.5*(z[row][o] + 1.0)
    wn = (u/ubar)[row][o]
    core = np.abs(xb - 0.5) < 0.1
    slope = np.polyfit(xb[core], wn[core], 1)[0]
    slope_paper = -R/(2*2*Ha)*4*G                 # -Gr/(2 Re Ha_paper)/Q, Q = 1/(4G)

    print(f"\n{title}, time {t.name}, {n} fluid cells")
    print(f"Paper convention: Ha_paper = {2*Ha:g}, Gr/Re = {R:g}\n")
    ok = report([("pressure gradient G a^2/(rho nu ubar)", gp, G),
                 ("u_max/ubar (hot-wall jet)", (u/ubar).max(), ua.max()),
                 ("u_min/ubar (cold side)", (u/ubar).min(), ua.min())], tol)
    for name, val, lim in [("velocity error (max, rel. to |u|max)", e_u, tol),
                           ("velocity error (L2, relative)", e_l2, None),
                           ("potential error (max, rel.)", e_p, tol),
                           ("temperature error (max, rel. to dT/2)", e_T, tol),
                           ("change over the last write interval", drift, 0.1*tol)]:
        flag = "  <-- above tolerance" if lim is not None and val > lim else ""
        print(f"{name:42s} {val:12.3e}{flag}")
        if lim is not None:
            ok = ok and val <= lim

    wall = None
    if conducting:
        Cs = read_internal_field(t/"solid"/"C")
        Ps = read_internal_field(t/"solid"/"PotE", len(Cs))
        ps = Ps - (P*V).sum()/V.sum()
        ys, zs = Cs[:, 1], Cs[:, 2]
        ymid = ys[ys > 0][np.argmin(np.abs(ys[ys > 0] - (1 + tw/2)))]
        sel = np.isclose(ys, ymid)
        pw = Ha*(u1.wall_potential(ys[sel], zs[sel]) + F*uz.wall_potential(ys[sel], zs[sel]))
        e_w = np.abs(ps[sel] - pw).max()/np.abs(pw).max()
        flag = "  <-- above tolerance" if e_w > tol else ""
        print(f"{'wall potential error (max, rel.)':42s} {e_w:12.3e}{flag}")
        ok = ok and e_w <= tol
        wall = (zs[sel], ps[sel], pw, ymid)

    print("\nComparison with the paper's statements (information only):")
    if conducting:
        c_paper = s*tw/2.0                        # c_hw = sigma_w delta_hw/(sigma l), delta = t/2
        slope_paper /= 1.0 + c_paper*2.0*Ha       # eq. (5.26): 1/(1 + c Ha), thin-wall, high Ha
        label = f"-Gr/(2 Re Ha (1 + c Ha)), c = {c_paper:g}"
    else:
        label = "-Gr/(2 Re Ha), eq. (3.22)"
    print(f"  core gradient dwbar/dxbar at ybar = 0.5: OpenFOAM {slope:.4g}, "
          f"paper scaling {label} = {slope_paper:.4g}")
    rev = (u < 0).any()
    print(f"  reverse flow present: {'yes' if rev else 'no'} "
          f"(inertialess core estimate: reverse flow for Gr/Re > 2)")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        s_ = np.linspace(-1, 1, 801)
        fig, ax = plt.subplots(2, 2, figsize=(11, 8.5))
        ax[0, 0].plot(0.5*(s_ + 1), u1.velocity(np.full_like(s_, y0), s_) + F*uz.velocity(np.full_like(s_, y0), s_),
                      "k-", label="exact (superposition)")
        ax[0, 0].plot(xb, wn, "o", ms=3, mfc="none", label="OpenFOAM")
        ax[0, 0].set(xlabel="xbar (hot at 0)", ylabel="wbar = u/ubar", title="along ybar = 0.5 (cf. paper fig. 5)")
        ax[0, 0].legend()
        col = np.isclose(z, np.abs(z).min())
        oc = np.argsort(y[col])
        ax[0, 1].plot(0.5*(s_ + 1), u1.velocity(s_, np.full_like(s_, np.abs(z).min()))
                      + F*uz.velocity(s_, np.full_like(s_, np.abs(z).min())), "k-")
        ax[0, 1].plot(0.5*(y[col][oc] + 1), (u/ubar)[col][oc], "o", ms=3, mfc="none")
        ax[0, 1].set(xlabel="ybar (along B)", ylabel="wbar", title="along xbar = 0.5")
        if wall is not None:
            zw, pnw, paw, ym = wall
            ow = np.argsort(zw)
            ax[1, 0].plot(0.5*(zw[ow] + 1), paw[ow], "k-")
            ax[1, 0].plot(0.5*(zw[ow] + 1), pnw[ow], "o", ms=3, mfc="none")
            ax[1, 0].set(xlabel="xbar", ylabel="PotE", title=f"wall potential, y = {ym:.3g}")
        else:
            ax[1, 0].plot(0.5*(s_ + 1), Ha*(u1.potential(np.full_like(s_, y0), s_) + F*uz.potential(np.full_like(s_, y0), s_)), "k-")
            ax[1, 0].plot(xb, pn[row][o], "o", ms=3, mfc="none")
            ax[1, 0].set(xlabel="xbar", ylabel="PotE", title="potential along ybar = 0.5")
        sc = ax[1, 1].tricontourf(0.5*(z + 1), 0.5*(y + 1), u/ubar, 30, cmap="RdBu_r")
        fig.colorbar(sc, ax=ax[1, 1], label="wbar")
        ax[1, 1].set(xlabel="xbar", ylabel="ybar (B)", title="OpenFOAM velocity", aspect="equal")
        fig.suptitle(title + f", Ha_paper = {2*Ha:g}, Gr/Re = {R:g}")
        fig.tight_layout()
        fig.savefig(case/"validation.png", dpi=130)
        print("\nWrote validation.png")
    except ImportError:
        pass
    print("\nPASS" if ok else "\nFAIL")
    sys.exit(0 if ok else 1)
