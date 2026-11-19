#!/usr/bin/env python3
"""Compare the core of the Hadley cavity with the exact parallel-flow solution
(Garandet, Alboussiere & Moreau 1992).

Far from the end walls the flow is parallel, u = u(y), and T = gamma x + theta(y).
These satisfy the full equations exactly, with
    u'' - Ha^2 u = Gr gamma y,   u(+-1) = 0      ->  u = -Gr gamma (y - sinh(Ha y)/sinh Ha)/Ha^2
    kappa theta'' = gamma u,     theta'(+-1) = 0 (adiabatic top and bottom)
theta is the temperature change produced by the flow carrying heat. The core
gradient gamma depends on the end regions, so it is measured from the solution
(a least-squares fit of the section-averaged T around x = L/2). u and theta
are then compared given that gamma.

Exits with status 1 if an error is above --tol, or the run is not steady.
"""
import argparse
import os
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, os.environ.get("MHD_COMMON", str(HERE.parents[1]/"common")))
from foamCase import read_internal_field, read_parameters, report, time_dirs  # noqa: E402
from mhdAnalytic import buoyant_channel_profile  # noqa: E402

p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("--tol", type=float, default=0.02, help="relative tolerance")
args = p.parse_args()

par = read_parameters(HERE/"system"/"caseParameters")
Ha, Gr, Pr, Lx = par["Ha"], par["Gr"], par["Pr"], par["length"]
kappa = 1.0/Pr
times = [t for t in time_dirs(HERE) if t.name != "0"]
t, tp = times[-1], (times[-2] if len(times) > 1 else None)
C = read_internal_field(t/"C")
n = len(C)
U = read_internal_field(t/"U", n)
T = read_internal_field(t/"T", n)
x, y = C[:, 0], C[:, 1]
xs = np.unique(np.round(x, 10))


def column(xc):
    """Cells of the column closest to xc, sorted by y."""
    x0 = xs[np.argmin(np.abs(xs - xc))]
    sel = np.where(np.isclose(x, x0))[0]
    return sel[np.argsort(y[sel])], x0


# --- core gradient gamma: fit the y-averaged temperature over L/2 +- 1
cols = [column(xc) for xc in xs[np.abs(xs - Lx/2) <= 1.0]]
xm = np.array([c[1] for c in cols])
Tm = np.array([np.trapz(T[c[0]], y[c[0]])/np.ptp(y[c[0]]) for c in cols])
gamma = np.polyfit(xm, Tm, 1)[0]

# --- reference profiles for this gamma
mid, x0 = column(Lx/2)
ym = y[mid]
ua = -Gr*gamma*buoyant_channel_profile(ym, Ha)
yy = np.linspace(-1.0, 1.0, 20001)
uu = -Gr*gamma*buoyant_channel_profile(yy, Ha)
dth = np.concatenate([[0.0], np.cumsum(0.5*(uu[1:] + uu[:-1])*np.diff(yy))])*gamma/kappa   # theta'
th = np.concatenate([[0.0], np.cumsum(0.5*(dth[1:] + dth[:-1])*np.diff(yy))])             # theta
th -= np.interp(0.0, yy, th)
tha = np.interp(ym, yy, th)

um = U[mid, 0]
Tcol = T[mid]
thn = Tcol - np.interp(0.0, ym, Tcol)            # theta(y) - theta(0) from OpenFOAM
e_u = np.abs(um - ua).max()/np.abs(ua).max()
e_th = np.abs(thn - tha).max()/np.abs(tha).max()
e_v = np.abs(U[mid, 1]).max()/np.abs(ua).max()    # parallel flow: v = 0 in the core
q, xq = column(Lx/2 - 2.0)
e_x = np.abs(U[q, 0] - um).max()/np.abs(ua).max()  # x-uniformity of the core (info)
drift = 0.0
if tp is not None:
    Up = read_internal_field(tp/"U", n)
    Tp = read_internal_field(tp/"T", n)
    drift = max(np.abs(U[:, 0] - Up[:, 0]).max()/np.abs(ua).max(),
                np.abs(T - Tp).max()/np.abs(tha).max())

print(f"\nHadley cavity, Ha = {Ha:g}, Gr = {Gr:g}, Pr = {Pr:g}, L/a = {Lx:g}, time {t.name}, {n} cells")
print(f"core gradient gamma = {gamma:.6g} (pure conduction would give {-1.0/Lx:.6g}); "
      f"Peclet-type number Gr |gamma| Pr / Ha^2 = {Gr*abs(gamma)*Pr/Ha**2:.3g}\n")
ok = report([("u_max at x = L/2 (top, towards cold)", um.max(), ua.max()),
             ("theta(+a) - theta(0)", thn[-1], tha[-1])], args.tol)
for name, val, lim in [(f"velocity error at x = {x0:.3g} (max, rel.)", e_u, args.tol),
                       ("advected temperature theta error (max, rel.)", e_th, args.tol),
                       ("vertical velocity in the core (rel.)", e_v, args.tol),
                       ("change over the last write interval", drift, 0.1*args.tol),
                       (f"u(x = {xq:.3g}) - u(x = {x0:.3g}) (info)", e_x, None)]:
    flag = "  <-- above tolerance" if lim is not None and val > lim else ""
    print(f"{name:46s} {val:12.3e}{flag}")
    if lim is not None:
        ok = ok and val <= lim

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.3))
    ax[0].plot(uu, yy, "k-", label="exact core (measured gamma)")
    ax[0].plot(um, ym, "o", ms=3, mfc="none", label=f"epotFoamBuo, x = {x0:.3g}")
    ax[0].plot(U[q, 0], y[q], "x", ms=3, label=f"epotFoamBuo, x = {xq:.3g}")
    ax[0].set(xlabel="u", ylabel="y / a", title=f"Ha = {Ha:g}, Gr = {Gr:g}")
    ax[0].legend(fontsize=8)
    ax[1].plot(th, yy, "k-", label="exact: kappa theta'' = gamma u")
    ax[1].plot(thn, ym, "o", ms=3, mfc="none", label="epotFoamBuo")
    ax[1].set(xlabel="theta(y) - theta(0)", ylabel="y / a", title="temperature carried by the flow")
    ax[1].legend(fontsize=8)
    sc = ax[2].tricontourf(x, y, T, 30, cmap="RdBu_r")
    fig.colorbar(sc, ax=ax[2], label="T")
    ax[2].set(xlabel="x / a", ylabel="y / a", title="temperature")
    fig.tight_layout()
    fig.savefig(HERE/"validation.png", dpi=130)
    print("\nWrote validation.png")
except ImportError:
    pass

print("\nPASS" if ok else "\nFAIL")
sys.exit(0 if ok else 1)
