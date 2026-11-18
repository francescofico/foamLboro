#!/usr/bin/env python3
"""Compare the Hartmann channel solution with Hartmann's (1937) analytical solution.

Reads the latest time directory (U, J and the cell centres C from
"postProcess -func writeCellCentres") and the meanVelocityForce output in the
solver log. Prints the errors and writes validation.png (if matplotlib is
available). Exits with status 1 if any error is above --tol.
"""
import argparse
import os
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, os.environ.get("MHD_COMMON", str(HERE.parents[1]/"common")))
from foamCase import last_pressure_gradient, latest_time, read_internal_field, read_parameters, report  # noqa: E402
from mhdAnalytic import hartmann_gradp, hartmann_profile  # noqa: E402

p = argparse.ArgumentParser(description=__doc__)
p.add_argument("--log", default=str(HERE/"log.epotFoam_isoTh"))
p.add_argument("--tol", type=float, default=0.01, help="relative tolerance")
args = p.parse_args()

par = read_parameters(HERE/"system"/"caseParameters")
Ha, electrodes = par["Ha"], par["electrodes"]
t = latest_time(HERE)
C = read_internal_field(t/"C")
n = len(C)
U = read_internal_field(t/"U", n)
J = read_internal_field(t/"J", n) if (t/"J").exists() else None
y = C[:, 1]
order = np.argsort(y)
y, u = y[order], U[order, 0]

# Cell-centre values against the analytical profile, both normalised by ubar
ubar = np.trapz(u, y)/2.0
ua = hartmann_profile(y, Ha)
e_inf = np.abs(u/ubar - ua).max()/ua.max()
e_l2 = np.sqrt(np.trapz((u/ubar - ua)**2, y)/np.trapz(ua**2, y))
gp = last_pressure_gradient(args.log)
gp_ref = hartmann_gradp(Ha, electrodes)

print(f"\nHartmann channel, Ha = {Ha:g}, {electrodes} electrodes, time {t.name}, {n} cells\n")
rows = [
    ("pressure gradient G a^2/(rho nu ubar)", gp, gp_ref),
    ("centreline velocity u(0)/ubar", np.interp(0.0, y, u)/ubar, float(hartmann_profile(0.0, Ha))),
]
ok = report(rows, args.tol)
print(f"{'profile error (max, rel. to u_max)':38s} {e_inf:14.3e}")
print(f"{'profile error (L2, relative)':38s} {e_l2:14.3e}")
ok = ok and e_inf <= args.tol
if J is not None:
    jz = np.trapz(J[order, 2], y)
    print(f"{'net spanwise current / (sigma ubar B a)':38s} {jz/(2*Ha):14.3e}"
          f"   (insulating: 0, shorted: 1)")

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    yy = np.linspace(-1, 1, 2001)
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    ax[0].plot(hartmann_profile(yy, Ha), yy, "k-", label="Hartmann (1937)")
    ax[0].plot(u/ubar, y, "o", ms=3, mfc="none", label="epotFoam_isoTh")
    ax[0].set(xlabel="u / ubar", ylabel="y / a", title=f"Ha = {Ha:g}")
    ax[0].legend()
    d = 1.0 - np.abs(y)
    ax[1].semilogx(Ha*d, u/ubar, "o", ms=3, mfc="none")
    dd = np.logspace(-3, 0, 400)
    ax[1].semilogx(Ha*dd, hartmann_profile(1 - dd, Ha), "k-")
    ax[1].set(xlabel="Ha (1 - |y|/a)", ylabel="u / ubar", title="Hartmann layer")
    fig.tight_layout()
    fig.savefig(HERE/"validation.png", dpi=130)
    print("\nWrote validation.png")
except ImportError:
    pass

print("\nPASS" if ok else "\nFAIL")
sys.exit(0 if ok else 1)
