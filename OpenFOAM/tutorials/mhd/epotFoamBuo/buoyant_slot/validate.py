#!/usr/bin/env python3
"""Compare the mid-height profile of the buoyant slot with the exact fully
developed solution.

Far from the end walls the flow is parallel and the temperature is pure
conduction, T = -y/2. The momentum balance is then
    u'' - Ha^2 u + F y = 0,   F = -Gr/2,
so u = F (y - sinh(Ha y)/sinh(Ha))/Ha^2. The net flow is zero, so the
spanwise current also integrates to zero and the periodic z boundary does not
need a potential jump.

The script checks the velocity and temperature at mid-height (x = H/2), the
change between the last two writes (steadiness), and prints the profile
error at x = H/4 to show how far the end effects reach. Exits with status 1 if
an error is above --tol.
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

p = argparse.ArgumentParser(description=__doc__)
p.add_argument("--tol", type=float, default=0.01, help="relative tolerance")
args = p.parse_args()

par = read_parameters(HERE/"system"/"caseParameters")
Ha, Gr, H = par["Ha"], par["Gr"], par["height"]
F = -0.5*Gr
times = [t for t in time_dirs(HERE) if t.name != "0"]
t, t_prev = times[-1], times[-2] if len(times) > 1 else None
C = read_internal_field(t/"C")
n = len(C)
U = read_internal_field(t/"U", n)
T = read_internal_field(t/"T", n)
x, y = C[:, 0], C[:, 1]


def row(xc):
    """Cells of the x-row closest to xc, sorted by y."""
    xs = x[np.argmin(np.abs(x - xc))]
    sel = np.where(np.isclose(x, xs))[0]
    return sel[np.argsort(y[sel])], xs


mid, xm = row(H/2)
ua = F*buoyant_channel_profile(y[mid], Ha)
um = U[mid, 0]
e_u = np.abs(um - ua).max()/np.abs(ua).max()
e_T = np.abs(T[mid] + 0.5*y[mid]).max()/0.5
q, xq = row(H/4)
e_q = np.abs(U[q, 0] - F*buoyant_channel_profile(y[q], Ha)).max()/np.abs(ua).max()
drift = 0.0
if t_prev is not None:
    Up = read_internal_field(t_prev/"U", n)
    drift = np.abs(U[mid, 0] - Up[mid, 0]).max()/np.abs(ua).max()

print(f"\nBuoyant slot, Ha = {Ha:g}, Gr = {Gr:g}, H/a = {H:g}, time {t.name}, {n} cells\n")
i_max = np.argmax(np.abs(ua))
ok = report([("u_max at mid-height (hot side, upward)", um.max(), ua.max()),
             ("u_min at mid-height (cold side)", um.min(), ua.min())], args.tol)
for name, val, tested in [(f"velocity error at x = {xm:.3g} (max, rel.)", e_u, True),
                          (f"temperature error at x = {xm:.3g} (max, rel.)", e_T, True),
                          ("change over the last write interval", drift, True),
                          (f"velocity error at x = {xq:.3g} (info)", e_q, False)]:
    lim = 0.1*args.tol if name.startswith("change") else args.tol
    flag = "  <-- above tolerance" if tested and val > lim else ""
    print(f"{name:42s} {val:12.3e}{flag}")
    if tested:
        ok = ok and val <= lim

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    yy = np.linspace(-1, 1, 801)
    fig, ax = plt.subplots(1, 3, figsize=(14, 4.2))
    ax[0].plot(F*buoyant_channel_profile(yy, Ha), yy, "k-", label="exact (fully developed)")
    ax[0].plot(um, y[mid], "o", ms=3, mfc="none", label=f"epotFoamBuo, x = {xm:.3g}")
    ax[0].plot(U[q, 0], y[q], "x", ms=3, label=f"epotFoamBuo, x = {xq:.3g}")
    ax[0].set(xlabel="u / (nu/a)", ylabel="y / a", title=f"Ha = {Ha:g}, Gr = {Gr:g}")
    ax[0].legend(fontsize=8)
    ax[1].plot(-0.5*yy, yy, "k-")
    ax[1].plot(T[mid], y[mid], "o", ms=3, mfc="none")
    ax[1].set(xlabel="T", ylabel="y / a", title="temperature at mid-height")
    sc = ax[2].tricontourf(y, x, U[:, 0], 30, cmap="RdBu_r")
    fig.colorbar(sc, ax=ax[2], label="u")
    ax[2].set(xlabel="y / a", ylabel="x / a (height)", title="vertical velocity")
    fig.tight_layout()
    fig.savefig(HERE/"validation.png", dpi=130)
    print("\nWrote validation.png")
except ImportError:
    pass

print("\nPASS" if ok else "\nFAIL")
sys.exit(0 if ok else 1)
