#!/usr/bin/env python3
"""Post-process the Mistrangelo et al. (2025) benchmark run.

Reads the latest time directory (time-averaged fields UMean, TMean, p_rghMean,
PotEMean if present, else the instantaneous ones) and the cell centres and
volumes (postProcess -region fluid -funcs "(writeCellCentres writeCellVolumes)"
-latestTime). It writes:
  benchmark/axial.csv     centreline u, p, T; side-wall potential; bulk temperature T_b(x);
                          heated-wall temperature T_w(x); Nusselt number Nu(x)
  benchmark/z_<x>.csv     u, T and phi along z at y = 0 for the stations of the paper
  benchmark/*.png         plots (if matplotlib is available)
and prints a few headline numbers next to approximate values read off the
paper's figures. Those are NOT benchmark data tables. The paper provides only
figures, so the comparison is indicative.

    python3 post.py [--nu-length a|2a] [--time <t>]
"""
import argparse
import os
import re
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, os.environ.get("MHD_COMMON", str(HERE.parents[1]/"common")))
from foamCase import read_internal_field, read_parameters, time_dirs  # noqa: E402

p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("--nu-length", choices=("a", "2a"), default="a",
               help="length L in Nu = h L/k; the paper's eq. (7) does not define L in the text")
p.add_argument("--time", help="time directory (default: latest)")
args = p.parse_args()

par = read_parameters(HERE/"system"/"caseParameters")
case = int(par["benchmarkCase"])
a, u0, B0, q = par["a"], par["u0"], par["B0"], par["qWall"]
rho, sig, kf, T0 = par["rho"], par["sigmaFluid"], par["kFluid"], par["T0"]
t = (HERE/args.time) if args.time else [d for d in time_dirs(HERE) if d.name != "0"][-1]
d = t/"fluid"


def field(name):
    for n in (name + "Mean", name):
        if (d/n).exists():
            return n, read_internal_field(d/n, N)
    raise FileNotFoundError(f"{name} not found in {d}")


def patch_values(path, patch):
    """Nonuniform 'value' list of one patch in an ASCII field file."""
    text = Path(path).read_text()
    m = re.search(patch + r"\s*\{(.*?)\n    \}", text, re.S)
    if not m:
        raise ValueError(f"patch {patch} not in {path}")
    blk = m.group(1)
    mv = re.search(r"value\s+nonuniform\s+List<(\w+)>\s*(\d+)\s*\(", blk)
    if not mv:
        mu = re.search(r"value\s+uniform\s+([^;]+);", blk)
        return np.array([float(v) for v in mu.group(1).strip("()").split()])
    body = blk[mv.end():]
    n = int(mv.group(2))
    if mv.group(1) == "scalar":
        return np.array(body.split(None, n)[:n], dtype=float)
    return np.array([[float(x) for x in v.split()] for v in re.findall(r"\(([^()]*)\)", body)[:n]])


C = read_internal_field(d/"C")
N = len(C)
V = read_internal_field(d/"V", N)
nU, U = field("U")
nT, T = field("T")
nP, P = field("p_rgh")
nE, E = field("PotE")
print(f"Benchmark case {case}, time {t.name}; fields {nU}, {nT}, {nP}, {nE}; {N} fluid cells")

x, y, z, u = C[:, 0], C[:, 1], C[:, 2], U[:, 0]
xs = np.unique(np.round(x, 9))
xi = np.searchsorted(xs, np.round(x, 9))
w = np.bincount(xi, u*V)
Tb = np.bincount(xi, u*T*V)/w                                      # bulk temperature, eq. (6)

# Centreline (cells next to y = z = 0) and the side wall z = -a (first fluid cell)
yc, zc = np.abs(y).min(), np.abs(z).min()
cl = np.isclose(np.abs(y), yc) & np.isclose(np.abs(z), zc)
x_, ucl_, Tcl_, Pcl_ = [], [], [], []
for xv in xs:                                                      # average the 4 cells around the axis
    s = cl & np.isclose(x, xv)
    x_.append(xv); ucl_.append(u[s].mean()); Tcl_.append(T[s].mean()); Pcl_.append(P[s].mean())
xcl, ucl, Tcl, Pcl = map(np.array, (x_, ucl_, Tcl_, Pcl_))
zw = z.min()
sw = np.isclose(np.abs(y), yc) & np.isclose(z, zw)
Esw = np.array([E[sw & np.isclose(x, xv)].mean() for xv in xs])
Eref = np.array([E[cl & np.isclose(x, xs[0])].mean()])[0]          # potential referenced to the inlet centre

# Heated-wall temperature at the interface, y = 0, z = -a (patch face values)
Cf = patch_values(d/"C", "fluid_to_solid")
Tf = patch_values(d/nT, "fluid_to_solid")
on = np.isclose(Cf[:, 2], Cf[:, 2].min()) & np.isclose(np.abs(Cf[:, 1]), np.abs(Cf[:, 1]).min())
Tw = np.interp(xs, *zip(*sorted(zip(Cf[on, 0], np.ravel(Tf[on]) if Tf.ndim else Tf[on]))))
L = a if args.nu_length == "a" else 2*a
heated = (xs >= -0.3) & (xs <= 0.3)
Nu = np.full_like(xs, np.nan)
Nu[heated] = q*L/(kf*(Tw[heated] - Tb[heated]))                    # eq. (7)

out = HERE/"benchmark"
out.mkdir(exist_ok=True)
np.savetxt(out/"axial.csv", np.column_stack([xs, ucl, rho*Pcl, Tcl - 273.15, (Esw - Eref)*1e6,
                                              Tb - 273.15, Tw - 273.15, Nu]),
           delimiter=",", header="x[m],u_axis[m/s],p_axis[Pa],T_axis[C],phi_sidewall_zminus[uV],"
                                 "T_bulk[C],T_wall_heated[C],Nu", comments="")
stations = [-0.3, 0.0, 0.15, 0.3, 0.6, 1.0]
zprof = {}
for xsx in stations:
    xv = xs[np.argmin(np.abs(xs - xsx))]
    s = np.isclose(x, xv) & np.isclose(np.abs(y), yc)
    zz = np.unique(np.round(z[s], 10))
    rows = []
    for zv in zz:
        ss = s & np.isclose(z, zv)
        rows.append((zv, u[ss].mean(), T[ss].mean() - 273.15, (E[ss].mean() - Eref)*1e6))
    zprof[xsx] = np.array(rows)
    np.savetxt(out/f"z_x{xsx:+.2f}.csv", zprof[xsx], delimiter=",",
               header=f"z[m] at x={xv:.4f} m y~0,u[m/s],T[C],phi[uV]", comments="")

# Headline numbers
core = (xcl > -0.3) & (xcl < 0.3)
dpdx = np.polyfit(xcl[core], rho*Pcl[core], 1)[0]/(sig*u0*B0**2)
z0 = zprof[0.0]
zm = zprof[-0.3]
paper = {1: dict(dpdx="-0.0855 (text; fully developed analytical value -0.0863)", ucore="~0.025 (figs 7, 13)",
                 jet="~0.08 (figs 12, 13)", Tb="~310 (fig 17)", phi="~-340 (fig 15)", Nu="-> 2.45 at x = 0.3 (text)"),
         2: dict(dpdx="(fig 30)", ucore="~0.0085 (fig 26)", jet="~0.028 (fig 26)", Tb="~376 (fig 34)",
                 phi="~-115 (fig 31)", Nu="(fig 37)")}[case]
print("\nquantity                                   this run          paper (approximate, from figures)")
rows = [("dp/dx / (sigma u0 B0^2), |x| < 0.3 m", f"{dpdx:.4g}", paper["dpdx"]),
        ("axial velocity at x = -0.3 m, centre [m/s]", f"{np.interp(0, zm[:, 0], zm[:, 1]):.4g}", paper["ucore"]),
        ("side-jet peak at x = -0.3 m [m/s]", f"{zm[:, 1].max():.4g}", paper["jet"]),
        ("bulk temperature at x = 1 m [C]", f"{np.interp(1.0, xs, Tb) - 273.15:.4g}", paper["Tb"]),
        ("side-wall potential at x = 0, z = -a [uV]", f"{(np.interp(0.0, xs, Esw) - Eref)*1e6:.4g}", paper["phi"]),
        (f"Nu at x = 0.3 m (L = {args.nu_length})", f"{np.interp(0.3, xs[heated], Nu[heated]):.4g}", paper["Nu"])]
for n_, v_, p_ in rows:
    print(f"{n_:42s} {v_:>12s}      {p_}")

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(2, 3, figsize=(15, 8))
    ax[0, 0].plot(xcl, ucl); ax[0, 0].set(xlabel="x [m]", ylabel="u [m/s]", title="axis velocity (fig 7/23)")
    ax[0, 1].plot(xcl, rho*Pcl); ax[0, 1].set(xlabel="x [m]", ylabel="p_rgh [Pa]", title="axis pressure (fig 8/30)")
    ax[0, 2].plot(xs, Tb - 273.15, label="T_b"); ax[0, 2].plot(xs, Tw - 273.15, label="T_w heated")
    ax[0, 2].set(xlabel="x [m]", ylabel="T [C]", title="bulk / wall temperature (fig 17/38)"); ax[0, 2].legend()
    for xsx in (-0.3, 0.0):
        ax[1, 0].plot(zprof[xsx][:, 0], zprof[xsx][:, 1], label=f"x = {xsx}")
        ax[1, 1].plot(zprof[xsx][:, 0], zprof[xsx][:, 2], label=f"x = {xsx}")
    ax[1, 0].set(xlabel="z [m]", ylabel="u [m/s]", title="velocity at y = 0 (figs 12, 13)"); ax[1, 0].legend()
    ax[1, 1].set(xlabel="z [m]", ylabel="T [C]", title="temperature at y = 0 (figs 9, 10)"); ax[1, 1].legend()
    ax[1, 2].plot(xs[heated], Nu[heated]); ax[1, 2].set(xlabel="x [m]", ylabel="Nu", title=f"Nusselt (L = {args.nu_length}) (fig 18/37)")
    fig.tight_layout()
    fig.savefig(out/"benchmark.png", dpi=120)
    print(f"\nWrote {out}/axial.csv, z_*.csv and benchmark.png")
except ImportError:
    print(f"\nWrote {out}/axial.csv and z_*.csv")
