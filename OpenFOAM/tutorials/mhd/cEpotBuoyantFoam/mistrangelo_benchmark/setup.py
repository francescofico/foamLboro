#!/usr/bin/env python3
"""Set up the magneto-convection benchmark of Mistrangelo et al. (2025),
"A code-to-code benchmark for magneto-convection in a horizontal duct",
Nucl. Fusion 65, 116006.

PbLi flow in a horizontal square duct (half-width a = 0.023 m, steel walls
t_w = 2 mm) through a fringing magnetic field B_y(x), heated from below over
-0.3 m <= x <= 0.3 m. Gravity is along -z. Dimensional (SI) units.

Writes system/caseParameters, system/{fluid,solid}/blockMeshDict and the
decomposeParDicts. The B0 field is written after meshing by setB0.py.

    python3 setup.py --case 1                      # Re = 3220, Gr = 1.42e7 (steady, laminar)
    python3 setup.py --case 2 --nproc 32           # Re = 1073, Gr = 3.55e7 (time-dependent)
    python3 setup.py --case 1 --cells 28 --axial-cells 300 --wall-cells 8
"""
import argparse
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, os.environ.get("MHD_COMMON", str(HERE.parents[1]/"common")))
from foamCase import HEADER, expansion, read_parameters, symmetric_grading, write_parameters  # noqa: E402

PARAMS = HERE/"system"/"caseParameters"
old = read_parameters(PARAMS)
p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("--case", type=int, choices=(1, 2), default=int(old.get("benchmarkCase", 1)))
p.add_argument("--cells", type=int, default=int(old.get("cellsHalf", 28)),
               help="fluid cells per half-width (y and z)")
p.add_argument("--wall-cells", type=int, default=int(old.get("wallCells", 8)), help="cells across the wall")
p.add_argument("--axial-cells", type=int, default=int(old.get("axialCells", 300)), help="cells along x")
p.add_argument("--layer-cells", type=float, default=old.get("layerCells", 6.0),
               help="first cell at the Hartmann walls = (a/Ha)/layer-cells")
p.add_argument("--side-cells", type=float, default=old.get("sideCells", 12.0),
               help="first cell at the side walls = (a/sqrt(Ha))/side-cells")
p.add_argument("--x-in", type=float, default=old.get("xIn", -1.2), help="inlet position [m]")
p.add_argument("--x-out", type=float, default=old.get("xOut", 1.5), help="outlet position [m]")
p.add_argument("--nproc", type=int, default=int(old.get("nProc", 1)), help="processors (split along x)")
p.add_argument("--delta-t", type=float, default=old.get("deltaT", 0.01), help="time step [s]")
a = p.parse_args()

# ----------------------------------------------------------------------------
# Benchmark specification (paper, sections 2-3, tables 1 and 4, figure 2)
# ----------------------------------------------------------------------------
A = 0.023          # half-width along B [m]
TW = 2e-3          # wall thickness [m]
B0 = 0.5           # [T] (figure 2 caption)
HA = 235.0
CASES = {1: dict(Re=3220.0, Gr=1.42e7, u0=0.03, q=4e4, end=300.0, avg=150.0),    # q'' = 4 W/cm2
         2: dict(Re=1073.0, Gr=3.55e7, u0=0.01, q=1e5, end=600.0, avg=150.0)}    # q'' = 10 W/cm2 (from Gr ratio)
c = CASES[a.case]
T0 = 573.15        # 300 C reference and inlet temperature [K]
G = 9.81

# Derived exactly from Re and Ha (Re = u0 a/nu, Ha = B0 a sqrt(sigma/(rho nu)))
nu = c["u0"]*A/c["Re"]
# ASSUMED material data (the paper takes them from refs [10] PbLi and [11] steel
# at 300 C but does not list values). Replace with the reference values when
# available; sigma and beta below are then re-derived from Ha and Gr.
rho, cp, k = 9838.0, 189.8, 13.1           # PbLi (assumed)
sigma = (HA/(B0*A))**2*rho*nu              # S/m, makes Ha exact for the assumed rho
beta = c["Gr"]*nu**2*k/(G*c["q"]*A**4)     # 1/K, makes Gr exact (dT = q'' a/k)
kappa = k/(rho*cp)
sigma_w, k_w, rho_w, cp_w = 1.2e6, 17.5, 7900.0, 530.0   # stainless steel (assumed)
kappa_w = k_w/(rho_w*cp_w)
dt_lorentz = rho/(sigma*B0**2)

# ----------------------------------------------------------------------------
# Mesh: axis-aligned blocks. x is split at the heated zone (-0.3, 0.3);
# y and z at -a-t, -a, a, a+t. Fluid: the centre column; solid: the 8 blocks around it.
# ----------------------------------------------------------------------------
X = [a.x_in, -0.3, 0.3, a.x_out]
Lx = a.x_out - a.x_in
nxs = [max(1, round(a.axial_cells*(X[i + 1] - X[i])/Lx)) for i in range(3)]
YZ = [-A - TW, -A, A, A + TW]
ny, gy = symmetric_grading(a.cells, (A/HA)/a.layer_cells, A)            # along B (Hartmann layers)
nz, gz = symmetric_grading(a.cells, (A/HA**0.5)/a.side_cells, A)        # normal to B (side layers)
first = (A/HA)/a.layer_cells
ew = expansion(TW, a.wall_cells, min(2.0*first, TW/a.wall_cells))       # fine at the interface


def vid(ix, iy, iz):
    return ix*16 + iy*4 + iz


verts = "\n".join(f"    ({X[ix]:.10g} {YZ[iy]:.10g} {YZ[iz]:.10g})  // {vid(ix, iy, iz)}"
                  for ix in range(4) for iy in range(4) for iz in range(4))


def hexv(ix, iy, iz):
    v = [vid(ix, iy, iz), vid(ix + 1, iy, iz), vid(ix + 1, iy + 1, iz), vid(ix, iy + 1, iz),
         vid(ix, iy, iz + 1), vid(ix + 1, iy, iz + 1), vid(ix + 1, iy + 1, iz + 1), vid(ix, iy + 1, iz + 1)]
    return v


def face(v, which):
    """Outward-oriented block face: xmin, xmax, ymin, ymax, zmin, zmax."""
    idx = dict(xmin=(0, 4, 7, 3), xmax=(1, 2, 6, 5), ymin=(0, 1, 5, 4),
               ymax=(3, 7, 6, 2), zmin=(0, 3, 2, 1), zmax=(4, 5, 6, 7))[which]
    return "(" + " ".join(str(v[i]) for i in idx) + ")"


def cells_grading(i, along):
    """(cells, grading) for index i of the y or z breakpoints: 0/2 are walls, 1 is the fluid range."""
    if i == 1:
        return (ny, gy) if along == "y" else (nz, gz)
    return a.wall_cells, f"{(1.0/ew if i == 0 else ew):.10g}"       # fine towards the fluid


def block(ix, iy, iz):
    v = hexv(ix, iy, iz)
    cy, gyy = cells_grading(iy, "y")
    cz, gzz = cells_grading(iz, "z")
    return v, f"    hex ({' '.join(map(str, v))}) ({nxs[ix]} {cy} {cz}) simpleGrading (1 {gyy} {gzz})"


fluid_blocks, interface_f, inlet_f, outlet_f = [], [], [], []
for ix in range(3):
    v, txt = block(ix, 1, 1)
    fluid_blocks.append(txt)
    interface_f += [face(v, w) for w in ("ymin", "ymax", "zmin", "zmax")]
    if ix == 0:
        inlet_f.append(face(v, "xmin"))
    if ix == 2:
        outlet_f.append(face(v, "xmax"))

solid_blocks, interface_s, heated, outer, ends_in, ends_out = [], [], [], [], [], []
for ix in range(3):
    blk = {}
    for iy in range(3):
        for iz in range(3):
            if (iy, iz) == (1, 1):
                continue
            v, txt = block(ix, iy, iz)
            solid_blocks.append(txt)
            blk[(iy, iz)] = v
            if ix == 0:
                ends_in.append(face(v, "xmin"))
            if ix == 2:
                ends_out.append(face(v, "xmax"))
            if iy == 0:
                outer.append(face(v, "ymin"))
            if iy == 2:
                outer.append(face(v, "ymax"))
            if iz == 2:
                outer.append(face(v, "zmax"))
            if iz == 0:
                (heated if ix == 1 else outer).append(face(v, "zmin"))
    # interface faces in the same order as the fluid patch (ymin, ymax, zmin, zmax)
    interface_s += [face(blk[(0, 1)], "ymax"), face(blk[(2, 1)], "ymin"),
                    face(blk[(1, 0)], "zmax"), face(blk[(1, 2)], "zmin")]


def mapped(name, region, sample, faces):
    return f"""    {name}
    {{
        type            mappedWall;
        sampleMode      nearestPatchFace;
        sampleRegion    {region};
        samplePatch     {sample};
        offsetMode      uniform;
        offset          (0 0 0);
        faces
        (
            {chr(10).join("            " + f for f in faces).strip()}
        );
    }}"""


def patch(name, ptype, faces, note=""):
    return f"""    {name}{"                    // " + note if note else ""}
    {{
        type {ptype};
        faces
        (
            {chr(10).join("            " + f for f in faces).strip()}
        );
    }}"""


def mesh_dict(blocks, patches, note):
    return HEADER.format(script="setup.py") + f"""
// {note}
scale 1;

vertices
(
{verts}
);

blocks
(
{chr(10).join(blocks)}
);

boundary
(
{chr(10).join(patches)}
);

// ************************************************************************* //
"""


(HERE/"system"/"fluid").mkdir(parents=True, exist_ok=True)
(HERE/"system"/"solid").mkdir(parents=True, exist_ok=True)
(HERE/"system"/"fluid"/"blockMeshDict").write_text(mesh_dict(
    fluid_blocks,
    [patch("inlet", "patch", inlet_f, f"x = {a.x_in:g} m"),
     patch("outlet", "patch", outlet_f, f"x = {a.x_out:g} m"),
     mapped("fluid_to_solid", "solid", "solid_to_fluid", interface_f)],
    "Fluid: -a < y, z < a, three x blocks split at the heated zone."))
(HERE/"system"/"solid"/"blockMeshDict").write_text(mesh_dict(
    solid_blocks,
    [mapped("solid_to_fluid", "fluid", "fluid_to_solid", interface_s),
     patch("heatedWall", "wall", heated, "outer bottom surface, -0.3 m <= x <= 0.3 m"),
     patch("outerWall", "wall", outer, "insulating, adiabatic"),
     patch("wallInlet", "wall", ends_in), patch("wallOutlet", "wall", ends_out)],
    "Solid: steel walls t_w = 2 mm around the duct, eight blocks per x segment."))

for d in ("system", "system/fluid", "system/solid"):
    (HERE/d/"decomposeParDict").write_text(HEADER.replace("blockMeshDict", "decomposeParDict").format(script="setup.py") + f"""
// Split along x only, with the same partition in both regions, so that the
// interface faces of the assembled (coupled) potential and temperature solves
// stay on the same processor.
numberOfSubdomains {a.nproc};

method          hierarchical;

coeffs
{{
    n           ({a.nproc} 1 1);
    order       xyz;
}}

// ************************************************************************* //
""")

write_parameters(PARAMS, "setup.py", [
    ("benchmarkCase", float(a.case), "1: Re 3220, Gr 1.42e7; 2: Re 1073, Gr 3.55e7"),
    ("Ha", HA, "B0 a sqrt(sigma/(rho nu)), a = 0.023 m"),
    ("Re", c["Re"], "u0 a / nu"),
    ("Gr", c["Gr"], "g beta dT a^3/nu^2, dT = q'' a/k"),
    ("a", A, "duct half-width along B [m]"),
    ("tWall", TW, "wall thickness [m]"),
    ("B0", B0, "field in the magnet centre [T]; B_y(x) from setB0.py, eq. (1)"),
    ("u0", c["u0"], "mean inlet velocity [m/s]"),
    ("qWall", c["q"], "heat flux on the outer bottom wall, -0.3 <= x <= 0.3 m [W/m2]"),
    ("qGradient", c["q"]/k_w, "q''/k_w: fixedGradient of T on heatedWall [K/m]"),
    ("T0", T0, "reference and inlet temperature [K] (300 C)"),
    ("nu", nu, "derived from Re [m2/s]"),
    ("rho", rho, "PbLi density [kg/m3] (ASSUMED, see readme)"),
    ("cp", cp, "PbLi heat capacity [J/kg/K] (ASSUMED)"),
    ("kFluid", k, "PbLi thermal conductivity [W/m/K] (ASSUMED)"),
    ("kappaFluid", kappa, "k/(rho cp) [m2/s]"),
    ("Pr", nu/kappa, "nu/kappa"),
    ("sigmaFluid", sigma, "derived from Ha with the assumed rho [S/m]"),
    ("beta", beta, "derived from Gr with the assumed k [1/K]"),
    ("sigmaWall", sigma_w, "steel [S/m] (ASSUMED)"),
    ("kWall", k_w, "steel [W/m/K] (ASSUMED)"),
    ("kappaWall", kappa_w, "steel k/(rho cp) [m2/s] (ASSUMED rho 7900, cp 530)"),
    ("cWall", sigma_w*TW/(sigma*A), "wall conductance ratio sigma_w t_w/(sigma a)"),
    ("deltaT", a.delta_t, f"explicit Lorentz force needs deltaT < rho/(sigma B0^2) = {dt_lorentz:.4g} s"),
    ("endTime", c["end"], "[s]"),
    ("averageStart", c["avg"], "start of time averaging [s]"),
    ("writeInterval", 50.0, "[s]"),
    ("xIn", a.x_in, "[m]"), ("xOut", a.x_out, "[m]"),
    ("cellsHalf", float(a.cells), ""), ("wallCells", float(a.wall_cells), ""),
    ("axialCells", float(sum(nxs)), ""), ("layerCells", float(a.layer_cells), ""),
    ("sideCells", float(a.side_cells), ""), ("nProc", float(a.nproc), ""),
], [f"Mistrangelo et al. (2025) benchmark, case {a.case}. Regenerate with: python3 setup.py --case <1|2>"])

nf = sum(nxs)*ny*nz
ns = sum(nxs)*((ny + 2*a.wall_cells)*(nz + 2*a.wall_cells) - ny*nz)
print(f"Case {a.case}: Re = {c['Re']:g}, Ha = {HA:g}, Gr = {c['Gr']:.3g}; nu = {nu:.4g}, sigma = {sigma:.4g}, "
      f"beta = {beta:.4g}, Pr = {nu/kappa:.3g}, c_w = {sigma_w*TW/(sigma*A):.3g}")
print(f"mesh: fluid {sum(nxs)} x {ny} x {nz} = {nf/1e6:.2f} M cells, solid {ns/1e6:.2f} M cells; "
      f"deltaT = {a.delta_t:g} s (Lorentz limit {dt_lorentz:.3g} s), {c['end']/a.delta_t:.0f} steps")
if a.delta_t > 0.5*dt_lorentz:
    print("WARNING: deltaT is close to the explicit Lorentz-force limit")
