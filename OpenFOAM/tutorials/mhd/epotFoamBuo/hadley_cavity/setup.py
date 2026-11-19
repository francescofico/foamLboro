#!/usr/bin/env python3
"""Set up the Hadley cavity: a long horizontal 2D cavity, 0 < x < L, -a < y < a,
heated at x = 0, cooled at x = L, adiabatic top and bottom, gravity along -y and
a vertical magnetic field (along y). z is periodic with one cell.

Writes system/blockMeshDict and system/caseParameters.

    python3 setup.py                          # reuse system/caseParameters
    python3 setup.py --ha 20 --gr 4e5 --pr 0.05 --length 12
"""
import argparse
import math
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, os.environ.get("MHD_COMMON", str(HERE.parents[1]/"common")))
from foamCase import HEADER, read_parameters, symmetric_grading, write_parameters  # noqa: E402

PARAMS = HERE/"system"/"caseParameters"
old = read_parameters(PARAMS)
p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("--ha", type=float, default=old.get("Ha", 20.0), help="Hartmann number (half-height based)")
p.add_argument("--gr", type=float, default=old.get("Gr", 4e5), help="Grashof number g beta dT a^3/nu^2")
p.add_argument("--pr", type=float, default=old.get("Pr", 0.05), help="Prandtl number")
p.add_argument("--length", type=float, default=old.get("length", 12.0), help="cavity length L/a")
p.add_argument("--cells", type=int, default=int(old.get("cellsHalf", 24)), help="cells per half-height")
p.add_argument("--length-cells", type=int, default=int(old.get("lengthCells", 96)), help="cells along the length")
p.add_argument("--layer-cells", type=float, default=old.get("layerCells", 6.0),
               help="first cell at the top/bottom walls = (1/Ha)/layer-cells")
a = p.parse_args()

Ha, Lx, Lz = a.ha, a.length, 0.05
kappa = 1.0/a.pr
dt = 0.25/Ha**2
# T must diffuse along the cavity: slowest conduction mode decays as L^2/(pi^2 kappa)
end = max(1.0, round(2.0*Lx**2/(math.pi**2*kappa), 2))
ny, gy = symmetric_grading(a.cells, 1.0/(Ha*a.layer_cells))
nx, gx = symmetric_grading(a.length_cells//2, 0.02, Lx/2)

write_parameters(PARAMS, "setup.py", [
    ("Ha", Ha, "Hartmann number B a sqrt(sigma/(rho nu)); a = half-height; a = rho = nu = sigma = 1"),
    ("B", Ha, "applied field magnitude (= Ha in these units)"),
    ("Gr", a.gr, "Grashof number g beta dT a^3/nu^2 (dT = end-wall difference, beta = 1)"),
    ("gravityY", -a.gr, "g = (0 gravityY 0)"),
    ("Pr", a.pr, "Prandtl number; kappa = 1/Pr"),
    ("length", Lx, "cavity length L/a"),
    ("deltaT", dt, "0.25/Ha^2 (explicit Lorentz force)"),
    ("endTime", end, "2 L^2/(pi^2 kappa): two decay times of the slowest conduction mode"),
    ("writeInterval", end/10, "the last two writes are kept for the steadiness check"),
    ("cellsHalf", float(a.cells), ""),
    ("lengthCells", float(nx), ""),
    ("layerCells", float(a.layer_cells), ""),
], ["Hadley cavity parameters. Regenerate with: python3 setup.py --ha <Ha> --gr <Gr>"])

(HERE/"system"/"blockMeshDict").write_text(HEADER.format(script="setup.py") + f"""
scale 1;

// x horizontal (0 < x < {Lx:g}), y vertical (-1 < y < 1, B and g along y),
// one periodic cell in z.
vertices
(
    (0  -1 0) ({Lx:g} -1 0) ({Lx:g} 1 0) (0 1 0)
    (0  -1 {Lz}) ({Lx:g} -1 {Lz}) ({Lx:g} 1 {Lz}) (0 1 {Lz})
);

blocks
(
    hex (0 1 2 3 4 5 6 7) ({nx} {ny} 1) simpleGrading ({gx} {gy} 1)
);

boundary
(
    hotWall                     // x = 0
    {{
        type wall;
        faces ((0 4 7 3));
    }}
    coldWall                    // x = L
    {{
        type wall;
        faces ((1 2 6 5));
    }}
    horizontalWalls             // y = -1 (bottom) and y = +1 (top)
    {{
        type wall;
        faces ((0 1 5 4) (3 7 6 2));
    }}
    zPeriodic0
    {{
        type cyclic;
        neighbourPatch zPeriodic1;
        faces ((0 3 2 1));
    }}
    zPeriodic1
    {{
        type cyclic;
        neighbourPatch zPeriodic0;
        faces ((4 5 6 7));
    }}
);

// ************************************************************************* //
""")
print(f"Ha = {Ha:g}, Gr = {a.gr:g}, Pr = {a.pr:g}, L/a = {Lx:g}; {nx} x {ny} cells; "
      f"deltaT = {dt:.3g}, endTime = {end:g} ({end/dt:.0f} steps)")
