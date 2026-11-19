#!/usr/bin/env python3
"""Set up the buoyant slot: a tall, closed vertical slot, -a < y < a, 0 < x < H,
heated at y = -a and cooled at y = +a, with B along y. z is periodic with one
cell (an infinitely wide slot).

Writes system/blockMeshDict and system/caseParameters.

    python3 setup.py                         # reuse system/caseParameters
    python3 setup.py --ha 20 --gr 4000 --height 8
"""
import argparse
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, os.environ.get("MHD_COMMON", str(HERE.parents[1]/"common")))
from foamCase import HEADER, read_parameters, symmetric_grading, write_parameters  # noqa: E402

PARAMS = HERE/"system"/"caseParameters"
old = read_parameters(PARAMS)
p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("--ha", type=float, default=old.get("Ha", 20.0), help="Hartmann number (half-width based)")
p.add_argument("--gr", type=float, default=old.get("Gr", 4000.0), help="Grashof number g beta dT a^3/nu^2")
p.add_argument("--pr", type=float, default=old.get("Pr", 0.01), help="Prandtl number")
p.add_argument("--height", type=float, default=old.get("height", 8.0), help="slot height H/a")
p.add_argument("--cells", type=int, default=int(old.get("cellsHalf", 30)), help="cells per half-width")
p.add_argument("--height-cells", type=int, default=int(old.get("heightCells", 32)), help="cells along the height")
p.add_argument("--layer-cells", type=float, default=old.get("layerCells", 6.0),
               help="first wall cell = (1/Ha)/layer-cells")
a = p.parse_args()

Ha, H, L = a.ha, a.height, 0.05
dt = 0.25/Ha**2
end = 1.0                                   # viscous times a^2/nu; steady (to ~1e-10) well before
ny, gy = symmetric_grading(a.cells, 1.0/(Ha*a.layer_cells))
nx, gx = symmetric_grading(a.height_cells//2, 0.05, H/2)

write_parameters(PARAMS, "setup.py", [
    ("Ha", Ha, "Hartmann number B a sqrt(sigma/(rho nu)); a = rho = nu = sigma = 1"),
    ("B", Ha, "applied field magnitude (= Ha in these units)"),
    ("Gr", a.gr, "Grashof number g beta dT a^3/nu^2 (beta = dT = 1, so |g| = Gr)"),
    ("gravityX", -a.gr, "g = (gravityX 0 0)"),
    ("Pr", a.pr, "Prandtl number nu/kappa"),
    ("height", H, "slot height H/a"),
    ("deltaT", dt, "0.25/Ha^2 (explicit Lorentz force)"),
    ("endTime", end, ""),
    ("writeInterval", end/10, "the last two writes are kept for the steadiness check"),
    ("cellsHalf", float(a.cells), ""),
    ("heightCells", float(nx), ""),
    ("layerCells", float(a.layer_cells), ""),
], ["Buoyant slot parameters. Regenerate with: python3 setup.py --ha <Ha> --gr <Gr>"])

(HERE/"system"/"blockMeshDict").write_text(HEADER.format(script="setup.py") + f"""
scale 1;

// x vertical (0 < x < {H:g}), -1 < y < 1 across the slot (B along y),
// one periodic cell in z.
vertices
(
    (0  -1 0) ({H:g} -1 0) ({H:g} 1 0) (0 1 0)
    (0  -1 {L}) ({H:g} -1 {L}) ({H:g} 1 {L}) (0 1 {L})
);

blocks
(
    hex (0 1 2 3 4 5 6 7) ({nx} {ny} 1) simpleGrading ({gx} {gy} 1)
);

boundary
(
    hotWall                     // y = -1
    {{
        type wall;
        faces ((0 1 5 4));
    }}
    coldWall                    // y = +1
    {{
        type wall;
        faces ((3 7 6 2));
    }}
    endWalls                    // x = 0 (bottom) and x = H (top)
    {{
        type wall;
        faces ((0 4 7 3) (1 2 6 5));
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
print(f"Ha = {Ha:g}, Gr = {a.gr:g}, Pr = {a.pr:g}, H/a = {H:g}; {nx} x {ny} cells; "
      f"deltaT = {dt:.3g}, endTime = {end:g} ({end/dt:.0f} steps)")
