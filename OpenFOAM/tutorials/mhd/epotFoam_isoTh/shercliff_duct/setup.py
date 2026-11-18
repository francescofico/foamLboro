#!/usr/bin/env python3
"""Set up the Shercliff duct for a given Hartmann number.

Writes system/blockMeshDict and system/caseParameters. The cross-section
(-1 < y, z < 1) is graded towards all four walls:
  y (along B): first cell = (1/Ha)/layer-cells, to resolve the Hartmann layers
  z (normal to B): first cell = (1/sqrt(Ha))/side-cells, to resolve the side layers
The streamwise direction has one periodic cell, so the flow is fully developed.

    python3 setup.py                 # reuse Ha from system/caseParameters
    python3 setup.py --ha 500 --cells 60
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
p.add_argument("--ha", type=float, default=old.get("Ha", 100.0), help="Hartmann number")
p.add_argument("--cells", type=int, default=int(old.get("cellsHalf", 40)),
               help="cells per half-width, in both y and z")
p.add_argument("--layer-cells", type=float, default=old.get("layerCells", 6.0),
               help="first y cell = (1/Ha)/layer-cells")
p.add_argument("--side-cells", type=float, default=old.get("sideCells", 10.0),
               help="first z cell = (1/sqrt(Ha))/side-cells")
args = p.parse_args()

Ha = args.ha
L = 0.05
dt = 0.25/Ha**2                       # explicit Lorentz force: dt < rho/(sigma B^2)
end = 5.0/Ha                          # steady after ~3 side-layer times a^2/(nu Ha)
ny, gy = symmetric_grading(args.cells, 1.0/(Ha*args.layer_cells))
nz, gz = symmetric_grading(args.cells, 1.0/(Ha**0.5*args.side_cells))

write_parameters(PARAMS, "setup.py", [
    ("Ha", Ha, "Hartmann number B a sqrt(sigma/(rho nu)); a = rho = nu = sigma = 1"),
    ("B", Ha, "applied field magnitude (= Ha in these units)"),
    ("deltaT", dt, "0.25/Ha^2"),
    ("endTime", end, "5/Ha"),
    ("cellsHalf", float(args.cells), ""),
    ("layerCells", float(args.layer_cells), ""),
    ("sideCells", float(args.side_cells), ""),
], ["Shercliff duct parameters. Regenerate with: python3 setup.py --ha <Ha>"])

(HERE/"system"/"blockMeshDict").write_text(HEADER.format(script="setup.py") + f"""
scale 1;

// Square duct -1 < y, z < 1 (B along y), one periodic cell in x.
vertices
(
    (0  -1 -1) ({L} -1 -1) ({L} 1 -1) (0 1 -1)
    (0  -1  1) ({L} -1  1) ({L} 1  1) (0 1  1)
);

blocks
(
    hex (0 1 2 3 4 5 6 7) (1 {ny} {nz}) simpleGrading (1 {gy} {gz})
);

boundary
(
    hartmannWalls               // y = -1, +1  (normal to B)
    {{
        type wall;
        faces ((0 1 5 4) (3 7 6 2));
    }}
    sideWalls                   // z = -1, +1  (parallel to B)
    {{
        type wall;
        faces ((0 3 2 1) (4 5 6 7));
    }}
    xPeriodic0
    {{
        type cyclic;
        neighbourPatch xPeriodic1;
        faces ((0 4 7 3));
    }}
    xPeriodic1
    {{
        type cyclic;
        neighbourPatch xPeriodic0;
        faces ((1 2 6 5));
    }}
);

// ************************************************************************* //
""")
print(f"Ha = {Ha:g}, {ny} x {nz} cells, deltaT = {dt:.3g}, endTime = {end:.3g} ({end/dt:.0f} steps)")
