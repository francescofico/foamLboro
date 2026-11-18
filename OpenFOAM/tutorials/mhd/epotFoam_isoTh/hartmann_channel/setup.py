#!/usr/bin/env python3
"""Set up the Hartmann channel for a given Hartmann number.

Writes system/blockMeshDict and system/caseParameters. Every Ha-dependent
value (B, time step, end time, potential jump) is derived here, so the
dictionaries themselves never need editing.

    python3 setup.py                     # reuse Ha from system/caseParameters
    python3 setup.py --ha 100
    python3 setup.py --ha 100 --electrodes shorted
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
p.add_argument("--ha", type=float, default=old.get("Ha", 50.0), help="Hartmann number")
p.add_argument("--electrodes", choices=("insulating", "shorted"),
               default=old.get("electrodes", "insulating"),
               help="spanwise electrical condition (see 0/PotE)")
p.add_argument("--cells", type=int, default=int(old.get("cellsHalf", 40)),
               help="cells per half channel")
p.add_argument("--layer-cells", type=float, default=old.get("layerCells", 6.0),
               help="first-cell size = (1/Ha)/layer-cells")
args = p.parse_args()

Ha = args.ha
L = 0.05                                    # x and z extent: one periodic cell
ubar = 1.0
dt = 0.25/Ha**2                             # explicit Lorentz force: dt < rho/(sigma B^2)
end = 400*dt                                # 100 magnetic damping times rho/(sigma B^2)
jump = ubar*Ha*L if args.electrodes == "insulating" else 0.0

ny, gy = symmetric_grading(args.cells, 1.0/(Ha*args.layer_cells))

write_parameters(PARAMS, "setup.py", [
    ("Ha", Ha, "Hartmann number B a sqrt(sigma/(rho nu)); a = rho = nu = sigma = 1"),
    ("B", Ha, "applied field magnitude (= Ha in these units)"),
    ("electrodes", args.electrodes, "insulating | shorted"),
    ("potentialJump", jump, "ubar B Lz for insulating, 0 for shorted"),
    ("deltaT", dt, "0.25/Ha^2"),
    ("endTime", end, ""),
    ("cellsHalf", float(args.cells), ""),
    ("layerCells", float(args.layer_cells), ""),
], ["Hartmann channel parameters. Regenerate with: python3 setup.py --ha <Ha>"])

(HERE/"system"/"blockMeshDict").write_text(HEADER.format(script="setup.py") + f"""
scale 1;

// -1 < y < 1 graded towards both Hartmann walls; one periodic cell in x and z.
vertices
(
    (0  -1 0) ({L} -1 0) ({L} 1 0) (0 1 0)
    (0  -1 {L}) ({L} -1 {L}) ({L} 1 {L}) (0 1 {L})
);

blocks
(
    hex (0 1 2 3 4 5 6 7) (1 {ny} 1) simpleGrading (1 {gy} 1)
);

boundary
(
    hartmannWalls
    {{
        type wall;
        faces ((0 1 5 4) (3 7 6 2));
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
print(f"Ha = {Ha:g} ({args.electrodes}), {ny} cells across, deltaT = {dt:.3g}, endTime = {end:.3g}")
