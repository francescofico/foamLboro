#!/usr/bin/env python3
"""Set up the buoyant vertical duct with conducting Hartmann walls (solid slabs
1 < |y| < 1 + t_w, as in cEpotFoam/hunt_duct).

Writes system/caseParameters, both blockMeshDicts, 0/fluid/T and 0/solid/T (the
heated faces depend on the orientation) and constant/fluid/fvOptions (the
Boussinesq body force).

    python3 setup.py                         # reuse system/caseParameters
    python3 setup.py --heating parallel      # grad T along B: heat crosses the walls
    python3 setup.py --heating normal --thickness 0.1 --sigma-ratio 1 --k-ratio 2
"""
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, os.environ.get("MHD_COMMON", str(HERE.parents[1]/"common")))
import buoyantDuct  # noqa: E402

DEFAULTS = dict(Ha=100.0, Gr=20000.0, Pr=0.01, heating="normal",
                cellsHalf=40, layerCells=6.0, sideCells=10.0,
                wallThickness=0.1, sigmaRatio=1.0, kRatio=2.0, wallCells=10)
args = buoyantDuct.parse(HERE, __doc__, DEFAULTS, conducting=True)
buoyantDuct.write(HERE, args, "Buoyant duct (conducting)")
