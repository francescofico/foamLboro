#!/usr/bin/env python3
"""Set up the buoyant vertical duct with insulating walls.

Writes system/blockMeshDict, system/caseParameters, 0/T (the heated walls depend
on the orientation) and constant/fvOptions (the Boussinesq body force).

    python3 setup.py                         # reuse system/caseParameters
    python3 setup.py --heating parallel      # grad T along B (Hartmann walls heated)
    python3 setup.py --heating normal --ha 100 --gr 20000
"""
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, os.environ.get("MHD_COMMON", str(HERE.parents[1]/"common")))
import buoyantDuct  # noqa: E402

DEFAULTS = dict(Ha=100.0, Gr=20000.0, Pr=0.01, heating="normal",
                cellsHalf=40, layerCells=6.0, sideCells=10.0)
args = buoyantDuct.parse(HERE, __doc__, DEFAULTS, conducting=False)
buoyantDuct.write(HERE, args, "Buoyant duct (insulating)")
