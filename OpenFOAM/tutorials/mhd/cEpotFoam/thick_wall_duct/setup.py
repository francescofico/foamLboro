#!/usr/bin/env python3
"""Set up the thick-wall duct (thick conducting Hartmann walls, insulating side walls).

Writes system/caseParameters and the fluid and solid blockMeshDicts. The
Hartmann walls are solid slabs of finite thickness t_w and conductivity
sigma_w. The default (t_w/a = 0.5, sigma_w/sigma = 0.2) has the same thin-wall
conductance ratio c = 0.1 as hunt_duct, so the two cases differ only through
the wall thickness.

    python3 setup.py                               # reuse system/caseParameters
    python3 setup.py --ha 100 --thickness 0.2 --c 0.1
    python3 setup.py --ha 100 --thickness 0.5 --sigma-ratio 1
"""
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, os.environ.get("MHD_COMMON", str(HERE.parents[1]/"common")))
import conductingDuct  # noqa: E402

DEFAULTS = dict(Ha=100.0, wallThickness=0.5, sigmaRatio=0.2,     # c = 0.1
                cellsHalf=40, wallCells=20, layerCells=6.0, sideCells=10.0)

args = conductingDuct.parse(HERE, __doc__, DEFAULTS)
conductingDuct.write(HERE, args, "Thick-wall duct")
