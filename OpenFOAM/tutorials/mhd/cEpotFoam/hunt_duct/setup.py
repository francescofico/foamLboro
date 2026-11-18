#!/usr/bin/env python3
"""Set up the Hunt duct (thin conducting Hartmann walls, insulating side walls).

Writes system/caseParameters and the fluid and solid blockMeshDicts. The
Hartmann walls are resolved solid slabs of thickness t_w. Keep t_w/a small so
that only the conductance ratio c = sigma_w t_w/(sigma a) matters, as in
Hunt's thin-wall model. validate.py reports the remaining thickness effect.

    python3 setup.py                         # reuse system/caseParameters
    python3 setup.py --ha 500 --c 0.05
    python3 setup.py --ha 100 --c 0.1 --thickness 0.005
"""
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, os.environ.get("MHD_COMMON", str(HERE.parents[1]/"common")))
import conductingDuct  # noqa: E402

DEFAULTS = dict(Ha=100.0, wallThickness=0.01, sigmaRatio=10.0,   # c = 0.1
                cellsHalf=40, wallCells=4, layerCells=6.0, sideCells=10.0)

args = conductingDuct.parse(HERE, __doc__, DEFAULTS)
conductingDuct.write(HERE, args, "Hunt duct")
