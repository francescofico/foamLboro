#!/usr/bin/env python3
"""Set up the conducting pipe (wall 1 < r < 1 + t_w, conductivity sigma_w).

Writes system/caseParameters and the fluid (O-grid) and solid (annulus)
blockMeshDicts. Their interface faces match, for the assembled potential solve.

    python3 setup.py                               # reuse system/caseParameters
    python3 setup.py --ha 100 --thickness 0.1 --c 0.1
    python3 setup.py --ha 100 --thickness 0.02 --sigma-ratio 5
"""
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, os.environ.get("MHD_COMMON", str(HERE.parents[1]/"common")))
import pipeCase  # noqa: E402

DEFAULTS = dict(Ha=100.0, wallThickness=0.1, sigmaRatio=1.0,     # c = 0.1
                cellsQuarter=40, radialCells=40, layerCells=6.0, wallCells=10)
args = pipeCase.parse(HERE, __doc__, DEFAULTS, conducting=True)
pipeCase.write(HERE, args, "Conducting pipe")
