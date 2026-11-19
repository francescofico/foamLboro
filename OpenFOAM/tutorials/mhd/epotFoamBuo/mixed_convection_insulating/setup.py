#!/usr/bin/env python3
"""Set up MHD mixed convection in a vertical duct with insulating walls
(Zheng, Chen & Ni 2025, case I).

Side walls z = -a (hot) and z = +a (cold), B along y, ubar = 1 held by
meanVelocityForce, buoyancy (Gr/Re) G (T - TRef) added as a coded body force.
Writes system/blockMeshDict, system/caseParameters, 0/T, constant/fvOptions.

    python3 setup.py                         # reuse system/caseParameters
    python3 setup.py --gr-re 1 --ha 50       # Ha is half-width based: Ha_paper = 2 Ha
"""
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, os.environ.get("MHD_COMMON", str(HERE.parents[1]/"common")))
import zhengDuct  # noqa: E402

DEFAULTS = dict(Ha=50.0, grRe=5.0, Pr=0.01, cellsHalf=40, layerCells=6.0, sideCells=10.0)
args = zhengDuct.parse(HERE, __doc__, DEFAULTS, conducting=False)
zhengDuct.write(HERE, args, "Mixed convection (insulating)")
