#!/usr/bin/env python3
"""Set up MHD mixed convection in a vertical duct with thick conducting Hartmann
walls and insulating side walls (Zheng, Chen & Ni 2025, case II).

Defaults follow the paper's case II: sigma_w = sigma, Hartmann wall thickness
delta_hw = 0.05 l, i.e. t_w/a = 0.1. Writes system/caseParameters, both
blockMeshDicts, 0/fluid/T, 0/solid/T and constant/fluid/fvOptions.

    python3 setup.py                                   # reuse system/caseParameters
    python3 setup.py --gr-re 1 --ha 50 --thickness 0.1 --sigma-ratio 1
"""
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, os.environ.get("MHD_COMMON", str(HERE.parents[1]/"common")))
import zhengDuct  # noqa: E402

DEFAULTS = dict(Ha=50.0, grRe=5.0, Pr=0.01, cellsHalf=40, layerCells=6.0, sideCells=10.0,
                wallThickness=0.1, sigmaRatio=1.0, kRatio=1.0, wallCells=10)
args = zhengDuct.parse(HERE, __doc__, DEFAULTS, conducting=True)
zhengDuct.write(HERE, args, "Mixed convection (conducting Hartmann walls)")
