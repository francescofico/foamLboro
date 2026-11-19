#!/usr/bin/env python3
"""Set up the insulating pipe for a given Hartmann number.

Writes system/blockMeshDict (O-grid) and system/caseParameters. The radial
grading puts the first wall cell at (1/Ha)/layer-cells, to resolve the
Hartmann layers.

    python3 setup.py                     # reuse Ha from system/caseParameters
    python3 setup.py --ha 300 --cells 60 --radial-cells 50
"""
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, os.environ.get("MHD_COMMON", str(HERE.parents[1]/"common")))
import pipeCase  # noqa: E402

DEFAULTS = dict(Ha=100.0, cellsQuarter=40, radialCells=40, layerCells=6.0)
args = pipeCase.parse(HERE, __doc__, DEFAULTS, conducting=False)
pipeCase.write(HERE, args, "Insulating pipe")
