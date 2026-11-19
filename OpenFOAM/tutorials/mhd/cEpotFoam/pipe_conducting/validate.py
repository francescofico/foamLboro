#!/usr/bin/env python3
"""Compare the conducting pipe with the exact solution for a wall of finite thickness.

The reference is a spectral solution of the induced-field equations
(common/pipeAnalytic.py). The wall enters exactly, mode by mode: Laplace's
equation is solved for b in the annular wall. The thin-wall value at the
same c and the Chang & Lundgren (1961) core-flow estimate are printed for
context. Exits with status 1 if an error is above --tol.
"""
import argparse
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, os.environ.get("MHD_COMMON", str(HERE.parents[1]/"common")))
import pipeCase  # noqa: E402
from foamCase import read_parameters  # noqa: E402

p = argparse.ArgumentParser(description=__doc__)
p.add_argument("--log", default=str(HERE/"log.cEpotFoam"))
p.add_argument("--tol", type=float, default=0.02, help="relative tolerance")
args = p.parse_args()
par = read_parameters(HERE/"system"/"caseParameters")
title = (f"Conducting pipe, Ha = {par['Ha']:g}, t_w/a = {par['wallThickness']:g}, "
         f"sigma_w/sigma = {par['sigmaRatio']:g}")
pipeCase.validate(HERE, title, args.tol, args.log, conducting=True)
