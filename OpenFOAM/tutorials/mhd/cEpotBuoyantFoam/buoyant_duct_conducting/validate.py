#!/usr/bin/env python3
"""Compare the fully developed buoyant duct with the exact series solution
(common/mhdAnalytic.py, BuoyantDuctSolution): velocity, potential and
temperature, and the potential and temperature inside the walls.
Exits with status 1 if an error is above --tol."""
import argparse
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, os.environ.get("MHD_COMMON", str(HERE.parents[1]/"common")))
import buoyantDuct  # noqa: E402
from foamCase import read_parameters  # noqa: E402

p = argparse.ArgumentParser(description=__doc__)
p.add_argument("--tol", type=float, default=0.02, help="relative tolerance")
args = p.parse_args()
par = read_parameters(HERE/"system"/"caseParameters")
title = (f"Buoyant duct (conducting), Ha = {par['Ha']:g}, Gr = {par['Gr']:g}, "
         f"t_w/a = {par['wallThickness']:g}, sigma_w/sigma = {par['sigmaRatio']:g}, k_w/k = {par['kRatio']:g}")
buoyantDuct.validate(HERE, title, args.tol, conducting=True)
