#!/usr/bin/env python3
"""Compare the thick-wall duct with the exact solution for walls of finite thickness.

The reference solves Laplace's equation for the potential inside the wall
exactly, mode by mode (common/mhdAnalytic.py, wall="thick"). That is the same
electromagnetically coupled problem as Bluck & Wolfendale (2015). Fluid
velocity, fluid potential and the potential inside the wall are all checked.
Hunt's thin-wall value at the same c is printed to show the thickness effect.
Exits with status 1 if an error is above --tol.
"""
import argparse
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, os.environ.get("MHD_COMMON", str(HERE.parents[1]/"common")))
import conductingDuct  # noqa: E402
from foamCase import read_parameters  # noqa: E402

p = argparse.ArgumentParser(description=__doc__)
p.add_argument("--log", default=str(HERE/"log.cEpotFoam"))
p.add_argument("--tol", type=float, default=0.02, help="relative tolerance")
args = p.parse_args()

par = read_parameters(HERE/"system"/"caseParameters")
title = (f"Thick-wall duct, Ha = {par['Ha']:g}, t_w/a = {par['wallThickness']:g}, "
         f"sigma_w/sigma = {par['sigmaRatio']:g}")
conductingDuct.validate(HERE, title, args.tol, args.log, primary="thick")
