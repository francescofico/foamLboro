#!/usr/bin/env python3
"""Compare the Hunt duct solution with Hunt's (1965) thin-wall solution.

The reference is Hunt's Fourier-series solution for conductance ratio
c = sigma_w t_w/(sigma a) (common/mhdAnalytic.py). The script also prints the
exact solution for the actual wall thickness, which shows how much of the
difference comes from the thin-wall model. Exits with status 1 if an error
is above --tol.
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
title = f"Hunt duct, Ha = {par['Ha']:g}, c = {par['c']:g} (t_w/a = {par['wallThickness']:g})"
conductingDuct.validate(HERE, title, args.tol, args.log, primary="thin")
