#!/usr/bin/env python3
"""Compare the Shercliff duct solution with Shercliff's (1953) exact solution.

The exact solution is evaluated as a Fourier series (common/mhdAnalytic.py,
insulating walls). It needs the latest time directory with U, PotE, C and V
(postProcess -func writeCellCentres / writeCellVolumes) and the solver log.
Exits with status 1 if an error is above --tol.
"""
import argparse
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, os.environ.get("MHD_COMMON", str(HERE.parents[1]/"common")))
from foamCase import read_parameters  # noqa: E402
from mhdAnalytic import DuctSolution  # noqa: E402
from ductValidation import compare  # noqa: E402

p = argparse.ArgumentParser(description=__doc__)
p.add_argument("--log", default=str(HERE/"log.epotFoam_isoTh"))
p.add_argument("--tol", type=float, default=0.02, help="relative tolerance")
args = p.parse_args()

Ha = read_parameters(HERE/"system"/"caseParameters")["Ha"]
sol = DuctSolution(Ha, wall="insulating")
ok = compare(HERE, sol, args.log, f"Shercliff duct, Ha = {Ha:g}", args.tol)
print("\nPASS" if ok else "\nFAIL")
sys.exit(0 if ok else 1)
