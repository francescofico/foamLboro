#!/usr/bin/env python3
"""Compare the insulating pipe with the exact solution (the Gold 1962 problem).

The reference is a spectral solution of the induced-field equations
(common/pipeAnalytic.py). It is converged to about 10 digits, and agrees to
round-off with a closed-form Bessel-series solution of Gold's problem for
Ha <= 20. Exits with status 1 if an error is above --tol.
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
p.add_argument("--log", default=str(HERE/"log.epotFoam_isoTh"))
p.add_argument("--tol", type=float, default=0.02, help="relative tolerance")
args = p.parse_args()
Ha = read_parameters(HERE/"system"/"caseParameters")["Ha"]
pipeCase.validate(HERE, f"Insulating pipe, Ha = {Ha:g}", args.tol, args.log, conducting=False)
