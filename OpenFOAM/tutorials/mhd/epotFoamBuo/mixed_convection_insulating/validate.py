#!/usr/bin/env python3
"""Compare MHD mixed convection (insulating walls, case I of Zheng, Chen & Ni 2025) with the
exact superposition of pressure-driven and buoyancy-driven duct solutions
(common/zhengDuct.py). Also prints the paper's core-gradient scaling and
reverse-flow criterion for comparison. Exits with status 1 if an error is
above --tol."""
import argparse
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, os.environ.get("MHD_COMMON", str(HERE.parents[1]/"common")))
import zhengDuct  # noqa: E402

p = argparse.ArgumentParser(description=__doc__)
p.add_argument("--log", default=str(HERE/"log.epotFoamBuo"))
p.add_argument("--tol", type=float, default=0.02, help="relative tolerance")
args = p.parse_args()
zhengDuct.validate(HERE, "Mixed convection (insulating walls, case I)", args.tol, args.log, conducting=False)
