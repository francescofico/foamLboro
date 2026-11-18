# epotFoam_isoTh: Shercliff flow (square duct, insulating walls)

Validation case. Fully developed laminar flow in a square duct with **all walls electrically insulating**, in a uniform transverse field (Shercliff 1953). It is one of the fully developed cases of the Smolentsev et al. (2015) fusion V&V database. Compared with the Hartmann channel it adds:

- side (Shercliff) layers on the walls parallel to B, of thickness ~ a/√Ha
- a 2D current loop that closes through the Hartmann layers
- a mesh graded strongly in two directions at once

## Run

```sh
./Allrun              # Ha = 100 (system/caseParameters), about 1.5 min
HA=500 ./Allrun       # about 8 min
./Allclean
```

`Allrun` runs `setup.py`, `blockMesh`, `epotFoam_isoTh`, `postProcess` (cell centres and volumes) and `validate.py`. `validate.py` prints the errors, writes `validation.png`, and exits non-zero if an error is above 2 % (`--tol`) or the run is not steady.

## Problem definition

| Item | Value |
|---|---|
| Geometry | Square duct −a < y, z < a. One periodic cell in x |
| Magnetic field | B0 = (0, B, 0), uniform |
| Hartmann walls (y = ±a) | No-slip, insulating |
| Side walls (z = ±a) | No-slip, insulating |
| Driving | ū = 1 held by `meanVelocityForce` |
| Units | Non-dimensional: a = ρ = ν = σ = ū = 1, so B = Ha |
| Mesh | 80 × 80. The first y cell is (1/Ha)/6 and the first z cell is (1/√Ha)/10 (`setup.py`) |
| Time step | 0.25/Ha². endTime = 5/Ha, and the solution is steady from about 3/Ha |

## Analytical solution
Shercliff's Fourier-series solution, evaluated in `../../common/mhdAnalytic.py` (`DuctSolution(Ha, wall="insulating")`) for velocity, potential and pressure gradient. The implementation is checked against:
- the hydrodynamic duct limit (exact)
- a wide duct, which tends to the insulating Hartmann channel
- an independent 2D finite-difference solution, which converges to it at second order

## Results (v2606, default mesh)

| Ha | G: OpenFOAM | G: exact | error in G | max error in u/ū | max error in φ | run time |
|---|---|---|---|---|---|---|
| 100 | 109.794 | 110.444 | 0.59 % | 0.29 % | 0.05 % | 93 s |
| 500 | 516.951 | 520.841 | 0.75 % | 0.28 % | 0.03 % | 8 min |

## Files

| File | Purpose |
|---|---|
| `setup.py` | Writes `system/blockMeshDict` and `system/caseParameters` from Ha |
| `validate.py` | Comparison with the series solution (`../../common/ductValidation.py`) |

## References
- J.A. Shercliff (1953). *Steady motion of conducting fluids in pipes under transverse magnetic fields.* Proc. Cambridge Phil. Soc. 49, 136–144.
- S. Smolentsev et al. (2015). *An approach to verification and validation of MHD codes for fusion applications.* Fusion Eng. Des. 100.
- M.-J. Ni et al. (2007). J. Comput. Phys. 227, 174–204 (Shercliff and Hunt cases at Ha = 10²–10⁴).
