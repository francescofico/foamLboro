# epotFoam_isoTh: Hartmann channel flow

Validation case. Steady, fully developed laminar flow between two infinite parallel plates in a uniform magnetic field normal to the plates (Hartmann 1937). The solution is exact and closed-form. The case checks:

- the Lorentz-force coupling in the momentum equation
- the conservative current/force interpolation (Ni et al. 2007)
- how well the mesh resolves the Hartmann layers (thickness ~ a/Ha)

## Run

```sh
./Allrun                          # Ha = 50 (system/caseParameters)
HA=100 ./Allrun                   # other Hartmann number
HA=50 ELECTRODES=shorted ./Allrun # short-circuited variant
./Allclean
```

`Allrun` runs `setup.py`, `blockMesh`, `epotFoam_isoTh`, `postProcess -func writeCellCentres` and `validate.py`. `validate.py` prints the errors, writes `validation.png`, and exits non-zero if an error is above 1 % (`--tol`). The run takes about 1 s.

## Problem definition

| Item | Value |
|---|---|
| Geometry | Plates at y = ±a. One periodic cell in x (streamwise) and in z (spanwise) |
| Magnetic field | B0 = (0, B, 0), uniform |
| Walls (y = ±a) | No-slip, electrically insulating (∂φ/∂n = 0) |
| Driving | Bulk velocity ū = 1 held by the `meanVelocityForce` fvOption (`constant/fvOptions`) |
| Units | Non-dimensional: a = ρ = ν = σ = ū = 1, so B = Ha |
| Mesh | 80 cells across, graded so the first cell is (1/Ha)/6 (`setup.py --cells --layer-cells`) |
| Time step | 0.25/Ha², because the Lorentz force is explicit; 400 steps |

Fully developed parallel flow has no inertia, so the result does not depend on Re.

### Spanwise electrical condition
Periodic z boundaries need care. A plain `cyclic` for φ forces ⟨∂φ/∂z⟩ = 0. That is the **short-circuited** channel (j_z = σuB). In the **insulating** channel the net spanwise current is zero, which needs ∂φ/∂z = ūB. The case imposes this as a `fixedJump` of ūB·L_z on the z cyclic pair (`0/PotE`, `potentialJump` in `system/caseParameters`). `ELECTRODES=shorted` sets the jump to 0. The velocity profile shape is the same for both, but the pressure gradient differs by about Ha².

## Analytical solution
With η = y/a:

```
u/ū = (cosh Ha − cosh Ha η) / (cosh Ha − sinh Ha / Ha)

G a²/(ρ ν ū) = Ha² / (Ha coth Ha − 1)      insulating
             = Ha² / (1 − tanh Ha / Ha)     short-circuited
```

## Results (v2606, default mesh)

| Ha | electrodes | G: OpenFOAM | G: exact | error in G | max error in u/ū |
|---|---|---|---|---|---|
| 10 | insulating | 11.064 | 11.111 | 0.42 % | 0.31 % |
| 50 | insulating | 50.748 | 51.020 | 0.53 % | 0.30 % |
| 100 | insulating | 100.41 | 101.01 | 0.59 % | 0.29 % |
| 50 | shorted | 2550.75 | 2551.02 | 0.011 % | 0.30 % |

The error in G falls as the Hartmann layer is refined. At Ha = 50 it is 1.6 % with `--layer-cells 3`, 0.53 % with 6 and 0.26 % with 12. In the insulating case G is a small difference between large Lorentz terms, so it is the most sensitive quantity.

## Files

| File | Purpose |
|---|---|
| `setup.py` | Writes `system/blockMeshDict` and `system/caseParameters` from Ha |
| `system/caseParameters` | Ha, B, potential jump, deltaT, endTime. Included by `0/B0`, `0/PotE` and `controlDict` |
| `validate.py` | Comparison with the analytical solution (`../../common/mhdAnalytic.py`) |

## References
- J. Hartmann (1937). *Hg-dynamics I: Theory of the laminar flow of an electrically conductive liquid in a homogeneous magnetic field.* Kgl. Danske Videnskab. Selskab, Mat.-Fys. Medd. 15(6).
- M.-J. Ni, R. Munipalli, N.B. Morley, P. Huang, M.A. Abdou (2007). *A current density conservative scheme for incompressible MHD flows at a low magnetic Reynolds number. Part I.* J. Comput. Phys. 227, 174–204.
