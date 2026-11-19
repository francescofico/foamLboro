# cEpotFoam: Hunt flow (conducting Hartmann walls, insulating side walls)

Validation case. This is the same square duct as `epotFoam_isoTh/shercliff_duct`, but the **Hartmann walls conduct** (Hunt 1965). Current closes through the walls, which produces **side-layer jets**: velocity overshoots near the side walls. It is the second fully developed case of the Smolentsev et al. (2015) V&V database. The case tests:

- the fluid/solid potential coupling (`solidWallMixedPotentialCoupled` over `mappedWall`)
- side-layer resolution
- how the solution depends on the wall conductance ratio c

## Run

```sh
./Allrun                  # Ha = 100, c = 0.1 (system/caseParameters), about 2 min
HA=500 C=0.05 ./Allrun    # C sets sigma_w at the current wall thickness
./Allclean
```

`Allrun` runs `setup.py`, `blockMesh` for both regions, `cEpotFoam`, `postProcess` for both regions and `validate.py`. `validate.py` exits non-zero if an error is above 2 % (`--tol`) or the run is not steady.

## Problem definition

| Item | Value |
|---|---|
| Fluid | Square duct −a < y, z < a. One periodic cell in x |
| Solid | Two slabs a < \|y\| < a + t_w on the Hartmann walls only, conductivity σ_w |
| Magnetic field | B0 = (0, B, 0), uniform |
| Side walls (z = ±a) | No-slip, insulating, no solid region |
| Outer slab faces and slab ends | Insulating |
| Driving | ū = 1 held by `meanVelocityForce` (`constant/fluid/fvOptions`) |
| Units | Non-dimensional: a = ρ = ν = σ = ū = 1, so B = Ha |
| Default | Ha = 100, t_w/a = 0.01, σ_w/σ = 10, so c = σ_w t_w/(σ a) = 0.1 |
| Mesh | Fluid 80 × 80 as in shercliff_duct; 4 cells across each wall; matching faces at the interface |
| Coupling | `coupled true` in `system/fvSolution`: one assembled potential matrix for both regions |

**Thin-wall approximation.** In Hunt's solution only the product σ_w t_w (that is, c) enters. This case resolves a thin solid wall instead. `validate.py` also evaluates the exact solution for the actual t_w, so the thin-wall error is reported separately from the discretisation error. With the defaults it is 0.02 % in G.

**Coupling choice.** The segregated mapped coupling (`coupled false`) is also available, but here the interface converges very slowly. After 2000 steps it was still 38 % high in G and drifting. The assembled solve reaches steady state in about 3/Ha.

## Analytical solution
Hunt's Fourier-series solution for thin walls of conductance ratio c, evaluated in `../../common/mhdAnalytic.py` (`DuctSolution(Ha, wall="thin", c=c)`). It is checked against an independent 2D finite-difference solution, which converges to it at second order.

## Results (v2606, default mesh)

**G** is the mean streamwise pressure gradient −dp/dx needed to drive the bulk velocity ū. It is the value the `meanVelocityForce` fvOption logs as `pressure gradient = ...`. The tables give it in the non-dimensional form G a²/(ρ ν ū), with a = the duct half-width. With no magnetic field it would be 7.11, so the rise above that is the MHD drag.

| Ha | c | G: OpenFOAM | G: Hunt | error in G | jet peak: OpenFOAM / Hunt | max error in u/ū | max error in wall φ | run time |
|---|---|---|---|---|---|---|---|---|
| 100 | 0.1 | 821.01 | 822.53 | 0.18 % | 2.466 / 2.460 | 0.39 % | 0.09 % | 107 s |
| 500 | 0.05 | 9518.2 | 9535.1 | 0.18 % | 5.082 / 5.056 | 0.52 % | 0.10 % | 9 min |

## Files

| File | Purpose |
|---|---|
| `setup.py` | Writes both `blockMeshDict`s and `system/caseParameters` (`--ha`, `--c` or `--sigma-ratio`, `--thickness`) |
| `validate.py` | Comparison with Hunt (`../../common/conductingDuct.py`) |

## References
- J.C.R. Hunt (1965). *Magnetohydrodynamic flow in rectangular ducts.* J. Fluid Mech. 21(4), 577–590.
- S. Smolentsev et al. (2015). *An approach to verification and validation of MHD codes for fusion applications.* Fusion Eng. Des. 100.
- FreeMHD V&V: arXiv:2409.08950 (Hunt case, c = 0.05).
- M.-J. Ni et al. (2007). J. Comput. Phys. 227, 174–204.
