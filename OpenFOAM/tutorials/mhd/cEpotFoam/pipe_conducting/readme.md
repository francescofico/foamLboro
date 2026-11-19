# cEpotFoam: circular pipe, conducting wall

Validation case. This is the same pipe as `epotFoam_isoTh/pipe_insulating`, but the wall **conducts** and is resolved as a solid annulus, a < r < a + t_w. Current closes through the wall, which raises the pressure gradient by about an order of magnitude at Ha = 100. The case tests the conjugate potential coupling on a curved, non-orthogonal interface.

## Run

```sh
./Allrun                                   # Ha = 100, t_w/a = 0.1, sigma_w/sigma = 1, about 2 min
python3 setup.py --thickness 0.02 --c 0.1 && ./Allrun    # thinner wall, same c
./Allclean
```

`validate.py` exits non-zero if an error is above 2 % (`--tol`) or the run is not steady.

## Problem definition

| Item | Value |
|---|---|
| Fluid | Pipe r < a, O-grid as in pipe_insulating. One periodic cell in x |
| Solid | Annulus a < r < a + t_w, four blocks whose faces match the fluid at r = a. Conductivity σ_w |
| Outer wall surface | Insulating |
| Magnetic field | B0 = (0, B, 0), uniform |
| Driving | ū = 1 held by `meanVelocityForce` |
| Units | Non-dimensional: a = ρ = ν = σ = ū = 1, so B = Ha |
| Default | Ha = 100, t_w/a = 0.1, σ_w/σ = 1, so c = σ_w t_w/(σ a) = 0.1 |
| Coupling | `coupled true`: assembled potential matrix for fluid and solid (see `cEpotFoam/hunt_duct`) |
| Numerics | `leastSquares` gradients and the same time-step cap as pipe_insulating (see its readme) |

## Reference solution
The same spectral solver as pipe_insulating (`../../common/pipeAnalytic.py`), with the wall treated exactly. Laplace's equation for b in the annular wall, with b and (1/σ) ∂b/∂n continuous at r = a and b = 0 outside, gives a Robin condition for each Fourier mode:

```
∂b_m/∂r + κ_m b_m = 0 at r = a,   κ_m = (σ/σ_w) |m| coth(|m| ln(1 + t_w/a)),   κ_0 = (σ/σ_w)/ln(1 + t_w/a)
```

As t_w → 0 at fixed c this tends to the thin-wall condition ∂b/∂r + b/c = 0. The solver reproduces that limit (G = 1040.06 against 1040.10 at Ha = 100, c = 0.1). `validate.py` also prints, for context:
- the thin-wall value at the same c
- the core-flow estimate G ≈ Ha² c/(1 + c) (Chang & Lundgren 1961), valid for c Ha ≫ 1

## Results (v2606, default mesh)

**G** is the mean streamwise pressure gradient −dp/dx needed to drive the bulk velocity ū. It is the value the `meanVelocityForce` fvOption logs as `pressure gradient = ...`. The tables give it in the non-dimensional form G a²/(ρ ν ū), with a = the pipe radius. With no magnetic field it would be 8, so the rise above that is the MHD drag.

| Ha | t_w/a | σ_w/σ | G: OpenFOAM | G: reference | error in G | max error in u/ū | run time |
|---|---|---|---|---|---|---|---|
| 100 | 0.1 | 1 | 997.88 | 998.63 | 0.07 % | 1.45 % | 1.7 min |

For comparison at the same c = 0.1: the thin-wall value is 1040.1 (+4.2 %) and the core-flow estimate is 909.1.

## Files

| File | Purpose |
|---|---|
| `setup.py` | Writes both `blockMeshDict`s and `system/caseParameters` (`--ha`, `--thickness`, `--c` or `--sigma-ratio`) |
| `validate.py` | Comparison with the reference solution |

## References
- J.A. Shercliff (1956). *The flow of conducting fluids in circular pipes under transverse magnetic fields.* J. Fluid Mech. 1(6), 644–666.
- C.C. Chang, T.S. Lundgren (1961). *Duct flow in magnetohydrodynamics.* Z. Angew. Math. Phys. 12, 100–114.
- S. Vantieghem, X. Albets-Chico, B. Knaepen (2009). *The velocity profile of laminar MHD flows in circular conducting pipes.* Theor. Comput. Fluid Dyn. 23 (numerical data for conducting walls up to Ha ≈ 3000; useful for extending this case to higher Ha).
