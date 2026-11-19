# epotFoam_isoTh: circular pipe, insulating wall

Validation case. Fully developed laminar flow in a circular pipe with an **electrically insulating wall**, in a uniform transverse magnetic field (Gold 1962; asymptotic theory by Shercliff 1956). Compared with the square-duct cases it adds:

- a curved wall, meshed as an **O-grid**: non-orthogonal, with skewed cells at the corners of the central block
- Hartmann layers whose thickness varies around the wall, and Roberts layers where the wall is tangent to B (z = ±a)

## Run

```sh
./Allrun              # Ha = 100 (system/caseParameters), about 2 min
HA=20 ./Allrun        # about 5 min (time step limited by the mesh, see below)
./Allclean
```

`validate.py` prints the errors, writes `validation.png`, and exits non-zero if an error is above 2 % (`--tol`) or the run is not steady.

## Problem definition

| Item | Value |
|---|---|
| Geometry | Pipe of radius a, axis x. One periodic cell in x |
| Magnetic field | B0 = (0, B, 0), uniform |
| Wall (r = a) | No-slip, insulating (∂φ/∂n = 0) |
| Driving | ū = 1 held by `meanVelocityForce` |
| Units | Non-dimensional: a = ρ = ν = σ = ū = 1, so B = Ha (radius-based) |
| Mesh | O-grid with 8000 cells: central square (corners at ±0.35a) and four outer blocks. First wall cell (1/Ha)/6 (`setup.py --cells --radial-cells --layer-cells`) |
| Time step | min(0.25/Ha², 0.15 Δ²/ν), Δ = central-block cell size. endTime = 5/Ha |

### Numerical settings that matter on this mesh
- **`gradSchemes leastSquares`.** With `Gauss linear` the velocity had local errors of up to 13 % at the O-grid corners, spreading along B, although G was still within 0.5 %. The cell gradient enters the non-orthogonal correction of the face current, and in the core that error is amplified by about Ha, because −∇φ and u×B nearly cancel there. With `leastSquares` the maximum error is 1.4 %. Skew-corrected interpolation did not help.
- **Time-step cap.** With the Lorentz-force limit alone (0.25/Ha²), low-Ha runs developed a growing spurious cross-plane flow (Ha = 50) or diverged (Ha = 20). The finer mesh (`--cells 80`) diverged even at Ha = 100. All of these fit a diffusion-type limit of about 0.16 Δ²/ν in the skewed corner cells, which points to the explicit viscous terms there. The cap of 0.15 Δ²/ν is empirical.

## Reference solution
For an insulating wall the problem has an exact solution (Gold 1962). `../../common/pipeAnalytic.py` solves the induced-field equations

```
∇²u + Ha ∂b/∂y + 1 = 0,     ∇²b + Ha ∂u/∂y = 0,     u = b = 0 at r = a
```

with a Chebyshev–Fourier spectral method, converged to about 10 digits. It is checked against:
- a closed-form Bessel-series solution of the same problem: W = u + b solves ∇²W + Ha ∂W/∂y = −1. They agree to round-off at Ha = 5, 10 and 20; above Ha ≈ 20 the Bessel series loses precision.
- Poiseuille flow at Ha → 0: G = 8 exactly.
- the high-Ha trend: G/Ha = 1.2267 (Ha 100) and 1.1883 (Ha 500), approaching 3π/8 ≈ 1.178 with a 1/Ha correction.

## Results (v2606, default mesh)

| Ha | G: OpenFOAM | G: reference | error in G | max error in u/ū | L2 error | run time |
|---|---|---|---|---|---|---|
| 20 | 28.218 | 28.291 | 0.26 % | 0.12 % | 0.07 % | 4.7 min |
| 50 | 63.450 | 63.674 | 0.35 % | 0.58 % | 0.11 % | 1.9 min |
| 100 | 122.158 | 122.671 | 0.42 % | 1.45 % | 0.15 % | 1.7 min |

The potential is not compared for the pipe; the reference is in induced-field form.

## Files

| File | Purpose |
|---|---|
| `setup.py` | Writes the O-grid `system/blockMeshDict` and `system/caseParameters` (`../../common/pipeCase.py`) |
| `validate.py` | Comparison with the reference solution |

## References
- R.R. Gold (1962). *Magnetohydrodynamic pipe flow. Part 1.* J. Fluid Mech. 13(4), 505–512.
- J.A. Shercliff (1956). *The flow of conducting fluids in circular pipes under transverse magnetic fields.* J. Fluid Mech. 1(6), 644–666.
- M.-J. Ni et al. (2007). *A current density conservative scheme ... Part II: on an arbitrary collocated mesh.* J. Comput. Phys. 227, 205ff.
