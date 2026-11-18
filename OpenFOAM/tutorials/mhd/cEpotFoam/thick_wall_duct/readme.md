# cEpotFoam: duct flow with thick conducting Hartmann walls

Validation case. This extends `hunt_duct` beyond the thin-wall approximation. The Hartmann walls have **finite thickness**, and the potential inside them is solved as its own region. That is the problem the conjugate `cEpotFoam` solver is built for. The case tests:

- the fluid↔solid potential coupling with real wall thickness
- the potential distribution inside the wall
- the departure from Hunt's thin-wall result as t_w/a grows

## Run

```sh
./Allrun        # Ha = 100, t_w/a = 0.5, sigma_w/sigma = 0.2 (system/caseParameters), about 2 min
python3 setup.py --ha 100 --thickness 0.2 --c 0.1 && ./Allrun
./Allclean
```

`validate.py` exits non-zero if an error is above 2 % (`--tol`) or the run is not steady.

## Problem definition

| Item | Value |
|---|---|
| Fluid | Square duct −a < y, z < a. One periodic cell in x |
| Solid | Two slabs a < \|y\| < a + t_w on the Hartmann walls, conductivity σ_w. They span −a < z < a |
| Side walls (z = ±a) | No-slip, insulating, no solid region |
| Outer slab faces and slab ends | Insulating |
| Driving | ū = 1 held by `meanVelocityForce` |
| Units | Non-dimensional: a = ρ = ν = σ = ū = 1, so B = Ha |
| Default | Ha = 100, t_w/a = 0.5, σ_w/σ = 0.2. The thin-wall c = 0.1 is the same as hunt_duct |
| Mesh | Fluid 80 × 80. 20 cells across each wall, graded towards the interface |
| Coupling | `coupled true` (see hunt_duct) |

## Analytical solution
The potential in the wall satisfies Laplace's equation, with continuity of φ and of the normal current at the interface and insulating outer faces. For each Fourier mode sin(α_k z) it can be solved exactly. The result is a Robin condition at the fluid side of the interface:

```
∂φ/∂y + (σ_w/σ) α_k tanh(α_k t_w) φ = 0      at y = a,
```

As t_w → 0 at fixed c this reduces to Hunt's thin-wall condition, with c α_k² in place of (σ_w/σ) α_k tanh(α_k t_w). The fluid problem is then solved exactly, as in the Shercliff and Hunt cases (`../../common/mhdAnalytic.py`, `DuctSolution(Ha, wall="thick", sigma_ratio, thickness)`). This is the same electromagnetically coupled problem treated by Bluck & Wolfendale (2015). The wall configuration here, with slabs on the Hartmann walls only and insulating side walls and ends, is defined by this case. Check it against their paper before comparing numbers one to one.

At the default parameters the wall thickness lowers G by **18 %** relative to Hunt's thin-wall value at the same c (673.3 against 822.5). That makes this a sharp test of the conjugate coupling.

## Results (v2606, default mesh)

| Ha | t_w/a | σ_w/σ | G: OpenFOAM | G: exact | G: Hunt (same c) | error in G | max error in u/ū | max error in wall φ | run time |
|---|---|---|---|---|---|---|---|---|---|
| 100 | 0.5 | 0.2 | 671.57 | 673.34 | 822.53 | 0.26 % | 0.31 % | 0.10 % | 112 s |

## Files

| File | Purpose |
|---|---|
| `setup.py` | Writes both `blockMeshDict`s and `system/caseParameters` |
| `validate.py` | Comparison with the exact thick-wall solution, including the potential inside the wall |

## References
- M.J. Bluck, M.J. Wolfendale (2015). *An analytical solution to electromagnetically coupled duct flow in MHD.* J. Fluid Mech. 771. doi:10.1017/jfm.2015.202
- J.C.R. Hunt (1965). *Magnetohydrodynamic flow in rectangular ducts.* J. Fluid Mech. 21(4), 577–590 (thin-wall limit).
