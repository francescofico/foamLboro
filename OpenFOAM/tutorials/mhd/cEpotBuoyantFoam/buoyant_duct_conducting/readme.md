# cEpotBuoyantFoam: buoyant vertical duct, conducting Hartmann walls, two heating orientations

Validation case. This is the buoyant duct of `epotFoamBuo/buoyant_duct_insulating` with **conducting Hartmann walls**, resolved as solid slabs a < |y| < a + t_w as in `cEpotFoam/hunt_duct`. Fluid and wall are coupled both electrically and thermally. It checks:

- the conjugate potential coupling (`solidWallMixedPotentialCoupled`) driven by buoyancy
- the conjugate temperature coupling (`solidWallMixedTemperatureCoupled`, conductivity ratio k_w/k)
- the potential and temperature inside the walls

| `--heating` | Heated faces | Temperature | What it adds |
|---|---|---|---|
| `normal` (default) | fluid side walls z = ∓a **and** the slab ends z = ∓a, hot/cold; outer slab faces adiabatic | T = −z/2 in fluid and walls; no heat crosses the interface | Zero interface flux: the coupled condition must not disturb T |
| `parallel` | outer slab faces y = ∓(a + t_w), hot/cold; slab ends and fluid side walls adiabatic | piecewise linear in y; the fluid gradient is set by the series resistance of walls and fluid | Heat crosses the interface; the fluid gradient depends on k_w/k |

## Run

```sh
./Allrun                      # heating normal to B, Ha = 100, about 2.5 min
HEATING=parallel ./Allrun     # grad T along B, heat flows through the walls
python3 setup.py --thickness 0.2 --sigma-ratio 0.5 --k-ratio 5 && ./Allrun
./Allclean
```

`validate.py` exits non-zero if an error is above 2 % (`--tol`) or the run is not steady.

## Problem definition

| Item | Value |
|---|---|
| Fluid | Square duct −a < y, z < a, x vertical, one periodic cell in x |
| Solid | Slabs a < \|y\| < a + t_w spanning −a < z < a, electrical conductivity σ_w, thermal conductivity k_w |
| Side walls (z = ±a) | No-slip, electrically insulating, no solid region |
| Electrical | Outer slab faces and slab ends insulating |
| Buoyancy | Coded body force Gr (T − T_ref) e_x (`constant/fluid/fvOptions`); see the insulating case for why |
| Units | Non-dimensional: a = ν = ρ = σ = β = ΔT = 1, B = Ha, fluid k = 1; κ = 1/Pr in fluid and walls |
| Default | Ha = 100, Gr = 2 × 10⁴, Pr = 0.01, t_w/a = 0.1, σ_w/σ = 1, k_w/k = 2 |
| Coupling | `coupled true`: assembled matrices for both potential and temperature |
| Time step | 0.25/Ha², endTime = 5/Ha |

For `parallel` heating the fluid temperature gradient is −1/(2 + 2t_w k/k_w). That is −0.476 with the defaults, so the body force is F = −9524 y.

## Analytical solution
`BuoyantDuctSolution(Ha, forcing, wall="thick", sigma_ratio, thickness)` in `../../common/mhdAnalytic.py`. It uses the exact thick-wall Robin condition of `cEpotFoam/thick_wall_duct`, which holds for both z-bases. The wall potential follows from the same modes. The wall temperature is linear (`parallel`: through the wall thickness; `normal`: −z/2).

## Results (v2606, default mesh)

| heating | u_max: OpenFOAM | u_max: exact | error in u_max | max error in u | max error in φ (fluid / wall) | max error in T (fluid / wall) | run time |
|---|---|---|---|---|---|---|---|
| normal to B | 29.514 | 29.384 | 0.44 % | 0.46 % | 0.27 % / 0.23 % | 7 × 10⁻¹⁰ / 7 × 10⁻¹⁰ | 2.3 min |
| parallel to B | 7.635 | 7.591 | 0.58 % | 1.33 % | 0.16 % / 0.07 % | 1 × 10⁻⁸ / 1 × 10⁻⁹ | 2.5 min |

Compared with the insulating duct, the conducting walls cut the peak velocity from 74 to 29.5 (normal) and from 11.6 to 7.6 (parallel). With ∇T ∥ B the flow sits in side-wall jets (see `validation.png`).

## Files

| File | Purpose |
|---|---|
| `setup.py` | Writes both `blockMeshDict`s, `system/caseParameters`, `0/fluid/T`, `0/solid/T` and `constant/fluid/fvOptions` |
| `validate.py` | Compares velocity, potential and temperature in the fluid and in the walls |

## References
- G. Authié, T. Tagawa, R. Moreau (2003). Eur. J. Mech. B/Fluids 22 (fully established buoyant flow in vertical enclosures with a horizontal field).
- J.C.R. Hunt (1965). J. Fluid Mech. 21(4), 577–590 (thin conducting Hartmann walls).
- Zheng et al. (2025). *Magnetohydrodynamic mixed convection in a vertical rectangular duct with all walls of arbitrary conductivity.* J. Fluid Mech. (the mixed-convection extension; candidate for Tier 3d).
