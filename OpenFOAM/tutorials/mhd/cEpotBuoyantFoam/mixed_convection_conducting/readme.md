# cEpotBuoyantFoam: MHD mixed convection in a vertical duct, conducting Hartmann walls (Zheng et al. 2025, case II)

> ⚠️ **Testing phase.** This case is new and has not yet been tested beyond the single run reported below. Treat its setup and results as preliminary.

This is the mixed-convection duct of `epotFoamBuo/mixed_convection_insulating` with **thick conducting Hartmann walls** and insulating side walls: case II of Zheng, Chen & Ni (2025). The walls are resolved as solid slabs a < |y| < a + t_w, coupled to the fluid electrically (assembled potential solve). The defaults follow the paper's case II: σ_w = σ, and δ_hw = 0.05 l, which is t_w/a = 0.1, so c_hw = σ_w δ_hw/(σ l) = 0.05.

## Run

```sh
./Allrun                                        # Gr/Re = 5, Ha = 50 (Ha_paper = 100), about 1.5 min
GRRE=1 ./Allrun
python3 setup.py --gr-re 1 --thickness 0.2 && ./Allrun     # thicker walls (paper fig. 9 uses delta_hw up to 0.5)
./Allclean
```

See the insulating case for the mapping to the paper's variables and how the case is driven. The heated faces are the fluid side walls and the slab ends at z = ∓a, so T = −z/2 in the fluid and in the walls.

## Reference solution
The sum of the thick-wall pressure-driven duct solution (`cEpotFoam/thick_wall_duct`) and the thick-wall buoyant solution (`cEpotBuoyantFoam/buoyant_duct_conducting`). Both use the exact Robin condition from Laplace's equation in the wall. Checks against the paper:
- **Fig. 7 (Gr → 0, Ha_paper = 100, δ_hw = 0.05):** the jets peak at 1.74 and the core is 0.89, matching the figure (about 1.75 and 0.88; the paper compares it with Bluck & Wolfendale 2015).
- **Core gradient (eq. 5.26):** −Gr/(2 Re Ha (1 + c Ha)) is −4.09; the OpenFOAM run gives −4.02. Eq. 5.26 is a thin-wall, high-Ha estimate.

Case III of the paper (conducting side walls as well) uses an approximation for highly conducting Hartmann walls. The exact series used here does not cover it, so it is not included.

## Results (single run, v2606, default mesh 80 × 80, 10 cells across each wall)

| quantity | OpenFOAM | exact | error |
|---|---|---|---|
| G a²/(ρνū) | 245.16 | 245.66 | 0.20 % |
| u_max/ū (hot-side jet) | 5.716 | 5.702 | 0.25 % |
| u_min/ū (reverse flow, cold side) | −2.254 | −2.242 | 0.50 % |
| max velocity error | | | 0.36 % |
| max potential error (fluid / wall) | | | 0.14 % / 0.13 % |

Run time: 79 s serial (1000 steps).

## References
- X. Zheng, L. Chen, M.-J. Ni (2025). J. Fluid Mech. 1018, A48. doi:10.1017/jfm.2025.10515
- M.J. Bluck, M.J. Wolfendale (2015). *An analytical solution to electromagnetically coupled duct flow in MHD.* J. Fluid Mech. 771, 595–623.
