# epotFoamBuo: MHD mixed convection in a vertical duct, insulating walls (Zheng et al. 2025, case I)

> ⚠️ **Testing phase.** This case is new and has not yet been tested beyond the single run reported below. Treat its setup and results as preliminary.

Fully developed upward flow in a long vertical square duct. A pressure gradient drives the flow (bulk velocity fixed), buoyancy acts on it, the side walls are hot and cold, and the field is horizontal along y (∇T ⊥ B). This is case I of Zheng, Chen & Ni (2025). The mixed-convection parameter Gr/Re sets the balance:
- Gr/Re ≪ 1: symmetric Shercliff-type flow
- Gr/Re ~ 1: asymmetric flow with a constant core velocity gradient
- Gr/Re > 2: reverse flow near the cold wall

## Run

```sh
./Allrun                       # Gr/Re = 5, Ha = 50 (Ha_paper = 100), about 2.5 min
GRRE=1 ./Allrun                # mixed regime
GRRE=0 ./Allrun                # pure forced flow (Shercliff)
./Allclean
```

## Mapping to the paper
The paper uses the full side length l = 2a and the unit square 0 ≤ x̄, ȳ ≤ 1. Here the duct is −a < y, z < a.

| paper | here |
|---|---|
| x̄ (across the heated walls, ⊥ B) | (z + 1)/2 |
| ȳ (along B) | (y + 1)/2 |
| Ha (based on l) | 2 Ha (Ha here is based on a) |
| Gr/Re = gβΔTρ/(∂p/∂z) | same (no length scale) |
| w̄ with flow rate Q = 1 | u/ū with ū = 1 |
| T̄ = ½ − x̄ | T = −z/2 |

## How the case is driven
- `meanVelocityForce` holds ū = 1 (`constant/fvOptions`).
- The buoyancy is a coded body force g(T − T_ref), with g = (Gr/Re)·G_ref, written by `setup.py`. As in `buoyant_duct_insulating`, a periodic `p_rgh` cannot carry streamwise buoyancy, so `constant/g` = 0.
- The heated walls give T = −z/2 exactly: the fully developed flow does not carry heat across the duct, which matches the paper's inertialess assumption.

## Reference solution
The equations are linear, so the exact solution is the sum of the pressure-driven duct flow (Shercliff) and the buoyancy-driven flow for forcing linear in z (`BuoyantDuctSolution`, already verified in Tier 3b). The buoyant part carries no net flow, so G equals the pure-forced value for every Gr/Re. Checks against the paper:
- **Figs. 3a and 5a (Ha_paper = 100):** the Gr/Re = 0 plateau is 1.16; at Gr/Re = 10 the peak is 4.91 at x̄ = 0.135. Both match the figures.
- **Core gradient (eq. 3.22):** ∂w̄/∂x̄ = −Gr/(2 Re Ha), in the paper's normalisation, is reproduced. The OpenFOAM run gives −5.835 against −5.806 at Gr/Re = 5.
- **Reverse-flow threshold:** the exact solution first shows reverse flow at Gr/Re = 2.29 (Ha_paper = 100) and 2.12 (Ha_paper = 500). This tends to the paper's inertialess core estimate (Gr/Re)_c = 2.

## Results (single run, v2606, default mesh 80 × 80)

| quantity | OpenFOAM | exact | error |
|---|---|---|---|
| G a²/(ρνū) | 57.752 | 58.057 | 0.52 % |
| u_max/ū (hot-wall jet) | 3.0019 | 2.9945 | 0.25 % |
| u_min/ū (reverse flow, cold side) | −0.8980 | −0.8894 | 0.97 % |
| max velocity error | | | 0.30 % |
| max potential error | | | 0.19 % |

Run time: 134 s serial (2000 steps).

## References
- X. Zheng, L. Chen, M.-J. Ni (2025). *Magnetohydrodynamic mixed convection in a vertical rectangular duct with all walls of arbitrary conductivity.* J. Fluid Mech. 1018, A48. doi:10.1017/jfm.2025.10515
- T. Tagawa, G. Authié, R. Moreau (2002). *Buoyant flow in long vertical enclosures in the presence of a strong horizontal magnetic field. Part 1. Fully-established flow.* Eur. J. Mech. B/Fluids 21(4), 383–398.
