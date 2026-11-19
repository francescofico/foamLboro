# epotFoamBuo: Hadley cavity (horizontal temperature gradient, vertical field)

Validation case. A long horizontal 2D cavity is heated at one end wall and cooled at the other, with gravity and the magnetic field both vertical. The fluid rises at the hot end, flows along the top to the cold end, sinks and returns along the bottom. Away from the end walls this becomes a parallel "Hadley" flow, damped by Hartmann layers at the top and bottom walls. This is the configuration of Garandet, Alboussière & Moreau (1992).

Unlike the vertical-duct cases, **the flow carries heat here**. The core temperature is T = γx + θ(y), and the correction θ exists only because u moves heat along the cavity. The case checks:

- the solver's own Boussinesq buoyancy (`p_rgh`, real gravity, closed cavity)
- advection of T and its coupling to the flow (κ θ'' = γ u)
- the Lorentz force on a buoyancy-driven core flow with Hartmann layers

## Run

```sh
./Allrun                      # Ha = 20, Gr = 4e5, Pr = 0.05 (system/caseParameters), about 2 min
HA=30 GR=9e5 ./Allrun
./Allclean
```

`validate.py` exits non-zero if an error is above 2 % (`--tol`) or the last write interval changed u or T by more than 0.2 % (not steady).

## Problem definition

| Item | Value |
|---|---|
| Geometry | 0 < x < L (horizontal), −a < y < a (vertical), L = 12a (aspect ratio 6). One periodic cell in z |
| Gravity | g = (0, −Gr, 0) in `constant/g` |
| Magnetic field | B0 = (0, B, 0): vertical, normal to the long walls |
| End walls | x = 0: T = +0.5; x = L: T = −0.5; no-slip, insulating |
| Top and bottom (y = ±a) | Adiabatic, no-slip, insulating |
| Units | Non-dimensional: a = ν = ρ = σ = β = ΔT = 1, so B = Ha, \|g\| = Gr, κ = 1/Pr |
| Default | Ha = 20, Gr = 4 × 10⁵, Pr = 0.05; mesh 96 × 48, graded towards all walls |
| Time step | 0.25/Ha². endTime = 2 L²/(π² κ), two decay times of the slowest conduction mode along the cavity |

The net flow through every vertical section is zero, so a plain cyclic for the potential in z is already the insulating condition (as in `buoyant_slot`). After the start-up the maximum Courant number settles at about 1.6, near the end walls. It peaks at about 6 in the first few steps, because T starts at 0 against end walls at ±0.5. That transient does not affect the steady result.

## Analytical solution
Far from the end walls, a parallel flow u(y) with T = γx + θ(y) satisfies the full equations exactly:

```
u'' − Ha² u = Gr γ y,  u(±a) = 0      →  u = −Gr γ (y − sinh(Ha y)/sinh Ha) / Ha²
κ θ'' = γ u,           θ'(±a) = 0
```

The core gradient γ is set by the end regions and is not known in closed form. `validate.py` measures it from the solution (a linear fit of the section-averaged T over L/2 ± a), then compares u and θ for that γ. With the defaults the flow lowers γ to −0.0765, against −1/L = −0.0833 for pure conduction. The heat carried by the flow shows directly as θ, which spans 8 % of ΔT across the cavity.

## Results (v2606, default mesh)

| quantity | OpenFOAM | exact core | error |
|---|---|---|---|
| u_max at x = L/2 | 61.111 | 61.114 | 0.005 % |
| θ(+a) − θ(0) | 0.08352 | 0.08357 | 0.06 % |
| max velocity error at mid-length (rel. to u_max) | | | 0.36 % |
| max error in θ (rel. to its range) | | | 0.06 % |
| vertical velocity in the core (rel.) | | | 2 × 10⁻⁵ |
| change in u between x = L/2 − 2a and L/2 | | | 0.09 % |

Run time: about 2 min serial (2336 steps, 4608 cells).

## Files

| File | Purpose |
|---|---|
| `setup.py` | Writes `system/blockMeshDict` and `system/caseParameters` (`--ha`, `--gr`, `--pr`, `--length`) |
| `validate.py` | Measures γ and compares the core u(y) and θ(y) with the exact parallel-flow solution |

## References
- J.P. Garandet, T. Alboussière, R. Moreau (1992). *Buoyancy driven convection in a rectangular enclosure with a transverse magnetic field.* Int. J. Heat Mass Transfer 35(4). doi:10.1016/0017-9310(92)90242-K. The core solution here corresponds to their parallel-flow analysis. Check the field orientation and scalings against the paper before comparing numbers one to one.
- A. Alchaar, P. Vasseur, E. Bilgen (1995). Numerical study of the same problem, including where the parallel-flow approximation is valid.
