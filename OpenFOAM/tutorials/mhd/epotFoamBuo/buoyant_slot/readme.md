# epotFoamBuo: buoyant vertical slot (buoyant Hartmann flow)

Validation case. A closed, tall vertical slot is heated on one wall and cooled on the other, with the magnetic field normal to the heated walls (∇T ∥ B). Away from the end walls the flow is fully developed: hot fluid rises, cold fluid sinks, and the Lorentz force flattens the profile into Hartmann layers. The case uses the solver's own Boussinesq buoyancy (`p_rgh` and `constant/g`), with no extra source terms, so it checks:

- the buoyancy term (ρk = 1 − β(T − T_ref) with `p_rgh`) coupled to the Lorentz force
- the temperature equation (pure conduction in the fully developed region)
- that the end walls only disturb the flow locally

## Run

```sh
./Allrun                    # Ha = 20, Gr = 4000 (system/caseParameters), about 30 s
HA=40 GR=16000 ./Allrun
./Allclean
```

`validate.py` exits non-zero if an error is above 1 % (`--tol`) or the last write interval changed the velocity by more than 0.1 % (not steady).

## Problem definition

| Item | Value |
|---|---|
| Geometry | −a < y < a, 0 < x < H with H = 8a (x vertical). One periodic cell in z: an infinitely wide slot |
| Gravity | g = (−Gr, 0, 0) in `constant/g` |
| Magnetic field | B0 = (0, B, 0), normal to the heated walls |
| Hot wall (y = −a) | T = +0.5, no-slip, insulating |
| Cold wall (y = +a) | T = −0.5, no-slip, insulating |
| End walls (x = 0, H) | Adiabatic, no-slip, insulating |
| Units | Non-dimensional: a = ν = ρ = σ = β = ΔT = 1, so B = Ha and \|g\| = Gr = gβΔT a³/ν² |
| Default | Ha = 20, Gr = 4000, Pr = 0.01; mesh 32 × 60 cells |
| Time step | 0.25/Ha² (explicit Lorentz force); endTime = 1 viscous time |

**Spanwise electrical condition.** The net flow across any horizontal section is zero, so the net spanwise current σB∫u dy is zero too. A plain cyclic for the potential in z is therefore already the insulating condition, and no potential jump is needed (unlike `epotFoam_isoTh/hartmann_channel`).

## Analytical solution
In the fully developed region the temperature is pure conduction, T = −y/2, and the momentum balance reduces to

```
u'' − Ha² u + F y = 0,   u(±1) = 0,   F = −Gr/2
u = F (y − sinh(Ha y)/sinh Ha) / Ha²
```

(`../../common/mhdAnalytic.py`, `buoyant_channel_profile`). The peak velocity is about F/Ha², here 4.0, and it sits just outside the Hartmann layers.

## Results (v2606, default mesh)

There is no net flow, so there is no pressure gradient to compare. The headline numbers are the peak upward and downward velocities at mid-height.

| quantity | OpenFOAM | exact | error |
|---|---|---|---|
| u_max (hot side, upward) | 4.0009 | 4.0010 | 0.002 % |
| u_min (cold side, downward) | −4.0011 | −4.0010 | 0.002 % |
| max velocity error at mid-height (rel. to u_max) | | | 0.37 % |
| max temperature error at mid-height | | | 3 × 10⁻⁵ |
| velocity error at x = H/4 | | | 0.37 % |

The error is the same at x = H/4 as at mid-height. The end walls only affect about one half-width next to them.

## Files

| File | Purpose |
|---|---|
| `setup.py` | Writes `system/blockMeshDict` and `system/caseParameters` (`--ha`, `--gr`, `--pr`, `--height`) |
| `validate.py` | Compares the mid-height profile with the exact solution |

## References
- J. Hartmann (1937). Kgl. Danske Videnskab. Selskab, Mat.-Fys. Medd. 15(6) (the Hartmann-layer structure).
- G. Authié, T. Tagawa, R. Moreau (2003). *Buoyant flow in long vertical enclosures in the presence of a strong horizontal magnetic field. Part 1. Fully-established flow.* Eur. J. Mech. B/Fluids 22. The fully established regime of the duct version of this problem.
