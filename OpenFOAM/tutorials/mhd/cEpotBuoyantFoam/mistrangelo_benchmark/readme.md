# cEpotBuoyantFoam: magneto-convection benchmark in a horizontal duct (Mistrangelo et al. 2025)

> ⚠️ **Testing phase: set up only, not tested.** This case has not been run beyond a 5-step smoke test on a very coarse mesh, serial and on 2 processors, which only checks that the setup runs end to end. No benchmark results have been produced yet. Several material properties are assumptions (see below). Treat the case as work in progress.

Setup of the code-to-code benchmark of Mistrangelo, Bühler, Smolentsev et al. (2025), *A code-to-code benchmark for magneto-convection in a horizontal duct*, Nucl. Fusion 65, 116006. PbLi flows through a horizontal square duct with electrically and thermally conducting steel walls, inside a fringing magnetic field. It is heated from below over the central part of the magnet. Five groups (FLUENT, CFX, COMSOL, and two foam-extend OpenFOAM solvers) compared:
- time-averaged velocity, pressure, temperature and electric potential
- the bulk temperature, the Nusselt number and the transverse kinetic energy

This is a full 3D case with inlet and outlet. A production run is an **HPC job**.

## Run

```sh
CASE=1 NPROC=32 ./Allrun      # case 1: Re = 3220, Gr = 1.42e7 (laminar)
CASE=2 NPROC=32 ./Allrun      # case 2: Re = 1073, Gr = 3.55e7 (time-dependent)
python3 setup.py --case 1 --cells 40 --wall-cells 10 --axial-cells 600   # finer mesh
./Allclean
```

`Allrun` does the following:
1. `setup.py` writes the meshes, parameters and decomposition.
2. `restore0Dir` copies `0.orig` to `0`.
3. `blockMesh` and `checkMesh` run for both regions.
4. `setB0.py` writes B_y(x).
5. The solver runs, serial or with `decomposePar -allRegions` and `runParallel`.
6. `post.py` writes `benchmark/*.csv` and `benchmark/benchmark.png`.

**Cost (default mesh).** Fluid 300 × 56 × 56 = 0.94 M cells and solid 0.61 M cells. Δt = 0.01 s, with 30 000 steps for case 1 (300 s) and 60 000 for case 2 (600 s). For comparison, the participants used 0.3–2.9 M cells on 16–512 CPUs (paper table 3).

`checkMesh` reports high-aspect-ratio cells, because the Hartmann-layer cells are thin and the axial cells long. The mesh is orthogonal, and this is normal for MHD duct meshes.

## Problem definition

| Item | Value | Source |
|---|---|---|
| Duct | square, half-width a = 0.023 m (along B), steel walls t_w = 2 mm | paper §2 |
| Coordinates | x axial (magnet centre at x = 0), B along y, gravity along −z | paper fig. 1 |
| Domain | −1.2 m ≤ x ≤ 1.5 m (`--x-in`, `--x-out`); participants used 2.0–3.0 m | table 4 |
| Field | B_y(x) = ½B₀[1 − tanh((\|x\| − 0.7315)/0.1)], B₀ = 0.5 T; uniform for \|x\| < 0.4 m | eq. (1), fig. 2 |
| Heating | uniform flux q″ on the outer bottom surface (z = −a − t_w), −0.3 m ≤ x ≤ 0.3 m; q″ = 4 W/cm² (case 1), 10 W/cm² (case 2) | §2, fig. 2; case 2 from the Gr ratio |
| Other outer surfaces | electrically insulating, adiabatic | §3 |
| Inlet | uniform u₀ (0.03 / 0.01 m/s), T = 300 °C, ∂φ/∂x = 0 | table 4 (as CFX, OF_ORNL, COMSOL) |
| Outlet | p = 0, zero-gradient U and T (inletOutlet), ∂φ/∂x = 0 | table 4 |
| Groups | Ha = B₀a√(σ/ρν) = 235, Re = u₀a/ν, Gr = gβΔTa³/ν², ΔT = q″a/k | eq. (5), table 1 |
| Model | Boussinesq, inductionless, conjugate (φ, j_n, T and heat flux continuous at the wall) | §3 |

Cross-checks on the inputs:
- **Bulk temperature rise.** With the assumed ρc_p, the rise over the heated zone, q″(2a + 2t_w)(0.6 m)/(ρc_p u₀ 4a²), is 10.1 K (case 1) and 76 K (case 2). The paper's figs. 17 and 34 show 300 → 310 °C and 300 → 376 °C, which confirms q″ for case 2 and the heated width.
- **Time step.** The explicit Lorentz-force limit ρ/(σB₀²) = 0.045 s agrees with the Δt = 0.004–0.015 s the participants used.

## Material properties: what is exact and what is assumed
The paper takes PbLi properties from Martelli et al. (2019) [its ref. 10] and steel properties from the Stahl-Eisen-Werkstoffblätter [ref. 11], at 300 °C, but does not list the values. In `system/caseParameters`:

| quantity | value | status |
|---|---|---|
| ν (PbLi) | u₀a/Re = 2.143 × 10⁻⁷ m²/s | **exact** (from Re) |
| σ/ρ (PbLi) | (Ha/(B₀a))² ν | **exact** (from Ha) |
| ρ, c_p, k (PbLi) | 9838 kg/m³, 189.8 J/kg K, 13.1 W/m K | **assumed** (older correlations; ρc_p is consistent with the bulk-temperature rise above) |
| σ (PbLi) | 8.80 × 10⁵ S/m | follows from σ/ρ and the assumed ρ |
| β (PbLi) | 7.78 × 10⁻⁵ 1/K | chosen so that Gr is exact with the assumed k (only β/k is fixed by Gr) |
| σ_w, k_w, ρ_w, c_p,w (steel) | 1.2 × 10⁶ S/m, 17.5 W/m K, 7900 kg/m³, 530 J/kg K | **assumed**; c_w = σ_w t_w/(σa) = 0.119 |

**Before using results, replace the assumed values with those of refs. [10] and [11].** σ and β are then re-derived automatically. The wall conductance ratio and the Prandtl number (0.031 here) depend on these choices.

## Numerical setup
- `coupled true`: potential and temperature are solved as single matrices assembled across fluid and solid. `setup.py` builds matching interface faces, and decomposes both regions along x only (hierarchical, n = (N 1 1)) so that the interface pairs stay on the same processor. GAMG cannot be used for the assembled matrices, so PCG/PBiCGStab are used.
- Schemes: backward in time, linearUpwind for U, limitedLinear for T.
- Fixed Δt = 0.01 s. Time averaging starts at 150 s (`averageStart`); the paper's table 5 lists 50–450 s averaging windows.
- The inlet velocity is uniform. Several participants used a fully developed profile instead, and the paper notes that inlet and outlet choices affect the results near the ends.

## Outputs for comparison with the paper
- **Line samples (`postProcessing/benchmarkLines`):**
  - the axis
  - lines along both side walls and the jet lines at z = ±0.0215 m
  - z-profiles at y = 0 for x = −0.3, 0, 0.15, 0.3, 0.6 and 1 m (the stations of figs. 6–16, 23–33)
- **`post.py`:**
  - the axial pressure gradient in the uniform-field region, scaled by σu₀B₀²; the paper gives −0.0855 for case 1, against −0.0863 for fully developed flow
  - the side-wall potential, referenced to the inlet centre
  - the bulk temperature T_b(x) (eq. 6), the heated-wall temperature, and Nu(x) (eq. 7)

  It prints these next to approximate values read off the paper's figures. The paper provides no data tables. The length L in eq. (7) is not defined in the text, so `--nu-length a|2a` selects it.

## References
- C. Mistrangelo, L. Bühler, S. Smolentsev, A. Tassone, G.G. Fonfría, Y. Jiang, L. Melchiorri, G. Politis, D. Suarez, F.R. Urgorri (2025). *A code-to-code benchmark for magneto-convection in a horizontal duct.* Nucl. Fusion 65, 116006. doi:10.1088/1741-4326/ae0800 (open access, CC BY 4.0).
- D. Martelli, A. Venturini, M. Utili (2019). *Literature review of lead-lithium thermophysical properties.* Fusion Eng. Des. 138, 183–195.
- S. Smolentsev et al. (2020). *Code-to-code comparison for a PbLi mixed-convection MHD flow.* Fusion Sci. Technol. 76, 653–669.
