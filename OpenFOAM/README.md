# OpenFOAM.com solvers

The `of2606-port` branch targets OpenFOAM.com v2606. `v1912-stable`
preserves the original v1912 baseline. Use a separate checkout of that branch
when running the original version; uncommitted edits can follow a branch switch.

## Build

Source OpenFOAM.com v2606, then build the shared library and all five solvers:

```sh
./Allwmake
```

`Allwmake` first builds `src/mhdBoundaryConditions` with `wmake libso`, then
builds the solvers. The library is written to
`FOAM_USER_LIBBIN/libmhdBoundaryConditions.so`. This branch uses the v2606 and has no compatibility paths for older OpenFOAM releases. Use
`v1912-stable` for the original v1912 implementation.

## Shared boundary-condition library

`src/mhdBoundaryConditions` provides one implementation of each condition:

- `solidWallMixedPotentialCoupled`
- `solidWallMixedTemperatureCoupled`
- `solidWallHeatFluxTemperature`

All three `cEpot*` solvers link to this library. Their existing field type
names and dictionary entries work without case changes. The potential and
temperature conditions remain separate classes; their matrix interface helper
is part of the library. The region equation assembly remains in
`applications/solvers/include` and the solver equation files.

To build just the library, run this from the `OpenFOAM` directory:

```sh
wmake libso src/mhdBoundaryConditions
```

Other applications and utilities can load the conditions by adding this to
`system/controlDict`:

```foam
libs ("libmhdBoundaryConditions.so");
```

## Region interfaces

`solidWallMixedPotentialCoupled` and `solidWallMixedTemperatureCoupled` use
OpenFOAM.com's `mappedWall` for both matching and nonmatching
interfaces. Select the mapping in the mesh patch dictionary:

| Interface | `sampleMode` |
|---|---|
| Matching faces | `nearestPatchFace` |
| Nonmatching faces | `nearestPatchFaceAMI` |

For example, the fluid-side mesh patch configuration is:

```foam
fluid_to_solid
{
    type            mappedWall;
    sampleRegion    solid;
    samplePatch     solid_to_fluid;
    sampleMode      nearestPatchFace;
    offsetMode      uniform;
    offset          (0 0 0);
    // faces (...) in blockMeshDict, or nFaces/startFace in polyMesh/boundary
}
```

On the solid side, point `sampleRegion` back to `fluid` and `samplePatch` to
`fluid_to_solid`. For a nonmatching interface, select `nearestPatchFaceAMI` on
both sides. AMI requires overlapping interface surfaces. Matching geometry
does not imply matching face order or processor ownership; the mapper handles
both. Check
[v2606 mappedPatchBase API](https://api.openfoam.com/2606/classFoam_1_1mappedPatchBase.html).

Set the field conditions on **both** sides:

```foam
// In the electric potential field:
fluid_to_solid
{
    type            solidWallMixedPotentialCoupled;
    PotEnbr         PotE;       // neighbour region's potential field
    sigma           sigma;      // this region's conductivity field
    value           uniform 0;
}

// In the temperature field:
fluid_to_solid
{
    type            solidWallMixedTemperatureCoupled;
    Tnbr            T;          // neighbour region's temperature field
    K               K;          // this region's thermal conductivity field
    value           uniform 300;
}
```

The neighbour's condition supplies its own conductivity field name. The two
regions can therefore use different field names and material properties.
The interface value is weighted by each side's conductivity and inverse
cell-to-face distance. By default, the regions are solved separately.

## Assembled region solves

The three `cEpot*` solvers can assemble one potential matrix across all fluid
and solid regions. The two buoyant solvers can also assemble one temperature
matrix. Add these switches at the top level of `system/fvSolution`:

```foam
coupledPotentialField true;
coupledEnergyField    true; // cEpotBuoyantFoam / cEpotBuoyantBulkFoam
```

Alternatively, `coupled true;` enables both fields in the buoyant solvers and
potential in `cEpotFoam`. Each field-specific switch overrides `coupled`.
The switches are read at startup. When no switch is specified, a saved
`useImplicit true;` on a field's interface enables assembly for that field.
Explicit `false` switches override these restart flags.

The solver enables implicit coefficients on both sides of each selected
`solidWallMixed*Coupled` interface. Potential and temperature are separate
linear systems. They are solved after all regional equations have been
assembled, before fluid pressure corrections. Current density, Lorentz force
and buoyancy are then updated from the solved fields. Frozen-flow and
`energyCoupling` loops use the same assembly path.

Native v2606 assembly requires matching, coincident interface faces with
opposite normals, **the same local face order**, and corresponding faces on
**the same MPI processor**. Use reciprocal `mappedWall` patches with
`sampleMode nearestPatchFace`. Startup checks reject incompatible interfaces.
For nonmatching AMI, reordered faces, or independently decomposed interfaces,
use the separate-region mode (`coupled false;`), which retains mapped exchange.

Optional global solver dictionaries in `system/fvSolution` control the
assembled systems, for example:

```foam
solvers
{
    PotE
    {
        solver          PCG;
        preconditioner  DIC;
        tolerance       1e-9;
        relTol          0;
    }
    PotEFinal { $PotE; }
    T
    {
        solver          PBiCGStab;
        preconditioner  DILU;
        tolerance       1e-9;
        relTol          0;
    }
    TFinal { $T; }
}
```

Without global dictionaries, the first region's field solver settings are
used. Temperature assembly is asymmetric: the equations use thermal
diffusivity while interface continuity uses thermal conductivity. Use a
solver/preconditioner that supports an asymmetric matrix for `T`.

Set `PIMPLE/nNonOrthogonalCorrectors` in the **top-level** `fvSolution` for
additional assembled scalar corrections. Regional settings continue to
control pressure and separately solved solid equations. Each unanchored
connected group of electrical regions receives one potential reference,
using the first region's `PotE` reference settings (default: zero at its
first global cell). Regions are assumed to be internally connected.

## Outer-iteration residual control

The three `cEpot*` solvers support `PIMPLE/residualControl` in each region's
`system/<region>/fvSolution`. For example, a buoyant fluid region can use:

```foam
PIMPLE
{
    // Keep the region's other PIMPLE controls here.
    residualControl
    {
        U     { tolerance 1e-4; relTol 0; }
        p_rgh { tolerance 1e-4; relTol 0; }
        PotE  { tolerance 1e-6; relTol 0; }
        T     { tolerance 1e-6; relTol 0; }
    }
}
```


This residual control in the final-correction follows
[OpenFOAM-dev multi-region residual-control patch](https://github.com/OpenFOAM/OpenFOAM-dev/commit/a884a44e4729ae662fa1e4361ccf67e7c498ca04),
using v2606's stored solver performance and adding relative tolerances.

## Regression tests (Codex)

With the OpenFOAM environment sourced:

```sh
python3 tests/coupledBoundaryConditions/run.py
```

The runner builds the shared library in a temporary directory and loads it
through `controlDict` in a probe that does not link to it. It checks analytic
interface values, equal/opposite integrated heat and electric-current fluxes,
and the prescribed wall heat flux. Cases cover matching faces in different orders,
nonmatching AMI faces, unequal conductivities and cell distances, different
field names, serial and two-rank MPI execution, empty local interface patches,
and write/read restarts. Invalid mesh patches and sampling modes must fail.

Use `--serial` to skip MPI. `--work-dir /path/to/new-directory` keeps
cases, binaries and logs at a chosen location. No tutorial or installed solver
is overwritten by these tests.

After building the solvers, run the actual executable tests:

```sh
python3 tests/coupledSolvers/run.py
```

Use `--bin-dir /path/to/bin` for an isolated build, `--serial` to skip MPI,
and `--work-dir /path/to/new-directory` to retain the test cases and logs.
The tests check a single assembled solve against analytic electric and
transient thermal slab balances, including current conservation, unequal
material properties, independent field switches, floating potential,
different thermal/electrical connection graphs, restarts, correction loops,
and rejection of incompatible interfaces. They also check the separate-region
AMI solve.

For outer-iteration convergence, run:

```sh
python3 tests/coupledSolvers/residualControl.py
```

This uses the same binary/work-directory options, plus `--parallel-only` for
just the MPI checks. Two cells across each slab give a spatial potential
gradient and a nonvanishing residual normalization. The tests compare early
exit solutions with analytic electrical profiles and discrete transient heat
balances over two time steps. Distinct final linear solvers verify `Final`
selection. Other cases check all-region/field gating, missing or disabled
equations, vector components, absolute and relative tolerances, iteration
limits, and energy/nonorthogonal loops in serial and with two MPI ranks.

## Validation record (codex)

Validated in November 2026:

- The shared boundary library and all five solvers compiled with OpenFOAM.com
  v2606, API 2606, patch 0.
- 
- The complete boundary-condition suite passed on v2606 in serial and with
  two MPI ranks, including dynamic library loading, prescribed heat flux,
  restarts and invalid-configuration checks.
