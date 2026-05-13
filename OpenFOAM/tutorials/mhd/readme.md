# Tutorials & Examples (Work in Progress)

This directory contains example cases to demonstrate the setup of the inductionless MHD solvers. 

## ⚠️ Disclaimer
These cases are currently **functional examples** intended to demonstrate dictionary setups and boundary conditions. 

### Standard Potential Solvers

*   **`epotFoam_isoTh`**: Isothermal inductionless MHD solver. Focuses on the coupling between the Navier-Stokes equations and the electric potential ($\phi$) formulation for Lorentz force calculation.
*   **`epotFoamBuo`**: Buoyant inductionless MHD solver. Incorporates Boussinesq or coupled density-driven flows, specifically designed for magnetoconvection studies.

#### Conjugate Solvers (Arbitrary Conductivity)
*   **`cEpotFoam`**: A conjugate solver designed to handle **arbitrary electrical conductivity** across the fluid-solid interface.
*   **`cEpotBouyantFoam`**: The buoyant counterpart to `cEpotFoam`, used for complex magnetoconvection scenarios in geometries with conductive walls.
*   **`cEpotBouyantBulkFoam`** `[Experimental]`: Extension of `cEpotBouyantFoam` featuring variable bulk temperature formulations.

## Mesh

Used https://github.com/damogranlabs/classy_blocks/tree/development library for generating the meshes.