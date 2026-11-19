
# OpenFOAM Code & Utilities - FF @lboro
![alt text](https://img.shields.io/badge/bash-gray)
![alt text](https://img.shields.io/badge/c%2B%2B-red)
![alt text](https://img.shields.io/badge/python-blueviolet)

This repository contains custom solvers, libraries, and post-processing tools developed for Magnetohydrodynamics (MHD) and Heat Transfer research at Loughborough University.

> **OpenFOAM version:** the code builds and runs with **OpenFOAM.com v2606**. The original v1912 implementation is kept on the `v1912-stable` branch.
>
> **Tutorials:** many tutorial cases are available, including validation cases checked against analytical and literature solutions. See [OpenFOAM/tutorials/mhd/readme.md](OpenFOAM/tutorials/mhd/readme.md) for the list and details.

## 🛠 OpenFOAM Code Development

### Solvers
*   **Inductionless MHD Solvers**: Custom solvers based on the low-magnetic Reynolds number approximation ($Rm \ll 1$). 
    *   Implementation of the electric potential ($\phi$) formulation.
    *   Coupling of Lorentz forces with momentum equations.
    *   Support for Buoyancy-driven flows (Magnetoconvection).
*   **Conjugate approach**: conjugate solver for arbiitrary electrical conductivity.


Check this dissertation for more info: [High-fidelity numerical investigation for liquid metal magnetoconvection](https://repository.lboro.ac.uk/articles/thesis/High-fidelity_numerical_investigation_of_liquid_metal_magnetoconvection/27837024?file=50598696).

### Libraries & Routines
*   **Boundary Conditions**: Custom `patch` implementations for electrically conductive walls.
*   **Statistics**: Classes for calculating statistics on the fly.

## 🐍 Post-Processing Tools

Located in the `/post` directory, these tools are for handling of large openfoam datasets.

*   **Field Manipulations**: Python scripts for calculating derived quantities
*   **mmap Field Reader**: A high-efficiency reader using memory mapping to handle OpenFOAM fields without loading entire datasets into RAM—ideal for large LES datasets.
*   **Automated Reconstruction**: Bash wrappers for `reconstructPar` and field mapping across different mesh topologies.

### Libraries & Routines
*   **Boundary Conditions**: Custom `patch` implementations for electrically conductive walls.
*   **Statistics**: Classes for calculating statistics on the fly.

## 📖 External Resources
*   **OpenFOAM.com**: [Source Code](https://dl.openfoam.com/source/) | [Installation](http://www.openfoam.com/download/)
*   **Documentation**: [User Guide](http://www.openfoam.com/documentation)
*   **Community**: [Code Wiki](https://develop.openfoam.com/Development/openfoam/-/wikis/) | [CFD Online](https://www.cfd-online.com/Forums/openfoam/)

## 📂 Repository Structure - [ in progress]
*   **`branch`**: Each branch corresponds to a specific OpenFOAM version (`main`: v2606, `v1912-stable`: v1912).
*   **`solvers/`**: Source code for the custom applications.
*   **`src/`**: Shared libraries and classes.
*   **`tutorials/`**: Example and validation cases (see the [tutorials readme](OpenFOAM/tutorials/mhd/readme.md))
## 📑 Citation
If you use this code for your research, please cite the following publications:

> Fico, F., Langella, I., & Xia, H. (2023). **Large-eddy simulation of magnetohydrodynamics and heat transfer in annular pipe liquid metal flow.** *Physics of Fluids*, 35(5), 055106. [DOI: 10.1063/5.0143687](https://doi.org/10.1063/5.0143687)

> Fico, F., Xia, H., Langella, I., & Politis, G. (2024). **Numerical study of liquid metal magnetoconvection and heat transfer in an electrically conductive square duct.** *International Journal of Heat and Mass Transfer*, 231, 125857. [DOI: 10.1016/j.ijheatmasstransfer.2024.125857](https://doi.org/10.1016/j.ijheatmasstransfer.2024.125857)

## 🚧 Status & Roadmap
- [ ] Complete tutorial documentation
- [ ] Add automated validation scripts
- [x] Update the code to latest OpenFOAM version (v2606)
- [ ] Extend mmap reader to support polyhedral meshes.