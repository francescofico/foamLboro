# cEpotBuoyantFoam

## Description
This tutorial is a simple test case for the MHD equations with conductive walls and buoyancy forces.
The solid/fluid interface is dealt with an ad-hoc boundary condition based on mapped patch.
In this case only two wall are conductive - see mappedPatches in boundary file.

Solve the inductionless formulation for the MHD equations with arbitrary conductivity.
Temperature equation is solved and buoyancy is considered through the Boussinesq approximation. 
Tref is the reference temperature for buoyancy computation.

Only incompressible Navier-Stokes equations are solved.

Note that the transport quantities are random and the results might not physical.
Mesh is very coarse and might not be able to converge.

## Steps
After compiling the solver, mesh and run the solver.

## info
Fluid and solid are separated into different directories and the structure follows the original openfoam conjugate tutorial.

- B0 is the magnetic field and it is not a variable - required
- sigma is the electrical conductivity -required
- p is the pressure -required
- U is the velocity - required
- PotE is the electric potential - required
- K is the thermal conductivity - required
- alphat is the thermal expansion coefficient - required

- J is the current density (not a variable) - not required
- lorentz is the Lorentz force - not required

- TRef is the reference temperature for buoyancy computation - required
- Pr is the Prandtl number - required
- Prt is the Prandtl number for the thermal transport - required
- beta is the volumetric thermal expansion coefficient - required
- cp is the specific heat - required

## Note
alphat is extrapolated from Pr and nu, so make sure to use the correct value.
K is mostly used for the coupled boundary, however make sure to use the correct value to avoid any problems


## details

The coupling boundary is defined as:

```
front
{
        type            solidWallMixedPotentialCoupled;
        value           uniform 0;
        PotEnbr         PotE;
        sigma           sigma;
}
```

in the field dictionary and as:

```
front
    {
        type            mappedWall;
        inGroups        2(wall mappedPatch);
        nFaces          2200;
        startFace       288464;
        sampleMode      nearestPatchFace;
        sampleRegion    solid;
        samplePatch     back;
    }

``` 
 
in the boundary file. - This need to be done for each coupling boundary and in both the solid and fluid directories.
Similary can be done for the temperature boundary.

## extra info 
python files use classy blocks library to generate the mesh, very useful for generating complex geometries
(https://github.com/damogranlabs/classy_blocks)