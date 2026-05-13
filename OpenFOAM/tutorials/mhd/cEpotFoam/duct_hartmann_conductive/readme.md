# cEpotFoam

## Description
This tutorial is a simple test case for the MHD equations with conductive walls.
The solid/fluid interface is dealt with an ad-hoc boundary condition based on mapped patch.
In this case only two wall are conductive - see mappedPatches in boundary file.

Solve the inductionless formulation for the MHD equations with arbitrary conductivity.
Temperature equation is not solved.

Only incompressible Navier-Stokes equations are solved.
Note that the transport quantities are random and the results might not physical.

## Steps
After compiling the solver, run the following commands:

1. blockMesh
2. cEpotFoam

or 

- run

## info
Fluid and solid are separated into different directories and the structure follows the original openfoam conjugate tutorial.

- B0 is the magnetic field and it is not a variable - required
- sigma is the electrical conductivity -required
- p is the pressure -required
- U is the velocity - required
- PotE is the electric potential - required

- J is the current density (not a variable) - not required
- lorentz is the Lorentz force - not required

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

## extra info 
python files use classy blocks library to generate the mesh, very useful for generating complex geometries
(https://github.com/damogranlabs/classy_blocks)