# epotFoam_isoTh

## Description
This tutorial is a simple test case for the MHD equations.


Solve the inductionless formulation for the MHD equations in insulated case (perfectly conductive walls can also be used)
Temperature equation is not solved.

Only incompressible Navier-Stokes equations are solved.

Note that the transport quantities are random and the results might not physical.

## Steps

- 1. blockMesh
- 2.epotFoam_isoTh

-or run

# info
- B0 is the magnetic field and it is not a variable - required
- sigma is the electrical conductivity -required
- p is the pressure -required
- U is the velocity - required
- PotE is the electric potential - required
- J is the current density (not a variable) - not required
- lorentz is the Lorentz force - not required



## extra info 
duct.py file uses classy blocks library to generate the mesh, very useful for generating complex geometries
(https://github.com/damogranlabs/classy_blocks)