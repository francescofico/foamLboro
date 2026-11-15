//Computation of electric current density at cell faces
surfaceScalarField jnSolid = -fvc::interpolate(sigma) * fvc::snGrad(PotE)*mesh.magSf();	 

//Electric current density at face center (Ni2007:doi:10.1016/j.jcp.2007.07.023) 
surfaceVectorField jnvSolid = jnSolid * mesh.Cf();

//Interpolation of electric current density at cell center
J = fvc::surfaceIntegrate(jnvSolid) - (fvc::surfaceIntegrate(jnSolid) * mesh.C());	

//Update electric current density distribution and boundary condition
J.correctBoundaryConditions();
