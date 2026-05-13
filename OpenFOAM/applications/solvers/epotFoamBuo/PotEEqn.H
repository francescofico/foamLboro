  // Description: Solving Poisson equation for electric potential
  // sigma must be interpolated here, cause it is a vol scalar field
  
  //Interpolating cross product u x B over mesh faces
  surfaceScalarField psiub = fvc::interpolate(sigma) * ((fvc::interpolate(U)^fvc::interpolate(B0)) & mesh.Sf());

  //Poisson equation for electric potential
  fvScalarMatrix PotEEqn
          (
              fvm::laplacian(sigma, PotE) 
              == 
              fvc::div(psiub)
          );

  //Reference potential
  //PotEEqn.setReference(PotERefCell, getRefCellValue(PotE, PotERefCell));
  PotEEqn.setReference(PotERefCell, PotERefValue);
  //Relaxation of Poisson equation
  //PotEEqn.relax();
  //Solving Poisson equation
  PotEEqn.solve();
  // Update potential distribution and boundary condition
  PotE.correctBoundaryConditions();

  //Computation of current density at cell faces
  surfaceScalarField jn = - fvc::interpolate(sigma)  *(fvc::snGrad(PotE) * mesh.magSf()) + psiub;

  //Current density at face center  (Ni2007:doi:10.1016/j.jcp.2007.07.023)
  surfaceVectorField jnv = jn * mesh.Cf();

  //Interpolation of current density at cell center
  J = fvc::surfaceIntegrate(jnv) - (fvc::surfaceIntegrate(jn) * mesh.C());

  //Update current density distribution and boundary condition
  J.correctBoundaryConditions();
  
  //Current density at face center (Ni2007:doi:10.1016/j.jcp.2007.07.023)
  surfaceVectorField jnvb = jn * (fvc::interpolate(B0)^mesh.Cf());

  //Conservative expression of Lorentz force (Ni2007:doi:10.1016/j.jcp.2007.07.023)
  lorentz = - fvc::surfaceIntegrate(jnvb) - (mesh.C()^(fvc::surfaceIntegrate(jn*(fvc::interpolate(B0)))));

  //Non conservative expression of Lorentz force
  //lorentz = sigma* (J ^ B0);

  //Update Lorentz force distribution and boundary conditions
  lorentz.correctBoundaryConditions();
  
  // printing the error of div J
  #include "electricCurrentErr.H"
  //Calculating and storing the error of div J
  divJ = fvc::div(jn);
  