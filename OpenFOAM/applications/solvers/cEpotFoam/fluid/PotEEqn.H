    //Interpolating cross product (u x B) over mesh faces
    surfaceScalarField psiub = fvc::interpolate(sigma) * ((fvc::interpolate(U) ^ fvc::interpolate(B0)) & mesh.Sf());

    //Poisson equation for electric potential
    fvScalarMatrix PotEEqn
            (
                fvm::laplacian(sigma,PotE) 
                == 
                fvc::div(psiub)
            );

    //Reference potential
    PotEEqn.setReference(PotERefCell, PotERefValue);

    //Solving Poisson equation
    PotEEqn.solve();
    // Update electric potential distribution and boundary conditions
    PotE.correctBoundaryConditions();

    //Computation of current density path at cell faces
    surfaceScalarField jn = -fvc::interpolate(sigma) * (fvc::snGrad(PotE) * mesh.magSf()) + psiub;

    //Current density at face center (Ni2007:doi:10.1016/j.jcp.2007.07.023)
    surfaceVectorField jnv = jn * mesh.Cf();

    //Interpolation of current density at cell center
    J = fvc::surfaceIntegrate(jnv) - (fvc::surfaceIntegrate(jn) * mesh.C());

    //Update current density distribution and boundary conditions
    J.correctBoundaryConditions();

    //Current density at face center
    surfaceVectorField jnvb = jn * (fvc::interpolate(B0)^mesh.Cf());
   
    //Conservative expression of Lorentz force (Ni2007:doi:10.1016/j.jcp.2007.07.023)
    lorentz = - fvc::surfaceIntegrate(jnvb) - (mesh.C()^(fvc::surfaceIntegrate(jn*(fvc::interpolate(B0)))));
    
    //Update Lorentz force distribution and boundary conditions
    lorentz.correctBoundaryConditions();

    //Calculate the error of div J
    divJ = fvc::div(jn);

