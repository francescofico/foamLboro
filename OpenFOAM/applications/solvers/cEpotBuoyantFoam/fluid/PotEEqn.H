    //Interpolating cross product (u x B) over mesh faces
    surfaceScalarField psiub = fvc::interpolate(sigma) * ((fvc::interpolate(U) ^ fvc::interpolate(B0)) & mesh.Sf());

    //Poisson equation for electric potential
    fvScalarMatrix PotEEqn
            (
                fvm::laplacian(sigma,PotE) 
                == 
                fvc::div(psiub)
            );

    if (coupledPotential.enabled())
    {
        coupledPotential.add(PotEEqn);
    }
    else
    {
        PotEEqn.setReference(PotERefCell, PotERefValue);
        PotEEqn.solve();
        PotE.correctBoundaryConditions();
        #include "correctElectricFields.H"
    }
