fvScalarMatrix PotEEqn
    (
        fvm::laplacian(sigma,PotE) 
        
    );

    if (coupledPotential.enabled())
    {
        coupledPotential.add(PotEEqn);
    }
    else
    {
        PotEEqn.solve();
        PotE.correctBoundaryConditions();
        #include "correctElectricFields.H"
    }
