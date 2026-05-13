{
    volScalarField rAU("rAU", 1.0/UEqn.A());
    surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));
    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p_rgh));
 
    surfaceScalarField phig(-rAUf*ghf*fvc::snGrad(rhok)*mesh.magSf());
 
    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::flux(HbyA)
      + phig
    );
 
    //MRF.makeRelative(phiHbyA);
 
    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p_rgh, U, phiHbyA, rAUf);
 
    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix p_rghEqn
        (
            fvm::laplacian(rAUf, p_rgh) == fvc::div(phiHbyA)
        );
	
        p_rghEqn.setReference(pRefCell, getRefCellValue(p_rgh, pRefCell));

        p_rghEqn.solve
            (
                mesh.solver
                (
            
                    p_rgh.select
                    ( 
            
                        (
                                oCorr == nOuterCorr-1
                            && corr == nCorr-1
                            && nonOrth == nNonOrthCorr
                        )
                    )
                )
            );
        
        // old fashing way to do it
        // if (
        //     oCorr == nOuterCorr-1
        //     && corr == nCorr-1
        //     && nonOrth == nNonOrthCorr
        //    )
        //    {
        //        p_rghEqn.solve(mesh.solver("pFinal"));
        //    }
        //    else
        //    {
        //        p_rghEqn.solve();
        //    }

        //p_rghEqn.solve(mesh.solver(p_rgh.select(pimple.finalInnerIter()))); // pimple controls way to do it
 
        if (nonOrth == nNonOrthCorr)
        {
            // Calculate the conservative fluxes
            phi = phiHbyA - p_rghEqn.flux();

            U = HbyA + rAU*fvc::reconstruct((phig - p_rghEqn.flux())/rAUf);
            U.correctBoundaryConditions();
            fvOptions.correct(U);

        }
    
    }
    
   #include "newContinuityErrs.H" //modified continuityErrs to work with multiple fluids
    
    // Explicitly relax pressure for momentum corrector
    p_rgh.relax();

    // Correct the momentum source with the pressure gradient flux
    // calculated from the relaxed pressure
    // U = HbyA + rAU*fvc::reconstruct((phig - p_rghEqn.flux())/rAUf);
    // U.correctBoundaryConditions();
    // fvOptions.correct(U);

    p = p_rgh + rhok*gh;
 
    if (p_rgh.needReference())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
        p_rgh = p - rhok*gh;
    }
}