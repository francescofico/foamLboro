

volScalarField rAU(1.0/UEqn.A());
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
surfaceScalarField phiHbyA("phiHbyA", fvc::flux(HbyA));

  if (p.needReference())
    {
        adjustPhi(phiHbyA, U, p);   
    }
    tmp<volScalarField> rAtU(rAU);
      if (consistent)
     {
        rAtU = 1.0/max(1.0/rAU - UEqn.H1(), 0.1/rAU);
        phiHbyA +=
            fvc::interpolate(rAtU() - rAU)*fvc::snGrad(p)*mesh.magSf();
        HbyA -= (rAU - rAtU())*fvc::grad(p);
    }

    if (nCorr <= 1)
    {
        tUEqn.clear();
    }
    

    constrainPressure(p, U, phiHbyA, rAtU());

    // Non-orthogonal pressure corrector loop
    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        // Pressure corrector
        fvScalarMatrix pEqn
        (
            fvm::laplacian(rAtU(), p) == fvc::div(phiHbyA)
        );

        pEqn.setReference(pRefCell, pRefValue);
        
        if (
                finalIter
            && corr == nCorr-1
            && nonOrth == nNonOrthCorr
           )
            {
                pEqn.solve(mesh.solver("pFinal"));
            }
            else
            {
                pEqn.solve();
            }
           

       //pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));

        if (nonOrth == nNonOrthCorr)
        {
            phi = phiHbyA - pEqn.flux();
        }
    
    }

    #include "newContinuityErrs.H" //modified continuityErrs to work with multiple fluids

    // Explicitly relax pressure for momentum corrector
    p.relax();
   

    U = HbyA - rAtU*fvc::grad(p);
    U.correctBoundaryConditions();
    fvOptions.correct(U); // to be included

  

