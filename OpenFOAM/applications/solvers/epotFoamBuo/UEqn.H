tmp<fvVectorMatrix> tUEqn
(
        fvm::ddt(U)
        + fvm::div(phi, U)
        + turbulence->divDevReff(U)
        - (1.0 / rho) * lorentz //Lorentz term
        ==
        fvOptions(U)
);
 
fvVectorMatrix& UEqn = tUEqn.ref();
 
UEqn.relax();
fvOptions.constrain(UEqn);
 
if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                  - ghf*fvc::snGrad(rhok)
                  - fvc::snGrad(p_rgh)
                )*mesh.magSf()
            )
        );
 
        fvOptions.correct(U);
    }
