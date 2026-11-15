// Solve the Momentum equation


tmp<fvVectorMatrix>& tUEqn = fluidMomentumMatrices[i];
tUEqn =
(
fvm::ddt(U) + fvm::div(phi, U)
+ turbulence.divDevReff(U)
-(1.0/rho)*lorentz
==
fvOptions(U)
);

fvVectorMatrix& UEqn = tUEqn.ref();


UEqn.relax();

fvOptions.constrain(UEqn);

if (momentumPredictor)
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
