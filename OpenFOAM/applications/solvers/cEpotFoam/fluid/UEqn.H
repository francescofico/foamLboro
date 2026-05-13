// Solve the Momentum equation


tmp<fvVectorMatrix> tUEqn
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
    solve(UEqn == -fvc::grad(p));

    fvOptions.correct(U);
}


// tmp<fvVectorMatrix> tUEqn
// (
//     fvm::ddt(U)
//   + fvm::div(phi, U)
//   - fvm::laplacian(nu, U)
//   + turb.divDevReff(U)
// );

// fvVectorMatrix& UEqn = tUEqn.ref();
// UEqn.relax();



// volScalarField rUA = 1.0/UEqn.A();

// if (momentumPredictor)
// {
//     if (oCorr == nOuterCorr-1)
//     {
//         solve(UEqn == -fvc::grad(p), mesh.solver("UFinal"));
//     }
//     else
//     {
//         solve(UEqn == -fvc::grad(p));
//     }
// }
// else
// {
//     U = rUA*(UEqn.H() - fvc::grad(p));
//     U.correctBoundaryConditions();
// }
