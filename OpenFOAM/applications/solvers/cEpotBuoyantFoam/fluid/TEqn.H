{

    // note effective thermal diffusivity
    // laminar + turbulent
    // laminar is given as a field
    // thermal diffusivity
    //alphat = turbulence.nut()/Prt;
    //alphat.correctBoundaryConditions();

    // effective thermal diffusivity
    //volScalarField alphaEff("alphaEff", nu/Pr + alphat);
    volScalarField alphaEff("alphaEff", turbulence.nut()/Prt + alphat);

    // temperature matrix

    fvScalarMatrix TEqn
    (
        fvm::ddt(T)
      + fvm::div(phi, T)
      - fvm::laplacian(alphaEff,T)
      ==
      fvOptions(T)
    );

    fvOptions.constrain(TEqn);
    TEqn.relax();


    if (oCorr == nOuterCorr-1)
    {
        TEqn.solve(mesh.solver("TFinal"));
    }
    else
    {
        TEqn.solve();
    }


    fvOptions.correct(T);

    // update of variable density
    rhok = 1.0 - beta*(T - TRef);

    //Info<< "Min/max T:" << min(kappaEff) << ' ' << max(kappaEff) << endl;
    Info<< "Min/max T:" << min(T) << ' ' << max(T) << endl;

}
