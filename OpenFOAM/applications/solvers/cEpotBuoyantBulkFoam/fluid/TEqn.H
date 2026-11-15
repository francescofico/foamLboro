{

    // note effective thermal diffusivity
    // laminar + turbulent
    // laminar is given as a field
    // thermal diffusivity
    alphat = turbulence.nut()/Prt;
    alphat.correctBoundaryConditions();

    // effective thermal diffusivity
    //volScalarField alphaEff("alphaEff", nu/Pr + alphat);
    volScalarField alphaEff("alphaEff", nu/Pr + alphat);

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


    if (coupledEnergy.enabled())
    {
        coupledEnergy.add(TEqn);
    }
    else
    {
        TEqn.solve(mesh.solver(T.select(finalIter)));
        #include "correctThermalFields.H"
    }
}
