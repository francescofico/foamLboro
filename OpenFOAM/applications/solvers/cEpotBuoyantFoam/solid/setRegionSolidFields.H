    fvMesh& mesh = solidRegions[i];

    volScalarField& PotE = PotEs[i];
    volScalarField& T = Ts[i];
    volScalarField& sigma = sigmaSolid[i];
    volScalarField& K = KSolid[i];
    volVectorField& J = JSolid[i];
    volScalarField& alphat = alphatSolid[i];