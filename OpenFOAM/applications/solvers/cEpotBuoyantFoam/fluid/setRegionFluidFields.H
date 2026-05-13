    fvMesh& mesh = fluidRegions[i];
    //fvMesh& meshName = fluidRegions[i];


    //pimpleControl pimple(meshName);

    singlePhaseTransportModel& thermo = thermoFluid[i];
    volScalarField& PotE = PotEFluid[i];
    volScalarField& T = TFluid[i];
    volScalarField& sigma = sigmaFluid[i];
    volVectorField& U = UFluid[i];
    volVectorField& J = JFluid[i];
    volVectorField& B0 = B0Fluid[i];
    volVectorField& lorentz = lorentzFluid[i];
    volScalarField& divJ = divJFluid[i];
    surfaceScalarField& phi = phiFluid[i];
    incompressible::turbulenceModel& turbulence = turbulenceFluid[i];
    volScalarField& p = pFluid[i];
    volScalarField& p_rgh = p_rghFluid[i];
    volScalarField& alphat = alphatFluid[i];
    volScalarField& K = KFluid[i];
    volScalarField& rhok = rhokFluid[i];

    const volScalarField& gh = ghFluid[i];
    const surfaceScalarField& ghf = ghfFluid[i];
    
    
    dimensionedScalar& nu = nuFluid[i];
    dimensionedScalar& beta = betaFluid[i];
    dimensionedScalar& rho = rhoFluid[i];
    dimensionedScalar& Pr = PrFluid[i];
    dimensionedScalar& Prt = PrtFluid[i];
    dimensionedScalar& TRef = TRefFluid[i];
    fv::options& fvOptions = fluidFvOptions[i];

    const label pRefCell = pRefCellFluid[i];
    const scalar pRefValue = pRefValueFluid[i];

    const label PotERefCell = PotERefCellFluid[i];
    const scalar PotERefValue = PotERefValueFluid[i];;

    bool frozenFlow = frozenFlowFluid[i];
    