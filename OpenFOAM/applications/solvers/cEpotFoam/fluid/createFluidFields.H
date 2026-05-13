
    // Initialise fluid field pointer lists
    PtrList<singlePhaseTransportModel> thermoFluid(fluidRegions.size());
    PtrList<volVectorField> UFluid(fluidRegions.size());
    PtrList<volVectorField> JFluid(fluidRegions.size());
    PtrList<volVectorField> lorentzFluid(fluidRegions.size());
    PtrList<surfaceScalarField> phiFluid(fluidRegions.size());
	PtrList<volScalarField> pFluid(fluidRegions.size());
    PtrList<volScalarField> PotEFluid(fluidRegions.size());
    PtrList<volScalarField> sigmaFluid(fluidRegions.size());
    PtrList<volVectorField> B0Fluid(fluidRegions.size());
    PtrList<volScalarField> divJFluid(fluidRegions.size());

    PtrList<volScalarField> rhokFluid(fluidRegions.size());
    //PtrList<volScalarField> alphatFluid(fluidRegions.size());

   
	PtrList<dimensionedScalar> nuFluid(fluidRegions.size());
    PtrList<dimensionedScalar> rhoFluid(fluidRegions.size());
    

    PtrList<incompressible::turbulenceModel> turbulenceFluid(fluidRegions.size());

    PtrList<fvVectorMatrix> UEqFluid(fluidRegions.size());
    PtrList<fv::options> fluidFvOptions(fluidRegions.size());

    List<bool> frozenFlowFluid(fluidRegions.size(), false);

    List<label> pRefCellFluid(fluidRegions.size(),0);
    List<scalar> pRefValueFluid(fluidRegions.size(),0.0);

    List<label> PotERefCellFluid(fluidRegions.size(),0);
    List<scalar> PotERefValueFluid(fluidRegions.size(),0.0);


    // Populate fluid field pointer lists
    forAll(fluidRegions, i)
    {
        Info<< "*** Reading fluid mesh thermophysical properties for region "
            << fluidRegions[i].name() << nl << endl;

        Info<< "    Adding to sigmaFluid\n" << endl;
        sigmaFluid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "sigma",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                fluidRegions[i]
            )
        );

        Info<< "    Adding to PotEFluid\n" << endl;
        PotEFluid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "PotE",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                fluidRegions[i]
            )
        );


        Info<< "    Adding to pFluid\n" << endl;
        pFluid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "p",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                fluidRegions[i]
            )
        );
         fluidRegions[i].setFluxRequired(pFluid[i].name());

        Info<< "    Adding to UFluid\n" << endl;
        UFluid.set
        (
            i,
            new volVectorField
            (
                IOobject
                (
                    "U",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                fluidRegions[i]
            )
        );

        Info<< "    Adding to JFluid\n" << endl;
        JFluid.set
        (
            i,
            new volVectorField
            (
                IOobject
                (
                    "J",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                fluidRegions[i],
                dimensionedVector("J", dimensionSet(0,-2,0,0,0,1,0), vector::zero)
            )
        );

        Info<< "    Adding to lorentzFluid\n" << endl;
        lorentzFluid.set
        (
            i,
            new volVectorField
            (
                IOobject
                (
                    "lorentz",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                fluidRegions[i],
                dimensionedVector("lorentz",  dimensionSet(1, -2, -2, 0, 0, 0, 0), vector::zero )
            )
        );

        Info<< "    Adding to J divergence\n" << endl;
        divJFluid.set
        (
        i,
        new volScalarField
        (
        IOobject
        (
            "divJ",
            runTime.timeName(),
            fluidRegions[i],
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fluidRegions[i],
        dimensionedScalar("divJ", dimensionSet(0,-3,0,0,0,1,0), Foam::scalar(0))
        )
    );       

        Info<< "    Adding to B0Fluid\n" << endl;
        B0Fluid.set
        (
            i,
            new volVectorField
            (
                IOobject
                (
                    "B0",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                fluidRegions[i]
            )
        );


        Info<< "    Adding to phiFluid\n" << endl;
        phiFluid.set
        (
            i,
            new surfaceScalarField
            (
                IOobject
                (
                    "phi",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                linearInterpolate(UFluid[i]) & fluidRegions[i].Sf()
            )
        );

        Info<< "    Adding to thermoFluid\n" << endl;
        thermoFluid.set
        (
            i,
            //new singlePhaseTransportModel(UFluid[i], phiFluid[i])
            new singlePhaseTransportModel (UFluid[i], phiFluid[i])

        );

        Info<< "    Adding to turbulence\n" << endl;
        turbulenceFluid.set
        (
            i,
            autoPtr<incompressible::turbulenceModel>
            (
            	incompressible::turbulenceModel::New
                (
                    UFluid[i],
                    phiFluid[i],
                    thermoFluid[i]
                )
            ).ptr()
        );

        const dictionary& pimpleDict =
        fluidRegions[i].solutionDict().subDict("PIMPLE");
        pimpleDict.readIfPresent("frozenFlow", frozenFlowFluid[i]);


        setRefCell
        (
            pFluid[i],
            fluidRegions[i].solutionDict().subDict("PIMPLE"),
            pRefCellFluid[i],
            pRefValueFluid[i]
        );

         setRefCell
        (
            PotEFluid[i],
            fluidRegions[i].solutionDict().subDict("PotE"),
            PotERefCellFluid[i],
            PotERefValueFluid[i]
        );
      
        Info<< "    Adding fvOptions\n" << endl;
        fluidFvOptions.set
        (
            i,
            new fv::options(fluidRegions[i])
        );

        UEqFluid.set
        (
        i,
        new fvVectorMatrix(UFluid[i], dimForce)
        );

        turbulenceFluid[i].validate();

        // Kynematic viscosity
        nuFluid.set
        (
        	i,
        	new dimensionedScalar(thermoFluid[i].lookup("nu"))
        );



         // density 
        rhoFluid.set
        (
        	i,
        	new dimensionedScalar(thermoFluid[i].lookup("rho"))
        );
         
       

        

    }
