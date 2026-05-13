
    // Initialise fluid field pointer lists
    PtrList<singlePhaseTransportModel> thermoFluid(fluidRegions.size());
    PtrList<volVectorField> UFluid(fluidRegions.size());
    PtrList<volVectorField> JFluid(fluidRegions.size());
    PtrList<volVectorField> lorentzFluid(fluidRegions.size());
    PtrList<surfaceScalarField> phiFluid(fluidRegions.size());
	PtrList<volScalarField> pFluid(fluidRegions.size());
    PtrList<volScalarField> PotEFluid(fluidRegions.size());
    PtrList<volScalarField> TFluid(fluidRegions.size());
    PtrList<volScalarField> sigmaFluid(fluidRegions.size());
    PtrList<volVectorField> B0Fluid(fluidRegions.size());
    PtrList<volScalarField> divJFluid(fluidRegions.size());

    // fields for buoyant force
    PtrList<uniformDimensionedScalarField> hRefFluid(fluidRegions.size());
    PtrList<volScalarField> ghFluid(fluidRegions.size());
    PtrList<surfaceScalarField> ghfFluid(fluidRegions.size());
    PtrList<volScalarField> p_rghFluid(fluidRegions.size());

    PtrList<volScalarField> rhokFluid(fluidRegions.size());
    PtrList<volScalarField> alphatFluid(fluidRegions.size());

    PtrList<volScalarField> KFluid(fluidRegions.size());

   
	PtrList<dimensionedScalar> nuFluid(fluidRegions.size());
    PtrList<dimensionedScalar> rhoFluid(fluidRegions.size());
    PtrList<dimensionedScalar> betaFluid(fluidRegions.size());
    PtrList<dimensionedScalar> TRefFluid(fluidRegions.size());
    PtrList<dimensionedScalar> PrFluid(fluidRegions.size());
    PtrList<dimensionedScalar> PrtFluid(fluidRegions.size());

    // turbulence model initialization
    PtrList<incompressible::turbulenceModel> turbulenceFluid(fluidRegions.size());

    PtrList<fvVectorMatrix> UEqFluid(fluidRegions.size());
    PtrList<fv::options> fluidFvOptions(fluidRegions.size());

    List<bool> frozenFlowFluid(fluidRegions.size(), false);

    // reference cell for pressure
    List<label> pRefCellFluid(fluidRegions.size(),0);
    List<scalar> pRefValueFluid(fluidRegions.size(),0.0);

    List<label> PotERefCellFluid(fluidRegions.size(),0);
    List<scalar> PotERefValueFluid(fluidRegions.size(),0.0);

    // gravity
    const uniformDimensionedVectorField& g = meshObjects::gravity::New(runTime);

    // Populate fluid field pointer lists
    forAll(fluidRegions, i)
    {
        Info<< "*** Reading fluid mesh thermophysical properties for region "
            << fluidRegions[i].name() << nl << endl;



        Info<< "    Adding to alphaFluid\n" << endl;
            alphatFluid.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        "alphat",
                        runTime.timeName(),
                        fluidRegions[i],
                        IOobject::MUST_READ,
                        IOobject::AUTO_WRITE
                    ),
                    fluidRegions[i]
                )
            );

              Info<< "    Adding to alphaFluid\n" << endl;
            KFluid.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        "K",
                        runTime.timeName(),
                        fluidRegions[i],
                        IOobject::MUST_READ,
                        IOobject::AUTO_WRITE
                    ),
                    fluidRegions[i]
                )
            );


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
        Info<< "    Adding to TFluid\n" << endl;
        TFluid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "T",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                fluidRegions[i]
            )
        );

         // buoyancy related 
         //------------------------------
         Info<< "    Adding to hRefFluid\n" << endl;
         hRefFluid.set
         (
             i,
             new uniformDimensionedScalarField
             (
                 IOobject
                 (
                     "hRef",
                     runTime.constant(),
                     fluidRegions[i],
                     IOobject::READ_IF_PRESENT,
                     IOobject::NO_WRITE
                 ),
                 dimensionedScalar("hRef", dimLength, Zero) // uses name
             )
         );

         dimensionedScalar ghRef
         (
             mag(g.value()) > SMALL
           ? g & (cmptMag(g.value())/mag(g.value()))*hRefFluid[i]
           : dimensionedScalar("ghRef", g.dimensions()*dimLength, 0)
         );
     
         Info<< "    Adding to ghFluid\n" << endl;
         ghFluid.set
         (
             i,
             new volScalarField
             (
                 "gh",
                 (g & fluidRegions[i].C()) - ghRef
             )
         );
         //for pressure equation usage 
         Info<< "    Adding to ghfFluid\n" << endl;
         ghfFluid.set
         (
             i,
             new surfaceScalarField
             (
                 "ghf",
                 (g & fluidRegions[i].Cf()) - ghRef
             )
         );

                

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
        
        // turbulence validate
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
        
        // thermal expansion coefficient
        betaFluid.set
        (
        	i,
        	new dimensionedScalar(thermoFluid[i].lookup("beta"))
        );

        // thermal expansion coefficient
        TRefFluid.set
        (
        	i,
        	new dimensionedScalar(thermoFluid[i].lookup("TRef"))
        );


        Info<< "    Adding to rhokFluid\n" << endl;
         rhokFluid.set
         (
             i,
             new volScalarField
             (
                 IOobject
                 (
                     "rhok",
                     runTime.timeName(),
                     fluidRegions[i],
                     IOobject::NO_READ,
                     IOobject::AUTO_WRITE
                 ),
                 1.0 - betaFluid[i] * (TFluid[i]-TRefFluid[i])
             )
         );

         Info<< "    Adding to p_rghFluid\n" << endl;

         p_rghFluid.set
         (
             i,
             new volScalarField
             (
                 IOobject
                 (
                     "p_rgh",
                     runTime.timeName(),
                     fluidRegions[i],
                     IOobject::MUST_READ,
                     IOobject::AUTO_WRITE
                 ),
                 fluidRegions[i] 
             )
         );

         Info<< "    Adding to p_Fluid\n" << endl;
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
                     IOobject::NO_READ,     // mod based on chtMultiRegionFoam
                     IOobject::AUTO_WRITE
                 ),
                 p_rghFluid[i] + rhokFluid[i]*ghFluid[i] // is this right?
             )
         );
        
        Info<< " Pimple dict add entries\n" << endl;
        const dictionary& pimpleDict =
        fluidRegions[i].solutionDict().subDict("PIMPLE");
        pimpleDict.readIfPresent("frozenFlow", frozenFlowFluid[i]);

        Info<< "  set ref cells \n" << endl;
       
        
        fluidRegions[i].setFluxRequired(p_rghFluid[i].name());
        fluidRegions[i].setFluxRequired(PotEFluid[i].name());


        pRefCellFluid[i] = -1;
        pRefValueFluid[i] = 0.0;

       
        if (p_rghFluid[i].needReference())
        {
            setRefCell
            (
                pFluid[i],
                p_rghFluid[i],
                pimpleDict,
                pRefCellFluid[i],
                pRefValueFluid[i]
            );
        }
        // ref cell PotE
        setRefCell
        (
            PotEFluid[i],
            fluidRegions[i].solutionDict().subDict("PotE"),
            PotERefCellFluid[i],
            PotERefValueFluid[i]
        );


        // Laminar Prandtl number
        PrFluid.set
        (
        	i,
        	new dimensionedScalar(thermoFluid[i].lookup("Pr"))
        );
       
        // Turbulent Prandtl number
        PrtFluid.set
        (
        	i,
            new dimensionedScalar(thermoFluid[i].lookup("Prt"))
        );
        
    
        

    }
