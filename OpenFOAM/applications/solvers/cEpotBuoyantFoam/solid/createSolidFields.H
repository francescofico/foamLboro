    // Initialise solid field pointer lists
    PtrList<volScalarField> PotEs(solidRegions.size());
    PtrList<volScalarField> Ts(solidRegions.size());
    PtrList<volScalarField> sigmaSolid(solidRegions.size());
    PtrList<volVectorField> JSolid(solidRegions.size());
    PtrList<volScalarField> alphatSolid(solidRegions.size());
    PtrList<volScalarField> KSolid(solidRegions.size());
    // Populate solid field pointer lists
    forAll(solidRegions, i)
    {
        Info<< "*** Reading solid mesh thermophysical properties for region "
            << solidRegions[i].name() << nl << endl;


        Info<< "    Adding to sigmaSolid\n" << endl;
        sigmaSolid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "sigma",
                    runTime.timeName(),
                    solidRegions[i],
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                solidRegions[i]
            )
        );

         Info<< "    Adding to KSolid\n" << endl;
        KSolid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "K",
                    runTime.timeName(),
                    solidRegions[i],
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                solidRegions[i]
            )
        );
        

        Info<< " Adding to PotEs\n" << endl;
        PotEs.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "PotE",
                    runTime.timeName(),
                    solidRegions[i],
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                solidRegions[i]
            )
        );
        Info<< " Adding to T\n" << endl;
        Ts.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "T",
                    runTime.timeName(),
                    solidRegions[i],
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                solidRegions[i]
            )
        );


        JSolid.set
        (   
        i,
        new volVectorField
        (
        IOobject
        (
            "J",
            runTime.timeName(),
            solidRegions[i],
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        solidRegions[i],
        dimensionedVector("J", dimensionSet(0,-2,0,0,0,1,0), Foam::vector(0,0,0))
        )
        );

        Info << "    Adding to alphat \n" << endl;
        alphatSolid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "alphat",
                    runTime.timeName(),
                    solidRegions[i],
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                solidRegions[i]
            )
        );

        
    }
