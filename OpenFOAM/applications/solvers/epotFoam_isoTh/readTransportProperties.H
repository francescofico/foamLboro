singlePhaseTransportModel laminarTransport(U, phi);

Info<< "Reading transportProperties\n" << endl;

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);


//Kinematic viscosity
dimensionedScalar nu
(
    "nu",
    dimViscosity,
    transportProperties
);

//Density
dimensionedScalar rho
(
    "rho",
    dimDensity,
    transportProperties
);
