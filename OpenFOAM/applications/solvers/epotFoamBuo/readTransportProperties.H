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

// Thermal expansion coefficient [1/K]
dimensionedScalar beta
(
	"beta",
	dimless / dimTemperature,
	laminarTransport
);

// Reference temperature [K]
dimensionedScalar TRef("TRef", dimTemperature, laminarTransport);


// Laminar Prandtl number
dimensionedScalar Pr("Pr", dimless, laminarTransport);

// Turbulent Prandtl number
dimensionedScalar Prt("Prt", dimless, laminarTransport);
