
Info << "Reading field T\n" << endl;

//Pressure field

Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


//Velocity field
Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"

//Electric potential field
Info<< "Reading field PotE\n" << endl;
volScalarField PotE
(
    IOobject
    (
        "PotE",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Including transport properties
#include "readTransportProperties.H"

//Current density field
Info << "Preparing J\n" << endl;
volVectorField J
(
	IOobject
	(
		"J",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	),
	mesh,
	dimensionedVector("J", dimensionSet(0,-2,0,0,0,1,0), vector::zero)
);

volScalarField divJ
    (
        IOobject
        (
            "divJ",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("divJ", dimensionSet(0,-3,0,0,0,1,0), scalar(0))
    );

//Lorentz force - just for saving purpose
Info << "Preparing lorentz\n" << endl;
volVectorField lorentz
(
	IOobject
	(
		"lorentz",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	),
	mesh,
	dimensionedVector("lorentz",  dimensionSet(1, -2, -2, 0, 0, 0, 0), vector::zero )
);

//B0

Info << "Create B0\n" << endl;
volVectorField B0
(
	IOobject
	(
		"B0",
		runTime.timeName(),
		mesh,
		IOobject::MUST_READ,
		IOobject::AUTO_WRITE
	),
	mesh
);

Info << "Create sigma\n" << endl;
volScalarField sigma
(
	IOobject
	(
		"sigma",
		runTime.timeName(),
		mesh,
		IOobject::MUST_READ,
		IOobject::AUTO_WRITE
	),
	mesh
);


Info << "Creating turbulence model\n" << endl;
autoPtr<incompressible::turbulenceModel> turbulence
(
	incompressible::turbulenceModel::New(U, phi, laminarTransport)
);

// Kinematic density for buoyancy force

//Pressure reference point
label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, pimple.dict(), pRefCell, pRefValue);

mesh.setFluxRequired(p.name());

//Electric potential reference point
label PotERefCell = 0;
scalar PotERefValue = 0.0;
setRefCell(PotE, mesh.solutionDict().subDict("PotE"), PotERefCell, PotERefValue);
mesh.setFluxRequired(PotE.name());

#include "createFvOptions.H"
