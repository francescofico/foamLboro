#include "fvCFD.H"

// Check the analytic interface value and equal/opposite integrated fluxes.
// Both regions must be registered before evaluating either boundary condition.
int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"

    IOdictionary settings
    (
        IOobject
        (
            "couplingTestDict", runTime.system(), runTime,
            IOobject::MUST_READ, IOobject::NO_WRITE
        )
    );
    wordList regionNames(2);
    regionNames[0] = "fluid";
    regionNames[1] = "solid";
    PtrList<fvMesh> meshes(2);
    PtrList<volScalarField> temperatures(2), potentials(2);
    PtrList<volScalarField> heatFluxTemperatures(2);
    PtrList<volScalarField> conductivities(2), sigmas(2);

    forAll(regionNames, regioni)
    {
        meshes.set
        (
            regioni,
            new fvMesh
            (
                IOobject
                (
                    regionNames[regioni], runTime.timeName(), runTime,
                    IOobject::MUST_READ
                )
            )
        );
    }

    forAll(regionNames, regioni)
    {
        fvMesh& mesh = meshes[regioni];
        const dictionary& dict = settings.subDict(regionNames[regioni]);
        auto readField = [&](const word& key) -> volScalarField*
        {
            return new volScalarField
            (
                IOobject
                (
                    dict.get<word>(key), runTime.timeName(), mesh,
                    IOobject::MUST_READ, IOobject::AUTO_WRITE
                ),
                mesh
            );
        };
        conductivities.set(regioni, readField("conductivity"));
        sigmas.set(regioni, readField("sigma"));
        temperatures.set(regioni, readField("temperature"));
        potentials.set(regioni, readField("potential"));
        heatFluxTemperatures.set(regioni, readField("heatFluxTemperature"));
    }

    const int oldTag = UPstream::msgType();
    scalar heatFlux = 0, currentFlux = 0;
    forAll(regionNames, regioni)
    {
        temperatures[regioni].correctBoundaryConditions();
        potentials[regioni].correctBoundaryConditions();
        const fvMesh& mesh = meshes[regioni];
        const label patchi = mesh.boundaryMesh().findPatchID("interface");
        if (patchi < 0)
        {
            FatalErrorInFunction << "Missing interface patch" << exit(FatalError);
        }
        const fvPatch& patch = mesh.boundary()[patchi];

        auto checkField =
        [&](const volScalarField& field, const word& expectedKey)
        {
            const scalar expected = settings.get<scalar>(expectedKey);
            const scalar slope = settings.get<scalar>(expectedKey + "Slope");
            const fvPatchScalarField& pf = field.boundaryField()[patchi];
            scalar error = 0;
            forAll(pf, facei)
            {
                error = max
                (
                    error,
                    mag(pf[facei] - expected - slope*patch.Cf()[facei].y())
                );
            }
            reduce(error, maxOp<scalar>());
            if (error > 1e-8)
            {
                FatalErrorInFunction
                    << mesh.name() << ':' << field.name()
                    << " interface error " << error << exit(FatalError);
            }
        };
        checkField(temperatures[regioni], "expectedTemperature");
        checkField(potentials[regioni], "expectedPotential");

        heatFlux += gSum
        (
            conductivities[regioni].boundaryField()[patchi]*patch.magSf()
           *temperatures[regioni].boundaryField()[patchi].snGrad()
        );
        currentFlux += gSum
        (
            sigmas[regioni].boundaryField()[patchi]*patch.magSf()
           *potentials[regioni].boundaryField()[patchi].snGrad()
        );

        heatFluxTemperatures[regioni].correctBoundaryConditions();
        const label walli = mesh.boundaryMesh().findPatchID("walls");
        const scalar imposedHeatFlux = settings.get<scalar>("prescribedHeatFlux");
        const scalar fluxError = gMax
        (
            mag
            (
                conductivities[regioni].boundaryField()[walli]
               *heatFluxTemperatures[regioni].boundaryField()[walli].snGrad()
               -imposedHeatFlux
            )
        );
        if (fluxError > 1e-8)
        {
            FatalErrorInFunction << "Prescribed wall heat-flux error "
                << fluxError << exit(FatalError);
        }

        // The runner reads these files again to exercise mixed-BC restarts.
        temperatures[regioni].write();
        potentials[regioni].write();
        heatFluxTemperatures[regioni].write();
    }

    if
    (
        mag(heatFlux) > 1e-8 || mag(currentFlux) > 1e-8
     || UPstream::msgType() != oldTag
    )
    {
        FatalErrorInFunction
            << "Unbalanced heat/current flux or changed communication tag: "
            << heatFlux << ' ' << currentFlux << exit(FatalError);
    }
    Info<< "PASS: interface values, heat/current conservation, prescribed heat flux and message tag"
        << endl;
    return 0;
}
