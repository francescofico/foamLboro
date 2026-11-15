fvSolution couplingControls(runTime);
UPtrList<volScalarField> potentialFields(fluidRegions.size() + solidRegions.size());
forAll(fluidRegions, i) potentialFields.set(i, &PotEFluid[i]);
forAll(solidRegions, i) potentialFields.set(fluidRegions.size() + i, &PotEs[i]);
regionCoupledScalarSystem coupledPotential
(
    potentialFields, couplingControls, "coupledPotentialField",
    "solidWallMixedPotentialCoupled", true
);

if (couplingControls.lookupOrDefault<bool>("coupledEnergyField", false))
{
    FatalErrorInFunction << "cEpotFoam has no temperature equation. Use a buoyant solver."
        << exit(FatalError);
}
const bool coupledScalars = coupledPotential.enabled();

List<tmp<fvVectorMatrix>> fluidMomentumMatrices(fluidRegions.size());
