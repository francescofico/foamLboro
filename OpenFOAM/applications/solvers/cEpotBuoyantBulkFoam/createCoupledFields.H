fvSolution couplingControls(runTime);
UPtrList<volScalarField> potentialFields(fluidRegions.size() + solidRegions.size());
forAll(fluidRegions, i) potentialFields.set(i, &PotEFluid[i]);
forAll(solidRegions, i) potentialFields.set(fluidRegions.size() + i, &PotEs[i]);
regionCoupledScalarSystem coupledPotential
(
    potentialFields, couplingControls, "coupledPotentialField",
    "solidWallMixedPotentialCoupled", true
);

UPtrList<volScalarField> temperatureFields(fluidRegions.size() + solidRegions.size());
forAll(fluidRegions, i) temperatureFields.set(i, &TFluid[i]);
forAll(solidRegions, i) temperatureFields.set(fluidRegions.size() + i, &Ts[i]);
regionCoupledScalarSystem coupledEnergy
(
    temperatureFields, couplingControls, "coupledEnergyField",
    "solidWallMixedTemperatureCoupled", false
);
coupledPotential.shareAssemblyWith(coupledEnergy);
const bool coupledScalars = coupledPotential.enabled() || coupledEnergy.enabled();

List<tmp<fvVectorMatrix>> fluidMomentumMatrices(fluidRegions.size());
