if (coupledScalars)
{
    const label nCoupledNonOrthCorr =
        pimple.lookupOrDefault<label>("nNonOrthogonalCorrectors", 0);
    for (label coupledNonOrth=0; coupledNonOrth<=nCoupledNonOrthCorr; ++coupledNonOrth)
    {
        if (coupledNonOrth)
        {
            forAll(fluidRegions, i)
            {
                #include "setRegionFluidFields.H"
                #include "setFinalRegionIteration.H"
                if (coupledPotential.enabled())
                {
                    #include "fluid/PotEEqn.H"
                }
                #include "clearFinalRegionIteration.H"
            }
            forAll(solidRegions, i)
            {
                #include "setRegionSolidFields.H"
                #include "setFinalRegionIteration.H"
                if (coupledPotential.enabled())
                {
                    #include "solid/PotEEqn.H"
                }
                #include "clearFinalRegionIteration.H"
            }
        }
        coupledPotential.solve(finalIter);
    }

    forAll(fluidRegions, i)
    {
        #include "setRegionFluidFields.H"
        if (coupledPotential.enabled())
        {
            surfaceScalarField psiub = fvc::interpolate(sigma)
                *((fvc::interpolate(U)^fvc::interpolate(B0)) & mesh.Sf());
            #include "fluid/correctElectricFields.H"
        }
    }
    if (coupledPotential.enabled())
    {
        forAll(solidRegions, i)
        {
            #include "setRegionSolidFields.H"
            #include "solid/correctElectricFields.H"
        }
    }
}
