if (coupledScalars)
{
    forAll(fluidRegions, i)
    {
        #include "setRegionFluidFields.H"
        #include "readFluidMultiRegionPIMPLEControls.H"
        if (!frozenFlow)
        {
            #include "setFinalRegionIteration.H"
            tmp<fvVectorMatrix>& tUEqn = fluidMomentumMatrices[i];
            fvVectorMatrix& UEqn = tUEqn.ref();
            #include "correctFluidPressure.H"
            tUEqn.clear();
            #include "clearFinalRegionIteration.H"
        }
    }
}
