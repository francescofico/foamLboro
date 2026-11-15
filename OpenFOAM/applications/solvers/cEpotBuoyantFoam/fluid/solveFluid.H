#include "setFinalRegionIteration.H"

if (!frozenFlow)
{
    #include "UEqn.H"
}

{
    #include "PotEEqn.H"
}
#include "TEqn.H"

if (!frozenFlow && !coupledScalars)
{
    tmp<fvVectorMatrix>& tUEqn = fluidMomentumMatrices[i];
    fvVectorMatrix& UEqn = tUEqn.ref();
    #include "correctFluidPressure.H"
    tUEqn.clear();
}

#include "clearFinalRegionIteration.H"
