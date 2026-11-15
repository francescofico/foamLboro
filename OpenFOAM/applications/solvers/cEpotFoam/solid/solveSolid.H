#include "setFinalRegionIteration.H"

{
    for (int nonOrth=0; nonOrth<=nNonOrthCorr; ++nonOrth)
    {
        if (!coupledPotential.enabled() || nonOrth == 0)
        {
            #include "PotEEqn.H"
        }
    }
}

#include "clearFinalRegionIteration.H"
