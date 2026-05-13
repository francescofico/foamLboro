if (finalIter)
{
    mesh.data::add("finalIteration", true);
}

if (frozenFlow) //usually false
{
    #include "PotEEqn.H"
    #include "TEqn.H"
}
else
{


//Coefficients matrix for velocity
// Solve the Momentum equation
#include "UEqn.H"
// Solve the electric potential
#include "PotEEqn.H"
// Solve the energy equation
#include "TEqn.H"
//PISO Loop
for (int corr=0; corr<nCorr; corr++)
{
	#include "pEqn.H"

}
    turbulence.correct();

} //end of else


if (finalIter)
{
    mesh.data::remove("finalIteration");
}

