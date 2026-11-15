/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2022 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
This file is result of research done at: Loughborough University
under the project:
"Development of cooling strategies  and advanced numerical approaches 
for heat transfer in nuclear fusion reactor components under extreme 
heat loads" (EP/R513088/1).
[Link](https://gtr.ukri.org/projects?ref=studentship-2498033#/tabOverview) 

Please aknowledge/cite Loughborugh University when using this code.
How to cite:
>   Francesco Fico, Ivan Langella, Hao Xia; Large-eddy simulation of
    magnetohydrodynamics and heat transfer in annular pipe liquid metal flow.
    Physics of Fluids 1 May 2023; 35 (5): 055106.
    https://doi.org/10.1063/5.0143687

>   Francesco Fico, Hao Xia, Ivan Langella, Gerasimos Politis; Numerical study of 
    liquid metal magnetoconvection and heat transfer in an electrically conductive square duct.
    International Journal of Heat and Mass Transfer, 2024.
    10.1016/j.ijheatmasstransfer.2024.125857

License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    cEpotFoam

Description
     A multiregion heat transfer solver. The fluid equation is of a turbulent, 
    electrically conductive, incompressible flow.
    The mhd equations are solved with the inductionless formulation.
    Turbulence model is valid for LES, DES and RANS. 
    This solver is based on an incompressible variant of chtMultiRegionFoam.

Contact
    Francesco Fico <f.fico@lboro.ac.uk> 01/03/2022 

\*---------------------------------------------------------------------------*/


#include "fvCFD.H"
//#include "turbulentFluidThermoModel.H" //thermofluid
#include "fixedGradientFvPatchFields.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "turbulenceModel.H"
#include "regionProperties.H"
#include "icoCourantNo.H"
//#include "solidRegionDiffNo.H"
#include "fvOptions.H"
//#include "pimpleControl.H" //pimple control re-created
#include "loopControl.H" // for last energy coupling only
#include "regionCoupledScalarSystem.H"
#include "regionResidualControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #define NO_CONTROL
    #define CREATE_MESH createMeshesPostProcess.H
    #include "postProcess.H"
    #include "setRootCaseLists.H"

    #include "createTime.H"    
    #include "createMeshes.H"
    #include "createFields.H"
    #include "createCoupledFields.H"
   
    #include "initContinuityErrs.H"  
    #include "createTimeControls.H" 
    #include "readTimeControls.H"
    #include "readSolidTimeControls.H"
    #include "icoMultiRegionCourantNo.H"
    #include "setInitialMultiRegionDeltaT.H"

   
    forAll(fluidRegions,i)
    { 
        Info << "Initial estimation of Lorentz force" << nl << endl;
        #include "setRegionFluidFields.H"
        #include "readFluidMultiRegionPIMPLEControls.H"
        lorentz = sigma * (-fvc::grad(PotE) ^ B0) + sigma * ((U ^ B0) ^ B0);

    }
    
    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "readSolidTimeControls.H"
        #include "readPIMPLEControls.H"
        #include "icoMultiRegionCourantNo.H"
        #include "setMultiRegionDeltaT.H"

        ++runTime;
        Info<< "Time = " << runTime.timeName() << nl << endl;

        if (nOuterCorr != 1)
        {
            forAll(fluidRegions, i)
            {
                #include "storeOldFluidFields.H"
            }
        }

        regionResidualControl residualControl(fluidRegions, solidRegions);

        // --- PIMPLE loop
        for (int oCorr=0; oCorr<nOuterCorr; oCorr++)
        {
            const bool finalIter =
                (oCorr == nOuterCorr-1 || residualControl.converged());
            Info<< nl << "Pimple Outer Correction = " << oCorr << endl;

            forAll(fluidRegions, i)
            {
                Info<< "\nSolving for fluid region "
                        << fluidRegions[i].name() << endl;
                #include "setRegionFluidFields.H"
                #include "readFluidMultiRegionPIMPLEControls.H"
                #include "solveFluid.H"
            }

            forAll(solidRegions, i)
            {
                Info<< "\nSolving for solid region "
                        << solidRegions[i].name() << endl;
                #include "setRegionSolidFields.H"
                #include "readSolidMultiRegionPIMPLEControls.H"
                #include "solveSolid.H"
            }

            #include "solveCoupledScalars.H"
            #include "correctCoupledFluidPressure.H"

            // Additional loops for energy solution only
            if (!oCorr && nOuterCorr > 1)
            {
                loopControl looping(runTime, pimple, "energyCoupling");

                while (looping.loop())
                {
                    Info<< nl << looping << nl;

                    forAll(fluidRegions, i)
                    {
                        Info<< "\nSolving for fluid region "
                            << fluidRegions[i].name() << endl;
                    #include "setRegionFluidFields.H"
                    #include "readFluidMultiRegionPIMPLEControls.H"
                    frozenFlow = true;
                    #include "solveFluid.H"
                    }

                    forAll(solidRegions, i)
                    {
                        Info<< "\nSolving for solid region "
                            << solidRegions[i].name() << endl;
                        #include "setRegionSolidFields.H"
                        #include "readSolidMultiRegionPIMPLEControls.H"
                        #include "solveSolid.H"
                    }

                    #include "solveCoupledScalars.H"
                }
            }

            if (residualControl.check(oCorr, nOuterCorr)) break;
        }

            runTime.write();
            runTime.printExecutionTime(Info);
    }
    Info << "End\n" << endl;

return 0;
}
// ************************************************************************* //
