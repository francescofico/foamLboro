/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 Original author
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
    epotFoamBuoyant

Description
    Transient solver for incompressible, turbulent MHD flow of Newtonian fluids
    employing the electric potential formulation.
    Buoyancy is treated explicitely with the Boussinesq approximation.
    (suitable for cylcic boundary conditions)
    The Lorentz force term is
    treated in a CONSERVATIVE way according Ni et al. (doi:10.1016/j.jcp.2007.07.023).

Contact
    Francesco Fico <f.fico@lboro.ac.uk> 15/09/2021 

\*---------------------------------------------------------------------------*/
#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "fvOptions.H"
#include "pisoControl.H"
#include "pimpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
	(
		"Transient solver for turbulent flow,"
		" electrical conductive of incompressible fluids.\n"
		".\n"
		"No radiation model"
        "Using electrical potential formulation"
	);

    #include "addCheckCaseOptions.H"
	#include "setRootCaseLists.H"
	#include "createTime.H"
	#include "createMesh.H"
    #include "createControl.H"
    #include "postProcess.H"

	#include "createFields.H"
	#include "createTimeControls.H"
	#include "CourantNo.H"
	#include "setInitialDeltaT.H"
	#include "initContinuityErrs.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    turbulence->validate();

    Info<< "\nStarting time loop\n" << endl;

	//Lorentz force term initialization
	//volVectorField lorentz = sigma * (-fvc::grad(PotE) ^ B0) + sigma * ((U ^ B0) ^ B0);
    lorentz = sigma  * (-fvc::grad(PotE) ^ B0) + sigma * ((U ^ B0) ^ B0);

    while (runTime.run())
    {
        #include "readTimeControls.H"
		#include "CourantNo.H"
		#include "setDeltaT.H"

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        while (pimple.loop())
        {
            // including momentum equation
		    #include "UEqn.H"
            // including electric potential equation
            #include "PotEEqn.H"

		// ----- momentum equation modified ----- //
		// --- Pressure corrector loop (PISO LOOP)
			while (pimple.correct())
			{
				#include "pEqn.H"
			}
			if (pimple.turbCorr())
			{
				laminarTransport.correct();
				turbulence->correct();
			}
        
        }

      
        runTime.write(); //ensure writing of all IObjects with AUTO.WRITE
        runTime.printExecutionTime(Info);
     
    } //end of time loop

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //