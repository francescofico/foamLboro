/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2010-2010 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

\*---------------------------------------------------------------------------*/

#include "solidRegionDiffNo.H"
#include "fvc.H"

Foam::scalar Foam::solidRegionDiffNo
(
    const fvMesh& mesh,
    const Time& runTime,
    const volScalarField& Cp,
    const volScalarField& rho,
    const volScalarField& K
)
{
    /*
    scalar DiNum = 0.0;
    scalar meanDiNum = 0.0;
    */

    //- Can have fluid domains with 0 cells so do not test.
   
    surfaceScalarField KrhoCpbyDelta =
        mesh.surfaceInterpolation::deltaCoeffs()
        * fvc::interpolate(K)
        / fvc::interpolate(Cp*rho);

    const scalar DiNum = max(KrhoCpbyDelta).value()*runTime.deltaTValue();

    const scalar meanDiNum = average(KrhoCpbyDelta).value()*runTime.deltaTValue();
    

    Info<< "Region: " << mesh.name() << " Diffusion Number mean: " << meanDiNum
        << " max: " << DiNum << endl;

    return DiNum;
}


// ************************************************************************* //
