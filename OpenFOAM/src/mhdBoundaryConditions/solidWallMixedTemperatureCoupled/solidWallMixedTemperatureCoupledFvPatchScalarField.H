/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
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

Class
    solidWallMixedTemperatureCoupledFvPatchScalarField

Description
    Mixed boundary condition for temperature at a region interface.
    Both sides use a conductivity-weighted mix of zeroGradient and the
    neighbour region's internal value.

    Example usage:
        myInterfacePatchName
        {
            type                solidWallMixedTemperatureCoupled;
            Tnbr                T;
            K                   K;
            value               uniform 300;
        }

    Requires a mappedWall (or mappedPatchBase-derived patch) with sampleMode
    nearestPatchFace for matching faces or nearestPatchFaceAMI for nonmatching
    faces. Both modes exchange neighbour data across arbitrary decompositions.
    Generic cyclic/processor patches are not region interfaces for this BC.

SourceFiles
    solidWallMixedTemperatureCoupledFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef solidWallMixedTemperatureCoupledFvPatchScalarField_H
#define solidWallMixedTemperatureCoupledFvPatchScalarField_H

#include "mixedFvPatchFields.H"
//#include "mappedPatchFieldBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class solidWallMixedTemperatureCoupledFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class solidWallMixedTemperatureCoupledFvPatchScalarField
:
    public mixedFvPatchScalarField
{
    // Private data

        //- Name of field on the neighbour region
        const word TnbrName_;
        
        //- Name of thermal conductivity field
        const word KName_;

public:

    //- Runtime type information
    TypeName("solidWallMixedTemperatureCoupled");


    // Constructors

        //- Construct from patch and internal field
        solidWallMixedTemperatureCoupledFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        solidWallMixedTemperatureCoupledFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  solidWallMixedTemperatureCoupledFvPatchScalarField onto a new patch
        solidWallMixedTemperatureCoupledFvPatchScalarField
        (
            const solidWallMixedTemperatureCoupledFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );


        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new solidWallMixedTemperatureCoupledFvPatchScalarField
                (
                    *this,
                    this->internalField()
                )
            );
        }

        //- Construct as copy setting internal field reference
        solidWallMixedTemperatureCoupledFvPatchScalarField
        (
            const solidWallMixedTemperatureCoupledFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new solidWallMixedTemperatureCoupledFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

        //- Get corresponding K field
        const fvPatchScalarField& K() const;


        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Insert the neighbour-cell coefficients into the assembled matrix
        virtual void manipulateMatrix
        (
            fvMatrix<scalar>&,
            const label regioni,
            const direction cmpt
        ) override;

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
