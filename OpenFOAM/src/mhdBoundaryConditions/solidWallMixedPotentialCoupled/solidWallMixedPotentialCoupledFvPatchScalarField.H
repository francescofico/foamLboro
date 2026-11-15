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
    solidWallMixedPotentialCoupledFvPatchScalarField

Description
    Mixed boundary condition for electric potential at a region interface.
    Both sides use a conductivity-weighted mix of zeroGradient and the
    neighbour region's internal value.

    Example usage:
        myInterfacePatchName
        {
            type                solidWallMixedPotentialCoupled;
            PotEnbr             PotE;
            sigma               sigma;
            value               uniform 0;
        }

    Requires a mappedWall (or mappedPatchBase-derived patch) with sampleMode
    nearestPatchFace for matching faces or nearestPatchFaceAMI for nonmatching
    faces. Both modes exchange neighbour data across arbitrary decompositions.
    Generic cyclic/processor patches are not region interfaces for this BC.

SourceFiles
    solidWallMixedPotentialCoupledFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef solidWallMixedPotentialCoupledFvPatchScalarField_H
#define solidWallMixedPotentialCoupledFvPatchScalarField_H

#include "mixedFvPatchFields.H"
//#include "mappedPatchFieldBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class solidWallMixedPotentialCoupledFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class solidWallMixedPotentialCoupledFvPatchScalarField
:
    public mixedFvPatchScalarField
{
    // Private data

        //- Name of field on the neighbour region
        const word PotEnbrName_;
        
        //- Name of electrical conductivity field
        const word sigmaName_;

public:

    //- Runtime type information
    TypeName("solidWallMixedPotentialCoupled");


    // Constructors

        //- Construct from patch and internal field
        solidWallMixedPotentialCoupledFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        solidWallMixedPotentialCoupledFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  solidWallMixedPotentialCoupledFvPatchScalarField onto a new patch
        solidWallMixedPotentialCoupledFvPatchScalarField
        (
            const solidWallMixedPotentialCoupledFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );


        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new solidWallMixedPotentialCoupledFvPatchScalarField
                (
                    *this,
                    this->internalField()
                )
            );
        }

        //- Construct as copy setting internal field reference
        solidWallMixedPotentialCoupledFvPatchScalarField
        (
            const solidWallMixedPotentialCoupledFvPatchScalarField&,
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
                new solidWallMixedPotentialCoupledFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

        //- Get corresponding electrical conductivity field
        const fvPatchScalarField& sigma() const;


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
