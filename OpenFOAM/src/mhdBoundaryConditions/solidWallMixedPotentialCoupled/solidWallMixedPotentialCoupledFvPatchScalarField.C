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

\*---------------------------------------------------------------------------*/

#include "solidWallMixedPotentialCoupledFvPatchScalarField.H"
#include "mappedRegionMatrix.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "mappedPatchBase.H"

// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidWallMixedPotentialCoupledFvPatchScalarField::
solidWallMixedPotentialCoupledFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    PotEnbrName_("undefined-PotEnbr"),
    sigmaName_("undefined-sigma")
{
    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    this->valueFraction() = 1.0;
}

Foam::solidWallMixedPotentialCoupledFvPatchScalarField::
solidWallMixedPotentialCoupledFvPatchScalarField
(
    const solidWallMixedPotentialCoupledFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    PotEnbrName_(ptf.PotEnbrName_),
    sigmaName_(ptf.sigmaName_)
{}

Foam::solidWallMixedPotentialCoupledFvPatchScalarField::
solidWallMixedPotentialCoupledFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    PotEnbrName_(dict.lookupOrDefault<word>("PotEnbr", "PotE")),
    sigmaName_(dict.lookupOrDefault<word>("sigma", "sigma"))
{
    if (!isA<mappedPatchBase>(p.patch()))
    {
        FatalIOErrorInFunction(dict)
            << "\n    patch type '" << p.type()
            << "' does not provide neighbour mapping. Use mappedWall with"
            << " sampleMode nearestPatchFace (matching faces) or"
            << " nearestPatchFaceAMI (nonmatching faces)"
            << "\n    for patch " << p.name()
            << " of field " << internalField().name()
            << " in file " << internalField().objectPath()
            << exit(FatalIOError);
    }

    const mappedPatchBase& mpp = refCast<const mappedPatchBase>(p.patch());
    if
    (
        mpp.mode() != mappedPatchBase::NEARESTPATCHFACE
     && mpp.mode() != mappedPatchBase::NEARESTPATCHFACEAMI
    )
    {
        FatalIOErrorInFunction(dict)
            << "Patch " << p.name()
            << " requires sampleMode nearestPatchFace or nearestPatchFaceAMI"
            << exit(FatalIOError);
    }

    useImplicit(dict.lookupOrDefault<bool>("useImplicit", false));

    fvPatchScalarField::operator=(scalarField("value", dict, p.size()));

    if (dict.found("refValue"))
    {
        refValue() = scalarField("refValue", dict, p.size());
        refGrad() = scalarField("refGradient", dict, p.size());
        valueFraction() = scalarField("valueFraction", dict, p.size());
    }
    else
    {
        refValue() = *this;
        refGrad() = 0.0;
        valueFraction() = 1.0;
    }
}

Foam::solidWallMixedPotentialCoupledFvPatchScalarField::
solidWallMixedPotentialCoupledFvPatchScalarField
(
    const solidWallMixedPotentialCoupledFvPatchScalarField& wtcsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(wtcsf, iF),
    PotEnbrName_(wtcsf.PotEnbrName_),
    sigmaName_(wtcsf.sigmaName_)
{}

const Foam::fvPatchScalarField&
Foam::solidWallMixedPotentialCoupledFvPatchScalarField::sigma() const
{
    return this->patch().lookupPatchField<volScalarField, scalar>(sigmaName_);
}

void Foam::solidWallMixedPotentialCoupledFvPatchScalarField::manipulateMatrix
(
    fvMatrix<scalar>& matrix,
    const label regioni,
    const direction
)
{
    addMappedRegionMatrixCoeffs(patch(), matrix, regioni);
}


void Foam::solidWallMixedPotentialCoupledFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Keep mapping communication separate from ongoing field evaluations.
    const int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag + 1;

    const mappedPatchBase& mpp =
        refCast<const mappedPatchBase>(patch().patch());
    const fvMesh& nbrMesh = refCast<const fvMesh>(mpp.sampleMesh());
    const fvPatch& nbrPatch = nbrMesh.boundary()[mpp.samplePolyPatch().index()];

    const solidWallMixedPotentialCoupledFvPatchScalarField& nbrField =
        refCast<const solidWallMixedPotentialCoupledFvPatchScalarField>
        (
            nbrPatch.lookupPatchField<volScalarField, scalar>(PotEnbrName_)
        );

    // Start with neighbour-sized fields. distribute() handles face ordering,
    // processor exchange and, for nearestPatchFaceAMI, interpolation.
    scalarField nbrIntFld(nbrField.patchInternalField());
    scalarField nbrKDelta(nbrField.sigma()*nbrPatch.deltaCoeffs());
    mpp.distribute(nbrIntFld);
    mpp.distribute(nbrKDelta);

    const scalarField myKDelta(sigma()*patch().deltaCoeffs());
    this->refValue() = nbrIntFld;
    this->refGrad() = 0.0;
    this->valueFraction() = nbrKDelta/(nbrKDelta + myKDelta);

    mixedFvPatchScalarField::updateCoeffs();

    if (debug)
    {
        const scalar Q = gSum(sigma()*patch().magSf()*snGrad());
        Info<< patch().boundaryMesh().mesh().name() << ':'
            << patch().name() << ':'
            << this->internalField().name() << " <- "
            << mpp.sampleRegion() << ':'
            << mpp.samplePatch() << ':'
            << PotEnbrName_ << " :"
            << " electric current rate:" << Q
            << " boundary value"
            << " min:" << gMin(*this)
            << " max:" << gMax(*this)
            << " avg:" << gAverage(*this)
            << endl;
    }

    UPstream::msgType() = oldTag;
}

void Foam::solidWallMixedPotentialCoupledFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchScalarField::write(os);
    os.writeEntry("PotEnbr", PotEnbrName_);
    os.writeEntry("sigma", sigmaName_);
}

namespace Foam
{

makePatchTypeField
(
    fvPatchScalarField,
    solidWallMixedPotentialCoupledFvPatchScalarField
);

} // End namespace Foam

// ************************************************************************* //
