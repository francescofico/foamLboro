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

#include "solidWallMixedTemperatureCoupledFvPatchScalarField.H"
#include "mappedRegionMatrix.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "mappedPatchBase.H"

// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidWallMixedTemperatureCoupledFvPatchScalarField::
solidWallMixedTemperatureCoupledFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    TnbrName_("undefined-Tnbr"),
    KName_("undefined-K")
{
    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    this->valueFraction() = 1.0;
}

Foam::solidWallMixedTemperatureCoupledFvPatchScalarField::
solidWallMixedTemperatureCoupledFvPatchScalarField
(
    const solidWallMixedTemperatureCoupledFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    TnbrName_(ptf.TnbrName_),
    KName_(ptf.KName_)
{}

Foam::solidWallMixedTemperatureCoupledFvPatchScalarField::
solidWallMixedTemperatureCoupledFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    TnbrName_(dict.lookupOrDefault<word>("Tnbr", "T")),
    KName_(dict.lookupOrDefault<word>("K", "K"))
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

Foam::solidWallMixedTemperatureCoupledFvPatchScalarField::
solidWallMixedTemperatureCoupledFvPatchScalarField
(
    const solidWallMixedTemperatureCoupledFvPatchScalarField& wtcsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(wtcsf, iF),
    TnbrName_(wtcsf.TnbrName_),
    KName_(wtcsf.KName_)
{}

const Foam::fvPatchScalarField&
Foam::solidWallMixedTemperatureCoupledFvPatchScalarField::K() const
{
    return this->patch().lookupPatchField<volScalarField, scalar>(KName_);
}

void Foam::solidWallMixedTemperatureCoupledFvPatchScalarField::manipulateMatrix
(
    fvMatrix<scalar>& matrix,
    const label regioni,
    const direction
)
{
    addMappedRegionMatrixCoeffs(patch(), matrix, regioni);
}


void Foam::solidWallMixedTemperatureCoupledFvPatchScalarField::updateCoeffs()
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

    const solidWallMixedTemperatureCoupledFvPatchScalarField& nbrField =
        refCast<const solidWallMixedTemperatureCoupledFvPatchScalarField>
        (
            nbrPatch.lookupPatchField<volScalarField, scalar>(TnbrName_)
        );

    // Start with neighbour-sized fields. distribute() handles face ordering,
    // processor exchange and, for nearestPatchFaceAMI, interpolation.
    scalarField nbrIntFld(nbrField.patchInternalField());
    scalarField nbrKDelta(nbrField.K()*nbrPatch.deltaCoeffs());
    mpp.distribute(nbrIntFld);
    mpp.distribute(nbrKDelta);

    const scalarField myKDelta(K()*patch().deltaCoeffs());
    this->refValue() = nbrIntFld;
    this->refGrad() = 0.0;
    this->valueFraction() = nbrKDelta/(nbrKDelta + myKDelta);

    mixedFvPatchScalarField::updateCoeffs();

    if (debug)
    {
        const scalar Q = gSum(K()*patch().magSf()*snGrad());
        Info<< patch().boundaryMesh().mesh().name() << ':'
            << patch().name() << ':'
            << this->internalField().name() << " <- "
            << mpp.sampleRegion() << ':'
            << mpp.samplePatch() << ':'
            << TnbrName_ << " :"
            << " heat transfer rate:" << Q
            << " boundary value"
            << " min:" << gMin(*this)
            << " max:" << gMax(*this)
            << " avg:" << gAverage(*this)
            << endl;
    }

    UPstream::msgType() = oldTag;
}

void Foam::solidWallMixedTemperatureCoupledFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchScalarField::write(os);
    os.writeEntry("Tnbr", TnbrName_);
    os.writeEntry("K", KName_);
}

namespace Foam
{

makePatchTypeField
(
    fvPatchScalarField,
    solidWallMixedTemperatureCoupledFvPatchScalarField
);

} // End namespace Foam

// ************************************************************************* //
