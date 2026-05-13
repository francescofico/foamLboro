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
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
//#include "directMappedPatchBase.H"
//#include "mappedPatchFieldBase.H"
#include "mappedPatchBase.H"
//#include "regionProperties.H"
//#include "basicThermo.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

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
    if (!isA<mappedPatchBase>(this->patch().patch()))
    {
        FatalErrorInFunction
            << "\n    patch type '" << p.type()
            << "' not type '" << mappedPatchBase::typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << internalField().name()
            << " in file " << internalField().objectPath()
            << exit(FatalError);
    }

    fvPatchScalarField::operator=(scalarField("value", dict, p.size()));

    if (dict.found("refValue"))
    {
        // Full restart
        refValue() = scalarField("refValue", dict, p.size());
        refGrad() = scalarField("refGradient", dict, p.size());
        valueFraction() = scalarField("valueFraction", dict, p.size());
    }
    else
    {
        // Start from user entered data. Assume fixedValue.
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


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


const Foam::fvPatchScalarField&
Foam::solidWallMixedPotentialCoupledFvPatchScalarField::sigma() const
{
    return this->patch().lookupPatchField<volScalarField, scalar>(sigmaName_);
}


void Foam::solidWallMixedPotentialCoupledFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    //const polyMesh& mesh = patch().boundaryMesh().mesh();
    int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag +1;

    //const label patchi = patch().index();

    const mappedPatchBase& mpp =
        refCast<const mappedPatchBase>(patch().patch());
    const polyMesh& nbrMesh = mpp.sampleMesh();
    const label samplePatchi = mpp.samplePolyPatch().index();
    const fvPatch& nbrPatch =
        refCast<const fvMesh>(nbrMesh).boundary()[samplePatchi];


    scalarField Tc(patchInternalField());
    scalarField& Tp = *this;

    const solidWallMixedPotentialCoupledFvPatchScalarField&
        nbrField = refCast
            <const solidWallMixedPotentialCoupledFvPatchScalarField>
            (
                nbrPatch.lookupPatchField<volScalarField, scalar>(PotEnbrName_)
            );

     // Swap to obtain full local values of neighbour internal field
    tmp<scalarField> nbrIntFld(new scalarField(nbrField.size(), Zero));
    tmp<scalarField> nbrKDelta(new scalarField(nbrField.size(), Zero));

    nbrIntFld.ref() = nbrField.patchInternalField();
    nbrKDelta.ref() = nbrField.sigma()*nbrPatch.deltaCoeffs();

    mpp.distribute(nbrIntFld.ref());
    mpp.distribute(nbrKDelta.ref());

    tmp<scalarField> myKDelta = sigma()*patch().deltaCoeffs();

       // Both sides agree on
    // - potential value : (myKDelta*fld + nbrKDelta*nbrFld)/(myKDelta+nbrKDelta)
    // - gradient    : (temperature-fld)*delta
    // We've got a degree of freedom in how to implement this in a mixed bc.
    // (what gradient, what fixedValue and mixing coefficient)
    // Two reasonable choices:
    // 1. specify above potential on one side (preferentially the high side)
    //    and above gradient on the other. So this will switch between pure
    //    fixedvalue and pure fixedgradient
    // 2. specify gradient and potential such that the equations are the
    //    same on both sides. This leads to the choice of
    //    - refGradient = zero gradient
    //    - refValue = neighbour value
    //    - mixFraction = nbrKDelta / (nbrKDelta + myKDelta())

    this->refValue() = nbrIntFld();
    this->refGrad() = 0.0;
    this->valueFraction() = nbrKDelta()/(nbrKDelta() + myKDelta());

    mixedFvPatchScalarField::updateCoeffs();

   if (debug)
    {
        scalar Q = gSum(sigma()*patch().magSf()*snGrad());

        Info<< patch().boundaryMesh().mesh().name() << ':'
            << patch().name() << ':'
            << this->internalField().name() << " <- "
            << mpp.sampleRegion() << ':'
            << mpp.samplePatch() << ':'
            << this->internalField().name() << " :"
            << " heat transfer rate:" << Q
            << " walltemperature "
            << " min:" << gMin(Tp)
            << " max:" << gMax(Tp)
            << " avg:" << gAverage(Tp)
            << endl;
    }


    // Restore tag
    UPstream::msgType() = oldTag;

}


void Foam::solidWallMixedPotentialCoupledFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchScalarField::write(os);
    os.writeEntry("PotEnbr", PotEnbrName_);
    os.writeKeyword("PotEnbr")<< PotEnbrName_
        << token::END_STATEMENT << nl;
    os.writeKeyword("sigma") << sigmaName_ << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

makePatchTypeField
(
    fvPatchScalarField,
    solidWallMixedPotentialCoupledFvPatchScalarField
);

} // End namespace Foam

// ************************************************************************* //
