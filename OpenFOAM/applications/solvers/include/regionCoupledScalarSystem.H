/*---------------------------------------------------------------------------*\
  Optional assembled scalar solves across MHD regions. GPLv3 or later.
\*---------------------------------------------------------------------------*/
#ifndef regionCoupledScalarSystem_H
#define regionCoupledScalarSystem_H

#include "fvCFD.H"
#include "mappedPatchBase.H"
#include "findRefCell.H"

namespace Foam
{
class regionCoupledScalarSystem
{
    UPtrList<volScalarField> fields_;
    bool enabled_;
    bool potential_;
    bool refreshAssembly_;
    label added_;
    autoPtr<fvScalarMatrix> matrix_;
    labelList referenceCells_;
    scalarField referenceValues_;

public:
    regionCoupledScalarSystem
    (
        const UPtrList<volScalarField>& fields,
        const dictionary& controls,
        const word& option,
        const word& patchType,
        const bool potential
    )
    :
        fields_(fields),
        enabled_(controls.lookupOrDefault<bool>
        (
            option, controls.lookupOrDefault<bool>("coupled", false)
        )),
        potential_(potential),
        refreshAssembly_(false),
        added_(0),
        referenceCells_(fields.size(), -1),
        referenceValues_(fields.size(), Zero)
    {
        // An explicit solver switch takes precedence over saved restart flags.
        if (!controls.found(option) && !controls.found("coupled"))
        {
            forAll(fields_, i)
            {
                for (const auto& pf : fields_[i].boundaryField())
                {
                    enabled_ = enabled_ || pf.useImplicit();
                }
            }
        }

        forAll(fields_, i)
        {
            for (auto& pf : fields_[i].boundaryFieldRef())
            {
                if (pf.type() == patchType)
                {
                    pf.useImplicit(enabled_);
                }
                else if (pf.useImplicit())
                {
                    FatalErrorInFunction
                        << "Unsupported implicit condition " << pf.type()
                        << " for " << fields_[i].name() << exit(FatalError);
                }
            }
        }

        if (!enabled_) return;
        labelList component(fields_.size());
        forAll(component, i) component[i] = i;
        boolList anchored(fields_.size(), false);
        label nInterfaces = 0;

        forAll(fields_, i)
        {
            const auto& field = fields_[i];
            for (const auto& pf : field.boundaryField())
            {
                if (!pf.useImplicit())
                {
                    anchored[i] = anchored[i] || pf.fixesValue();
                    continue;
                }
                ++nInterfaces;
                const auto& mpp = refCast<const mappedPatchBase>(pf.patch().patch());
                if (!mpp.sameWorld() || mpp.mode() != mappedPatchBase::NEARESTPATCHFACE)
                {
                    FatalErrorInFunction
                        << "Implicit region coupling requires nearestPatchFace "
                        << "in the same simulation. Use separate-region solves "
                        << "for nearestPatchFaceAMI." << exit(FatalError);
                }
                label nbri = -1;
                forAll(fields_, j)
                {
                    if (fields_[j].mesh().name() == mpp.sampleRegion()) nbri = j;
                }
                if (nbri < 0)
                {
                    FatalErrorInFunction << "Neighbour region " << mpp.sampleRegion()
                        << " is absent from the assembled system" << exit(FatalError);
                }
                const fvPatch& nbrPatch = fields_[nbri].mesh().boundary()
                [
                    mpp.samplePolyPatch().index()
                ];
                const auto& nbr = fields_[nbri].boundaryField()[nbrPatch.index()];
                if (!nbr.useImplicit() || nbr.type() != patchType)
                {
                    FatalErrorInFunction << "Both sides of " << pf.patch().name()
                        << " must use " << patchType << exit(FatalError);
                }
                const auto& reverse = refCast<const mappedPatchBase>(nbrPatch.patch());
                if
                (
                    reverse.sampleRegion() != field.mesh().name()
                 || reverse.samplePatch() != pf.patch().name()
                )
                {
                    FatalErrorInFunction << "Non-reciprocal region interface "
                        << pf.patch().name() << exit(FatalError);
                }

                // Native mappedWall assembly pairs local face i with face i.
                bool aligned = pf.size() == nbr.size();
                if (aligned)
                {
                    forAll(pf, facei)
                    {
                        const scalar area = pf.patch().magSf()[facei];
                        aligned = aligned
                            && mag(pf.patch().Cf()[facei] - nbrPatch.Cf()[facei])
                               <= 1e-8*sqrt(area)
                            && mag(pf.patch().Sf()[facei] + nbrPatch.Sf()[facei])
                               <= 1e-8*area;
                    }
                }
                if (!returnReduce(aligned, andOp<bool>()))
                {
                    FatalErrorInFunction
                        << "Implicit interface " << field.mesh().name() << ':'
                        << pf.patch().name() << " requires matching faces in "
                        << "the same local order on the same processor. "
                        << "Align the patches and their decomposition, or disable "
                        << option << " to use mapped separate-region solves."
                        << exit(FatalError);
                }
                const label from = component[nbri], to = component[i];
                forAll(component, j)
                {
                    if (component[j] == from) component[j] = to;
                }
            }
        }
        if (!nInterfaces)
        {
            FatalErrorInFunction << option << " requires a mapped region interface"
                << exit(FatalError);
        }

        if (potential_)
        {
            // One gauge per connected group without an external fixed potential.
            forAll(fields_, i)
            {
                bool first = true, fixed = false;
                forAll(fields_, j)
                {
                    if (component[j] == component[i])
                    {
                        first = first && j >= i;
                        fixed = fixed || anchored[j];
                    }
                }
                fixed = returnReduce(fixed, orOp<bool>());
                if (first && !fixed)
                {
                    dictionary refs(fields_[i].mesh().solutionDict().subOrEmptyDict("PotE"));
                    const word name = fields_[i].name();
                    if (!refs.found(name + "RefValue"))
                    {
                        refs.add(word(name + "RefValue"), scalar(0));
                    }
                    if (refs.found(name + "RefCell") || refs.found(name + "RefPoint"))
                    {
                        setRefCell(fields_[i], refs, referenceCells_[i], referenceValues_[i], true);
                    }
                    else
                    {
                        const label rank = returnReduce
                        (
                            fields_[i].size() ? Pstream::myProcNo() : labelMax,
                            minOp<label>()
                        );
                        referenceCells_[i] = rank == Pstream::myProcNo() ? 0 : -1;
                        referenceValues_[i] = refs.get<scalar>(name + "RefValue");
                    }
                    Info<< "Coupled potential reference in region "
                        << fields_[i].mesh().name() << ": " << referenceValues_[i] << nl;
                }
            }
        }
        Info<< option << " = " << enabled_ << nl;
    }

    bool enabled() const { return enabled_; }

    void shareAssemblyWith(regionCoupledScalarSystem& other)
    {
        if (!enabled_ || !other.enabled_) return;
        // OpenFOAM names the cached assembly from the first implicit patch
        // indices, without the field name. Refresh if the two graphs differ.
        forAll(fields_, i)
        {
            forAll(fields_[i].boundaryField(), patchi)
            {
                if
                (
                    fields_[i].boundaryField()[patchi].useImplicit()
                 != other.fields_[i].boundaryField()[patchi].useImplicit()
                )
                {
                    refreshAssembly_ = other.refreshAssembly_ = true;
                }
            }
        }
    }

    void add(fvScalarMatrix& equation)
    {
        if (!enabled_ || added_ >= fields_.size() || &equation.psi() != &fields_[added_])
        {
            FatalErrorInFunction << "Missing, duplicate or unordered region equation"
                << exit(FatalError);
        }
        if (!potential_)
        {
            // Temperature equations use diffusivity, but the interface weights
            // use conductivity. Their two row coefficients can differ.
            equation.lower();
            equation.upper();
        }
        const label ref = referenceCells_[added_];
        if (ref >= 0)
        {
            scalar diagonal = equation.diag()[ref];
            forAll(equation.internalCoeffs(), patchi)
            {
                const labelUList& cells = equation.psi().mesh().boundary()[patchi].faceCells();
                forAll(cells, facei)
                {
                    if (cells[facei] == ref) diagonal += equation.internalCoeffs()[patchi][facei];
                }
            }
            equation.diag()[ref] += diagonal;
            equation.source()[ref] += diagonal*referenceValues_[added_];
        }
        if (!matrix_)
        {
            matrix_.reset(new fvScalarMatrix(equation.psi(), equation.dimensions()));
        }
        matrix_->addFvMatrix(equation);
        ++added_;
    }

    void solve(const bool finalIter)
    {
        if (!enabled_) return;
        if (added_ != fields_.size())
        {
            FatalErrorInFunction << "Expected " << fields_.size()
                << " region equations, received " << added_ << exit(FatalError);
        }
        const fvMesh& mesh = fields_[0].mesh();
        fvSolution controls(mesh.time());
        const word fieldName = fields_[0].name();
        const word solverName = fieldName + (finalIter ? "Final" : "");
        const dictionary* solver = nullptr;
        if (controls.found("solvers"))
        {
            const dictionary& solvers = controls.subDict("solvers");
            if (solvers.found(solverName)) solver = &solvers.subDict(solverName);
            else if (solvers.found(fieldName)) solver = &solvers.subDict(fieldName);
        }
        Info<< "Solving coupled " << fieldName << " across " << added_ << " regions" << nl;
        if (refreshAssembly_ && matrix_->lduMeshPtr())
        {
            // update() retains the old patch maps: a changed connectivity
            // needs a fresh registry-owned assembly, not a geometry update.
            matrix_->lduMeshPtr()->checkOut();
        }
        matrix_->solve(solver ? *solver : mesh.solver(solverName));
        // Refresh after all internal fields have been unpacked by the assembly.
        for (auto& field : fields_) field.correctBoundaryConditions();
        matrix_.clear();
        added_ = 0;
    }
};
}
#endif
