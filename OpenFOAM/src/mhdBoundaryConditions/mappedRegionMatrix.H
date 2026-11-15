/*---------------------------------------------------------------------------*\
  Implicit scalar interface coefficients for OpenFOAM.com v2606.
  Based on OpenCFD's mixedEnergyFvPatchScalarField matrix assembly (GPLv3).
\*---------------------------------------------------------------------------*/
#ifndef mappedRegionMatrix_H
#define mappedRegionMatrix_H

#include "fvMatrix.H"

namespace Foam
{
inline void addMappedRegionMatrixCoeffs
(
    const fvPatch& patch,
    fvScalarMatrix& matrix,
    const label regioni
)
{
    const auto& assembly = matrix.lduMeshAssembly();
    const label patchi = patch.index();
    const label nbrPatchi = patch.patch().neighbPolyPatchID();
    const label nbri = assembly.findNbrMeshId(patch.patch(), regioni);
    const label globalPatchi = assembly.patchLocalToGlobalMap()[regioni][patchi];
    const scalarField& coeffs = matrix.internalCoeffs()[globalPatchi];
    const labelList& faces = assembly.faceBoundMap()[regioni][patchi];
    const labelList& localCells = assembly.cellBoundMap()[nbri][nbrPatchi];
    const labelList& neighbourCells = assembly.cellBoundMap()[regioni][patchi];
    const labelList& patchFaces = assembly.facePatchFaceMap()[regioni][patchi];

    forAll(faces, i)
    {
        const scalar conductance = coeffs[patchFaces[i]];
        const label celli = localCells[i];
        const label nbrCelli = neighbourCells[i];
        matrix.diag()[celli] += conductance;
        if (celli < nbrCelli)
        {
            matrix.upper()[faces[i]] -= conductance;
        }
        else if (matrix.asymmetric())
        {
            matrix.lower()[faces[i]] -= conductance;
        }
    }
}
}
#endif
