/*---------------------------------------------------------------------------*\
  Residual control for the MHD multi-region PIMPLE loops. GPLv3 or later.

  Uses the final-correction approach of OpenFOAM-dev commit a884a44e4729,
  with v2606's recorded solver performance instead of per-equation storage.
\*---------------------------------------------------------------------------*/

#ifndef regionResidualControl_H
#define regionResidualControl_H

#include "solutionControl.H"
#include <cmath>

namespace Foam
{

class regionResidualControl
{
    const PtrList<fvMesh>& fluidRegions_;
    const PtrList<fvMesh>& solidRegions_;
    bool converged_ = false;

    bool checkRegions
    (
        const PtrList<fvMesh>& regions,
        const bool firstOuterIter,
        bool& used
    ) const
    {
        bool achieved = true;

        for (const fvMesh& mesh : regions)
        {
            const dictionary pimple =
                mesh.solutionDict().subOrEmptyDict("PIMPLE");
            const dictionary controls = pimple.subOrEmptyDict("residualControl");
            const dictionary& performances = mesh.data().solverPerformanceDict();

            for (const entry& control : controls)
            {
                used = true;
                const dictionary& tolerances = control.dict();
                const scalar absTol = tolerances.get<scalar>("tolerance");
                const scalar relTol = tolerances.getOrDefault<scalar>("relTol", 0);
                if
                (
                    !std::isfinite(absTol) || absTol < 0
                 || !std::isfinite(relTol) || relTol < 0
                )
                {
                    FatalIOErrorInFunction(tolerances)
                        << "Residual tolerances must be finite and non-negative"
                        << exit(FatalIOError);
                }

                bool checked = false;
                for (const entry& performance : performances)
                {
                    const word& fieldName = performance.keyword();
                    if (!control.keyword().match(fieldName)) continue;

                    // Do not reuse a flow residual if its solve was disabled
                    // by a runtime change to the regional PIMPLE controls.
                    if
                    (
                        (fieldName == "U"
                      && !pimple.getOrDefault("momentumPredictor", true))
                     || ((fieldName == "p" || fieldName == "p_rgh")
                      && pimple.getOrDefault<label>("nCorrectors", 1) < 1)
                    )
                    {
                        achieved = false;
                        continue;
                    }

                    const Pair<scalar> residuals =
                        solutionControl::maxResidual(mesh, performance);
                    // v2606 retains the first and latest initial residuals
                    // in this time step, including assembled scalar solves.
                    const bool absolute = residuals.last() < absTol;
                    const bool relative =
                        !firstOuterIter && relTol > 0
                     && residuals.last()/(residuals.first() + ROOTVSMALL) < relTol;
                    checked = true;
                    achieved = achieved
                     && std::isfinite(residuals.first())
                     && std::isfinite(residuals.last())
                     && (absolute || relative);
                }

                if (!checked)
                {
                    achieved = false;
                    if (firstOuterIter)
                    {
                        WarningInFunction
                            << "No solver residual for " << control.keyword()
                            << " in region " << mesh.name()
                            << "; this residual criterion cannot be satisfied. "
                            << "Check the field name and whether its equation "
                            << "is enabled." << nl;
                    }
                }
            }
        }

        return achieved;
    }

public:

    // Construct anew each time step to reset the convergence decision.
    regionResidualControl
    (
        const PtrList<fvMesh>& fluidRegions,
        const PtrList<fvMesh>& solidRegions
    )
    :
        fluidRegions_(fluidRegions),
        solidRegions_(solidRegions)
    {}

    bool converged() const { return converged_; }

    // Called after all solves in an outer correction. Return true to stop.
    bool check(const label outerCorr, const label nOuterCorr)
    {
        if (converged_)
        {
            Info<< "PIMPLE: completed final correction after " << outerCorr + 1
                << " outer correctors" << nl;
            return true;
        }

        bool used = false;
        bool achieved = checkRegions(fluidRegions_, outerCorr == 0, used);
        achieved = checkRegions(solidRegions_, outerCorr == 0, used) && achieved;
        if (!used) return false;

        if (outerCorr == nOuterCorr - 1)
        {
            Info<< "PIMPLE: residual criteria "
                << (achieved ? "satisfied" : "not satisfied") << " after "
                << nOuterCorr << " outer correctors (iteration limit)" << nl;
        }
        else if (achieved)
        {
            converged_ = true;
            Info<< "PIMPLE: residual criteria satisfied after " << outerCorr + 1
                << " outer correctors; performing final correction" << nl;
        }

        return false;
    }
};

} // End namespace Foam

#endif
