#!/usr/bin/env python3
"""Write the fringing magnetic field of the benchmark into 0/fluid/B0:

    B_y(x) = B0/2 [1 - tanh((|x| - 0.7315)/0.1)]          (paper, eq. 1)

Needs the cell centres in 0/fluid/C (postProcess -region fluid -func
writeCellCentres -time 0). Only B_y is set, as in the benchmark definition
(the field is not curl-free; the inductionless model takes it as given).
"""
import os
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, os.environ.get("MHD_COMMON", str(HERE.parents[1]/"common")))
from foamCase import read_internal_field, read_parameters  # noqa: E402

par = read_parameters(HERE/"system"/"caseParameters")
B0 = par["B0"]
C = read_internal_field(HERE/"0"/"fluid"/"C")
By = 0.5*B0*(1.0 - np.tanh((np.abs(C[:, 0]) - 0.7315)/0.1))
body = "\n".join(f"(0 {b:.10g} 0)" for b in By)
(HERE/"0"/"fluid"/"B0").write_text(f"""/*--------------------------------*- C++ -*----------------------------------*\\
| Written by setB0.py: B_y(x) = B0/2 [1 - tanh((|x| - 0.7315)/0.1)]          |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       volVectorField;
    location    "0/fluid";
    object      B0;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [1 0 -2 0 0 -1 0];

internalField   nonuniform List<vector>
{len(By)}
(
{body}
)
;

boundaryField
{{
    ".*"
    {{
        type            zeroGradient;
    }}
}}

// ************************************************************************* //
""")
print(f"B0: {len(By)} cells, B_y from {By.min():.4g} to {By.max():.4g} T")
