"""Setup and validation for a fully developed duct with conducting Hartmann walls
(cEpotFoam, fluid + solid regions). Shared by hunt_duct and thick_wall_duct.

Geometry (lengths scaled by the half-width a):
    fluid   -1 < y < 1,        -1 < z < 1        (B along y)
    solid   1 < |y| < 1 + t,   -1 < z < 1        (two slabs on the Hartmann walls)
x is periodic with one cell. The side walls z = +-1 are insulating and carry
no solid region. The outer face and the ends of each slab are insulating.
"""
import argparse
from pathlib import Path

import numpy as np

from foamCase import (HEADER, expansion, latest_time, read_internal_field,
                      read_parameters, symmetric_grading, write_parameters)

L = 0.05          # streamwise extent (one periodic cell)


def parse(case, description, defaults):
    params = case/"system"/"caseParameters"
    old = read_parameters(params)
    get = lambda k: old.get(k, defaults[k])
    p = argparse.ArgumentParser(description=description,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--ha", type=float, default=get("Ha"), help="Hartmann number")
    p.add_argument("--thickness", type=float, default=get("wallThickness"),
                   help="Hartmann wall thickness t_w/a")
    g = p.add_mutually_exclusive_group()
    g.add_argument("--c", type=float, help="wall conductance ratio sigma_w t_w/(sigma a); sets sigma_w")
    g.add_argument("--sigma-ratio", type=float, help="sigma_w/sigma")
    p.add_argument("--cells", type=int, default=int(get("cellsHalf")), help="fluid cells per half-width")
    p.add_argument("--wall-cells", type=int, default=int(get("wallCells")), help="solid cells across the wall")
    p.add_argument("--layer-cells", type=float, default=get("layerCells"),
                   help="first y cell = (1/Ha)/layer-cells")
    p.add_argument("--side-cells", type=float, default=get("sideCells"),
                   help="first z cell = (1/sqrt(Ha))/side-cells")
    a = p.parse_args()
    if a.c is not None:
        a.sigma_ratio = a.c/a.thickness
    elif a.sigma_ratio is None:
        a.sigma_ratio = get("sigmaRatio")
    a.c = a.sigma_ratio*a.thickness
    return a


def write(case, a, title):
    Ha, t = a.ha, a.thickness
    dt = 0.25/Ha**2
    end = 5.0/Ha                       # steady after ~3 side-layer times a^2/(nu Ha)
    ny, gy = symmetric_grading(a.cells, 1.0/(Ha*a.layer_cells))
    nz, gz = symmetric_grading(a.cells, 1.0/(Ha**0.5*a.side_cells))
    # Solid: graded towards the interface with the fluid's first-cell size
    first = min(1.0/(Ha*a.layer_cells), t/a.wall_cells)
    ew = expansion(t, a.wall_cells, first)

    write_parameters(case/"system"/"caseParameters", "setup.py", [
        ("Ha", Ha, "Hartmann number B a sqrt(sigma/(rho nu)); a = rho = nu = sigma = 1"),
        ("B", Ha, "applied field magnitude (= Ha in these units)"),
        ("wallThickness", t, "t_w / a"),
        ("sigmaRatio", a.sigma_ratio, "sigma_w / sigma"),
        ("sigmaWall", a.sigma_ratio, "sigma_w (sigma = 1)"),
        ("c", a.c, "wall conductance ratio sigma_w t_w / (sigma a)"),
        ("deltaT", dt, "0.25/Ha^2"),
        ("endTime", end, "5/Ha"),
        ("cellsHalf", float(a.cells), ""),
        ("wallCells", float(a.wall_cells), ""),
        ("layerCells", float(a.layer_cells), ""),
        ("sideCells", float(a.side_cells), ""),
    ], [f"{title} parameters. Regenerate with: python3 setup.py --ha <Ha> --c <c>"])

    (case/"system"/"fluid"/"blockMeshDict").write_text(HEADER.format(script="setup.py") + f"""
scale 1;

// Fluid: square duct -1 < y, z < 1 (B along y), one periodic cell in x.
vertices
(
    (0  -1 -1) ({L} -1 -1) ({L} 1 -1) (0 1 -1)
    (0  -1  1) ({L} -1  1) ({L} 1  1) (0 1  1)
);

blocks
(
    hex (0 1 2 3 4 5 6 7) (1 {ny} {nz}) simpleGrading (1 {gy} {gz})
);

boundary
(
    fluid_to_solid_bottom       // y = -1, conducting Hartmann wall
    {{
        type            mappedWall;
        sampleMode      nearestPatchFace;
        sampleRegion    solid;
        samplePatch     solid_to_fluid_bottom;
        offsetMode      uniform;
        offset          (0 0 0);
        faces           ((0 1 5 4));
    }}
    fluid_to_solid_top          // y = +1, conducting Hartmann wall
    {{
        type            mappedWall;
        sampleMode      nearestPatchFace;
        sampleRegion    solid;
        samplePatch     solid_to_fluid_top;
        offsetMode      uniform;
        offset          (0 0 0);
        faces           ((3 7 6 2));
    }}
    sideWalls                   // z = -1, +1, insulating
    {{
        type wall;
        faces ((0 3 2 1) (4 5 6 7));
    }}
    xPeriodic0
    {{
        type cyclic;
        neighbourPatch xPeriodic1;
        faces ((0 4 7 3));
    }}
    xPeriodic1
    {{
        type cyclic;
        neighbourPatch xPeriodic0;
        faces ((1 2 6 5));
    }}
);

// ************************************************************************* //
""")

    T = 1.0 + t
    (case/"system"/"solid"/"blockMeshDict").write_text(HEADER.format(script="setup.py") + f"""
scale 1;

// Solid: two slabs, 1 < |y| < {T:g}, -1 < z < 1, one periodic cell in x.
// z grading matches the fluid so the interface faces coincide.
vertices
(
    // bottom slab  -{T:g} < y < -1
    (0 {-T} -1) ({L} {-T} -1) ({L} -1 -1) (0 -1 -1)
    (0 {-T}  1) ({L} {-T}  1) ({L} -1  1) (0 -1  1)
    // top slab      1 < y < {T:g}
    (0 1 -1) ({L} 1 -1) ({L} {T} -1) (0 {T} -1)
    (0 1  1) ({L} 1  1) ({L} {T}  1) (0 {T}  1)
);

blocks
(
    hex (0 1 2 3 4 5 6 7)         (1 {a.wall_cells} {nz}) simpleGrading (1 {1.0/ew:.10g} {gz})
    hex (8 9 10 11 12 13 14 15)   (1 {a.wall_cells} {nz}) simpleGrading (1 {ew:.10g} {gz})
);

boundary
(
    solid_to_fluid_bottom
    {{
        type            mappedWall;
        sampleMode      nearestPatchFace;
        sampleRegion    fluid;
        samplePatch     fluid_to_solid_bottom;
        offsetMode      uniform;
        offset          (0 0 0);
        faces           ((3 7 6 2));
    }}
    solid_to_fluid_top
    {{
        type            mappedWall;
        sampleMode      nearestPatchFace;
        sampleRegion    fluid;
        samplePatch     fluid_to_solid_top;
        offsetMode      uniform;
        offset          (0 0 0);
        faces           ((8 9 13 12));
    }}
    outerWalls                  // y = -{T:g}, +{T:g}, insulating
    {{
        type wall;
        faces ((0 1 5 4) (11 15 14 10));
    }}
    wallEnds                    // z = -1, +1 ends of the slabs, insulating
    {{
        type wall;
        faces ((0 3 2 1) (4 5 6 7) (8 11 10 9) (12 13 14 15));
    }}
    xPeriodic0
    {{
        type cyclic;
        neighbourPatch xPeriodic1;
        faces ((0 4 7 3) (8 12 15 11));
    }}
    xPeriodic1
    {{
        type cyclic;
        neighbourPatch xPeriodic0;
        faces ((1 2 6 5) (9 10 14 13));
    }}
);

// ************************************************************************* //
""")
    print(f"Ha = {Ha:g}, t_w/a = {t:g}, sigma_w/sigma = {a.sigma_ratio:g}, c = {a.c:g}; "
          f"fluid {ny} x {nz}, wall {a.wall_cells} cells; deltaT = {dt:.3g}, "
          f"endTime = {end:.3g} ({end/dt:.0f} steps)")


def validate(case, title, tol, log, primary):
    """primary: 'thick' (exact for the resolved wall) or 'thin' (Hunt)."""
    import sys
    from mhdAnalytic import DuctSolution
    from ductValidation import compare

    par = read_parameters(case/"system"/"caseParameters")
    Ha, t, s = par["Ha"], par["wallThickness"], par["sigmaRatio"]
    c = s*t
    exact = DuctSolution(Ha, wall="thick", sigma_ratio=s, thickness=t)
    hunt = DuctSolution(Ha, wall="thin", c=c)
    ref = exact if primary == "thick" else hunt

    tdir = latest_time(case)

    # Potential inside the solid, along the wall mid-plane
    Cs = read_internal_field(tdir/"solid"/"C")
    Ps = read_internal_field(tdir/"solid"/"PotE", len(Cs))
    Cf = read_internal_field(tdir/"fluid"/"C")
    Vf = read_internal_field(tdir/"fluid"/"V", len(Cf))
    Pf = read_internal_field(tdir/"fluid"/"PotE", len(Cf))
    Uf = read_internal_field(tdir/"fluid"/"U", len(Cf))
    ubar = (Uf[:, 0]*Vf).sum()/Vf.sum()
    gauge = (Pf*Vf).sum()/Vf.sum()
    ps = (Ps - gauge)/(ubar*Ha)
    ys, zs = Cs[:, 1], Cs[:, 2]
    top = ys > 0
    ymid = ys[top][np.argmin(np.abs(ys[top] - (1 + t/2)))]
    sel = top & np.isclose(ys, ymid)
    order = np.argsort(zs[sel])
    pa_wall = exact.wall_potential(ys[sel], zs[sel])
    e_wall = np.abs(ps[sel] - pa_wall).max()/np.abs(pa_wall).max()

    def plot(ax):
        zz = np.linspace(-1, 1, 801)
        ax.plot(zz, exact.wall_potential(np.full_like(zz, ymid), zz), "k-", label="exact (thick wall)")
        ax.plot(zs[sel][order], ps[sel][order], "o", ms=3, mfc="none", label="cEpotFoam (solid)")
        ax.set(xlabel="z / a", ylabel="phi_w / (ubar B a)",
               title=f"wall potential, y = {ymid:.4g}")
        ax.legend()

    print(f"\nFor context: G a^2/(rho nu ubar) = {hunt.gradp:.6g} for Hunt's thin wall (c = {c:g}), "
          f"{exact.gradp:.6g} for the exact wall of thickness {t:g} "
          f"({100*(exact.gradp/hunt.gradp - 1):+.2f}%)."
          f" Tested against: {'Hunt' if primary == 'thin' else 'exact thick wall'}.")

    ok = compare(case, ref, log, title, tol, region="fluid", B=Ha, extra_plot=plot)
    flag = "  <-- above tolerance" if e_wall > tol else ""
    print(f"{'wall potential error (max, rel.)':38s} {e_wall:14.3e}{flag}")
    ok = ok and e_wall <= tol
    print("\nPASS" if ok else "\nFAIL")
    sys.exit(0 if ok else 1)
