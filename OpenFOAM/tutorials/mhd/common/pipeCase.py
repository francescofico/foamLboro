"""Setup and validation for fully developed flow in a circular pipe in a
transverse field. Shared by epotFoam_isoTh/pipe_insulating and
cEpotFoam/pipe_conducting.

Geometry (lengths scaled by the radius a): fluid r < 1 in the y-z plane, B
along y, x periodic with one cell. The fluid mesh is an O-grid: a central
square block plus four outer blocks graded towards the wall. For the
conducting case the wall 1 < r < 1 + t is a solid region with matching faces.
"""
import argparse
import sys

import numpy as np

from foamCase import (HEADER, expansion, latest_time, pressure_gradient_history,
                      read_internal_field, read_parameters, report, write_parameters)

L = 0.05            # streamwise extent (one periodic cell)
S = 0.35            # half-size of the central square block (corners at y, z = +-S)
Q = 1.0/np.sqrt(2.0)


# ----------------------------------------------------------------------------
# Parameters
# ----------------------------------------------------------------------------
def parse(case, description, defaults, conducting):
    old = read_parameters(case/"system"/"caseParameters")
    get = lambda k: old.get(k, defaults[k])
    p = argparse.ArgumentParser(description=description,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--ha", type=float, default=get("Ha"), help="Hartmann number (radius-based)")
    p.add_argument("--cells", type=int, default=int(get("cellsQuarter")),
                   help="cells per quarter circumference (and across the central block)")
    p.add_argument("--radial-cells", type=int, default=int(get("radialCells")),
                   help="cells from the central block to the wall")
    p.add_argument("--layer-cells", type=float, default=get("layerCells"),
                   help="first wall cell = (1/Ha)/layer-cells")
    if conducting:
        p.add_argument("--thickness", type=float, default=get("wallThickness"), help="t_w/a")
        g = p.add_mutually_exclusive_group()
        g.add_argument("--c", type=float, help="thin-wall conductance ratio sigma_w t_w/(sigma a); sets sigma_w")
        g.add_argument("--sigma-ratio", type=float, help="sigma_w/sigma")
        p.add_argument("--wall-cells", type=int, default=int(get("wallCells")), help="cells across the wall")
    a = p.parse_args()
    a.conducting = conducting
    if conducting:
        if a.c is not None:
            a.sigma_ratio = a.c/a.thickness
        elif a.sigma_ratio is None:
            a.sigma_ratio = get("sigmaRatio")
        a.c = a.sigma_ratio*a.thickness
    return a


# ----------------------------------------------------------------------------
# blockMesh
# ----------------------------------------------------------------------------
def _ring(rad):
    """Vertices (y, z) of a ring at radius rad: A', B', C', D' at 45 degrees."""
    q = rad*Q
    return [(-q, -q), (q, -q), (q, q), (-q, q)]


def _vertices(pts):
    out = []
    for x in (0.0, L):
        for y, z in pts:
            out.append(f"    ({x:g} {y:.12g} {z:.12g})")
    return "\n".join(out)


# Outer blocks in order z<0, y>0, z>0, y<0. Each entry gives the ring indices
# (p0, p1) of the inner edge (dir1 runs p0 -> p1) so the hex is right-handed
# with dir2 pointing outwards and dir3 along x.
QUARTERS = [(1, 0), (2, 1), (3, 2), (0, 3)]
MIDPOINTS = [(0.0, -1.0), (1.0, 0.0), (0.0, 1.0), (-1.0, 0.0)]   # (y, z) on the unit circle


def _hex(inner, outer, i0, i1, nv):
    """Hex between ring `inner` and ring `outer` (vertex offsets) for edge i0->i1."""
    a, b = inner + i0, inner + i1
    c, d = outer + i1, outer + i0
    return f"{a} {b} {c} {d} {a + nv} {b + nv} {c + nv} {d + nv}"


def fluid_blockmesh(a, wall_patch):
    Ha = a.ha
    n, nr = a.cells, a.radial_cells
    first = 1.0/(Ha*a.layer_cells)
    e = expansion(1.0 - S, nr, first)          # fine at the wall
    pts = [(-S, -S), (S, -S), (S, S), (-S, S)] + _ring(1.0)
    nv = len(pts)
    blocks = [f"    hex (0 1 2 3 {nv} {nv + 1} {nv + 2} {nv + 3}) ({n} {n} 1) simpleGrading (1 1 1)"]
    for i0, i1 in QUARTERS:
        blocks.append(f"    hex ({_hex(0, 4, i0, i1, nv)}) ({n} {nr} 1) simpleGrading (1 {1.0/e:.10g} 1)")
    edges = []
    for k, (i0, i1) in enumerate(QUARTERS):
        y, z = MIDPOINTS[k]
        for off, x in ((0, 0.0), (nv, L)):
            edges.append(f"    arc {4 + i0 + off} {4 + i1 + off} ({x:g} {y:g} {z:g})")
    wall_faces, x0, x1 = [], [], []
    for blk in blocks:
        p = [int(v) for v in blk.split("(")[1].split(")")[0].split()]
        x0.append(f"({p[0]} {p[3]} {p[2]} {p[1]})")
        x1.append(f"({p[4]} {p[5]} {p[6]} {p[7]})")
    for blk in blocks[1:]:
        p = [int(v) for v in blk.split("(")[1].split(")")[0].split()]
        wall_faces.append(f"({p[3]} {p[7]} {p[6]} {p[2]})")
    return HEADER.format(script="setup.py") + f"""
scale 1;

// Fluid: pipe r < 1 (B along y), O-grid, one periodic cell in x.
vertices
(
{_vertices(pts)}
);

blocks
(
{chr(10).join(blocks)}
);

edges
(
{chr(10).join(edges)}
);

boundary
(
{wall_patch(" ".join(wall_faces))}
    xPeriodic0
    {{
        type cyclic;
        neighbourPatch xPeriodic1;
        faces ({" ".join(x0)});
    }}
    xPeriodic1
    {{
        type cyclic;
        neighbourPatch xPeriodic0;
        faces ({" ".join(x1)});
    }}
);

// ************************************************************************* //
""", (n, nr, first, e)


def solid_blockmesh(a):
    n, nw, t = a.cells, a.wall_cells, a.thickness
    first = min(1.0/(a.ha*a.layer_cells), t/nw)
    e = expansion(t, nw, first)                 # fine at the fluid interface
    pts = _ring(1.0) + _ring(1.0 + t)
    nv = len(pts)
    blocks, edges, inner, outer, x0, x1 = [], [], [], [], [], []
    for k, (i0, i1) in enumerate(QUARTERS):
        h = _hex(0, 4, i0, i1, nv)
        blocks.append(f"    hex ({h}) ({n} {nw} 1) simpleGrading (1 {e:.10g} 1)")
        p = [int(v) for v in h.split()]
        inner.append(f"({p[0]} {p[1]} {p[5]} {p[4]})")
        outer.append(f"({p[3]} {p[7]} {p[6]} {p[2]})")
        x0.append(f"({p[0]} {p[3]} {p[2]} {p[1]})")
        x1.append(f"({p[4]} {p[5]} {p[6]} {p[7]})")
        y, z = MIDPOINTS[k]
        for ring, rad in ((0, 1.0), (4, 1.0 + t)):
            for off, x in ((0, 0.0), (nv, L)):
                edges.append(f"    arc {ring + i0 + off} {ring + i1 + off} ({x:g} {rad*y:.12g} {rad*z:.12g})")
    return HEADER.format(script="setup.py") + f"""
scale 1;

// Solid: pipe wall 1 < r < {1 + t:g}, four blocks matching the fluid O-grid.
vertices
(
{_vertices(pts)}
);

blocks
(
{chr(10).join(blocks)}
);

edges
(
{chr(10).join(edges)}
);

boundary
(
    solid_to_fluid
    {{
        type            mappedWall;
        sampleMode      nearestPatchFace;
        sampleRegion    fluid;
        samplePatch     fluid_to_solid;
        offsetMode      uniform;
        offset          (0 0 0);
        faces           ({" ".join(inner)});
    }}
    outerWall                   // r = {1 + t:g}, insulating
    {{
        type wall;
        faces ({" ".join(outer)});
    }}
    xPeriodic0
    {{
        type cyclic;
        neighbourPatch xPeriodic1;
        faces ({" ".join(x0)});
    }}
    xPeriodic1
    {{
        type cyclic;
        neighbourPatch xPeriodic0;
        faces ({" ".join(x1)});
    }}
);

// ************************************************************************* //
"""


def write(case, a, title):
    Ha = a.ha
    # Two limits on the time step:
    #  - the Lorentz force is explicit: dt < rho/(sigma B^2) = 1/Ha^2
    #  - with the O-grid, the explicit (non-orthogonal) viscous terms in the
    #    skewed corner cells became unstable above about 0.16 Delta^2/nu,
    #    Delta = 2S/n being the central-block cell size (found empirically:
    #    stable at 2.5e-5 and 6.25e-5, unstable at 1e-4 for n = 40)
    dt = min(0.25/Ha**2, 0.15*(2*S/a.cells)**2)
    end = 5.0/Ha
    if a.conducting:
        def wall(faces):
            return f"""    fluid_to_solid              // r = 1, coupled to the conducting wall
    {{
        type            mappedWall;
        sampleMode      nearestPatchFace;
        sampleRegion    solid;
        samplePatch     solid_to_fluid;
        offsetMode      uniform;
        offset          (0 0 0);
        faces           ({faces});
    }}"""
        fluid_dict = case/"system"/"fluid"/"blockMeshDict"
    else:
        def wall(faces):
            return f"""    wall                        // r = 1, insulating
    {{
        type wall;
        faces ({faces});
    }}"""
        fluid_dict = case/"system"/"blockMeshDict"
    text, (n, nr, first, e) = fluid_blockmesh(a, wall)
    fluid_dict.write_text(text)
    entries = [
        ("Ha", Ha, "Hartmann number B a sqrt(sigma/(rho nu)), a = radius; a = rho = nu = sigma = 1"),
        ("B", Ha, "applied field magnitude (= Ha in these units)"),
    ]
    if a.conducting:
        (case/"system"/"solid"/"blockMeshDict").write_text(solid_blockmesh(a))
        entries += [
            ("wallThickness", a.thickness, "t_w / a"),
            ("sigmaRatio", a.sigma_ratio, "sigma_w / sigma"),
            ("sigmaWall", a.sigma_ratio, "sigma_w (sigma = 1)"),
            ("c", a.c, "thin-wall conductance ratio sigma_w t_w / (sigma a)"),
            ("wallCells", float(a.wall_cells), ""),
        ]
    entries += [
        ("deltaT", dt, "min(0.25/Ha^2, 0.15 (2S/n)^2), see pipeCase.write"),
        ("endTime", end, "5/Ha"),
        ("cellsQuarter", float(n), ""),
        ("radialCells", float(nr), ""),
        ("layerCells", float(a.layer_cells), ""),
    ]
    write_parameters(case/"system"/"caseParameters", "setup.py", entries,
                     [f"{title} parameters. Regenerate with: python3 setup.py --ha <Ha>"])
    ncell = n*n + 4*n*nr
    extra = f", t_w/a = {a.thickness:g}, sigma_w/sigma = {a.sigma_ratio:g}, c = {a.c:g}" if a.conducting else ""
    print(f"Ha = {Ha:g}{extra}; {ncell} fluid cells, first wall cell {first:.3g}; "
          f"deltaT = {dt:.3g}, endTime = {end:.3g} ({end/dt:.0f} steps)")


# ----------------------------------------------------------------------------
# Validation
# ----------------------------------------------------------------------------
def validate(case, title, tol, log, conducting):
    from pipeAnalytic import PipeSolution

    par = read_parameters(case/"system"/"caseParameters")
    Ha = par["Ha"]
    if conducting:
        t, s = par["wallThickness"], par["sigmaRatio"]
        ref = PipeSolution(Ha, "thick", sigma_ratio=s, thickness=t)
        thin = PipeSolution(Ha, "thin", c=s*t)
        print(f"\nFor context: G a^2/(rho nu ubar) = {thin.gradp:.6g} for a thin wall of the same "
              f"c = {s*t:g}, {ref.gradp:.6g} for the actual wall ({100*(ref.gradp/thin.gradp - 1):+.2f}%); "
              f"core-flow estimate Ha^2 c/(1+c) = {Ha**2*s*t/(1 + s*t):.6g} (Chang & Lundgren 1961).")
    else:
        ref = PipeSolution(Ha, "insulating")
        print(f"\nFor context: G/Ha = {ref.gradp/Ha:.5f}; the high-Ha limit is 3 pi/8 = {3*np.pi/8:.5f}.")

    tdir = latest_time(case)
    d = tdir/"fluid" if conducting else tdir
    C = read_internal_field(d/"C")
    n = len(C)
    V = read_internal_field(d/"V", n)
    U = read_internal_field(d/"U", n)
    y, z, u = C[:, 1], C[:, 2], U[:, 0]
    ubar = (u*V).sum()/V.sum()
    un = u/ubar
    ua = ref.velocity(y, z)
    e_inf = np.abs(un - ua).max()/ua.max()
    e_l2 = np.sqrt(((un - ua)**2*V).sum()/(ua**2*V).sum())

    hist = pressure_gradient_history(log)
    gp = hist[-1]
    tail = hist[int(0.9*len(hist)):]
    drift = (tail.max() - tail.min())/abs(gp)

    print(f"\n{title}, time {tdir.name}, {n} fluid cells\n")
    ok = report([("pressure gradient G a^2/(rho nu ubar)", gp, ref.gradp),
                 ("u_max / ubar", un.max(), ua.max())], tol)
    for name, val in [("velocity error (max, rel. to u_max)", e_inf),
                      ("velocity error (L2, relative)", e_l2),
                      ("steady-state drift of G, last 10%", drift)]:
        flag = "  <-- above tolerance" if val > tol else ""
        print(f"{name:38s} {val:14.3e}{flag}")
    ok = ok and e_inf <= tol and drift <= 0.1*tol

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        s_ = np.linspace(-1, 1, 1201)
        fig, ax = plt.subplots(2, 2, figsize=(11, 9))
        for k, (along, other, lab) in enumerate([(z, y, "z"), (y, z, "y")]):
            # cells whose centre lies within ~half a circumferential cell of the axis
            nq = par["cellsQuarter"]
            sel = np.abs(other) < 0.6*np.maximum(2*S/nq, 0.5*np.pi*np.hypot(y, z)/nq)
            yy, zz = (np.zeros_like(s_), s_) if lab == "z" else (s_, np.zeros_like(s_))
            ax[0, k].plot(s_, ref.velocity(yy, zz), "k-", label="reference")
            ax[0, k].plot(along[sel], un[sel], "o", ms=3, mfc="none", label="OpenFOAM")
            ax[0, k].set(xlabel=f"{lab} / a", ylabel="u / ubar",
                         title=("normal to B (side / Roberts layers)" if lab == "z"
                                else "along B (Hartmann layers)"))
        ax[0, 0].legend()
        err = (un - ua)/ua.max()
        m = np.abs(err).max()
        sc = ax[1, 0].tricontourf(z, y, err, np.linspace(-m, m, 31), cmap="RdBu_r")
        fig.colorbar(sc, ax=ax[1, 0], label="(u - u_ref)/u_max")
        ax[1, 0].set(xlabel="z / a", ylabel="y / a", title="velocity error", aspect="equal")
        sc = ax[1, 1].tricontourf(z, y, un, 30, cmap="viridis")
        fig.colorbar(sc, ax=ax[1, 1], label="u / ubar")
        ax[1, 1].set(xlabel="z / a", ylabel="y / a", title="OpenFOAM velocity (B along y)", aspect="equal")
        fig.suptitle(title)
        fig.tight_layout()
        fig.savefig(case/"validation.png", dpi=130)
        print("\nWrote validation.png")
    except ImportError:
        pass

    print("\nPASS" if ok else "\nFAIL")
    sys.exit(0 if ok else 1)
