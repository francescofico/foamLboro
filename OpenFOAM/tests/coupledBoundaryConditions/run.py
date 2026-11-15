#!/usr/bin/env python3
"""Build and dynamically load the shared boundary library in scratch cases."""

import argparse
import os
from pathlib import Path
import re
import resource
import shutil
import subprocess
import tempfile


HERE = Path(__file__).resolve().parent
OPENFOAM = HERE.parents[1]


def header(kind, name):
    return f"FoamFile {{ version 2.0; format ascii; class {kind}; object {name}; }}\n"


def run(command, log, env, expect_error=None):
    with log.open("w") as output:
        result = subprocess.run(command, stdout=output, stderr=output, env=env, timeout=120)
    text = log.read_text()
    if expect_error:
        if result.returncode == 0 or expect_error not in text:
            raise RuntimeError(f"Expected rejection missing in {log}:\n{text[-6000:]}")
    elif result.returncode:
        raise RuntimeError(f"Command failed: {command}\n{text[-6000:]}\nLog: {log}")
    return text


def build(work, env):
    source = OPENFOAM / "src" / "mhdBoundaryConditions"
    library = work / "library"
    application = work / "probe"
    shutil.copytree(source, library, ignore=shutil.ignore_patterns("lnInclude", "linux*"))
    for directory in (application / "Make", work / "lib", work / "bin"):
        directory.mkdir(parents=True)
    env["FOAM_USER_LIBBIN"] = str(work / "lib")
    env["FOAM_USER_APPBIN"] = str(work / "bin")
    lib = work / "lib" / "libmhdBoundaryConditions.so"
    includes = "-I$(LIB_SRC)/finiteVolume/lnInclude -I$(LIB_SRC)/meshTools/lnInclude"
    shutil.copy2(HERE / "probe.C", application / "probe.C")
    (application / "Make" / "files").write_text(f"probe.C\nEXE = {work / 'bin' / 'probe'}\n")
    # No direct library link: controlDict must load and register all three BCs.
    (application / "Make" / "options").write_text(
        f"EXE_INC = {includes}\n"
        "EXE_LIBS = -lfiniteVolume -lmeshTools\n"
    )
    run(["wmake", "libso", str(library)], work / "build-library.log", env)
    run(["wmake", str(application)], work / "build-probe.log", env)
    return str(lib), str(work / "bin" / "probe")


def make_case(case, mode, lib, renamed=False, invalid=None, empty_patches=False):
    nonmatching = mode == "nearestPatchFaceAMI"
    t_slope, v_slope = (0, 0) if nonmatching else (20, 2)
    names = {
        "fluid": ("T", "PotE", "K", "sigma"),
        "solid": ("Ts" if renamed else "T", "Vs" if renamed else "PotE", "Ks", "sigmaS"),
    }
    (case / "system").mkdir(parents=True)
    (case / "constant").mkdir()
    (case / "system" / "controlDict").write_text(
        header("dictionary", "controlDict")
        + "application probe; startFrom startTime; startTime 0; stopAt endTime; "
        + "endTime 1; deltaT 1; writeControl timeStep; writeInterval 1; writePrecision 16;\n"
        + f'libs ("{lib}");\n'
    )
    (case / "constant" / "regionProperties").write_text(
        header("dictionary", "regionProperties") + "regions (fluid (fluid) solid (solid));\n"
    )
    # Fluid/solid cell-centre distances are 0.5/1.0, conductivities 2/5
    # and electrical conductivities 3/7: analytic series-resistance solution.
    settings = (
        header("dictionary", "couplingTestDict")
        + f"expectedTemperature {4200 / 9:.16g}; expectedTemperatureSlope {t_slope};\n"
        + f"expectedPotential {270 / 13:.16g}; expectedPotentialSlope {v_slope};\n"
        + "prescribedHeatFlux 12;\n"
    )
    for region, fields in names.items():
        settings += (
            f"{region} {{ temperature {fields[0]}; potential {fields[1]}; "
            f"conductivity {fields[2]}; sigma {fields[3]}; heatFluxTemperature Tflux; }}\n"
        )
    (case / "system" / "couplingTestDict").write_text(settings)

    for region in ("fluid", "solid"):
        fluid = region == "fluid"
        neighbour = "solid" if fluid else "fluid"
        ny = 4 if fluid or not nonmatching else 6
        nx = 2 if empty_patches else 1
        x0, x1 = (0, 1) if fluid else (1, 3)
        system, constant, zero = (case / d / region for d in ("system", "constant", "0"))
        for directory in (system, constant, zero):
            directory.mkdir(parents=True)
        # The solid block reverses its tangential axes so local face ordering
        # differs even for matching interfaces. Mapping must reorder values.
        points = (
            [(x0, 0, 0), (x1, 0, 0), (x1, 1, 0), (x0, 1, 0),
             (x0, 0, 1), (x1, 0, 1), (x1, 1, 1), (x0, 1, 1)]
            if fluid else
            [(x0, 1, 1), (x1, 1, 1), (x1, 0, 1), (x0, 0, 1),
             (x0, 1, 0), (x1, 1, 0), (x1, 0, 0), (x0, 0, 0)]
        )
        vertices = " ".join(f"({x} {y} {z})" for x, y, z in points)
        interface = "(1 2 6 5)" if fluid else "(0 4 7 3)"
        outer = "(0 4 7 3)" if fluid else "(1 2 6 5)"
        patch_type = (
            "type wall;" if invalid == "wall" else
            f"type mappedWall; sampleRegion {neighbour}; samplePatch interface; "
            f"sampleMode {'nearestCell' if invalid == 'mode' else mode}; "
            "offsetMode uniform; offset (0 0 0);"
        )
        partner = ""
        wall_outer = outer
        if invalid == "cyclic":
            offset = x0 - x1 if fluid else x1 - x0
            patch_type = (
                "type cyclic; neighbourPatch periodicPartner; transform translational; "
                f"separationVector ({offset} 0 0);"
            )
            partner = (
                "periodicPartner { type cyclic; neighbourPatch interface; "
                f"transform translational; separationVector ({-offset} 0 0); faces ({outer}); }}\n"
            )
            wall_outer = ""
        (system / "blockMeshDict").write_text(
            header("dictionary", "blockMeshDict")
            + f"scale 1; vertices ({vertices});\n"
            + f"blocks (hex (0 1 2 3 4 5 6 7) ({nx} {ny} 1) simpleGrading (1 1 1));\n"
            + "edges (); boundary (\n"
            + f"interface {{ {patch_type} faces ({interface}); }}\n"
            + partner
            + f"walls {{ type wall; faces ({wall_outer} (0 1 5 4) (3 7 6 2) "
            + "(0 3 2 1) (4 5 6 7)); }\n);\n"
        )
        (system / "fvSchemes").write_text(
            header("dictionary", "fvSchemes")
            + "ddtSchemes { default steadyState; } gradSchemes { default Gauss linear; } "
            + "divSchemes { default none; } laplacianSchemes { default Gauss linear corrected; } "
            + "interpolationSchemes { default linear; } snGradSchemes { default corrected; }\n"
        )
        (system / "fvSolution").write_text(header("dictionary", "fvSolution") + "solvers {}\n")
        decomposition = (
            header("dictionary", "decomposeParDict")
            + "numberOfSubdomains 2; method manual; "
            + "manualCoeffs { dataFile cellDecomposition; }\n"
        )
        (system / "decomposeParDict").write_text(decomposition)
        (case / "system" / "decomposeParDict").write_text(decomposition)
        # Equal index partitions, reversed solid geometry: corresponding faces
        # live on different ranks, with unequal counts for the AMI case.
        ranks = ([1, 0] * ny if empty_patches else
                 [0 if j < ny // 2 else 1 for j in range(ny)])
        (constant / "cellDecomposition").write_text(
            header("labelList", "cellDecomposition")
            + f"{nx * ny}\n(" + " ".join(map(str, ranks)) + ")\n"
        )
        fields = names[region]
        for i, field in enumerate(fields):
            body = header("volScalarField", field)
            if i < 2:
                base = (300 if fluid else 600) if i == 0 else (10 if fluid else 30)
                slope = t_slope if i == 0 else v_slope
                y = [(j + 0.5) / ny for j in range(ny)]
                if not fluid:
                    y = [1 - value for value in y]
                values = " ".join(f"{base + slope * value:.16g}" for value in y for _ in range(nx))
                dimensions = "0 0 0 1 0 0 0" if i == 0 else "1 2 -3 0 0 -1 0"
                kind = "Temperature" if i == 0 else "Potential"
                key, material = ("Tnbr", "K") if i == 0 else ("PotEnbr", "sigma")
                body += f"dimensions [{dimensions}]; internalField nonuniform List<scalar> {nx * ny} ({values});\n"
                body += (
                    "boundaryField { interface { "
                    f"type solidWallMixed{kind}Coupled; {key} {names[neighbour][i]}; "
                    f"{material} {fields[i + 2]}; value uniform {base}; "
                    + ("patchType cyclic; " if invalid == "cyclic" else "")
                    + "} walls { type zeroGradient; } "
                    + ("periodicPartner { type cyclic; } " if invalid == "cyclic" else "")
                    + "}\n"
                )
            else:
                value = (2 if fluid else 5) if i == 2 else (3 if fluid else 7)
                dimensions = "1 1 -3 -1 0 0 0" if i == 2 else "-1 -3 3 0 0 2 0"
                body += f"dimensions [{dimensions}]; internalField uniform {value};\n"
                body += (
                    "boundaryField { interface { "
                    + ("type cyclic; " if invalid == "cyclic" else f"type fixedValue; value uniform {value}; ")
                    + "} walls { type zeroGradient; } "
                    + ("periodicPartner { type cyclic; } " if invalid == "cyclic" else "")
                    + "}\n"
                )
            (zero / field).write_text(body)

        (zero / "Tflux").write_text(
            header("volScalarField", "Tflux")
            + "dimensions [0 0 0 1 0 0 0]; internalField uniform 300;\n"
            + "boundaryField { interface { "
            + ("type cyclic; " if invalid == "cyclic" else "type zeroGradient; ")
            + "} walls { type solidWallHeatFluxTemperature; "
            + f"K {fields[2]}; q uniform 12; value uniform 300; gradient uniform 0; }} "
            + ("periodicPartner { type cyclic; } " if invalid == "cyclic" else "")
            + "}\n"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--serial", action="store_true", help="Skip the two-rank MPI checks")
    parser.add_argument("--work-dir", type=Path, help="New directory for cases, binaries and logs")
    args = parser.parse_args()
    if not os.environ.get("WM_PROJECT_DIR") or not shutil.which("wmake"):
        parser.error("Source the OpenFOAM.com etc/bashrc first")
    resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
    work = args.work_dir.resolve() if args.work_dir else Path(tempfile.mkdtemp(prefix="foam-coupling-test-"))
    if args.work_dir:
        work.mkdir(parents=True, exist_ok=False)
    print(f"OpenFOAM {os.environ.get('WM_PROJECT_VERSION')}; logs and cases: {work}", flush=True)
    env = os.environ.copy()
    env.pop("FOAM_ABORT", None)
    env["LD_LIBRARY_PATH"] = str(work / "lib") + ":" + env.get("LD_LIBRARY_PATH", "")
    lib, probe = build(work, env)

    for name, mode, renamed, empty_patches in (
        ("matching", "nearestPatchFace", False, False),
        ("nonmatching", "nearestPatchFaceAMI", False, False),
        ("renamed", "nearestPatchFace", True, False),
        ("empty-patches", "nearestPatchFaceAMI", False, True),
    ):
        case = work / name
        make_case(case, mode, lib, renamed, empty_patches=empty_patches)
        for region in ("fluid", "solid"):
            run(["blockMesh", "-case", str(case), "-region", region], case / f"mesh-{region}.log", env)
        command = [probe, "-case", str(case)]
        for stage in ("serial", "restart"):
            text = run(command, case / f"{stage}.log", env)
            if "PASS:" not in text:
                raise RuntimeError(f"Missing test result in {case / (stage + '.log')}")
        for field in (case / "0").glob("*/*"):
            for key in ("Tnbr", "PotEnbr"):
                if len(re.findall(rf"\b{key}\s", field.read_text())) > 1:
                    raise RuntimeError(f"Duplicate {key} entry in restart file {field}")
        if not args.serial:
            run(["decomposePar", "-case", str(case), "-allRegions"], case / "decompose.log", env)
            for stage in ("parallel", "parallel-restart"):
                run(["mpirun", "-np", "2", probe, "-case", str(case), "-parallel"],
                    case / f"{stage}.log", env)
        print(f"PASS {name}: serial/restart" + (" + MPI/restart" if not args.serial else ""), flush=True)

    for invalid, message in (("wall", "does not provide neighbour mapping"),
                             ("cyclic", "does not provide neighbour mapping"),
                             ("mode", "requires sampleMode nearestPatchFace")):
        case = work / ("invalid-" + invalid)
        make_case(case, "nearestPatchFace", lib, invalid=invalid)
        for region in ("fluid", "solid"):
            run(["blockMesh", "-case", str(case), "-region", region], case / f"mesh-{region}.log", env)
        run([probe, "-case", str(case)], case / "rejection.log", env, expect_error=message)
        print(f"PASS invalid {invalid}: rejected with configuration guidance", flush=True)


if __name__ == "__main__":
    main()
