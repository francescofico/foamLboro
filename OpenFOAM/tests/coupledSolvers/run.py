#!/usr/bin/env python3
"""Exercise the actual cEpot executables on analytic two-region slab cases."""
import argparse
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile

SOLVERS = ('cEpotFoam', 'cEpotBuoyantFoam', 'cEpotBuoyantBulkFoam')


def header(kind, name):
    return f'FoamFile {{ version 2.0; format ascii; class {kind}; object {name}; }}\n'


def write(path, text):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)


def run(command, log, error=None):
    with log.open('w') as out:
        result = subprocess.run(command, stdout=out, stderr=out, timeout=120)
    text = log.read_text()
    if (error and (not result.returncode or error not in text)) or (not error and result.returncode):
        raise RuntimeError(f'{command}\n{text[-5000:]}\nLog: {log}')
    return text


def make_case(case, solver, switches='coupled true;', floating=False, invalid=None,
              outer=1, extra='', frozen=False):
    write(case/'system/controlDict', header('dictionary', 'controlDict') +
          f'application {solver}; startFrom startTime; startTime 0; stopAt endTime; '
          'endTime 1; deltaT 1; adjustTimeStep no; writeControl timeStep; writeInterval 1; '
          'writeFormat ascii; writePrecision 16; timePrecision 12; runTimeModifiable false;\n')
    solve = 'solver PBiCGStab; preconditioner diagonal; tolerance 1e-13; relTol 0; maxIter 1000;'
    write(case/'system/fvSolution', header('dictionary', 'fvSolution') + switches +
          f'\nPIMPLE {{ nOuterCorrectors {outer}; {extra} }}\n' +
          f'solvers {{ "(PotE|T).*" {{ {solve} }} }}\n')
    write(case/'constant/regionProperties', header('dictionary', 'regionProperties') +
          'regions (fluid (fluid) solid (solid));\n')
    write(case/'constant/g', header('uniformDimensionedVectorField', 'g') +
          'dimensions [0 1 -2 0 0 0 0]; value (0 0 0);\n')
    decomposition = header('dictionary', 'decomposeParDict') + (
        'numberOfSubdomains 2; method simple; simpleCoeffs { n (1 2 1); delta 0.001; }\n')
    write(case/'system/decomposeParDict', decomposition)
    for region in ('fluid', 'solid'):
        fluid = region == 'fluid'
        neighbour = 'solid' if fluid else 'fluid'
        x0, x1 = (0, 1) if fluid else (1, 3)
        points = [(x0,0,0),(x1,0,0),(x1,1,0),(x0,1,0),
                  (x0,0,1),(x1,0,1),(x1,1,1),(x0,1,1)]
        if invalid == 'order' and not fluid:
            points = [(x,1-y,1-z) for x,y,z in points]
        vertices = ' '.join(f'({x} {y} {z})' for x,y,z in points)
        interface = '(1 2 6 5)' if fluid else '(0 4 7 3)'
        end = '(0 4 7 3)' if fluid else '(1 2 6 5)'
        mode = 'nearestPatchFaceAMI' if invalid == 'AMI' else 'nearestPatchFace'
        write(case/'system'/region/'blockMeshDict', header('dictionary','blockMeshDict') +
              f'scale 1; vertices ({vertices}); '
              'blocks (hex (0 1 2 3 4 5 6 7) (1 4 1) simpleGrading (1 1 1)); edges ();\n'
              f'boundary ( interface {{ type mappedWall; sampleRegion {neighbour}; '
              f'samplePatch interface; sampleMode {mode}; offsetMode uniform; offset (0 0 0); '
              f'faces ({interface}); }} end {{ type wall; faces ({end}); }} '
              'walls { type wall; faces ((0 1 5 4) (3 7 6 2) (0 3 2 1) (4 5 6 7)); } );\n')
        write(case/'system'/region/'decomposeParDict', decomposition)
        write(case/'system'/region/'fvSchemes', header('dictionary','fvSchemes') +
              'ddtSchemes { default Euler; } gradSchemes { default Gauss linear; } '
              'divSchemes { default Gauss linear; } laplacianSchemes { default Gauss linear corrected; } '
              'interpolationSchemes { default linear; } snGradSchemes { default corrected; } '
              'fluxRequired { default no; p; p_rgh; PotE; }\n')
        write(case/'system'/region/'fvSolution', header('dictionary','fvSolution') +
              f'solvers {{ ".*" {{ {solve} }} }} '
              f'PIMPLE {{ nCorrectors 2; nNonOrthogonalCorrectors 1; pRefCell 0; pRefValue 0; '
              f'frozenFlow {str(frozen).lower()}; }} '
              'PotE { PotERefCell 0; PotERefValue 2.5; }\n')
        def field(name, dims, value, interface_bc=None, end_bc=None, vector=False):
            default = f'type fixedValue; value uniform {value};' if vector else 'type zeroGradient;'
            write(case/'0'/region/name, header('volVectorField' if vector else 'volScalarField',name) +
                  f'dimensions [{dims}]; internalField uniform {value};\n'
                  f'boundaryField {{ interface {{ {interface_bc or default} }} '
                  f'end {{ {end_bc or default} }} walls {{ {default} }} }}\n')
        field('sigma','-1 -3 3 0 0 2 0',3 if fluid else 7)
        field('PotE','1 2 -3 0 0 -1 0',-40 if fluid else 70,
              'type solidWallMixedPotentialCoupled; value uniform 0;',
              'type zeroGradient;' if floating else
              f'type fixedValue; value uniform {0 if fluid else 10};')
        if solver != 'cEpotFoam':
            field('K','1 1 -3 -1 0 0 0',2 if fluid else 5)
            field('alphat','0 2 -1 0 0 0 0',0.1 if fluid else 0.3)
            field('T','0 0 0 1 0 0 0',300 if fluid else 600,
                  'type solidWallMixedTemperatureCoupled; value uniform 0;',
                  f'type fixedValue; value uniform {300 if fluid else 600};')
        if fluid:
            field('U','0 1 -1 0 0 0 0','(0 0 0)',vector=True)
            field('B0','1 0 -2 0 0 -1 0','(0 0 0)',vector=True)
            field('p' if solver=='cEpotFoam' else 'p_rgh','0 2 -2 0 0 0 0',0)
            if solver == 'cEpotBuoyantBulkFoam':
                field('Tb','0 0 0 1 0 0 0',300)
            write(case/'constant/fluid/turbulenceProperties',
                  header('dictionary','turbulenceProperties')+'simulationType laminar;\n')
            write(case/'constant/fluid/transportProperties', header('dictionary','transportProperties')+
                  'transportModel Newtonian; nu nu [0 2 -1 0 0 0 0] 0.1; '
                  'rho rho [1 -3 0 0 0 0 0] 1; beta beta [0 0 0 -1 0 0 0] 0.0001; '
                  'TRef TRef [0 0 0 1 0 0 0] 300; Pr Pr [0 0 0 0 0 0 0] 1; '
                  'Prt Prt [0 0 0 0 0 0 0] 1;\n')


def values(path, key='internalField', vector=False):
    text = path.read_text()
    if key != 'internalField':
        text = re.search(r'\b'+key+r'\s*\{([^}]+)\}', text, re.S)[1]
        key = 'value'
    value = re.search(r'\b'+key+r'\s+uniform\s+([^;]+);', text)
    if value:
        data = value[1]
    else:
        value = re.search(r'\b'+key+r'\s+nonuniform\s+List<\w+>\s+\d+\s*\((.*?)\)\s*;', text, re.S)
        if not value:
            raise AssertionError(f'Missing {key} in {path}')
        data = value[1]
    numbers = [float(x) for x in re.findall(r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?', data)]
    return list(zip(*(iter(numbers),)*3)) if vector else numbers


def close(actual, expected, label, tolerance=1e-8):
    error = max(abs(x-expected) for x in actual)
    if error > tolerance:
        raise AssertionError(f'{label}: expected {expected}, got {actual}, error {error}')


def thermal_step(previous=(300,600)):
    # One cell across each slab. Distances 1/2, K=2/5, alpha=0.1/0.3.
    # These are the exact backward-Euler finite-volume equations (dt=1).
    gf, gs = 1/9, 2/15
    a, d = 1+0.2+gf, 2+0.3+gs
    b, c = previous[0]+0.2*300, 2*previous[1]+0.3*600
    return ((b*d+gf*c)/(a*d-gf*gs), (a*c+gs*b)/(a*d-gf*gs))


def check(case, solver, time='1', potential=True, thermal=True, floating=False, expected_t=None):
    expected_t = expected_t or thermal_step()
    for i, region in enumerate(('fluid','solid')):
        root = case/time/region
        if potential:
            close(values(root/'PotE'), 2.5 if floating else (35/13,100/13)[i], f'{region} potential')
            j = 0 if floating else -210/13
            for actual in values(root/'J',vector=True):
                for axis, expected in zip(actual,(j,0,0)):
                    close([axis], expected, f'{region} current')
        if solver != 'cEpotFoam' and thermal:
            close(values(root/'T'),expected_t[i],f'{region} temperature')
            interface_t = (4*expected_t[0]+5*expected_t[1])/9
            close(values(root/'T','interface'),interface_t,f'{region} interface temperature')
            if i == 0:
                close(values(root/'rhok'),1-0.0001*(expected_t[0]-300),'buoyancy')
    if potential:
        close(values(case/time/'fluid/divJ'),0,'divergence of current')


def execute(case, solver, bin_dir, parallel=False, error=None):
    for region in ('fluid','solid'):
        run(['blockMesh','-case',str(case),'-region',region],case/f'mesh-{region}.log')
    command = [str(bin_dir/solver),'-case',str(case),'-fileHandler','uncollated']
    if parallel:
        run(['decomposePar','-case',str(case),'-allRegions'],case/'decompose.log')
        command = ['mpirun','-np','2']+command+['-parallel']
    text = run(command,case/'solve.log',error)
    if parallel and not error:
        run(['reconstructPar','-case',str(case),'-allRegions'],case/'reconstruct.log')
    return text


def restart(case, solver, bin_dir, switches):
    control = case/'system/controlDict'
    control.write_text(control.read_text().replace('startTime 0;', 'startTime 1;').replace('endTime 1;', 'endTime 2;'))
    solution = case/'system/fvSolution'
    solution.write_text(solution.read_text().replace('coupled true;', switches))
    return run([str(bin_dir/solver), '-case', str(case), '-fileHandler', 'uncollated'], case/'restart.log')


def different_topologies(case, bin_dir):
    # Add a third solid that is thermally connected but electrically isolated.
    # PotE and T now have different graphs with the same first implicit patch.
    solver = 'cEpotBuoyantFoam'
    make_case(case,solver,outer=2)
    for directory in ('0', 'system'):
        shutil.copytree(case/directory/'solid',case/directory/'solid2',dirs_exist_ok=True)
    mesh = case/'system/solid2/blockMeshDict'
    data = mesh.read_text()
    vertices = re.search(r'vertices \((.*?)\);',data)[1]
    vertices = re.sub(r'\((\d+) ',lambda m: f'({int(m[1])+2} ',vertices)
    data = re.sub(r'vertices \(.*?\);',f'vertices ({vertices});',data)
    mesh.write_text(data.replace('sampleRegion fluid;', 'sampleRegion solid;').replace('samplePatch interface;', 'samplePatch end;'))
    mesh = case/'system/solid/blockMeshDict'
    mesh.write_text(mesh.read_text().replace('end { type wall;',
        'end { type mappedWall; sampleRegion solid2; samplePatch interface; '
        'sampleMode nearestPatchFace; offsetMode uniform; offset (0 0 0);'))
    for region, patch in (('solid','end'),('solid2','interface'),('solid2','end')):
        path = case/'0'/region/'PotE'
        path.write_text(re.sub(patch+r' \{[^}]+\}',patch+' { type zeroGradient; }',path.read_text()))
    path = case/'0/solid/T'
    path.write_text(path.read_text().replace('end { type fixedValue; value uniform 600; }',
        'end { type solidWallMixedTemperatureCoupled; value uniform 600; }'))
    path = case/'constant/regionProperties'
    path.write_text(path.read_text().replace('solid (solid)', 'solid (solid solid2)'))
    for region in ('fluid','solid','solid2'):
        run(['blockMesh','-case',str(case),'-region',region],case/f'mesh-{region}.log')
    run([str(bin_dir/solver),'-case',str(case)],case/'solve.log')
    # Solve the three analytic backward-Euler slab balances by elimination.
    a = [[1+0.2+1/9,-1/9,0,360],
         [-2/15,2+2/15+0.15,-0.15,1200],
         [0,-0.15,2+0.15+0.3,1380]]
    for i in range(3):
        divisor = a[i][i]
        a[i] = [x/divisor for x in a[i]]
        for j in range(3):
            if i != j:
                factor = a[j][i]
                a[j] = [x-factor*y for x,y in zip(a[j],a[i])]
    for i, region in enumerate(('fluid','solid','solid2')):
        close(values(case/'1'/region/'T'),a[i][3],f'{region} three-region temperature')
        close(values(case/'1'/region/'PotE'),2.5 if i==2 else 0,f'{region} disconnected potential')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bin-dir',type=Path,default=Path(os.environ.get('FOAM_USER_APPBIN','.')))
    parser.add_argument('--work-dir',type=Path)
    parser.add_argument('--serial',action='store_true')
    args = parser.parse_args()
    work = args.work_dir or Path(tempfile.mkdtemp(prefix='coupled-solvers-'))
    work.mkdir(parents=True,exist_ok=True)
    print(f'Test results: {work}',flush=True)
    for solver in SOLVERS:
        for parallel in ((False,) if args.serial else (False,True)):
            case = work/(solver+('-mpi' if parallel else '-serial'))
            make_case(case,solver)
            log = execute(case,solver,args.bin_dir,parallel)
            assert 'Solving coupled PotE across 2 regions' in log
            if solver != 'cEpotFoam':
                assert 'Solving coupled T across 2 regions' in log
            check(case,solver)
            print(f'PASS {case.name}: single assembled solve, analytic fields and fluxes',flush=True)
    for floating, switches, label in ((True,'coupled true;','floating'),
            (False,'coupledPotentialField true; coupledEnergyField false;','potential-only'),
            (False,'coupledPotentialField false; coupledEnergyField true;','temperature-only')):
        solver = 'cEpotBuoyantFoam'
        case = work/label
        make_case(case,solver,switches,floating)
        log = execute(case,solver,args.bin_dir)
        check(case,solver,potential=label!='temperature-only',thermal=label!='potential-only',floating=floating)
        assert ('Solving coupled PotE' in log) == (label!='temperature-only')
        assert ('Solving coupled T' in log) == (label!='potential-only')
        print(f'PASS {label}',flush=True)
    if not args.serial:
        case = work/'floating-mpi'
        make_case(case,'cEpotBuoyantFoam',floating=True)
        execute(case,'cEpotBuoyantFoam',args.bin_dir,parallel=True)
        check(case,'cEpotBuoyantFoam',floating=True)
        print('PASS floating MPI reference',flush=True)
        case = work/'reject-processor-pairing'
        make_case(case,'cEpotFoam')
        for region in ('fluid','solid'):
            write(case/'system'/region/'decomposeParDict',header('dictionary','decomposeParDict')+
                  'numberOfSubdomains 2; method manual; manualCoeffs { dataFile cellDecomposition; }')
            ranks = '0 0 1 1' if region=='fluid' else '1 1 0 0'
            write(case/'constant'/region/'cellDecomposition',header('labelList','cellDecomposition')+f'4({ranks})')
        execute(case,'cEpotFoam',args.bin_dir,parallel=True,error='requires matching faces')
        print('PASS reject mismatched processor pairing',flush=True)
    different_topologies(work/'different-topologies',args.bin_dir)
    print('PASS different thermal/electrical graphs and isolated reference',flush=True)
    # Repeated PIMPLE/energy and non-orthogonal passes must not add a region twice.
    for frozen in (False, True):
        case = work/('correction-loops-frozen' if frozen else 'correction-loops')
        make_case(case,'cEpotBuoyantFoam',outer=2, frozen=frozen,
                  extra='nNonOrthogonalCorrectors 1; energyCoupling { iterations 2; }')
        execute(case,'cEpotBuoyantFoam',args.bin_dir)
        check(case,'cEpotBuoyantFoam')
        print(f'PASS {case.name}',flush=True)
    for switches, label in (('', 'restart-auto'), ('coupled false;', 'restart-disable')):
        case = work/label
        make_case(case,'cEpotBuoyantFoam',outer=60 if switches else 1)
        execute(case,'cEpotBuoyantFoam',args.bin_dir)
        check(case,'cEpotBuoyantFoam')
        log = restart(case,'cEpotBuoyantFoam',args.bin_dir,switches)
        check(case,'cEpotBuoyantFoam',time='2',expected_t=thermal_step(thermal_step()))
        assert ('Solving coupled PotE' in log) == (not switches)
        assert ('Solving coupled T' in log) == (not switches)
        print(f'PASS {label}',flush=True)
    case = work/'segregated-AMI'
    make_case(case,'cEpotBuoyantFoam','coupled false;',outer=60,invalid='AMI')
    log = execute(case,'cEpotBuoyantFoam',args.bin_dir)
    check(case,'cEpotBuoyantFoam')
    assert 'Solving coupled' not in log
    print('PASS segregated AMI solution',flush=True)
    for invalid, error in (('order','requires matching faces'),('AMI','Implicit region coupling requires nearestPatchFace')):
        case = work/('reject-'+invalid)
        make_case(case,'cEpotFoam',invalid=invalid)
        execute(case,'cEpotFoam',args.bin_dir,error=error)
        print(f'PASS reject {invalid}',flush=True)


if __name__ == '__main__':
    main()
