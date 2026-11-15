#!/usr/bin/env python3
"""Check residual-controlled outer loops using the actual cEpot executables."""
import argparse
import os
from pathlib import Path
import re
import shutil
import tempfile

from run import SOLVERS, check, close, execute, make_case, values


def controls(case, region, fields):
    entries = ' '.join(f'{field} {{ tolerance {absolute}; relTol {relative}; }}'
                       for field, (absolute, relative) in fields.items())
    path = case/'system'/region/'fvSolution'
    path.write_text(path.read_text().replace(
        'PIMPLE {', f'PIMPLE {{ residualControl {{ {entries} }}', 1))


def final_solvers(case):
    # A distinct solver name makes Final selection observable even at zero residual.
    final = ('solver smoothSolver; smoother symGaussSeidel; '
             'tolerance 1e-13; relTol 0; maxIter 1000;')
    entries = ' '.join(f'{name}Final {{ {final} }}' for name in ('U', 'p', 'p_rgh', 'PotE', 'T'))
    for region in ('', 'fluid', 'solid'):
        path = case/'system'/region/'fvSolution'
        path.write_text(path.read_text().replace('solvers {', f'solvers {{ {entries}', 1))


def two_steps(case):
    path = case/'system/controlDict'
    path.write_text(path.read_text().replace('endTime 1;', 'endTime 2;'))


def resolved_slabs(case):
    # Uniform potential inside a one-cell-wide slab has a vanishing residual
    # normalization near convergence. Resolve the physical gradient instead.
    for region in ('fluid', 'solid'):
        path = case/'system'/region/'blockMeshDict'
        path.write_text(path.read_text().replace('(1 4 1)', '(2 4 1)'))


def slab_step(previous):
    # Four independent backward-Euler balances, dx=0.5/1, alpha=0.1/0.3.
    # Interface conductances use K=2/5 and half-cell distances 0.25/0.5.
    a = [[0.5+0.4+0.2, -0.2, 0, 0, 0.5*previous[0]+0.4*300],
         [-0.2, 0.5+0.2+2/9, -2/9, 0, 0.5*previous[1]],
         [0, -4/15, 1+4/15+0.3, -0.3, previous[2]],
         [0, 0, -0.3, 1+0.3+0.6, previous[3]+0.6*600]]
    for i in range(4):
        divisor = a[i][i]
        a[i] = [value/divisor for value in a[i]]
        for j in range(4):
            if i != j:
                factor = a[j][i]
                a[j] = [x-factor*y for x, y in zip(a[j], a[i])]
    return [row[4] for row in a]


def check_slabs(case, solver, time, temperature):
    for i, region in enumerate(('fluid', 'solid')):
        root = case/time/region
        potential = ((17.5/13, 52.5/13), (85/13, 115/13))[i]
        for cell, value in enumerate(values(root/'PotE')):
            close([value], potential[cell % 2], f'{region} potential cell {cell}')
        for current in values(root/'J', vector=True):
            for actual, expected in zip(current, (-210/13, 0, 0)):
                close([actual], expected, f'{region} current')
        if solver != 'cEpotFoam':
            for cell, value in enumerate(values(root/'T')):
                close([value], temperature[2*i+cell % 2], f'{region} temperature cell {cell}')
            close(values(root/'T', 'interface'), (8*temperature[1]+10*temperature[2])/18,
                  f'{region} interface temperature')
            if i == 0:
                for cell, value in enumerate(values(root/'rhok')):
                    close([value], 1-0.0001*(temperature[cell % 2]-300), 'buoyancy')
    close(values(case/time/'fluid/divJ'), 0, 'divergence of current')


def steps(log):
    parts = re.split(r'^Time = (\S+)\s*$', log, flags=re.M)
    return list(zip(parts[1::2], parts[2::2]))


def corrections(log):
    parts = re.split(r'^Pimple Outer Correction = (\d+)\s*$', log, flags=re.M)
    indices = [int(x) for x in parts[1::2]]
    assert indices == list(range(len(indices))), indices
    return parts[2::2]


def residuals(log, field):
    return [float(value) for value in re.findall(
        r'Solving for '+re.escape(field)+r', Initial residual = ([^,]+)', log)]


def assert_early(log, maximum):
    outer = corrections(log)
    assert 2 <= len(outer) < maximum, len(outer)
    assert log.count('performing final correction') == 1, log[-3000:]
    assert log.count('completed final correction') == 1, log[-3000:]
    assert 'performing final correction' in outer[-2]
    assert 'completed final correction' in outer[-1]
    return outer


def basic(case, solver, bin_dir, coupled, parallel):
    maximum = 50
    make_case(case, solver, f'coupled {str(coupled).lower()};', outer=maximum)
    resolved_slabs(case)
    pressure = 'p' if solver == 'cEpotFoam' else 'p_rgh'
    scalar_fields = ['PotE'] + ([] if solver == 'cEpotFoam' else ['T'])
    controls(case, 'fluid', {name: (1e-10, 0) for name in ['U', pressure]+scalar_fields})
    controls(case, 'solid', {name: (1e-10, 0) for name in scalar_fields})
    final_solvers(case)
    two_steps(case)
    log = execute(case, solver, bin_dir, parallel)
    assert len(steps(log)) == 2
    counts = []
    temperature = [300, 300, 600, 600]
    for time, step in steps(log):
        outer = assert_early(step, maximum)
        counts.append(len(outer))
        # Ordinary solvers must be restored at the next time step.
        assert 'smoothSolver: ' not in outer[0]
        for field in ['Ux', 'Uy', 'Uz', pressure]+scalar_fields:
            assert re.search(r'smoothSolver:\s+Solving for '+field+r',', outer[-1]), field
        if not coupled:
            solid = outer[-1].split('Solving for solid region solid')[1]
            for field in scalar_fields:
                assert re.search(r'smoothSolver:\s+Solving for '+field+r',', solid), field
        temperature = slab_step(temperature)
        check_slabs(case, solver, time, temperature)
    print(f'PASS {case.name}: outer counts {counts}, Final selection, analytic fields', flush=True)


def edge_cases(work, bin_dir, parallel):
    suffix = '-mpi' if parallel else '-serial'

    def setup(name, solver='cEpotBuoyantFoam', **kwargs):
        case = work/(name+suffix)
        make_case(case, solver, **kwargs)
        return case

    def solve(case, solver='cEpotBuoyantFoam', **kwargs):
        return execute(case, solver, bin_dir, parallel, **kwargs)

    # Every configured region/field must pass, including the solid and T.
    for blocked in ('PotE', 'T'):
        case = setup('blocked-solid-'+blocked, outer=4)
        controls(case, 'fluid', {'U': (1, 0), 'PotE': (1, 0), 'T': (1, 0)})
        controls(case, 'solid', {'PotE': (1, 0), 'T': (1, 0), blocked: (0, 0)})
        log = solve(case)
        assert len(corrections(log)) == 4
        assert 'performing final correction' not in log
        assert 'residual criteria not satisfied after 4' in log
        print(f'PASS {case.name}: all fields and regions gate early exit', flush=True)

    # No criteria: the fixed iteration count still applies, even with a global block.
    case = setup('no-regional-controls', outer=4,
                 extra='residualControl { U { tolerance 1; relTol 0; } }')
    log = solve(case)
    assert len(corrections(log)) == 4
    assert 'PIMPLE: residual criteria' not in log
    print(f'PASS {case.name}', flush=True)

    case = setup('single-correction', outer=1)
    controls(case, 'fluid', {'U': (1, 0)})
    log = solve(case)
    assert len(corrections(log)) == 1
    assert 'iteration limit' in log
    assert 'performing final correction' not in log
    print(f'PASS {case.name}', flush=True)

    # A region without controls is allowed, but a missing configured field is not.
    for name, region, field, frozen, predictor in (
            ('unknown-field', 'fluid', 'PotETypo', False, True),
            ('solid-flow-field', 'solid', 'U', False, True),
            ('no-momentum-predictor', 'fluid', 'U', False, False),
            ('frozen-pressure', 'fluid', 'p_rgh', True, True)):
        case = setup(name, outer=3, frozen=frozen)
        controls(case, region, {field: (1, 0)})
        if not predictor:
            path = case/'system/fluid/fvSolution'
            path.write_text(path.read_text().replace('PIMPLE {', 'PIMPLE { momentumPredictor false;', 1))
        two_steps(case)
        log = solve(case)
        assert f'No solver residual for {field} in region {region}' in log
        assert len(steps(log)) == 2
        for _, step in steps(log):
            assert len(corrections(step)) == 3
            assert 'performing final correction' not in step
        print(f'PASS {case.name}: unavailable residual cannot pass', flush=True)

    # Independent oracle from reported fluid residuals, with a different baseline
    # on the second time step. Absolute tolerance is disabled throughout.
    case = setup('relative-only', switches='coupled false;', outer=30)
    resolved_slabs(case)
    ratio = 0.2
    controls(case, 'fluid', {'PotE': (0, ratio)})
    two_steps(case)
    log = solve(case)
    baselines = []
    for _, step in steps(log):
        outer = assert_early(step, 30)
        residual = [residuals(part, 'PotE')[0] for part in outer]
        baselines.append(residual[0])
        assert residual[0] > 0
        crossing = next(i for i in range(1, len(residual)) if residual[i]/residual[0] < ratio)
        assert len(outer) == crossing+2, (residual, crossing)
        assert 'performing final correction' not in outer[0]
    assert len(baselines) == 2 and abs(baselines[0]-baselines[1]) > 1e-3
    print(f'PASS {case.name}: relative stopping matches logged residual ratios {baselines}', flush=True)

    case = setup('zero-relative-baseline', outer=6)
    controls(case, 'fluid', {'U': (0, 0.1)})
    log = solve(case)
    assert len(assert_early(log, 6)) == 3
    assert all(value == 0 for value in residuals(log, 'Uy'))
    print(f'PASS {case.name}', flush=True)

    # Ux is zero while Uy has a nonzero initial residual; checking only Ux fails.
    case = setup('vector-components', outer=4)
    path = case/'0/fluid/U'
    path.write_text(path.read_text().replace('internalField uniform (0 0 0);',
                                            'internalField uniform (0 1 0);'))
    controls(case, 'fluid', {'U': (0.5, 0)})
    log = solve(case)
    first = corrections(log)[0]
    assert residuals(first, 'Ux')[0] < 0.5 < residuals(first, 'Uy')[0]
    assert 'performing final correction' not in first
    print(f'PASS {case.name}: non-x component blocks convergence', flush=True)

    # Additional energy and nonorthogonal solves finish before checking controls.
    case = setup('energy-and-nonorthogonal', outer=8, frozen=True,
                 extra='nNonOrthogonalCorrectors 1; energyCoupling { iterations 2; }')
    controls(case, 'solid', {'"(PotE|T)"': (1e-10, 0)})
    log = solve(case)
    outer = assert_early(log, 8)
    assert len(outer) == 2
    assert outer[0].count('energyCoupling:') == 2
    assert outer[0].index('energyCoupling: 2/2') < outer[0].index('performing final correction')
    check(case, 'cEpotBuoyantFoam')
    print(f'PASS {case.name}: scalar-only frozen flow and field patterns', flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bin-dir', type=Path, default=Path(os.environ.get('FOAM_USER_APPBIN', '.')))
    parser.add_argument('--work-dir', type=Path)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--serial', action='store_true')
    mode.add_argument('--parallel-only', action='store_true')
    args = parser.parse_args()
    required = ['blockMesh'] + ([] if args.serial else ['mpirun', 'decomposePar', 'reconstructPar'])
    missing = [name for name in required if not shutil.which(name)]
    missing += [str(args.bin_dir/name) for name in SOLVERS if not (args.bin_dir/name).is_file()]
    if missing:
        parser.error('Missing executables: '+', '.join(missing))
    work = args.work_dir or Path(tempfile.mkdtemp(prefix='residual-control-'))
    work.mkdir(parents=True, exist_ok=True)
    print(f'Test results: {work}', flush=True)
    modes = (False,) if args.serial else ((True,) if args.parallel_only else (False, True))
    for parallel in modes:
        for solver in SOLVERS:
            for coupled in (False, True):
                name = solver+('-assembled' if coupled else '-separate')+('-mpi' if parallel else '-serial')
                basic(work/name, solver, args.bin_dir, coupled, parallel)
        edge_cases(work, args.bin_dir, parallel)


if __name__ == '__main__':
    main()
