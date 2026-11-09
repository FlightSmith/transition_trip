"""Portable menu or command-line entry point; run setup_environment.py first."""
from __future__ import annotations

from pathlib import Path
import os
import sys

ROOT = Path(__file__).resolve().parent


def entered_path(prompt):
    value = input(prompt).strip().strip('"').strip("'")
    return Path(value).expanduser() if value else None


def section_menu(mode):
    from section_workflow import default_settings, write_json
    folder = entered_path(f'Section/curve folder [{ROOT / "synthetic_Cp"}]: ') or ROOT / 'synthetic_Cp'
    folder = folder.resolve()
    settings = default_settings(folder)
    settings['name'] = f'{folder.name} ({mode})'
    settings['mode'] = mode
    for key, prompt in (('mach', 'Mach'), ('reynolds_mac', 'Reynolds number based on MAC'),
                        ('reference_mac_m', 'MAC in metres')):
        settings[key] = float(input(f'{prompt} [{settings[key]:g}]: ').strip() or settings[key])
    angles = input('Angles for states 1, 2, 3 in degrees [0, 1, 2]: ').strip() or '0,1,2'
    values = [float(v.strip()) for v in angles.split(',')]
    if len(values) != 3:
        raise ValueError('Enter three comma-separated angles; use a workflow JSON for other state IDs/counts')
    settings['states'] = [dict(state=i + 1, alpha_deg=value) for i, value in enumerate(values)]
    positions = input('Trip x/c values [0.07, 0.10, 0.12]: ').strip() or '.07,.10,.12'
    settings['positions_xc'] = [float(v.strip()) for v in positions.split(',')]
    if mode == 'curves':
        kind = input('Values: cp, pressure (absolute Pa), or gauge_pressure [cp]: ').strip().lower() or 'cp'
        settings['value_kind'] = kind
        if kind in ('pressure', 'absolute_pressure', 'gauge_pressure'):
            settings['pressure_reference'] = dict(
                dynamic_pressure_pa=float(input('Dynamic pressure q in Pa: ')),
                static_pressure_pa=float(input('Reference static pressure in Pa (same absolute/gauge convention): ')))
        demonstration = folder == (ROOT / 'synthetic_Cp').resolve()
        default_init = 'synthetic_le_zero' if demonstration else 'stagnation_ramp'
        print('Initial thickness: stagnation_ramp approximates the upstream flow; provided uses your theta0.')
        print('synthetic_le_zero is only the supplied synthetic demonstration assumption.')
        method = input(f'Initialization method [{default_init}]: ').strip() or default_init
        settings['cp_initialization'] = dict(method=method)
        if method == 'provided':
            settings['cp_initialization']['theta0_m'] = float(input('Initial momentum thickness theta0 in metres: '))
    else:
        executable = input('XFOIL executable path [XFOIL_EXE environment variable or xfoil on PATH]: ').strip().strip('"')
        if executable:
            settings['xfoil']['executable'] = str(Path(executable).expanduser().resolve()) if Path(executable).is_file() else executable
    default_output = ROOT / f'results/section_{mode}'
    output = (entered_path(f'Output folder [{default_output}]: ') or default_output).resolve()
    output.mkdir(parents=True, exist_ok=True)
    try:
        settings['input_directory'] = Path(os.path.relpath(folder, output)).as_posix()
    except ValueError:  # Different Windows drives need an absolute input path.
        settings['input_directory'] = str(folder)
    config = output / 'workflow.json'
    write_json(config, settings)
    print(f'Saved reusable settings: {config}')
    job = ['--section-config', str(config), '--output', str(output)]
    if mode == 'xfoil' and executable:
        job.extend(['--xfoil-exe', settings['xfoil']['executable']])
    return [job]


def menu():
    print('\nTransition-trip sizing')
    print('1. M6 synthetic example using Cp')
    print('2. Left wing using idealized Blasius profiles (7%, 10%, 12%)')
    print('3. Left wing Cp-input demonstration (Cp=0; 7%, 10%, 12%)')
    print('4. My Cp CSV, using the left-wing geometry')
    print('5. My JSON configuration (Cp or profiles)')
    print('6. Section files with Cp or pressure curves (states and Y stations)')
    print('7. Section files with XFOIL Cp and laminar boundary layer')
    print('0. Exit')
    while True:
        choice = input('Choose 0-7: ').strip()
        if choice in {'0', '1', '2', '3', '4', '5', '6', '7'}:
            break
        print('Enter a number from 0 through 7.')
    if choice == '0':
        return []
    if choice in ('6', '7'):
        return section_menu('curves' if choice == '6' else 'xfoil')
    if choice == '1':
        return [['--config', str(ROOT / 'examples/m6_synthetic/config.json'),
                 '--output', str(ROOT / 'results/m6_cp')]]
    if choice in {'2', '3'}:
        folder = 'left_wing_flatplate' if choice == '2' else 'left_wing_cp'
        if choice == '3':
            print('Cp=0 is a format demonstration, not a prediction for this wing.')
        return [['--config', str(ROOT / f'examples/{folder}/config_xc{tag}.json'),
                 '--output', str(ROOT / f'results/{folder}/xc{tag}')] for tag in ('07', '10', '12')]
    if choice == '4':
        print('Uses the supplied left-wing geometry: MAC=4.629 m, Mach=0.2, Re_MAC=8 million.')
        print('The template AoA=0 is a placeholder; use option 5 to set actual angle-of-attack conditions.')
        print('Use case_01, eta stations 0 to 1 in steps of 0.05, and upper/lower rows.')
        print('For other stations or conditions, choose option 5 and edit a copied configuration.')
        cp = entered_path('Path to your Cp CSV: ')
        if cp is None:
            raise ValueError('A Cp CSV path is required')
        raw_position = input('Trip x/c [0.10]: ').strip() or '0.10'
        position = float(raw_position)
        if not 0 < position < 1:
            raise ValueError('Trip x/c must be between 0 and 1')
        output = entered_path(f'Output folder [{ROOT / "results/my_cp_run"}]: ') or ROOT / 'results/my_cp_run'
        return [['--config', str(ROOT / 'examples/left_wing_cp/config_xc10.json'),
                 '--cp-file', str(cp.resolve()), '--position-xc', str(position), '--output', str(output.resolve())]]
    config = entered_path('Path to your JSON configuration: ')
    if config is None:
        raise ValueError('A configuration path is required')
    output = entered_path(f'Output folder [{ROOT / "results/my_case"}]: ') or ROOT / 'results/my_case'
    return [['--config', str(config.resolve()), '--output', str(output.resolve())]]


def main(argv=None):
    arguments = list(sys.argv[1:] if argv is None else argv)
    if arguments and ('--help' in arguments or '-h' in arguments):
        # Help works before dependencies have been installed.
        print('Run without arguments for the menu. Or use:\n'
              '  python run_trip.py --config CONFIG.json --output OUTPUT_DIR\n'
              '  python run_trip.py --config CONFIG.json --cp-file CP.csv --position-xc 0.10 --output OUTPUT_DIR\n'
              '  python run_trip.py --section-config examples/section_curves/workflow.json --output results/section_curves\n'
              '  python run_trip.py --section-config WORKFLOW.json --mode xfoil --xfoil-exe PATH --output OUTPUT_DIR\n'
              'Section options: --value-kind cp|pressure|gauge_pressure, --dynamic-pressure-pa Q,\n'
              '  --static-pressure-pa P_REF, --prepare-only, --no-plot\n'
              'Optional: --profiles-file PROFILES.csv, --no-plot\n'
              'CLI paths are relative to your current folder; paths inside JSON are relative to that JSON file.\n'
              'First install dependencies with: python setup_environment.py')
        return 0
    try:
        from trip_design import main as calculate
    except ModuleNotFoundError as exc:
        print(f'Missing dependency: {exc.name}. Run setup_environment.py, then use the .venv Python.', file=sys.stderr)
        return 1
    try:
        jobs = [arguments] if arguments else menu()
        codes = []
        for job in jobs:
            if '--section-config' in job:
                from section_workflow import main as section_calculate
                job = ['--config' if item == '--section-config' else item for item in job]
                codes.append(section_calculate(job))
            else:
                codes.append(calculate(job))
        return max(codes, default=0)
    except (EOFError, KeyboardInterrupt):
        print('\nCancelled.')
        return 130
    except (ValueError, OSError) as exc:
        print(f'Input error: {exc}', file=sys.stderr)
        return 2


if __name__ == '__main__':
    raise SystemExit(main())
