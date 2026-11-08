"""Portable menu or command-line entry point; run setup_environment.py first."""
from __future__ import annotations

from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent


def entered_path(prompt):
    value = input(prompt).strip().strip('"').strip("'")
    return Path(value).expanduser() if value else None


def menu():
    print('\nTransition-trip sizing')
    print('1. M6 synthetic example using Cp')
    print('2. Left wing using idealized Blasius profiles (7%, 10%, 12%)')
    print('3. Left wing Cp-input demonstration (Cp=0; 7%, 10%, 12%)')
    print('4. My Cp CSV, using the left-wing geometry')
    print('5. My JSON configuration (Cp or profiles)')
    print('0. Exit')
    while True:
        choice = input('Choose 0-5: ').strip()
        if choice in {'0', '1', '2', '3', '4', '5'}:
            break
        print('Enter a number from 0 through 5.')
    if choice == '0':
        return []
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
        codes = [calculate(job) for job in jobs]
        return max(codes, default=0)
    except (EOFError, KeyboardInterrupt):
        print('\nCancelled.')
        return 130
    except (ValueError, OSError) as exc:
        print(f'Input error: {exc}', file=sys.stderr)
        return 2


if __name__ == '__main__':
    raise SystemExit(main())
