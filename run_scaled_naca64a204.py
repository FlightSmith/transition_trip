"""Generate and run the original left wing at 1:8 scale with NACA 64A204.

Re_MAC is specified on the MODEL MAC and is not divided by the scale factor.
Run once for all nine Reynolds/angle conditions; derive per-Re envelopes from
the same solver records. No previous example inputs are modified.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

from curve_io import load_section
from generate_left_wing import (ROOT_Y_M, TIP_Y_M, ROOT_X_LE_M, TIP_X_LE_M,
                               ROOT_CHORD_M, TIP_CHORD_M, REFERENCE_MAC_M,
                               REFERENCE_AREA_M2, geometry_metadata)
from section_workflow import default_settings, execute, write_json

ROOT = Path(__file__).resolve().parent
CASE = ROOT / 'examples/left_wing_naca64a204_scale8'
SCALE = 1 / 8
REYNOLDS = (3e6, 5e6, 8e6)
ANGLES = (0., 1., 2.)


def write_section(path, coordinates):
    # Seven header lines, numeric data starting at line eight, as in user files.
    header = ['*KEYWORD', '*DEFINE_CURVE_TITLE', 'NACA 64A204 model section; X Z in metres',
              '$ LCID SIDR SFA SFO OFFA OFFO DATTYP', '1 0 1 1 0 0 0',
              '$ Model scale 1:8; dimensions already scaled', '$ ABSCISSA ORDINATE']
    text = '\n'.join(header + [f'{x:.15g} {z:.15g}' for x, z in coordinates] + ['*END', ''])
    path.write_text(text, encoding='ascii', newline='\n')


def prepare(case=CASE, airfoil=None):
    case = Path(case).resolve()
    airfoil = Path(airfoil or CASE / 'naca64a204.dat').resolve()
    xy = np.loadtxt(airfoil, skiprows=1)
    if xy.ndim != 2 or xy.shape[1] != 2 or not np.all(np.isfinite(xy)):
        raise ValueError('Expected a labeled normalized NACA 64A204 coordinate file')
    case.mkdir(parents=True, exist_ok=True)
    sections = case / 'sections'
    sections.mkdir(exist_ok=True)
    # Preserve the nominal NACA chord orientation. The geometric minimum-X
    # point differs slightly from the nominal origin on a cambered section.
    # Scale its true LE-to-mid-TE chord; let the importer rotate for XFOIL and
    # correct local alpha through alpha_reference='global_x'.
    leading = xy[np.argmin(xy[:, 0])]
    trailing = (xy[0] + xy[-1]) / 2
    chord_vector = trailing - leading
    reference_chord = np.linalg.norm(chord_vector)
    normalized = (xy - leading) / reference_chord
    stations = []
    for eta in np.linspace(0, 1, 21):
        chord = (ROOT_CHORD_M + eta * (TIP_CHORD_M - ROOT_CHORD_M)) * SCALE
        x_le = (ROOT_X_LE_M + eta * (TIP_X_LE_M - ROOT_X_LE_M)) * SCALE
        y = (ROOT_Y_M + eta * (TIP_Y_M - ROOT_Y_M)) * SCALE
        path = sections / f'section_state1_station{y:.8f}'
        write_section(path, normalized * chord + [x_le, 0.])
        actual = load_section(path)
        if not np.isclose(actual['chord_m'], chord, rtol=1e-8, atol=1e-10):
            raise ValueError('Generated section chord does not match the requested model scale')
        stations.append(dict(eta=float(eta), y_m=float(y), x_le_m=float(x_le), chord_m=float(chord)))
    origin = ('Original user left-wing planform scaled uniformly 1:8 (all lengths /8 and area /64); '
              'model MAC=0.578625 m; root/tip chords=0.695/0.125 m. '
              'NACA 64A204 from the recorded PDAS source is used at every section with zero twist. '
              'Angles refer to the nominal NACA chord/global X; the geometric-LE chord rotation is corrected in XFOIL local alpha. '
              'Mach=0.2; model-MAC Reynolds numbers 3e6, 5e6, 8e6; angles 0,1,2 degrees. '
              'These Reynolds numbers are not divided by 8. Independent 2D XFOIL strips do not solve 3D sweep flow.')
    settings = default_settings('sections')
    settings.update(name='NACA 64A204 left wing, model scale 1:8', mode='xfoil',
                    alpha_reference='global_x',
                    reference_mac_m=REFERENCE_MAC_M * SCALE,
                    reference_area_m2=REFERENCE_AREA_M2 * SCALE**2, data_origin=origin,
                    states=[dict(state=i + 1, alpha_deg=alpha, reynolds_mac=reynolds)
                            for i, (reynolds, alpha) in enumerate((re, a) for re in REYNOLDS for a in ANGLES)])
    settings['xfoil']['panels'] = 280
    settings['xfoil']['retry_panels'] = [280, 240, 200]
    settings['xfoil']['retry_iterations'] = 600
    write_json(case / 'workflow.json', settings)
    write_json(case / 'model_geometry.json', dict(
        scale=SCALE, full_scale=geometry_metadata(), model_reference_mac_m=REFERENCE_MAC_M*SCALE,
        model_reference_area_m2=REFERENCE_AREA_M2*SCALE**2, model_stations=stations,
        airfoil_file=airfoil.name, airfoil_sha256=hashlib.sha256(airfoil.read_bytes()).hexdigest(),
        coordinate_chord_normalization=dict(leading_edge=leading.tolist(), trailing_edge=trailing.tolist(),
                                            chord=float(reference_chord)),
        notes=origin))
    return case / 'workflow.json'


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--xfoil-exe', help='External XFOIL executable; alternatively set XFOIL_EXE')
    parser.add_argument('--output', type=Path, default=ROOT / 'results/left_wing_naca64a204_scale8')
    parser.add_argument('--prepare-only', action='store_true')
    parser.add_argument('--summarize-only', action='store_true', help='Rebuild reports from existing solver records')
    parser.add_argument('--no-plot', action='store_true')
    args = parser.parse_args(argv)
    try:
        if args.prepare_only and args.summarize_only:
            raise ValueError('Choose either prepare-only or summarize-only')
        if not args.summarize_only:
            config = prepare()
            code = execute(config, args.output / 'all_conditions', xfoil_executable=args.xfoil_exe,
                           prepare_only=args.prepare_only, no_plot=args.no_plot)
            if args.prepare_only:
                return code
        from scaled_wing_report import summarize
        summarize(args.output, plot=not args.no_plot)
        print(f'Comparison report: {(args.output / "report.html").resolve()}')
        # Reports can be complete while the requested aerodynamic envelope is
        # infeasible. Preserve the same exit-status convention as other routes.
        states = [json.loads(path.read_text(encoding='utf-8'))['status']
                  for path in (args.output / 'all_conditions').glob('xc_*/summary.json')]
        return 0 if states and all(s == 'feasible' for s in states) else 2
    except (OSError, ValueError, KeyError, TypeError) as exc:
        print(f'Scaled-wing run error: {exc}', file=sys.stderr)
        return 2


if __name__ == '__main__':
    raise SystemExit(main())
