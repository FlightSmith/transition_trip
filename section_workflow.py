"""Import ordered section curves or run optional XFOIL from dimensional sections."""
from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import sys

import numpy as np

from aero import edge_state, tunnel_state
from curve_io import load_section, map_curve_to_section, parse_curve_name, read_curve
from trip_design import positive, run_design, write_outputs

CP_FIELDS = ['condition_id', 'eta', 'surface', 'x_c', 's_m', 'cp', 'theta0_m']
BL_FIELDS = ['condition_id', 'eta', 'surface', 'x_c', 'status', 'reason', 'theta_m',
             'shape_factor', 'nu_m2_s', 'ue_m_s', 'edge_mach', 'laminar_confirmed',
             'transition_xc', 'cf', 'delta_star_m', 'h_mass', 'solver_method']


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n', encoding='utf-8', newline='\n')


def write_table(path, fields, rows):
    with Path(path).open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def default_settings(input_directory):
    """Explicit demo defaults; physical geometry always comes from section files."""
    return dict(schema_version=1, name='Section-curve example', input_directory=str(input_directory),
                data_start_line=8, section_state=1, mode='curves', value_kind='cp',
                mach=.2, reynolds_mac=8e6, reference_mac_m=4.629, temperature_k=288.15,
                states=[dict(state=i + 1, alpha_deg=float(i)) for i in range(3)],
                positions_xc=[.07, .10, .12], alpha_reference='section_chord',
                cp_initialization=dict(method='synthetic_le_zero'),
                trip=dict(type='grit', height_margin=1.1, max_k_delta99=.8, available_heights_mm=[]),
                xfoil=dict(executable='xfoil', ncrit=9., iterations=200, timeout_s=60., panels=200))


def discover(input_dir):
    records = {'cp': {}, 'section': {}}
    for path in sorted(Path(input_dir).iterdir()):
        if not path.is_file() or not path.name.lower().startswith(('cp_dist_', 'section_')):
            continue
        metadata = parse_curve_name(path)
        identity = (metadata['state'], round(metadata['station_y_m'], 9))
        bucket = records[metadata['kind']]
        if identity in bucket:
            raise ValueError(f'duplicate state/station file: {path.name}')
        bucket[identity] = path
    if not records['section']:
        raise ValueError('No section_stateN_stationY files were found')
    return records


def validate_settings(settings):
    if settings.get('schema_version') != 1:
        raise ValueError('workflow schema_version must be 1')
    if settings.get('mode') not in ('curves', 'xfoil'):
        raise ValueError('mode must be curves or xfoil')
    if settings.get('value_kind', 'cp') not in ('cp', 'pressure', 'absolute_pressure', 'gauge_pressure'):
        raise ValueError('value_kind must be cp, pressure, absolute_pressure, or gauge_pressure')
    positive(settings['reference_mac_m'], 'reference_mac_m')
    for name, default in (('section_state', 1), ('data_start_line', 8)):
        value = positive(settings.get(name, default), name)
        if not value.is_integer():
            raise ValueError(f'{name} must be a positive integer')
    states = settings['states']
    ids = [int(s['state']) for s in states]
    if not states or len(ids) != len(set(ids)) or any(i <= 0 for i in ids):
        raise ValueError('states must contain unique positive integer state IDs')
    if any(float(s['state']) != int(s['state']) for s in states):
        raise ValueError('state IDs must be integers')
    for state in states:
        if not math.isfinite(float(state['alpha_deg'])):
            raise ValueError('Every state requires a finite alpha_deg')
        mach = positive(state.get('mach', settings['mach']), 'mach')
        if mach > .3:
            raise ValueError('This trip-sizing workflow is restricted to freestream Mach <=0.30')
        positive(state.get('reynolds_mac', settings['reynolds_mac']), 'reynolds_mac')
        positive(state.get('temperature_k', settings.get('temperature_k', 288.15)), 'temperature_k')
    positions = [float(x) for x in settings['positions_xc']]
    if not positions or len(set(positions)) != len(positions) or any(not 0 < x < 1 for x in positions):
        raise ValueError('positions_xc must be unique fractions in (0,1)')
    if settings.get('alpha_reference', 'section_chord') not in ('section_chord', 'global_x'):
        raise ValueError('alpha_reference must be section_chord or global_x')
    method = settings.get('cp_initialization', {}).get('method')
    if settings['mode'] == 'curves' and method not in ('provided', 'stagnation_ramp', 'synthetic_le_zero'):
        raise ValueError('Cp initialization must explicitly choose provided, stagnation_ramp, or synthetic_le_zero')
    return positions


def prepare_geometry(settings, records, output):
    section_state = int(settings.get('section_state', 1))
    ys = sorted({y for _, y in records['section']}, key=abs)
    if len(ys) < 2:
        raise ValueError('At least two sampled Y stations are required for this wing workflow')
    if ys[0] * ys[-1] < 0:
        raise ValueError('Process left and right semispans separately')
    span = abs(ys[-1] - ys[0])
    sections, stations, provenance = {}, [], []
    for y in ys:
        source = records['section'].get((section_state, y))
        if source is None:
            raise ValueError(f'missing reference section_state{section_state} at Y={y:g}')
        section = load_section(source, start_line=int(settings.get('data_start_line', 8)))
        # One static wing geometry is required for an envelope over states.
        for state in settings['states']:
            other_path = records['section'].get((int(state['state']), y))
            if other_path is not None and other_path != source:
                other = load_section(other_path, start_line=int(settings.get('data_start_line', 8)))
                if (len(other['x_m']) != len(section['x_m']) or
                    not np.allclose(other['x_m'], section['x_m'], rtol=0, atol=1e-8) or
                    not np.allclose(other['z_m'], section['z_m'], rtol=0, atol=1e-8)):
                    raise ValueError('State-varying section geometry requires separate workflows')
        eta = abs(y - ys[0]) / span
        station = dict(eta=eta, y_m=y, x_le_m=section['x_le_m'], chord_m=section['chord_m'],
                       streamwise_chord_m=section['streamwise_chord_m'])
        sections[y] = section
        stations.append(station)
        section_dir = output / 'airfoils'
        section_dir.mkdir(exist_ok=True)
        airfoil = section_dir / f'station_{len(stations):02d}.dat'
        with airfoil.open('w', encoding='ascii', newline='\n') as handle:
            handle.write(f'Section_at_Y_{y:.6f}\n')
            for x, z in zip(section['normalized_x'], section['normalized_z']):
                handle.write(f'{x:.12g} {z:.12g}\n')
        section['airfoil_path'] = airfoil
        provenance.append(dict(file=source.name, sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                               y_m=y, eta=eta, chord_m=section['chord_m'],
                               streamwise_chord_m=section['streamwise_chord_m'],
                               x_le_m=section['x_le_m'], z_le_m=section['z_le_m'],
                               incidence_deg=section['incidence_deg'], normalized_file=airfoil.relative_to(output).as_posix()))
    sweep = math.degrees(math.atan2(stations[-1]['x_le_m'] - stations[0]['x_le_m'], span))
    wing = dict(semispan_m=span, sweep_le_deg=sweep, reference_length_m=float(settings['reference_mac_m']), stations=stations)
    if 'reference_area_m2' in settings:
        wing['reference_area_m2'] = positive(settings['reference_area_m2'], 'reference_area_m2')
    return sections, wing, provenance


def initialized_cp_rows(surface, settings, state, flow):
    """Keep real arc distance and make any upstream initialization assumption explicit."""
    x, cp, s = (np.asarray(surface[k]) for k in ('x_c', 'cp', 's_m'))
    init = state.get('cp_initialization', settings['cp_initialization'])
    method = init['method']
    if method == 'synthetic_le_zero':
        if abs(s[0]) > 1e-10 or abs(x[0]) > 1e-10:
            raise ValueError('synthetic_le_zero requires data at the actual geometric leading edge')
        start, theta0 = 0, 0.
    else:
        valid = np.flatnonzero((s > 0) & (x > 0))
        if not len(valid):
            raise ValueError('no nonstagnant downstream Cp point available for initialization')
        start = int(valid[0])
        if method == 'provided':
            theta0 = positive(init['theta0_m'], 'provided theta0_m')
        elif method == 'stagnation_ramp':
            ue = float(edge_state(np.array([cp[start]]), flow, flow['mach'])['U_e_m_s'][0])
            if ue <= 0:
                raise ValueError('stagnation-ramp initialization needs positive edge velocity')
            theta0 = math.sqrt(.075 * flow['nu_inf_m2_s'] * s[start] / ue)
        else:
            raise ValueError('unsupported Cp initialization method')
    return x[start:], cp[start:], s[start:], theta0, method


def execute(settings_path, output, *, mode=None, value_kind=None, dynamic_pressure_pa=None,
            static_pressure_pa=None, xfoil_executable=None, prepare_only=False, no_plot=False):
    settings_path = Path(settings_path).resolve()
    settings = json.loads(settings_path.read_text(encoding='utf-8-sig'))
    if mode is not None:
        settings['mode'] = mode
    if value_kind is not None:
        settings['value_kind'] = value_kind
    # Explicit CLI reference pressures override per-state as well as global values.
    for name, value in (('dynamic_pressure_pa', dynamic_pressure_pa), ('static_pressure_pa', static_pressure_pa)):
        if value is not None:
            settings.setdefault('pressure_reference', {})[name] = value
            for state in settings['states']:
                state[name] = value
    positions = validate_settings(settings)
    input_dir = (settings_path.parent / settings['input_directory']).resolve()
    records = discover(input_dir)
    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=True)
    sections, wing, geometry = prepare_geometry(settings, records, output)
    if settings['mode'] == 'curves':
        state_ids = {int(state['state']) for state in settings['states']}
        missing_geometry = sorted({y for state_id, y in records['cp'] if state_id in state_ids and y not in sections})
        if missing_geometry:
            raise ValueError(f'Curve stations lack section geometry: {missing_geometry}')
    conditions, curve_rows, integral_rows, solver_cp_rows = [], [], [], []
    diagnostics, input_files, calls = [], [], []
    for state in settings['states']:
        state_id = int(state['state'])
        condition = dict(id=f'state{state_id}', alpha_deg=float(state['alpha_deg']),
                         mach=float(state.get('mach', settings['mach'])),
                         reynolds=float(state.get('reynolds_mac', settings['reynolds_mac'])),
                         temperature_k=float(state.get('temperature_k', settings.get('temperature_k', 288.15))))
        conditions.append(condition)
        flow = tunnel_state(condition['mach'], condition['reynolds'], wing['reference_length_m'], condition['temperature_k'])
        for station in wing['stations']:
            y, eta = station['y_m'], station['eta']
            section = sections[y]
            common = dict(condition_id=condition['id'], eta=eta)
            if settings['mode'] == 'curves':
                path = records['cp'].get((state_id, y))
                if path is None:
                    raise ValueError(f'missing curve for state{state_id}, Y={y:g}; all configured cases are required')
                reference = settings.get('pressure_reference', {})
                q = state.get('dynamic_pressure_pa', reference.get('dynamic_pressure_pa'))
                p_ref = state.get('static_pressure_pa', reference.get('static_pressure_pa'))
                mapped = map_curve_to_section(read_curve(path, int(settings.get('data_start_line', 8))), section,
                                              settings.get('value_kind', 'cp'), q, p_ref)
                input_files.append(dict(file=path.name, sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                                        state=state_id, y_m=y, value_kind=settings.get('value_kind', 'cp'),
                                        dynamic_pressure_pa=q, static_pressure_pa=p_ref))
                for surface_name, surface in mapped.items():
                    x, cp, s, theta0, method = initialized_cp_rows(surface, settings, state, flow)
                    for i, (xx, pp, ss) in enumerate(zip(x, cp, s)):
                        curve_rows.append(dict(common, surface=surface_name, x_c=float(xx), cp=float(pp), s_m=float(ss),
                                               theta0_m=theta0 if i == 0 else ''))
                    diagnostics.append(dict(state=state_id, y_m=y, surface=surface_name, initialization=method,
                                            initial_xc=float(x[0]), initial_s_m=float(s[0]), theta0_m=theta0))
            else:
                from xfoil_backend import run_xfoil_with_retries, sample_laminar_station
                options = settings.get('xfoil', {})
                executable = xfoil_executable or os.environ.get('XFOIL_EXE') or options.get('executable', 'xfoil')
                local_re = condition['reynolds'] * section['chord_m'] / wing['reference_length_m']
                local_alpha = condition['alpha_deg']
                if settings.get('alpha_reference', 'section_chord') == 'global_x':
                    local_alpha -= section['incidence_deg']
                plan = dict(state=state_id, y_m=y, reynolds_local=local_re, mach=condition['mach'],
                            alpha_local_deg=local_alpha, section_incidence_deg=section['incidence_deg'],
                            airfoil_file=section['airfoil_path'].relative_to(output).as_posix())
                calls.append(plan)
                if prepare_only:
                    continue
                try:
                    solved = run_xfoil_with_retries(executable, section['airfoil_path'], output / 'xfoil_runs', local_re,
                                       condition['mach'], local_alpha, ncrit=float(options.get('ncrit', 9.)),
                                       iterations=int(options.get('iterations', 200)),
                                       timeout_s=float(options.get('timeout_s', 60.)), panels=int(options.get('panels', 200)),
                                       retry_panels=options.get('retry_panels', ()),
                                       retry_iterations=options.get('retry_iterations', 600))
                    audit = {k: solved[k] for k in ('converged', 'reason', 'polar', 'ncrit', 'panels',
                             'iterations', 'executable_sha256', 'airfoil_sha256', 'transition_model',
                             'dump_note', 'returncode', 'solver_warnings') if k in solved}
                    for k in ('workdir', 'stdout_path', 'commands_path'):
                        if solved.get(k):
                            audit[k] = Path(solved[k]).relative_to(output).as_posix()
                    audit['attempts'] = copy.deepcopy(solved.get('attempts', []))
                    for attempt in audit['attempts']:
                        for k in ('workdir', 'stdout_path', 'commands_path'):
                            if attempt.get(k):
                                attempt[k] = Path(attempt[k]).relative_to(output).as_posix()
                    if not solved.get('converged'):
                        failure = solved.get('reason') or 'XFOIL did not converge'
                        diagnostics.append(dict(plan, status='failed', reason=failure, solver=audit))
                        solved = None
                    else:
                        audit['h_column_interpretation'] = solved['boundary_layer']['h_column_interpretation']
                        diagnostics.append(dict(plan, status='converged', solver=audit))
                        failure = None
                    if solved is None:
                        raise ValueError(failure)
                    cp_data = solved['cp']
                    for i, (xx, pp) in enumerate(zip(cp_data['x_c'], cp_data['cp'])):
                        solver_cp_rows.append(dict(common, contour_index=i, x_c=float(xx), cp=float(pp)))
                    failure = None
                except (ValueError, OSError, RuntimeError) as exc:
                    solved, failure = None, str(exc)
                    if not diagnostics or diagnostics[-1].get('state') != state_id or diagnostics[-1].get('y_m') != y:
                        diagnostics.append(dict(plan, status='failed', reason=failure))
                for xc in positions:
                    for surface in ('upper', 'lower'):
                        record = dict(common, surface=surface, x_c=xc, status='invalid', reason=failure or '')
                        try:
                            if solved is None:
                                raise ValueError(failure)
                            sampled = sample_laminar_station(solved, surface, xc)
                            record.update(status='valid', reason='', theta_m=sampled['theta_c'] * section['chord_m'],
                                          delta_star_m=sampled['delta_star_c'] * section['chord_m'],
                                          shape_factor=sampled['h_kinematic'], h_mass=sampled['h_mass'], cf=sampled['cf'],
                                          ue_m_s=sampled['ue_over_uinf'] * flow['U_inf_m_s'],
                                          nu_m2_s=flow['nu_inf_m2_s'], edge_mach=sampled['edge_mach'],
                                          transition_xc=sampled['transition_xc'], laminar_confirmed='true',
                                          solver_method='XFOIL integral BL; Hk-matched quartic reconstruction for u(y)')
                        except (ValueError, KeyError) as exc:
                            record.update(reason=str(exc), laminar_confirmed='false')
                        integral_rows.append(record)
    notes = [
        'Geometry comes from supplied sections; coordinates and reference lengths must already be at wind-tunnel model scale.',
        'Eta spans only the supplied Y stations; no unprovided wing root, tip or full planform is inferred.',
        'All coordinates are interpreted in metres. Upper/lower curve order must match the section contour order.',
        'Angle mapping is explicit in the workflow settings. XFOIL is a 2D section method without crossflow.',
    ]
    if settings.get('data_origin'):
        notes.insert(0, str(settings['data_origin']))
    if settings['mode'] == 'curves':
        notes.append('Cp integration uses true section arc length and the explicitly selected initial-thickness assumption.')
        notes.append('synthetic_le_zero is a demonstration assumption; Cp=1 at the geometric LE is not exact compressible stagnation.')
        filename, source, fields, rows = 'cp.csv', 'cp', CP_FIELDS, curve_rows
    else:
        notes.extend(['XFOIL computes its own Cp from geometry; imported Cp is not prescribed to its viscous solver.',
                      'The undisturbed calculation allows natural transition; the requested trip is not imposed in XFOIL.',
                      'DUMP supplies integral quantities, not resolved u(y). Low-Mach Hk/theta reconstruction is approximate.',
                      'Re_k reconstruction uses constant freestream kinematic viscosity; DUMP does not provide nu(y) or wall temperature.',
                      'Laminar confirmation is predicted by the solver transition model, not established experimentally.'])
        filename, source, fields, rows = 'boundary_layer.csv', 'integral', BL_FIELDS, integral_rows
    metadata = dict(settings=settings, geometry=geometry, imported_files=input_files, solver_plan=calls,
                    diagnostics=diagnostics, notes=notes, status='prepared' if prepare_only else 'processed')
    # Retain portable source identification; no executable machine path is published.
    metadata['settings'] = copy.deepcopy(settings)
    metadata['settings']['input_directory'] = input_dir.name
    if 'xfoil' in metadata['settings']:
        metadata['settings']['xfoil']['executable'] = Path(str(metadata['settings']['xfoil'].get('executable', 'xfoil'))).name
    write_json(output / 'workflow_summary.json', metadata)
    if prepare_only:
        print(f'Prepared normalized section files and solver plan in {output}')
        return 0
    write_table(output / filename, fields, rows)
    if solver_cp_rows:
        write_table(output / 'xfoil_cp.csv', ['condition_id', 'eta', 'contour_index', 'x_c', 'cp'], solver_cp_rows)
    statuses = []
    for xc in positions:
        label = f'xc_{xc}'
        cfg = dict(schema_version=1, name=f"{settings.get('name', 'Section curves')}: {100*xc:g}% chord ({settings['mode']})",
                   data_origin=' '.join(notes), wing=wing, conditions=conditions,
                   trip=dict(settings.get('trip', {}), position_xc={'upper': xc, 'lower': xc}),
                   analysis={'source': source, 'cp_file' if source == 'cp' else 'integral_file': filename})
        if settings.get('airfoil_name'):
            cfg['airfoil_name'] = str(settings['airfoil_name'])
        cfg['trip'].setdefault('type', 'grit')
        cfg['trip'].setdefault('height_margin', 1.1)
        cfg['trip'].setdefault('max_k_delta99', .8)
        config_path = output / f'config_{label}.json'
        write_json(config_path, cfg)
        final_cfg, summary, condition_rows = run_design(config_path)
        write_outputs(output / label, final_cfg, summary, condition_rows, not no_plot)
        statuses.append(summary['status'])
        n = sum(d['status'] == 'feasible' for d in summary['design'])
        print(f'{100*xc:g}% chord: {n}/{len(summary["design"])} station/surface designs feasible; '
              f'uniform-height option {summary["uniform_height"]["status"]}.')
        print(f'Report: {output / label / "report.html"}')
    return 0 if all(status == 'feasible' for status in statuses) else 2


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, required=True, help='Section workflow JSON, not the ordinary sizing JSON')
    parser.add_argument('--output', type=Path, default=Path('results/section_workflow'))
    parser.add_argument('--mode', choices=('curves', 'xfoil'))
    parser.add_argument('--value-kind', choices=('cp', 'pressure', 'absolute_pressure', 'gauge_pressure'))
    parser.add_argument('--dynamic-pressure-pa', type=float)
    parser.add_argument('--static-pressure-pa', type=float)
    parser.add_argument('--xfoil-exe')
    parser.add_argument('--prepare-only', action='store_true')
    parser.add_argument('--no-plot', action='store_true')
    args = parser.parse_args(argv)
    try:
        return execute(args.config, args.output, mode=args.mode, value_kind=args.value_kind,
                       dynamic_pressure_pa=args.dynamic_pressure_pa, static_pressure_pa=args.static_pressure_pa,
                       xfoil_executable=args.xfoil_exe, prepare_only=args.prepare_only, no_plot=args.no_plot)
    except (OSError, ValueError, KeyError, TypeError) as exc:
        print(f'Section workflow error: {exc}', file=sys.stderr)
        return 2


if __name__ == '__main__':
    raise SystemExit(main())
