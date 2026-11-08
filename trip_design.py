"""Preliminary fixed-position wing transition-trip sizing. See README.md."""
from __future__ import annotations

import argparse
import csv
import hashlib
import html
import json
import math
import os
from pathlib import Path
import sys

import numpy as np

from aero import (tunnel_state, edge_state, thwaites, build_surrogate_profile,
                  profile_metrics, required_height, roughness_re)


def positive(value, label, allow_zero=False):
    value = float(value)
    if not math.isfinite(value) or (value < 0 if allow_zero else value <= 0):
        raise ValueError(f"{label} must be finite and {'nonnegative' if allow_zero else 'positive'}")
    return value


def key(condition, eta, surface):
    return str(condition), round(float(eta), 10), str(surface)


def laminar_limit(value, label):
    value = positive(value, label)
    if value > 1:
        raise ValueError(f'{label} must be in (0,1]')
    return value


def read_csv(path, columns):
    with path.open(newline='', encoding='utf-8-sig') as handle:
        reader = csv.DictReader(handle)
        missing = set(columns) - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path.name}: missing columns {sorted(missing)}")
        rows = list(reader)
    if not rows:
        raise ValueError(f"{path.name}: empty data")
    return rows


def numeric(rows, column):
    values = np.array([float(r[column]) for r in rows], dtype=float)
    if not np.all(np.isfinite(values)):
        raise ValueError(f"nonfinite values in {column}")
    return values


def validate_config(cfg):
    if cfg.get('schema_version') != 1:
        raise ValueError('schema_version must be 1')
    wing = cfg['wing']
    for field in ('semispan_m', 'reference_length_m'):
        positive(wing[field], field)
    sweep = float(wing['sweep_le_deg'])
    if not math.isfinite(sweep) or abs(sweep) >= 80:
        raise ValueError('sweep_le_deg must be finite and between -80 and 80')
    stations = wing['stations']
    etas = [float(s['eta']) for s in stations]
    if not etas or any(not 0 <= e <= 1 for e in etas) or len(set(etas)) != len(etas):
        raise ValueError('wing.stations must have unique eta values in [0,1]')
    for station in stations:
        positive(station['chord_m'], 'chord_m')
        if 'y_m' in station and not math.isfinite(float(station['y_m'])):
            raise ValueError('station y_m must be finite')
    conditions = cfg['conditions']
    ids = [str(c['id']) for c in conditions]
    if not ids or len(set(ids)) != len(ids):
        raise ValueError('conditions must have unique IDs')
    for c in conditions:
        for field in ('mach', 'reynolds'):
            positive(c[field], field)
        positive(c.get('temperature_k', 288.15), 'temperature_k')
        if not math.isfinite(float(c['alpha_deg'])):
            raise ValueError('alpha_deg must be finite')
    trip = cfg['trip']
    if trip['type'] not in ('grit', 'zigzag', 'wire', 'custom'):
        raise ValueError('trip.type must be grit, zigzag, wire, or custom')
    positions = trip['position_xc']
    if not positions or not set(positions).issubset({'upper', 'lower'}):
        raise ValueError('position_xc must map upper and/or lower to x/c')
    for surface, value in positions.items():
        if not 0 < float(value) < 1:
            raise ValueError(f'{surface} position must be in (0,1)')
    margin = positive(trip.get('height_margin', 1.1), 'height_margin')
    if margin < 1:
        raise ValueError('height_margin must be >=1')
    cap = positive(trip.get('max_k_delta99', 0.8), 'max_k_delta99')
    if cap > 1:
        raise ValueError('max_k_delta99 must be <=1 for this tool')
    sizes = trip.get('available_heights_mm', [])
    if sizes:
        for size in sizes:
            positive(size, 'available_heights_mm')
        if len(set(sizes)) != len(sizes):
            raise ValueError('available_heights_mm contains duplicates')
    calibration = trip.get('calibration')
    if calibration is None:
        if trip['type'] != 'grit':
            raise ValueError('zigzag, wire, and custom trips require trip.calibration')
        calibration = dict(critical_re_k=600., min_re_x=1e5, max_edge_mach=.95,
                           notes='Subsonic equiaxed distributed grit; NASA TN D-3579. Mach cap 0.95 is a tool guard.')
    else:
        for field in ('critical_re_k', 'min_re_x', 'max_edge_mach'):
            positive(calibration[field], f'calibration.{field}')
        if not str(calibration.get('notes', '')).strip():
            raise ValueError('calibration.notes must identify the applicability/source')
    analysis = cfg['analysis']
    if analysis['source'] not in ('cp', 'profiles'):
        raise ValueError('analysis.source must be cp or profiles')
    if 'laminar_until_xc' in analysis:
        laminar_limit(analysis['laminar_until_xc'], 'analysis.laminar_until_xc')
    return calibration


def cp_profile(rows, chord, xc, fs, condition, analysis):
    if float(condition['mach']) > .3 + 1e-10:
        raise ValueError('Cp fallback exceeds freestream Mach 0.30; supply laminar profiles')
    rows = sorted(rows, key=lambda r: float(r['x_c']))
    x = numeric(rows, 'x_c')
    if np.any(np.diff(x) <= 0) or np.any(x < 0) or np.any(x > 1):
        raise ValueError('Cp x_c must be unique and increasing in [0,1]')
    if not x[0] <= xc <= x[-1]:
        raise ValueError('Cp data do not bracket the requested x/c; no extrapolation')
    cp = numeric(rows, 'cp')
    explicit_s = all(r.get('s_m', '').strip() for r in rows)
    if any(r.get('s_m', '').strip() for r in rows) and not explicit_s:
        raise ValueError('s_m must be provided for every Cp row or omitted entirely')
    s = numeric(rows, 's_m') if explicit_s else x * chord
    if s[0] < 0 or np.any(np.diff(s) <= 0):
        raise ValueError('surface distance s_m must be nonnegative and strictly increasing')
    # Only upstream history is used; include an exact interpolated trip station.
    before = x < xc - 1e-12
    x_eval = np.r_[x[before], xc]
    s_eval = np.r_[s[before], np.interp(xc, x, s)]
    cp_eval = np.r_[cp[before], np.interp(xc, x, cp)]
    if len(s_eval) < 3:
        raise ValueError('at least three Cp stations through the trip are required')
    theta_text = rows[0].get('theta0_m', '').strip()
    if s_eval[0] > 0 and not theta_text:
        raise ValueError('Cp history starts downstream of s=0: first row requires theta0_m')
    theta0 = positive(theta_text, 'theta0_m', allow_zero=True) if theta_text else 0.
    if s_eval[0] > 0 and theta0 <= 0:
        raise ValueError('positive starting s_m requires positive theta0_m')
    edge = edge_state(cp_eval, fs, float(condition['mach']))
    if np.max(edge['M_e']) > .3 + 1e-10:
        raise ValueError('Cp fallback exceeds local Mach 0.30 upstream of trip; supply laminar profiles')
    bl = thwaites(s_eval, edge['U_e_m_s'], fs['nu_inf_m2_s'], theta0)
    if bool(np.asarray(bl['separated'])[-1]):
        raise ValueError('Thwaites predicts laminar separation at or upstream of trip')
    if 'laminar_until_xc' in analysis and xc >= float(analysis['laminar_until_xc']):
        raise ValueError('trip is at/downstream of supplied laminar extent')
    known_limits = [laminar_limit(r['known_transition_xc'], 'known_transition_xc')
                    for r in rows if r.get('known_transition_xc', '').strip()]
    if known_limits:
        if not np.allclose(known_limits, known_limits[0], rtol=0, atol=1e-12):
            raise ValueError('known_transition_xc must be consistent within a Cp series')
        if xc >= known_limits[0]:
            raise ValueError('trip is at/downstream of known transition')
    profile = build_surrogate_profile(float(bl['theta_m'][-1]), float(bl['lambda_thwaites'][-1]),
                                      float(edge['U_e_m_s'][-1]), fs['nu_inf_m2_s'])
    profile['edge_mach'] = float(edge['M_e'][-1])
    profile['ue_m_s'] = float(edge['U_e_m_s'][-1])
    profile['lambda_thwaites'] = float(bl['lambda_thwaites'][-1])
    profile['s_m'] = float(s_eval[-1])
    profile['distance_basis'] = 'supplied s_m' if explicit_s else 's=x: chordwise approximation'
    return profile


def supplied_profile(rows, xc, analysis):
    rows = [r for r in rows if abs(float(r['x_c']) - xc) < 1e-9]
    if not rows:
        raise ValueError('no supplied profile at requested x/c (spatial interpolation is not performed)')
    rows.sort(key=lambda r: float(r['y_m']))
    if any(r['laminar_confirmed'].strip().lower() not in ('true', '1', 'yes') for r in rows):
        raise ValueError('every supplied profile row must confirm undisturbed laminar flow')
    y, u, nu = (numeric(rows, c) for c in ('y_m', 'u_m_s', 'nu_m2_s'))
    ue = numeric(rows, 'ue_m_s')
    mach = numeric(rows, 'edge_mach')
    if not np.allclose(ue, ue[0], rtol=1e-8) or not np.allclose(mach, mach[0], rtol=1e-8):
        raise ValueError('ue_m_s and edge_mach must be constant within one profile')
    positive(ue[0], 'ue_m_s')
    positive(mach[0], 'edge_mach', allow_zero=True)
    if 'laminar_until_xc' in analysis and xc >= float(analysis['laminar_until_xc']):
        raise ValueError('trip is at/downstream of supplied laminar extent')
    metrics = profile_metrics(y, u, nu, float(ue[0]))
    return dict(metrics, y_m=y, u_m_s=u, nu_m2_s=nu, ue_m_s=float(ue[0]),
                edge_mach=float(mach[0]), lambda_thwaites=None, method='supplied laminar profile',
                distance_basis='profile at requested x/c')


def pick_height(required_m, hardware_mm):
    if not hardware_mm:
        return required_m
    return next((float(h) / 1000 for h in sorted(hardware_mm)
                 if float(h) / 1000 >= required_m * (1 - 1e-12)), None)


def select_common_height(profiles, required_m, hardware_mm, cap, target):
    """Smallest admissible common height, including nonmonotone Re_k from nu(y)."""
    upper = min(min(p['delta99_m'] * cap, p['y_m'][-1]) for p in profiles)
    if required_m > upper * (1 + 1e-12):
        return None
    def feasible(k):
        return k <= upper * (1 + 1e-12) and all(
            k <= p['y_m'][-1] and roughness_re(p, k) >= target * (1 - 1e-10) for p in profiles)
    if hardware_mm:
        candidates = sorted(float(h) / 1000 for h in hardware_mm
                            if required_m * (1 - 1e-12) <= float(h) / 1000 <= upper * (1 + 1e-12))
    elif feasible(required_m):
        return required_m
    else:
        # All changes in admissibility occur at roots of u(k)k-target*nu(k).
        # Both u and nu are linear on each profile interval; use local t in
        # [0,1] to keep the quadratic coefficients well scaled.
        candidates = [required_m, upper]
        for p in profiles:
            y, u, nu = (np.asarray(p[name]) for name in ('y_m', 'u_m_s', 'nu_m2_s'))
            for j in range(len(y) - 1):
                if y[j + 1] < required_m or y[j] > upper:
                    continue
                h, du, dn = y[j + 1] - y[j], u[j + 1] - u[j], nu[j + 1] - nu[j]
                a = du * h
                b = u[j] * h + du * y[j] - target * dn
                c = u[j] * y[j] - target * nu[j]
                roots = np.roots([a, b, c]) if a != 0 else ([-c / b] if b != 0 else [])
                for t in roots:
                    if np.isreal(t) and 0 <= float(np.real(t)) <= 1:
                        k = float(y[j] + float(np.real(t)) * h)
                        if required_m <= k <= upper:
                            candidates.append(k)
        candidates.sort()
    return next((k for k in candidates if feasible(k)), None)


def run_design(config_path, *, cp_file=None, profiles_file=None, position_xc=None):
    config_path = Path(config_path).resolve()
    cfg = json.loads(config_path.read_text(encoding='utf-8-sig'))
    if cp_file is not None and profiles_file is not None:
        raise ValueError('Choose either cp_file or profiles_file, not both')
    overrides = {}
    for source, supplied in (('cp', cp_file), ('profiles', profiles_file)):
        if supplied is not None:
            data_override = Path(supplied).expanduser().resolve()
            field = 'cp_file' if source == 'cp' else 'profiles_file'
            cfg['analysis']['source'] = source
            cfg['analysis'][field] = str(data_override)
            overrides.update(source=source, input_filename=data_override.name)
            cfg['name'] = f"User {source} input with configured wing geometry"
            cfg['data_origin'] = (f'User-selected {source} input: {data_override.name}. '
                                  'Geometry and tunnel conditions come from the configuration. '
                                  'Angle-of-attack labels and their limitations are retained from that configuration. '
                                  'The template provenance does not establish the origin or validity of this input.')
    if position_xc is not None:
        cfg['trip']['position_xc'] = {surface: float(position_xc) for surface in cfg['trip']['position_xc']}
        overrides['position_xc'] = float(position_xc)
        cfg['name'] = f"Trip sizing at {100*float(position_xc):g}% chord ({cfg['analysis']['source']} input)"
    calibration = validate_config(cfg)
    wing, trip, analysis = cfg['wing'], cfg['trip'], cfg['analysis']
    source = analysis['source']
    data_path = (config_path.parent / analysis['cp_file' if source == 'cp' else 'profiles_file']).resolve()
    columns = ['condition_id', 'eta', 'surface', 'x_c']
    columns += ['cp'] if source == 'cp' else ['y_m', 'u_m_s', 'nu_m2_s', 'ue_m_s', 'edge_mach', 'laminar_confirmed']
    data = read_csv(data_path, columns)
    groups = {}
    ids = {str(c['id']) for c in cfg['conditions']}
    etas = {round(float(s['eta']), 10) for s in wing['stations']}
    for row in data:
        k = key(row['condition_id'], row['eta'], row['surface'])
        if k[0] not in ids or k[1] not in etas or k[2] not in ('upper', 'lower'):
            raise ValueError(f'data row refers to unknown condition/station/surface: {k}')
        if not math.isfinite(float(row['x_c'])) or not 0 <= float(row['x_c']) <= 1:
            raise ValueError('data x_c must be finite and in [0,1]')
        groups.setdefault(k, []).append(row)
    margin = float(trip.get('height_margin', 1.1))
    cap = float(trip.get('max_k_delta99', .8))
    target = float(calibration['critical_re_k'])
    condition_rows, design_rows, profiles = [], [], {}
    for station in sorted(wing['stations'], key=lambda s: s['eta']):
        eta, chord = float(station['eta']), float(station['chord_m'])
        outboard_distance = eta * float(wing['semispan_m'])
        y_span = float(station.get('y_m', outboard_distance))
        x_le = float(station.get('x_le_m', outboard_distance * math.tan(math.radians(float(wing['sweep_le_deg'])))))
        if not math.isfinite(x_le):
            raise ValueError('x_le_m must be finite')
        for surface, xc_value in trip['position_xc'].items():
            xc = float(xc_value)
            local_rows = []
            for condition in cfg['conditions']:
                cid = str(condition['id'])
                record = dict(condition_id=cid, alpha_deg=condition['alpha_deg'],
                              mach=condition['mach'], reynolds=condition['reynolds'], eta=eta,
                              surface=surface, x_c=xc, status='invalid', reason='')
                if condition.get('alpha_label'):
                    record['alpha_label'] = str(condition['alpha_label'])
                try:
                    fs = tunnel_state(float(condition['mach']), float(condition['reynolds']),
                                      float(wing['reference_length_m']), float(condition.get('temperature_k', 288.15)))
                    re_x = fs['U_inf_m_s'] * xc * chord / fs['nu_inf_m2_s']
                    record.update(re_x_inf=re_x, tunnel_pressure_pa=fs['p_inf_pa'])
                    rows = groups.get(key(cid, eta, surface), [])
                    if not rows:
                        raise ValueError('missing data for condition/station/surface')
                    if re_x < float(calibration['min_re_x']):
                        raise ValueError('station Re_x is below calibration minimum; downstream position or calibrated criterion needed')
                    p = (cp_profile(rows, chord, xc, fs, condition, analysis) if source == 'cp'
                         else supplied_profile(rows, xc, analysis))
                    if p['edge_mach'] > float(calibration['max_edge_mach']):
                        raise ValueError('local edge Mach exceeds trip calibration range')
                    critical = required_height(p, target, float(p['y_m'][-1]))
                    record.update(delta99_mm=float(p['delta99_m']) * 1000,
                                  theta_mm=float(p['theta_m']) * 1000, ue_m_s=p['ue_m_s'],
                                  edge_mach=p['edge_mach'], lambda_thwaites=p['lambda_thwaites'],
                                  profile_method=p['method'], distance_basis=p['distance_basis'])
                    if critical is None:
                        record.update(status='infeasible', reason='target Re_k exceeds supplied profile height coverage')
                    else:
                        margin_height = critical * margin
                        required = select_common_height([p], margin_height, [], cap, target)
                        record.update(k_critical_mm=critical * 1000, k_margin_mm=margin_height * 1000,
                                      k_required_mm=required * 1000 if required is not None else margin_height * 1000,
                                      status='feasible' if required is not None else 'infeasible',
                                      reason='' if required is not None else 'no height above margin meets Re_k within k/delta99 cap')
                    profiles[key(cid, eta, surface)] = p
                except (ValueError, KeyError, ArithmeticError) as exc:
                    record.update(status='invalid', reason=str(exc))
                local_rows.append(record)
                condition_rows.append(record)
            design = dict(eta=eta, surface=surface, x_c=xc, chord_m=chord,
                          y_span_m=y_span, x_global_m=x_le + chord * xc,
                          status='infeasible', reason='', height_mm=None)
            bad = [r for r in local_rows if r['status'] != 'feasible']
            if bad:
                design['status'] = 'invalid' if any(r['status'] == 'invalid' for r in bad) else 'infeasible'
                design['reason'] = '; '.join(f"{r['condition_id']}: {r['reason']}" for r in bad)
            else:
                governing = max(local_rows, key=lambda r: r['k_required_mm'])
                local_profiles = [profiles[key(r['condition_id'], eta, surface)] for r in local_rows]
                selected = select_common_height(local_profiles, governing['k_required_mm'] / 1000,
                                                trip.get('available_heights_mm', []), cap, target)
                if selected is None:
                    design['reason'] = 'no common height satisfies requirement, k/delta99 cap, profile coverage and hardware choices'
                else:
                    ratios = []
                    rek = []
                    for r in local_rows:
                        p = profiles[key(r['condition_id'], eta, surface)]
                        ratios.append(selected / p['delta99_m'])
                        # Re-evaluate the actual common height, including hardware rounding.
                        rek.append(roughness_re(p, selected) if selected <= p['y_m'][-1] else None)
                    if max(ratios) > cap * (1 + 1e-10) or any(v is None or v < target * (1 - 1e-9) for v in rek):
                        design['reason'] = 'common/rounded height violates a condition cap or profile coverage'
                    else:
                        design.update(status='feasible', height_mm=selected * 1000,
                                      governing_condition=governing['condition_id'],
                                      max_k_delta99=max(ratios), min_re_k=min(rek))
                        for r, ratio, reynolds in zip(local_rows, ratios, rek):
                            r.update(selected_height_mm=selected * 1000, selected_k_delta99=ratio,
                                     selected_re_k=reynolds)
            design_rows.append(design)
    all_feasible = all(r['status'] == 'feasible' for r in design_rows)
    uniform = dict(status='unavailable', height_mm=None)
    if all_feasible:
        uniform_m = select_common_height(list(profiles.values()), max(r['height_mm'] for r in design_rows) / 1000,
                                         trip.get('available_heights_mm', []), cap, target)
        uniform_ok = uniform_m is not None
        uniform = dict(status='feasible' if uniform_ok else 'infeasible',
                       height_mm=uniform_m * 1000 if uniform_ok else None,
                       reason='' if uniform_ok else 'one height for the entire wing violates a local cap')
    warnings = [
        'Preliminary trip sizing, not a prediction or guarantee of transition or downstream turbulent development.',
        'Only the supplied discrete conditions are covered; no guarantee between sampled angles of attack or Reynolds numbers.',
        'Sweep sets placement geometry only. Crossflow, attachment-line transition, and 3D streamline histories are not solved.',
        'Height is exposed protrusion above the model surface, including mounting effects; verify actual installed height.',
        'Height margin and k/delta99 cap are user design policies, not universal transition correlations.',
    ]
    warnings.extend(f"Condition {condition['id']}: {condition['alpha_label']}"
                    for condition in cfg['conditions'] if condition.get('alpha_label'))
    if source == 'cp':
        warnings += ['Cp mode assumes the upstream boundary layer remains laminar unless a known limit is supplied.',
                     'Cp mode uses incompressible Thwaites and a reconstructed profile; freestream and upstream local Mach must stay <=0.30.',
                     'Angle of attack labels select supplied Cp cases; they do not generate aerodynamic loading.']
    else:
        warnings += ['Profile mode trusts the supplied laminar confirmation and local kinematic viscosity.',
                     'Scalar profile velocity must follow the direction applicable to the calibrated trip criterion; crossflow is not modeled.']
    summary = dict(name=cfg.get('name', config_path.stem), data_origin=cfg.get('data_origin', 'user supplied'),
                   source=source, status='feasible' if all_feasible else 'not_feasible', overrides=overrides,
                   calibration=calibration, height_margin=margin, max_k_delta99=cap,
                   design=design_rows, uniform_height=uniform, warnings=warnings,
                   inputs=[dict(path=p.relative_to(config_path.parent).as_posix()
                                if p.is_relative_to(config_path.parent) else p.name,
                                path_base='configuration directory' if p.is_relative_to(config_path.parent)
                                else 'filename only (outside configuration directory)',
                                sha256=hashlib.sha256(p.read_bytes()).hexdigest())
                           for p in (config_path, data_path)])
    return cfg, summary, condition_rows


def write_csv(path, rows):
    fields = list(dict.fromkeys(k for row in rows for k in row))
    with path.open('w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def make_plot(output, cfg, summary, rows):
    # Keep matplotlib's first-run font cache inside the writable output folder.
    os.environ.setdefault('MPLCONFIGDIR', str((output / '.plot_cache').resolve()))
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        return False
    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False})
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.7), layout='constrained')
    fig.suptitle(f"{summary['name']} | fixed-position trip sizing", fontsize=16, fontweight='bold')
    colors = {'upper': '#146a8a', 'lower': '#c66527'}
    wing = cfg['wing']
    stations = sorted(wing['stations'], key=lambda s: s['eta'])
    outboard = np.array([float(s['eta']) * float(wing['semispan_m']) for s in stations])
    span = np.array([float(s.get('y_m', distance)) for s, distance in zip(stations, outboard)])
    le = np.array([float(s.get('x_le_m', distance * math.tan(math.radians(float(wing['sweep_le_deg'])))))
                   for s, distance in zip(stations, outboard)])
    chord = np.array([s['chord_m'] for s in stations])
    axes[0].fill_between(span, le, le + chord, color='#e6ebef')
    for surface, color in colors.items():
        design = [r for r in summary['design'] if r['surface'] == surface]
        if not design:
            continue
        axes[0].plot([r['y_span_m'] for r in design], [r['x_global_m'] for r in design],
                     'o-' if surface == 'upper' else 'x--', color=color, label=surface)
        good = [r for r in design if r['height_mm'] is not None]
        if good:
            axes[1].plot([r['eta'] for r in good], [r['height_mm'] for r in good], 'o-', color=color, label=surface)
        for cid in [c['id'] for c in cfg['conditions']]:
            rr = sorted([r for r in rows if r['surface'] == surface and r['condition_id'] == cid and
                         'k_required_mm' in r], key=lambda r: r['eta'])
            if rr:
                axes[1].plot([r['eta'] for r in rr], [r['k_required_mm'] for r in rr], color=color, alpha=.2)
        rr = [r for r in rows if r['surface'] == surface and 'selected_k_delta99' in r]
        if rr:
            axes[2].scatter([r['eta'] for r in rr], [r['selected_k_delta99'] for r in rr], color=color, s=18, alpha=.65, label=surface)
    axes[0].invert_yaxis()
    axes[0].set(xlabel='Spanwise position y (m)', ylabel='Global chordwise x (m)', title='Trip location at sampled stations')
    axes[1].set(xlabel='Semispan fraction', ylabel='Exposed height k (mm)', title='Selected heights; faint lines: case requirements')
    axes[2].axhline(summary['max_k_delta99'], color='#ad2639', linestyle='--', label='User cap')
    axes[2].set(xlabel='Semispan fraction', ylabel='k / delta99', title='All sampled conditions')
    for ax in axes:
        ax.grid(alpha=.15)
        if ax.get_legend_handles_labels()[0]:
            ax.legend(fontsize=8)
    fig.text(.5, -.015, 'Preliminary sizing | supplied discrete conditions only | see report for input provenance', ha='center', fontsize=9)
    fig.savefig(output / 'overview.png', dpi=170, bbox_inches='tight')
    fig.savefig(output / 'overview.svg', bbox_inches='tight')
    plt.close(fig)
    return True


def write_outputs(output, cfg, summary, rows, plot=True):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    write_csv(output / 'design.csv', summary['design'])
    write_csv(output / 'conditions.csv', rows)
    plotted = make_plot(output, cfg, summary, rows) if plot else False
    summary['plot_created'] = plotted
    (output / 'summary.json').write_text(json.dumps(summary, indent=2, allow_nan=False), encoding='utf-8')
    def esc(value):
        return html.escape(str(value))
    headers = ('surface', 'eta', 'x_c', 'height_mm', 'status', 'reason')
    table = '<tr>' + ''.join(f'<th>{esc(h)}</th>' for h in headers) + '</tr>'
    for row in summary['design']:
        table += '<tr>' + ''.join(f'<td>{esc(row.get(h, ""))}</td>' for h in headers) + '</tr>'
    body = f'''<!doctype html><html lang="en"><meta charset="utf-8"><title>{esc(summary['name'])}</title>
    <style>body{{font:16px/1.5 system-ui;max-width:1150px;margin:40px auto;padding:0 20px;color:#173447}}
    table{{border-collapse:collapse;width:100%;font-size:14px}}td,th{{padding:9px;border-bottom:1px solid #ccd6dc;text-align:left}}
    img{{width:100%}}.note{{background:#eef4f6;padding:16px;border-radius:8px}}code{{background:#eef4f6}}</style>
    <h1>{esc(summary['name'])}</h1><p class="note"><strong>Status: {esc(summary['status'])}.</strong>
    {esc(summary['data_origin'])}</p><p>Source: {esc(summary['source'])}; critical Re_k: {summary['calibration']['critical_re_k']};
    height multiplier: {summary['height_margin']}; k/delta99 cap: {summary['max_k_delta99']}.</p>
    <p>Calibration: {esc(summary['calibration']['notes'])}</p>
    {('<img src="overview.png" alt="Trip position, height and boundary-layer ratio plots">' if plotted else '')}
    <h2>Common height across supplied conditions at each station</h2><table>{table}</table>
    <p>Whole-wing uniform-height alternative: {esc(summary['uniform_height'])}</p>
    <h2>Interpretation</h2><ul>{''.join('<li>'+esc(w)+'</li>' for w in summary['warnings'])}</ul>
    <p>See <a href="conditions.csv">per-condition results</a>, <a href="design.csv">design table</a>,
    and <a href="summary.json">full results and input hashes</a>.</p></html>'''
    (output / 'report.html').write_text(body, encoding='utf-8')


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--output', type=Path, default=Path('results'))
    inputs = parser.add_mutually_exclusive_group()
    inputs.add_argument('--cp-file', type=Path, help='Use this Cp CSV (path relative to current directory); configuration is not modified')
    inputs.add_argument('--profiles-file', type=Path, help='Use this laminar-profile CSV instead of the configured input')
    parser.add_argument('--position-xc', type=float, help='Override x/c for all configured surfaces, e.g. 0.10')
    parser.add_argument('--no-plot', action='store_true', help='Skip optional matplotlib figures')
    args = parser.parse_args(argv)
    try:
        cfg, summary, rows = run_design(args.config, cp_file=args.cp_file,
                                        profiles_file=args.profiles_file, position_xc=args.position_xc)
        write_outputs(args.output, cfg, summary, rows, not args.no_plot)
    except (ValueError, KeyError, OSError, TypeError) as exc:
        print(f'Input error: {exc}', file=sys.stderr)
        return 2
    n = sum(r['status'] == 'feasible' for r in summary['design'])
    print(f"{n}/{len(summary['design'])} station/surface designs feasible across all supplied conditions.")
    uniform = summary['uniform_height']
    if uniform['status'] == 'feasible':
        print(f"Uniform-height option across all stations: {uniform['height_mm']:.6f} mm")
    else:
        print(f"Uniform-height option: {uniform['status']}. {uniform.get('reason', '')}")
    print(f"Report: {(args.output / 'report.html').resolve()}")
    if not summary['plot_created'] and not args.no_plot:
        print('Charts not created: install matplotlib; tables and report are complete.')
    return 0 if summary['status'] == 'feasible' else 2


if __name__ == '__main__':
    raise SystemExit(main())
