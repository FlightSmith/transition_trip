"""Reproduce the NACA 64A204 source airfoil used by the left-wing example.

The algorithm and mapping arrays are from Ralph Carmichael's public-domain
PDAS NACA456 program: https://www.pdas.com/packages/naca456.zip . This is a
float64 Python translation of ScaleFactor, SetSixDigitPoints, MeanLine6M,
CombineThicknessAndCamber, and (for verification only) FMMspline.

The bundled source_data.json preserves EPS7/PSI7, the published direct table,
and source provenance. No network access, Fortran compiler, or .runtime files
are needed. The 401 output points use the original NACA chord, including the
tiny negative-x nose point; no incidence rotation or geometry fitting occurs.
Only the airfoil and its source metadata are generated, not a solver workflow.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


DEFAULT_OUTPUT = Path(__file__).resolve().parent / 'examples' / 'left_wing_naca64a204_scale8'
SOURCE_DATA = DEFAULT_OUTPUT / 'source_data.json'
SOURCE_ZIP_SHA256 = '327f9dae5552803990d0b683dd0e1d6c777562e4e3695b50d8cee80eb412f2c9'


def load_source_data(path=SOURCE_DATA):
    """Read and validate the fixed, attributed source for this exact airfoil."""
    source = json.loads(Path(path).read_text(encoding='utf-8'))
    if (source.get('airfoil') != 'NACA 64A204' or source.get('format_version') != 1
            or source.get('source_zip_sha256') != SOURCE_ZIP_SHA256
            or source.get('parameters') != {'family': 7, 'nominal_thickness_chord': .04,
                'design_lift_coefficient': .2, 'modified_meanline_a': .8}):
        raise ValueError('Unexpected NACA 64A204 source identity or parameters')
    for key in ('eps7', 'psi7'):
        array = np.asarray(source.get(key), dtype=float)
        if array.shape != (201,) or not np.all(np.isfinite(array)):
            raise ValueError(f'{key} must contain 201 finite mapping values')
    table = np.asarray(source.get('published_direct_table'), dtype=float)
    if (table.shape != (26, 4) or not np.all(np.isfinite(table))
            or source.get('published_direct_table_units') != 'percent chord'):
        raise ValueError('Expected the 26-row PDAS direct table in percent chord')
    return source


def native_thickness(source_data=None):
    """Return the 201 native 64A thickness points at nominal t/c=.04."""
    source_data = load_source_data() if source_data is None else source_data
    tc = .04
    # PDAS family 7 = 64A, ScaleFactor and SetSixDigitPoints.
    coeff = [0., 8.2125018, .76855961, 1.4922345, 3.6130133]
    sf = np.polynomial.polynomial.polyval(tc, coeff)
    eps = sf * np.asarray(source_data['eps7'], dtype=float)
    psi = sf * np.asarray(source_data['psi7'], dtype=float)
    phi = (3.14159265 / 200.) * np.arange(201)
    z = np.exp(psi[0] + 1j * phi)
    zp = z * np.exp(psi - psi[0] - 1j * eps)
    zeta = zp + 1. / zp
    final = (zeta[0] - zeta) / abs(zeta[-1] - zeta[0])
    x, yt = final.real, -final.imag
    # Remove tiny nonzero endpoint due to PDAS's rounded pi constant.
    x[0], x[-1], yt[0], yt[-1] = 0., 1., 0., 0.
    assert np.all(np.diff(x) > 0) and np.all(yt >= 0)
    return x, yt


def meanline(x, cl=.2):
    # PDAS MeanLine6M: modified 6-series mean line, fixed a=.8.
    a, oma, pi = .8, .2, 3.141592654
    y, yp = np.zeros_like(x), np.zeros_like(x)
    g = -(a*a*(.5*np.log(a)-.25)+.25)/oma
    h = g + (.5*oma*oma*np.log(oma)-.25*oma*oma)/oma
    for i, xx in enumerate(x):
        if xx < 1e-7:
            continue
        if 1 - xx < 1e-7:
            yp[i] = -.24521 * cl
            continue
        amx, omx = a - xx, 1 - xx
        t1 = 0. if abs(amx) < 1e-7 else amx*amx*(2*np.log(abs(amx))-1)
        t1p = 0. if abs(amx) < 1e-7 else -amx*np.log(abs(amx))
        t2 = omx*omx*(1-2*np.log(omx))
        t2p = omx*np.log(omx)
        y[i] = .25*(t1+t2)/oma - xx*np.log(xx)+g-h*xx
        yp[i] = (t1p+t2p)/oma - 1-np.log(xx)-h
    factor = cl*.97948/(2*pi*(a+1))
    y, yp = y*factor, yp*factor
    aft = x > .86
    y[aft], yp[aft] = (-.24521*cl)*(x[aft]-1), -.24521*cl
    return y, yp


def combine(x, yt):
    y, yp = meanline(x)
    s, c = np.sin(np.arctan(yp)), np.cos(np.arctan(yp))
    return np.column_stack([x-yt*s, y+yt*c, x+yt*s, y-yt*c])


def fmm(x, y):
    # Literal index translation of PDAS splprocs.f90 FMMspline, n>3.
    n = len(x)
    dx, delta = np.diff(x), np.diff(y)/np.diff(x)
    dd = np.diff(delta)
    alpha, beta, sigma = np.zeros(n), np.zeros(n), np.zeros(n)
    alpha[0] = -dx[0]
    alpha[1:-1] = 2*(dx[:-1]+dx[1:])
    for i in range(1,n-1):
        alpha[i] -= dx[i-1]**2/alpha[i-1]
    alpha[-1] = -dx[-1]-dx[-1]**2/alpha[-2]
    beta[0] = (dd[1]/(x[3]-x[1])-dd[0]/(x[2]-x[0]))*dx[0]**2/(x[3]-x[0])
    beta[1:-1] = dd
    beta[-1] = -(dd[-1]/(x[-1]-x[-3])-dd[-2]/(x[-2]-x[-4]))*dx[-1]**2/(x[-1]-x[-4])
    for i in range(1,n):
        beta[i] -= dx[i-1]*beta[i-1]/alpha[i-1]
    sigma[-1] = beta[-1]/alpha[-1]
    for i in range(n-2,-1,-1):
        sigma[i] = (beta[i]-dx[i]*sigma[i+1])/alpha[i]
    yp = np.zeros(n)
    yp[:-1] = delta-dx*(2*sigma[:-1]+sigma[1:])
    yp[-1] = yp[-2]+3*dx[-1]*(sigma[-1]+sigma[-2])
    return yp


def hermite(x, y, yp, q):
    k = int(np.clip(np.searchsorted(x,q)-1, 0, len(x)-2))
    h = x[k+1]-x[k]
    t = (q-x[k])/h
    return (2*t**3-3*t*t+1)*y[k]+(t**3-2*t*t+t)*h*yp[k]+(-2*t**3+3*t*t)*y[k+1]+(t**3-t*t)*h*yp[k+1]


def thickness_at(x, yt, q):
    contour = np.vstack([np.column_stack([x[::-1],yt[::-1]]), np.column_stack([x[1:],-yt[1:]])])
    s = np.r_[0.,np.cumsum(np.linalg.norm(np.diff(contour,axis=0),axis=1))]
    dx, dy = fmm(s,contour[:,0]), fmm(s,contour[:,1])
    ans = []
    for xx in q:
        lo, hi = s[200], s[-1]
        for _ in range(60):
            mid = (lo+hi)/2
            if hermite(s,contour[:,0],dx,mid) < xx: lo = mid
            else: hi = mid
        ans.append(-hermite(s,contour[:,1],dy,(lo+hi)/2))
    return np.array(ans)



def coordinates_and_provenance(source_data=None):
    """Return 401 unit-chord contour points after checking the published table.

    Order is upper TE to nominal LE to lower TE, including both coincident
    trailing-edge endpoints. FMM interpolation is used only for verification;
    output geometry comes directly from the native conformal mapping points.
    """
    source = load_source_data() if source_data is None else source_data
    x, yt = native_thickness(source)
    combined = combine(x, yt)
    points = np.vstack([combined[::-1, :2], combined[1:, 2:]])
    published = np.asarray(source['published_direct_table'], dtype=float) / 100.
    query = (published[:, 0] + published[:, 2]) / 2.
    calculated = combine(query, thickness_at(x, yt, query))
    error = abs(published - calculated)
    maximum_error = float(error.max())
    # Printed ordinates have 1e-6 chord resolution. A measured 5.164e-7
    # maximum also includes the original source's interpolation tolerance.
    if not np.isfinite(maximum_error) or maximum_error >= 1e-6:
        raise ValueError(f'PDAS published-coordinate check failed: {maximum_error:.6g} chord')
    if (points.shape != (401, 2) or not np.all(np.isfinite(points))
            or not np.array_equal(points[0], [1., 0.])
            or not np.array_equal(points[-1], [1., 0.])):
        raise ValueError('Generated contour must be 401 finite points with a closed TE')
    provenance = {
        'airfoil': 'NACA 64A204',
        'attribution': source['attribution'],
        'source_url': source['source_zip_url'],
        'source_zip_sha256': source['source_zip_sha256'],
        'published_table_url': source['table_url'],
        'source_html_sha256': source['source_html_sha256'],
        'source_data_file': 'source_data.json',
        'source_description': 'PDAS family 7 (64A); nominal thickness .04; modified meanline cl=.2, a=.8.',
        'method': 'Float64 translation of the PDAS native 201-point conformal map and modified mean line; no coarse-table interpolation for output geometry.',
        'point_count': len(points),
        'nominal_tc': .04,
        'computed_max_native_normal_thickness': float(2 * yt.max()),
        'published_table_comparison_max_abs_chord': maximum_error,
        'published_table_comparison_by_column_xu_yu_xl_yl': error.max(axis=0).tolist(),
        'validation': 'Every coordinate agrees with the published direct table within 1e-6 chord using the original FMM thickness interpolation for the check.',
        'closed_te': True,
        'min_x_chord': float(points[:, 0].min()),
        'coordinate_order': 'upper TE to nominal LE to lower TE; shared nominal LE once',
        'incidence_reference': 'Original NACA chord from (0,0) to (1,0). Minimum-x surface point is not the nominal chord origin.',
        'output_file': 'naca64a204.dat',
        'limitation': 'Source-code translation, not execution of the Fortran compiler. Geometric coordinate agreement does not validate XFOIL aerodynamics.',
    }
    return points, provenance


def coordinate_bytes(points):
    """Stable ASCII/LF format matching the source airfoil used in the sweep."""
    content = 'NACA 64A204 PDAS native conformal map 401 points\n'
    content += ''.join(f'{x:.12f} {y:.12f}\n' for x, y in points)
    return content.encode('ascii')


def _write_if_changed(path, data):
    if not path.exists() or path.read_bytes() != data:
        path.write_bytes(data)


def generate(output=DEFAULT_OUTPUT):
    """Write only naca64a204.dat, source_data.json, and airfoil_source.json.

    ``output`` is a destination directory. Existing identical files are left
    untouched. The active model geometry, section files, and workflow settings
    are outside the scope of this helper and are never written.
    """
    source = load_source_data()
    points, provenance = coordinates_and_provenance(source)
    data = coordinate_bytes(points)
    provenance['output_sha256'] = hashlib.sha256(data).hexdigest()
    destination = Path(output).expanduser().resolve()
    destination.mkdir(parents=True, exist_ok=True)
    _write_if_changed(destination / 'naca64a204.dat', data)
    _write_if_changed(destination / 'source_data.json', SOURCE_DATA.read_bytes())
    metadata = (json.dumps(provenance, indent=2) + '\n').encode('utf-8')
    _write_if_changed(destination / 'airfoil_source.json', metadata)
    return provenance


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, default=DEFAULT_OUTPUT,
                        help='destination directory for source airfoil and provenance')
    args = parser.parse_args(argv)
    result = generate(args.out)
    print(f"Generated {result['output_file']}: {result['point_count']} points; "
          f"published-coordinate error {result['published_table_comparison_max_abs_chord']:.3g} chord")


if __name__ == '__main__':
    main()
