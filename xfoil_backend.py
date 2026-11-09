"""Optional, isolated XFOIL execution and conservative laminar-station sampling.

XFOIL is a separate executable, not a Python dependency.  Its DUMP output
contains integral boundary-layer quantities, NOT resolved wall-normal u(y).
All lengths below are normalized by the unit input chord.  The official 6.99
Windows binary reports mass H, while the downloadable 6.99 source writes Hk
in the DUMP column labelled ``H``. Both are recognized from the integral data
and the source's HKIN relation; the interpretation is retained in metadata.

References: https://web.mit.edu/drela/Public/web/xfoil/xfoil_doc.txt and the
official xfoil6.99.tgz, src/xoper.f (BLDUMP, CPDUMP and ALFA/PACC), plutil.f.
The source only adds a viscous PACC point when LVCONV is true.  A matching
polar row is required here; a process exit code or a DUMP file is insufficient.
"""
from __future__ import annotations

import hashlib
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile

import numpy as np


def _finite(value, name, *, positive=False):
    value = float(value)
    if not math.isfinite(value) or (positive and value <= 0):
        raise ValueError(f"{name} must be finite" + (" and positive" if positive else ""))
    return value


def _number(token):
    return float(token.replace("D", "E").replace("d", "e"))


def _table(path, minimum_columns, *, allowed_columns=None):
    """Read numeric text; reject malformed rows once the data have started."""
    rows = []
    for line_number, line in enumerate(Path(path).read_text(encoding="utf-8-sig", errors="replace").splitlines(), 1):
        line = line.strip()
        if not line or line.startswith(("#", "!")):
            continue
        fields = line.replace(",", " ").split()
        try:
            _number(fields[0])
        except ValueError:
            if rows or any(character.isdigit() for character in fields[0]):
                raise ValueError(f"{Path(path).name}:{line_number}: malformed numeric row")
            continue
        if len(fields) < minimum_columns or (allowed_columns and len(fields) not in allowed_columns):
            raise ValueError(f"{Path(path).name}:{line_number}: unexpected column count")
        try:
            values = [_number(token) for token in fields[:minimum_columns]]
        except ValueError as exc:
            raise ValueError(f"{Path(path).name}:{line_number}: malformed numeric value") from exc
        if not np.all(np.isfinite(values)):
            raise ValueError(f"{Path(path).name}:{line_number}: nonfinite numeric value")
        rows.append((values, fields))
    if not rows:
        raise ValueError(f"{Path(path).name}: no numeric data")
    return rows


def parse_cp(path):
    """Read the 2-column CPWR format or a 3-column x,y,Cp variant."""
    rows = _table(path, 2, allowed_columns={2, 3})
    width = len(rows[0][1])
    if any(len(fields) != width for _, fields in rows):
        raise ValueError("Cp file mixes two- and three-column formats")
    values = np.array([[_number(v) for v in fields] for _, fields in rows])
    if not np.all(np.isfinite(values)) or len(values) < 5:
        raise ValueError("Cp file requires at least five finite rows")
    return {"x_c": values[:, 0], "y_c": values[:, 1] if width == 3 else None,
            "cp": values[:, -1]}


def parse_dump(path, mach=0.):
    """Read first eight DUMP columns; extra columns and shorter wake rows coexist."""
    values = np.array([row for row, _ in _table(path, 8)])
    # The first wake point can share the TE's s after five-decimal printing.
    if len(values) < 5 or np.any(np.diff(values[:, 0]) < 0):
        raise ValueError("DUMP requires at least five rows with nondecreasing contour distance")
    if np.any(values[:, 4:6] < 0) or np.any(values[:, 7] <= 0):
        raise ValueError("DUMP has negative thickness or nonpositive kinematic shape factor")
    keys = ("s_c", "x_c", "y_c", "ue_over_uinf_signed", "delta_star_c",
            "theta_c", "cf", "h_reported")
    result = dict(zip(keys, values.T))
    result["h_mass"] = np.divide(values[:, 4], values[:, 5],
                                   out=np.full(len(values), np.nan), where=values[:, 5] > 0)
    mach = _finite(mach, "mach")
    speed2 = values[:, 3]**2
    temperature_ratio = 1 + .2 * mach**2 * (1 - speed2)
    if np.any(temperature_ratio <= 0):
        raise ValueError("DUMP velocities imply invalid local temperature")
    msq = mach**2 * speed2 / temperature_ratio
    reported = values[:, 7]
    # Compare in thickness units, accounting for DUMP's printed precision.
    # Thick downstream points distinguish H versus Hk without amplifying the
    # six-decimal Theta roundoff at a very thin leading-edge boundary layer.
    good = values[:, 5] > max(0., .2 * float(np.max(values[:, 5])))
    if not np.any(good):
        raise ValueError("DUMP contains no resolved viscous thickness")
    as_mass = reported
    as_kinematic_mass = reported * (1 + .113 * msq) + .29 * msq
    tolerance = .5e-6 * (1 + as_kinematic_mass) + .0005 * values[:, 5] + 1e-10
    def score(h):
        return float(np.median(np.abs(values[good, 4] - values[good, 5] * h[good]) / tolerance[good]))
    mass_score, kinematic_score = score(as_mass), score(as_kinematic_mass)
    if min(mass_score, kinematic_score) > 3:
        raise ValueError("DUMP H column is inconsistent with both supported H/Hk formats")
    if mach == 0:
        result["h_kinematic"] = reported.copy()
        interpretation = "H and Hk coincide at Mach zero"
    elif mass_score <= kinematic_score:
        result["h_kinematic"] = (reported - .29 * msq) / (1 + .113 * msq)
        interpretation = "mass H converted to Hk using XFOIL HKIN relation"
    else:
        result["h_kinematic"] = reported.copy()
        interpretation = "kinematic Hk as written by newer source"
    result["h_column_interpretation"] = interpretation
    return result


def parse_polar(path, alpha_deg):
    """Require exactly one finite matching converged point in the fresh polar."""
    lines = Path(path).read_text(encoding="utf-8-sig", errors="replace").splitlines()
    headers = None
    points = []
    conditions = None
    number = r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eEdD][-+]?\d+)?"
    for line in lines:
        match = re.search(rf"\bMach\s*=\s*({number})\s+Re\s*=\s*({number})"
                          rf"(?:\s+[eE]\s*([-+]?\d+))?\s+Ncrit\s*=\s*({number})", line)
        if match:
            mantissa = match[2]
            exponent = int(match[3] or 0)
            embedded = re.split("[eEdD]", mantissa)
            if len(embedded) > 1:
                exponent += int(embedded[1])
                mantissa = embedded[0]
            decimals = len(mantissa.split(".")[1]) if "." in mantissa else 0
            conditions = {"mach": _number(match[1]), "reynolds": _number(mantissa) * 10.**exponent,
                          "ncrit": _number(match[4]),
                          "reynolds_print_tolerance": .51 * 10.**(exponent - decimals)}
        fields = line.split()
        lowered = [f.lower() for f in fields]
        if "alpha" in lowered and "top_xtr" in lowered and "bot_xtr" in lowered:
            headers = lowered
            continue
        if headers is None or not fields or set(line.strip()) <= {"-", " "}:
            continue
        try:
            values = [_number(f) for f in fields]
        except ValueError as exc:
            raise ValueError("Malformed polar data after column header") from exc
        if len(values) != len(headers) or not np.all(np.isfinite(values)):
            raise ValueError("Malformed or nonfinite polar row")
        points.append(dict(zip(headers, values)))
    if headers is None or len(points) != 1:
        raise ValueError("No unique converged polar point; XFOIL may not have converged")
    if conditions is None or not all(math.isfinite(v) for v in conditions.values()):
        raise ValueError("Missing or invalid Mach/Re/Ncrit metadata in polar")
    point = points[0]
    if abs(point["alpha"] - float(alpha_deg)) > 0.00051:
        raise ValueError("Polar angle does not match requested angle")
    if not all(0 <= point[k] <= 1 for k in ("top_xtr", "bot_xtr")):
        raise ValueError("Polar transition positions must lie in [0,1]")
    return dict(point, **conditions)


def _normalized_coordinates(path):
    coordinates = []
    for line in Path(path).read_text(encoding="utf-8-sig").splitlines():
        fields = line.strip().split()
        if not fields or fields[0].startswith(("#", "!")):
            continue
        try:
            pair = [_number(v) for v in fields]
        except ValueError:
            if coordinates:
                raise ValueError("Airfoil has nonnumeric content after coordinate rows")
            continue
        if len(pair) != 2 or not np.all(np.isfinite(pair)):
            raise ValueError("Airfoil needs finite x/c,y/c pairs")
        coordinates.append(pair)
    xy = np.asarray(coordinates)
    if len(xy) < 10:
        raise ValueError("Airfoil requires at least ten normalized contour points")
    x, y = xy.T
    lead = int(np.argmin(x))
    if not (abs(float(x.min())) <= 1e-5 and abs(float(x.max()) - 1) <= 1e-5
            and x[0] > .95 and x[-1] > .95 and 2 <= lead <= len(x) - 3):
        raise ValueError("Airfoil must be unit chord, LE at x=0, ordered TE-LE-TE")
    if np.mean(y[:lead]) < np.mean(y[lead + 1:]):
        xy = xy[::-1]
    return xy


def xfoil_commands(reynolds_local, mach, alpha_deg, ncrit=9., iterations=200, panels=200):
    """Build a one-point, free-transition script with graphics disabled."""
    return "\n".join([
        "PLOP", "G", "", "LOAD foil.dat", "PPAR", f"N {panels}", "", "", "PANE",
        "OPER", "TYPE 1", f"VISC {reynolds_local:.12g}", f"MACH {mach:.12g}",
        "VPAR", f"N {ncrit:.12g}", "XTR 1 1", "", f"ITER {iterations}",
        "PACC", "polar.txt", "", f"ALFA {alpha_deg:.12g}",
        "CPWR cp.txt", "DUMP bl.txt", "PACC", "", "QUIT", "",
    ])


def run_xfoil(executable, airfoil_path, workdir, reynolds_local, mach, alpha_deg,
              ncrit=9., timeout_s=60., iterations=200, panels=200):
    """Run one section in a new retained subdirectory and return structured data.

    Invalid inputs raise ValueError. Execution/solver/output failures return
    ``converged=False`` and ``reason``; logs and commands remain available.
    XFOIL is never forced turbulent at the candidate trip. XTR=1,1 merely puts
    imposed trips at the trailing edge while the natural e^N model stays active.
    Mach is deliberately limited to 0.30 for this tool's low-speed surrogate.
    """
    reynolds_local = _finite(reynolds_local, "reynolds_local", positive=True)
    mach = _finite(mach, "mach")
    alpha_deg = _finite(alpha_deg, "alpha_deg")
    ncrit = _finite(ncrit, "ncrit", positive=True)
    timeout_s = _finite(timeout_s, "timeout_s", positive=True)
    if not 0 <= mach <= .3 or abs(alpha_deg) >= 90:
        raise ValueError("Runner requires 0 <= Mach <= 0.30 and |alpha| < 90 degrees")
    if isinstance(iterations, bool) or int(iterations) != iterations or not 1 <= iterations <= 10000:
        raise ValueError("iterations must be an integer from 1 to 10000")
    if isinstance(panels, bool) or int(panels) != panels or not 40 <= panels <= 280:
        raise ValueError("panels must be an integer from 40 to 280")
    executable = Path(shutil.which(str(executable)) or executable).expanduser().resolve()
    if not executable.is_file():
        raise ValueError(f"XFOIL executable not found: {executable}")
    coordinates = _normalized_coordinates(airfoil_path)
    parent = Path(workdir).expanduser().resolve()
    parent.mkdir(parents=True, exist_ok=True)
    directory = Path(tempfile.mkdtemp(prefix="xfoil_", dir=parent))
    (directory / "foil.dat").write_text("_trip_section\n" + "".join(
        f"{x:.12g} {y:.12g}\n" for x, y in coordinates), encoding="ascii")
    commands = xfoil_commands(reynolds_local, mach, alpha_deg, ncrit, int(iterations), int(panels))
    (directory / "commands.txt").write_text(commands, encoding="ascii")
    result = {"converged": False, "reason": "", "workdir": str(directory),
              "cp": None, "boundary_layer": None, "polar": None,
              "reynolds_local": reynolds_local, "mach": mach, "alpha_deg": alpha_deg,
              "ncrit": ncrit, "panels": int(panels), "iterations": int(iterations),
              "executable_sha256": hashlib.sha256(executable.read_bytes()).hexdigest(),
              "airfoil_sha256": hashlib.sha256(Path(airfoil_path).read_bytes()).hexdigest(),
              "transition_model": "natural envelope e^N; imposed trips at x/c=1,1",
              "dump_note": "integral quantities, not u(y); H/Hk format inferred and recorded",
              "returncode": None}
    log = ""
    try:
        completed = subprocess.run([str(executable)], input=commands, cwd=directory,
                                   stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                   text=True, encoding="ascii", errors="replace",
                                   timeout=timeout_s, check=False, shell=False,
                                   creationflags=subprocess.CREATE_NO_WINDOW if os.name == "nt" else 0)
        log = completed.stdout or ""
        result["returncode"] = completed.returncode
        if completed.returncode:
            raise ValueError(f"XFOIL exited with code {completed.returncode}")
        if re.search(r"convergence\s+failed|not\s+converged|floating[ -]point|segmentation", log, re.I):
            raise ValueError("XFOIL reported convergence or numerical failure")
        polar = parse_polar(directory / "polar.txt", alpha_deg)
        if (abs(polar["mach"] - mach) > .00051 or abs(polar["ncrit"] - ncrit) > .00051
                or abs(polar["reynolds"] - reynolds_local) > polar["reynolds_print_tolerance"]):
            raise ValueError("Polar Mach/Re/Ncrit do not match requested conditions")
        cp = parse_cp(directory / "cp.txt")
        bl = parse_dump(directory / "bl.txt", mach)
        count = len(cp["x_c"])
        if len(bl["x_c"]) < count or not np.allclose(cp["x_c"], bl["x_c"][:count], rtol=0, atol=1.1e-5):
            raise ValueError("Cp and DUMP surface nodes do not match")
        if np.any(np.diff(bl["s_c"][:count]) <= 0):
            raise ValueError("DUMP surface contour distance must strictly increase")
        lead = int(np.argmin(cp["x_c"]))
        if not 1 < lead < count - 2:
            raise ValueError("Cannot identify both geometric airfoil surfaces")
        lower_start = int(np.flatnonzero(np.isclose(cp["x_c"], cp["x_c"][lead],
                                                    rtol=0, atol=1e-12))[-1])
        result.update(converged=True, cp=cp, boundary_layer=bl, polar=polar,
                      surface_node_count=count, leading_edge_index=lead,
                      lower_surface_start_index=lower_start)
    except subprocess.TimeoutExpired as exc:
        raw = exc.stdout or ""
        log = raw.decode("ascii", errors="replace") if isinstance(raw, bytes) else raw
        result["reason"] = f"XFOIL timed out after {timeout_s:g} seconds"
    except (OSError, ValueError) as exc:
        result["reason"] = str(exc)
    (directory / "stdout.txt").write_text(log, encoding="utf-8")
    result["stdout_path"] = str(directory / "stdout.txt")
    result["commands_path"] = str(directory / "commands.txt")
    return result


def sample_laminar_station(result, surface, x_c):
    """Interpolate one attached laminar station without crossing transition.

    Both bracketing nodes must be upstream of the polar transition location
    (including half of the polar's 4-decimal printing interval). Upstream
    skin-friction history is checked along the signed-velocity flow branch,
    so a reattached point downstream of laminar separation is not accepted.
    Raises ValueError with a usable rejection reason. No extrapolation.
    """
    if not result.get("converged"):
        raise ValueError("A converged XFOIL result is required")
    if surface not in ("upper", "lower"):
        raise ValueError("surface must be upper or lower")
    xc = _finite(x_c, "x_c")
    if not 0 < xc < 1:
        raise ValueError("x_c must be in (0,1)")
    transition = float(result["polar"]["top_xtr" if surface == "upper" else "bot_xtr"])
    if xc >= transition - .00005:
        raise ValueError(f"Trip is at/downstream of predicted {surface} transition x/c={transition:g}")
    bl = result["boundary_layer"]
    count, lead = result["surface_node_count"], result["leading_edge_index"]
    lower_start = int(np.flatnonzero(np.isclose(bl["x_c"][:count], bl["x_c"][lead],
                                               rtol=0, atol=1e-12))[-1])
    indices = np.arange(lead + 1) if surface == "upper" else np.arange(lower_start, count)
    indices = indices[np.argsort(bl["x_c"][indices])]
    x = bl["x_c"][indices]
    if np.any(np.diff(x) <= 0) or not x[0] <= xc <= x[-1]:
        raise ValueError("Surface x coordinates must bracket the trip monotonically")
    right = int(np.searchsorted(x, xc, side="left"))
    if right < len(x) and math.isclose(float(x[right]), xc, rel_tol=0, abs_tol=1e-12):
        bracket = indices[right:right + 1]
    else:
        bracket = indices[right - 1:right + 1]
    if np.max(bl["x_c"][bracket]) >= transition - .00005:
        raise ValueError("Trip interpolation bracket reaches transition; refine the paneling")
    signed = bl["ue_over_uinf_signed"][:count]
    positive = np.flatnonzero(signed > 1e-7)
    negative = np.flatnonzero(signed < -1e-7)
    if not len(positive) or not len(negative) or positive[-1] >= negative[0]:
        raise ValueError("Cannot establish a single stagnation point and flow branches")
    sign = 1 if surface == "upper" else -1
    if np.any(sign * signed[bracket] <= 0):
        raise ValueError("Trip is on the other side of the stagnation point for this surface")
    history = (np.arange(int(np.min(bracket)), int(positive[-1]) + 1) if surface == "upper"
               else np.arange(int(negative[0]), int(np.max(bracket)) + 1))
    if (np.any(bl["cf"][history] <= 0) or np.any(bl["theta_c"][history] <= 0)
            or np.any(bl["delta_star_c"][history] <= 0)):
        raise ValueError("XFOIL indicates separated/unresolved boundary layer at or upstream of trip")
    sampled = {name: float(np.interp(xc, x, bl[name][indices])) for name in
               ("s_c", "y_c", "theta_c", "delta_star_c", "cf", "h_kinematic")}
    ratio = sign * float(np.interp(xc, x, signed[indices]))
    m_inf = float(result["mach"])
    temperature_ratio = 1 + .2 * m_inf**2 * (1 - ratio**2)
    if temperature_ratio <= 0:
        raise ValueError("Invalid edge state inferred from XFOIL velocity")
    edge_mach = m_inf * ratio / math.sqrt(temperature_ratio)
    if edge_mach > .3 + 1e-10:
        raise ValueError("Local edge Mach exceeds 0.30; low-speed profile surrogate not applicable")
    cp = float(np.interp(xc, x, result["cp"]["cp"][indices]))
    return dict(sampled, x_c=xc, surface=surface, ue_over_uinf=ratio, cp=cp,
                h_mass=sampled["delta_star_c"] / sampled["theta_c"],
                edge_mach=edge_mach, transition_xc=transition,
                interpolation_bracket_xc=[float(v) for v in bl["x_c"][bracket]],
                laminar_confirmed_by_model=True,
                method="XFOIL integral BL; wall-normal profile must be reconstructed")
