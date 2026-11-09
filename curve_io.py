"""Read ordered two-column section and pressure curves without mixing surfaces.

The supplied format has seven header lines, numeric data beginning on line 8,
and a terminal *END. Geometry is a single TE-to-LE-to-TE contour. Its first
branch may be either surface; upper/lower labels come from geometry. Pressure
curves must use that same branch order, even when sampled on a different grid.
"""

from __future__ import annotations

import math
from pathlib import Path
import re
from typing import Any

import numpy as np


_NUMBER = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eEdD][+-]?\d+)?"
_NAME = re.compile(rf"^(cp_dist|section)_state(\d+)_station({_NUMBER})(?:\.[A-Za-z][A-Za-z0-9_-]*)?$", re.I)


def _strip_comment(line: str) -> str:
    return re.split(r"[$#]", line, maxsplit=1)[0].strip()


def read_curve(path: str | Path, start_line: int = 8) -> np.ndarray:
    """Read finite two-column E/D-notation data through a required *END marker.

    Blank lines and $/# comments are ignored in the data block. Unexpected
    directives, malformed rows, and any noncomment content after *END are
    rejected, so a second curve block cannot silently be read as this curve.
    ``start_line`` is one-based and defaults to the supplied seven-line header.
    """
    path = Path(path)
    if isinstance(start_line, bool) or not isinstance(start_line, int) or start_line < 1:
        raise ValueError("start_line must be a positive one-based integer")
    lines = path.read_text(encoding="utf-8-sig").splitlines()
    if len(lines) < start_line:
        raise ValueError(f"{path.name}: file ends before numeric start line {start_line}")
    values = []
    ended = False
    for number, original in enumerate(lines[start_line - 1:], start=start_line):
        line = _strip_comment(original)
        if not line:
            continue
        if ended:
            raise ValueError(f"{path.name}:{number}: content after *END; multiple blocks are not supported")
        if line.upper() == "*END":
            ended = True
            continue
        if line.startswith("*"):
            raise ValueError(f"{path.name}:{number}: unexpected directive or multiple curve blocks")
        fields = line.split()
        if len(fields) != 2:
            raise ValueError(f"{path.name}:{number}: expected exactly two numeric columns")
        try:
            row = [float(field.replace("D", "E").replace("d", "e")) for field in fields]
        except ValueError as exc:
            raise ValueError(f"{path.name}:{number}: malformed numeric row") from exc
        if not all(math.isfinite(value) for value in row):
            raise ValueError(f"{path.name}:{number}: nonfinite numeric value")
        values.append(row)
    if not ended:
        raise ValueError(f"{path.name}: missing terminal *END")
    if len(values) < 3:
        raise ValueError(f"{path.name}: at least three numeric points are required")
    return np.asarray(values, dtype=float)


def parse_curve_name(path: str | Path) -> dict[str, Any]:
    """Return kind ('cp'/'section'), positive integer state and global station Y."""
    name = Path(path).name
    match = _NAME.fullmatch(name)
    if not match:
        raise ValueError(f"Unrecognized curve filename: {name}")
    prefix, state_text, station_text = match.groups()
    state = int(state_text)
    station = float(station_text.replace("D", "E").replace("d", "e"))
    if state < 1 or not math.isfinite(station):
        raise ValueError(f"Invalid state/station in curve filename: {name}")
    return {"kind": "cp" if prefix.lower() == "cp_dist" else "section",
            "state": state, "station_y_m": station}


def _minimum_block(x: np.ndarray, tolerance: float) -> tuple[int, int]:
    indices = np.flatnonzero(np.abs(x - np.min(x)) <= tolerance)
    if indices[0] == 0 or indices[-1] == len(x) - 1:
        raise ValueError("Contour must run TE-to-LE-to-TE; minimum X cannot be an endpoint")
    if np.any(np.diff(indices) != 1):
        raise ValueError("Contour has multiple separated minimum-X points")
    return int(indices[0]), int(indices[-1])


def _clean_section_branch(indices: np.ndarray, points: np.ndarray,
                          leading_edge: np.ndarray, tolerance: float) -> np.ndarray:
    kept = [int(indices[0])]
    trailing_edge = points[int(indices[-1])]
    for raw_index in indices[1:]:
        index = int(raw_index)
        previous = kept[-1]
        distance = float(np.linalg.norm(points[index] - points[previous]))
        at_edge = (np.linalg.norm(points[index] - leading_edge) <= tolerance
                   or np.linalg.norm(points[index] - trailing_edge) <= tolerance)
        if distance <= tolerance and at_edge:
            # Preserve the complete raw contour separately; redundant identical
            # edge samples are omitted only from the interpolation branch.
            continue
        if points[index, 0] - points[previous, 0] <= tolerance:
            raise ValueError("Section branch is nonmonotone in X or has ambiguous duplicate X")
        kept.append(index)
    if len(kept) < 3:
        raise ValueError("Each section surface requires at least three distinct points")
    return np.asarray(kept, dtype=int)


def load_section(path: str | Path, start_line: int = 8) -> dict[str, Any]:
    """Read a section, normalize its chord, and label both surfaces geometrically.

    LE is the unique minimum-X point (coincident repeats are allowed). TE is
    the midpoint of the first/last contour endpoints, including a finite gap.
    The physical chord includes their relative Z displacement; normalization
    translates LE, rotates by minus the chord incidence, and divides by chord.
    Every raw point/order is retained. Surface arrays run LE-to-TE and include
    their original contour indices. Unsupported folds/crossings are rejected.
    """
    path = Path(path)
    metadata = parse_curve_name(path)
    if metadata["kind"] != "section":
        raise ValueError("load_section requires a section_state... filename")
    points = read_curve(path, start_line=start_line)
    span_x = float(np.ptp(points[:, 0]))
    if span_x <= 0:
        raise ValueError("Section must have positive streamwise extent")
    tolerance = max(1e-12, 1e-9 * span_x)
    first_le, last_le = _minimum_block(points[:, 0], tolerance)
    le_group = points[first_le:last_le + 1]
    if np.max(np.linalg.norm(le_group - le_group[0], axis=1)) > tolerance:
        raise ValueError("Distinct minimum-X section points make the leading edge ambiguous")
    leading_edge = np.array((float(np.min(points[:, 0])), float(np.mean(le_group[:, 1]))))
    trailing_edge = .5 * (points[0] + points[-1])
    chord_vector = trailing_edge - leading_edge
    streamwise_chord = float(chord_vector[0])
    chord = float(np.linalg.norm(chord_vector))
    if streamwise_chord <= tolerance or chord <= tolerance:
        raise ValueError("Section leading/trailing edges do not define a positive chord")
    incidence = math.atan2(float(chord_vector[1]), float(chord_vector[0]))
    cosine, sine = math.cos(incidence), math.sin(incidence)
    shifted = points - leading_edge
    normalized = np.column_stack((shifted[:, 0] * cosine + shifted[:, 1] * sine,
                                  -shifted[:, 0] * sine + shifted[:, 1] * cosine)) / chord
    # Canonicalize edge roundoff so a rotated mathematical x/c=1 is not
    # serialized as 1+one ulp and rejected as out of range downstream.
    for edge in (0., 1.):
        normalized[np.abs(normalized[:, 0] - edge) <= tolerance / chord, 0] = edge
    branches = []
    for order, indices in enumerate((np.arange(first_le, -1, -1), np.arange(last_le, len(points)))):
        indices = _clean_section_branch(indices, points, leading_edge, tolerance)
        branch_xy = normalized[indices]
        if np.any(np.diff(branch_xy[:, 0]) <= tolerance / chord):
            raise ValueError("Section branch folds after chord normalization; XFOIL mapping is ambiguous")
        physical = points[indices]
        s = np.r_[0., np.cumsum(np.linalg.norm(np.diff(physical, axis=0), axis=1))]
        branches.append({
            "x_c": branch_xy[:, 0], "z_c": branch_xy[:, 1],
            "x_global_m": physical[:, 0], "z_global_m": physical[:, 1],
            "s_m": s, "contour_indices": indices, "contour_branch": order,
        })
    left = max(branches[0]["x_c"][0], branches[1]["x_c"][0])
    right = min(branches[0]["x_c"][-1], branches[1]["x_c"][-1])
    if right <= left:
        raise ValueError("Section surfaces have no overlapping chord interval")
    # The difference of the two piecewise-linear surfaces is linear between
    # their combined knots. Checking every interior knot detects even very
    # narrow crossings that a fixed coarse sampling grid could miss.
    samples = np.unique(np.r_[branches[0]["x_c"], branches[1]["x_c"]])
    samples = samples[(samples > left + tolerance / chord) & (samples < right - tolerance / chord)]
    if len(samples) == 0:
        raise ValueError("Section has no resolved interior thickness")
    difference = (np.interp(samples, branches[0]["x_c"], branches[0]["z_c"])
                  - np.interp(samples, branches[1]["x_c"], branches[1]["z_c"]))
    geometric_tolerance = tolerance / chord
    if np.max(np.abs(difference)) <= geometric_tolerance:
        raise ValueError("Coincident section branches cannot be labeled upper/lower")
    if np.any(difference > geometric_tolerance) and np.any(difference < -geometric_tolerance):
        raise ValueError("Section branches cross; upper/lower surfaces are ambiguous")
    if np.any(np.abs(difference) <= geometric_tolerance):
        raise ValueError("Section branches touch or duplicate points away from the leading/trailing edge")
    upper_index = 0 if float(np.mean(difference)) > 0 else 1
    surfaces = {"upper": branches[upper_index], "lower": branches[1 - upper_index]}
    return {
        **metadata, "path": str(path), "x_m": points[:, 0], "z_m": points[:, 1],
        "x_le_m": float(leading_edge[0]), "z_le_m": float(leading_edge[1]),
        "x_te_m": float(trailing_edge[0]), "z_te_m": float(trailing_edge[1]),
        "chord_m": chord, "streamwise_chord_m": streamwise_chord,
        "incidence_deg": math.degrees(incidence), "rotation_deg": -math.degrees(incidence),
        "normalized_x": normalized[:, 0], "normalized_z": normalized[:, 1],
        "leading_edge_indices": np.arange(first_le, last_le + 1),
        "surfaces": surfaces, "coordinate_tolerance_m": tolerance,
    }


def _pressure_coefficients(ordinate: np.ndarray, value_kind: str,
                           dynamic_pressure_pa: float | None,
                           static_pressure_pa: float | None) -> np.ndarray:
    if value_kind == "cp":
        return ordinate.copy()
    if value_kind not in ("pressure", "absolute_pressure", "gauge_pressure"):
        raise ValueError("value_kind must be cp, pressure, absolute_pressure, or gauge_pressure")
    if dynamic_pressure_pa is None or static_pressure_pa is None:
        raise ValueError("Pressure conversion requires both dynamic_pressure_pa and static_pressure_pa")
    dynamic = float(dynamic_pressure_pa)
    reference = float(static_pressure_pa)
    if not math.isfinite(dynamic) or dynamic <= 0 or not math.isfinite(reference):
        raise ValueError("Dynamic pressure must be finite/positive and reference static pressure finite")
    if value_kind in ("pressure", "absolute_pressure") and (reference <= 0 or np.any(ordinate <= 0)):
        raise ValueError("Absolute pressure samples and reference must be positive; use gauge_pressure for signed gauge values")
    cp = (ordinate - reference) / dynamic
    if not np.all(np.isfinite(cp)):
        raise ValueError("Pressure conversion produced nonfinite Cp")
    return cp


def map_curve_to_section(values: Any, section: dict[str, Any], value_kind: str = "cp",
                         dynamic_pressure_pa: float | None = None,
                         static_pressure_pa: float | None = None) -> dict[str, dict[str, Any]]:
    """Map Cp/pressure to LE-to-TE surfaces using contour order, never global sorting.

    Identical grids retain exact section indices. Otherwise each monotone
    TE-to-LE/LE-to-TE branch is mapped separately using section X to interpolate
    Z and arc length, with no extrapolation beyond a branch's geometry. Branch
    identity follows the section contour; scalar Cp cannot reveal a reversed
    contour convention, so the caller must ensure the conventions match.

    ``pressure`` and ``absolute_pressure`` require positive pressure values.
    ``gauge_pressure`` permits signed values with an explicit gauge reference.
    Both require Cp=(p-static_pressure_pa)/dynamic_pressure_pa; the supplied
    static reference must use the same absolute/gauge convention as the data.
    """
    data = read_curve(values) if isinstance(values, (str, Path)) else np.asarray(values, dtype=float)
    if data.ndim != 2 or data.shape[1] != 2 or len(data) < 3 or not np.all(np.isfinite(data)):
        raise ValueError("Curve data must be a finite Nx2 array with at least three rows")
    coefficients = _pressure_coefficients(data[:, 1], value_kind, dynamic_pressure_pa, static_pressure_pa)
    tolerance = float(section["coordinate_tolerance_m"])
    identical = (len(data) == len(section["x_m"])
                 and np.allclose(data[:, 0], section["x_m"], rtol=0, atol=tolerance))
    if not identical:
        first_le, last_le = _minimum_block(data[:, 0], tolerance)
        if last_le - first_le > 1:
            raise ValueError("Cp curve has more than two minimum-X samples; branch split is ambiguous")
        curve_branches = (np.arange(first_le, -1, -1), np.arange(last_le, len(data)))
        for indices in curve_branches:
            if len(indices) < 2 or np.any(np.diff(data[indices, 0]) <= tolerance):
                raise ValueError("Cp branch must be strictly monotone in X; duplicate X within a surface is ambiguous")
    result = {}
    angle = math.radians(section["incidence_deg"])
    for name, surface in section["surfaces"].items():
        if identical:
            indices = surface["contour_indices"]
            x_global = surface["x_global_m"].copy()
            z_global = surface["z_global_m"].copy()
            x_c = surface["x_c"].copy()
            s = surface["s_m"].copy()
            method = "matching contour indices"
        else:
            indices = curve_branches[surface["contour_branch"]]
            x = data[indices, 0]
            section_x = surface["x_global_m"]
            if x[0] < section_x[0] - tolerance or x[-1] > section_x[-1] + tolerance:
                raise ValueError(f"{name} Cp grid exceeds section coordinates; extrapolation is not allowed")
            # Tolerance is only for text roundoff at identical physical edges.
            x_global = np.clip(x, section_x[0], section_x[-1])
            z_global = np.interp(x_global, section_x, surface["z_global_m"])
            s = np.interp(x_global, section_x, surface["s_m"])
            x_c = ((x_global - section["x_le_m"]) * math.cos(angle)
                   + (z_global - section["z_le_m"]) * math.sin(angle)) / section["chord_m"]
            for edge in (0., 1.):
                x_c[np.abs(x_c - edge) <= tolerance / section["chord_m"]] = edge
            if np.any(np.diff(x_c) <= 0):
                raise ValueError(f"{name} Cp branch is nonmonotone after chord normalization")
            method = "separate contour-branch interpolation"
        result[name] = {
            "x_c": x_c, "cp": coefficients[indices], "s_m": s,
            "x_global_m": x_global, "z_global_m": z_global,
            "curve_indices": np.asarray(indices, dtype=int), "mapping_method": method,
        }
    return result
