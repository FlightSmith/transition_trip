"""Generate the user's left-wing flat-plate baseline using a Blasius profile.

The supplied geometry and Reynolds reference are preserved.  The undisturbed
boundary layer is an idealized zero-pressure-gradient, incompressible flat
plate, not measured or simulated flow on the actual swept wing.  Both surfaces
are identical by this explicit modeling assumption; physical AoA is unknown.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np

from aero import profile_metrics, tunnel_state


ROOT_Y_M = -1.333
TIP_Y_M = -6.250
ROOT_X_LE_M = 9.280
TIP_X_LE_M = 12.557
ROOT_CHORD_M = 5.56
TIP_CHORD_M = 1.0
REFERENCE_MAC_M = 4.629
REFERENCE_AREA_M2 = 49.0
MACH = .2
REYNOLDS_MAC = 8_000_000.0
TEMPERATURE_K = 288.15
CONDITION_ID = "flat_plate_m02_re8e6"
TRIP_POSITIONS = (.07, .10, .12)
BLASIUS_WALL_SECOND_DERIVATIVE = .3320573362151963


def blasius_solution(integration_step: float = .001, sample_step: float = .01,
                     eta_max: float = 10.) -> dict:
    """Integrate f''' + .5 f f''=0, with f(0)=f'(0)=0, f'(infinity)=1.

    Classical fourth-order Runge--Kutta uses the established shooting value
    f''(0)=0.3320573362151963.  The similarity coordinate is
    eta = y sqrt(U/(nu*x)); velocity is u/U=f'(eta).

    At eta=10 the asymptotic velocity differs negligibly from one.  Normalize
    the final velocity to exactly one ONLY if its discrepancy is <=1e-8; this
    numerical correction is recorded in the returned metadata.  It keeps the
    generated profile compatible with strict edge/no-overshoot validation.
    """
    if not all(math.isfinite(v) and v > 0 for v in (integration_step, sample_step, eta_max)):
        raise ValueError("Blasius steps and eta_max must be finite and positive")
    count = round(eta_max / integration_step)
    stride = round(sample_step / integration_step)
    if (count < 1 or stride < 1 or count % stride
            or not math.isclose(count * integration_step, eta_max, abs_tol=1e-12)
            or not math.isclose(stride * integration_step, sample_step, abs_tol=1e-12)):
        raise ValueError("eta_max and sample_step must be integer multiples of integration_step")

    def rhs(state):
        return np.array((state[1], state[2], -.5 * state[0] * state[2]))

    state = np.array((0., 0., BLASIUS_WALL_SECOND_DERIVATIVE))
    samples = [state.copy()]
    step = integration_step
    for i in range(1, count + 1):
        k1 = rhs(state)
        k2 = rhs(state + .5 * step * k1)
        k3 = rhs(state + .5 * step * k2)
        k4 = rhs(state + step * k3)
        state = state + step / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        if i % stride == 0:
            samples.append(state.copy())
    samples = np.asarray(samples)
    edge_velocity = float(samples[-1, 1])
    if abs(edge_velocity - 1) > 1e-8:
        raise ValueError("Blasius integration has not reached the asymptotic edge accurately enough")
    fraction = samples[:, 1] / edge_velocity
    fraction[0], fraction[-1] = 0., 1.
    if np.any(np.diff(fraction) < 0) or np.any(fraction > 1):
        raise ValueError("Blasius integration produced a nonmonotone velocity profile")
    eta = np.arange(len(fraction), dtype=float) * sample_step
    metrics = profile_metrics(eta, fraction, 1., 1.)
    return {
        "eta": eta, "u_over_U": fraction,
        "eta_delta99": metrics["delta99_m"],
        "theta_coefficient": metrics["theta_m"],
        "displacement_coefficient": metrics["delta_star_m"],
        "H": metrics["H"],
        "integration_step": integration_step, "sample_step": sample_step,
        "eta_max": eta_max, "raw_edge_velocity_fraction": edge_velocity,
        "edge_normalization_multiplier": 1. / edge_velocity,
        "wall_second_derivative": BLASIUS_WALL_SECOND_DERIVATIVE,
    }


def geometry_metadata() -> dict:
    panel_span = abs(TIP_Y_M - ROOT_Y_M)
    taper = TIP_CHORD_M / ROOT_CHORD_M
    panel_area = panel_span * (ROOT_CHORD_M + TIP_CHORD_M) / 2
    trapezoid_mac = 2 / 3 * ROOT_CHORD_M * (1 + taper + taper**2) / (1 + taper)
    sweep = math.degrees(math.atan2(TIP_X_LE_M - ROOT_X_LE_M, panel_span))
    return {
        "root": {"y_m": ROOT_Y_M, "x_le_m": ROOT_X_LE_M, "chord_m": ROOT_CHORD_M},
        "tip": {"y_m": TIP_Y_M, "x_le_m": TIP_X_LE_M, "chord_m": TIP_CHORD_M},
        "analyzed_panel_span_m": panel_span, "leading_edge_sweep_deg": sweep,
        "provided_reference_area_m2": REFERENCE_AREA_M2,
        "provided_reference_mac_m": REFERENCE_MAC_M,
        "linear_panel_area_m2": panel_area,
        "linear_panel_trapezoid_mac_m": trapezoid_mac,
        "twice_linear_panel_area_m2": 2 * panel_area,
        "reference_note": (
            "User-provided Sref=49 m^2 and MAC=4.629 m are retained as references; "
            "they are not replaced by values derived from this single trapezoidal "
            "panel. Their different values may reflect a different reference "
            "definition or geometry not included here; no reconciliation is assumed."
        ),
        "coordinate_note": (
            "Station eta is fractional distance from the supplied root to tip. "
            "Station y_m is global spanwise Y (negative on this left panel); "
            "profile CSV y_m is instead local wall-normal distance. Global XLE "
            "and chord are linearly interpolated between the supplied endpoints."
        ),
    }


def make_config(x_c: float) -> dict:
    if x_c not in TRIP_POSITIONS:
        raise ValueError("Expected one of the requested 7%, 10%, or 12% chord positions")
    geometry = geometry_metadata()
    stations = []
    for eta in np.linspace(0, 1, 21):
        eta = float(eta)
        stations.append({
            "eta": eta,
            "y_m": ROOT_Y_M + eta * (TIP_Y_M - ROOT_Y_M),
            "x_le_m": ROOT_X_LE_M + eta * (TIP_X_LE_M - ROOT_X_LE_M),
            "chord_m": ROOT_CHORD_M + eta * (TIP_CHORD_M - ROOT_CHORD_M),
        })
    return {
        "schema_version": 1,
        "name": f"Left wing: Blasius flat-plate grit estimate at {100*x_c:g}% chord",
        "data_origin": (
            "IDEALIZED FLAT-PLATE ESTIMATE using user-supplied wing dimensions. "
            "Profiles are numerical solutions of the incompressible Blasius "
            "zero-pressure-gradient equations, not measured wing data or wing "
            "CFD. Both surfaces are identical by assumption. Surface distance "
            "s=x=(x/c)*local chord; Ue=Uinf and nu=nuinf. Sweep affects coordinates "
            "only; three-dimensional flow and crossflow are not modeled. Physical "
            "angle of attack was not supplied; alpha_deg=0 labels this idealized "
            "zero-pressure-gradient baseline and does not assert test AoA."
        ),
        "wing": {
            "semispan_m": geometry["analyzed_panel_span_m"],
            "sweep_le_deg": geometry["leading_edge_sweep_deg"],
            "reference_length_m": REFERENCE_MAC_M,
            "reference_area_m2": REFERENCE_AREA_M2,
            "stations": stations,
        },
        "conditions": [{
            "id": CONDITION_ID, "mach": MACH, "reynolds": REYNOLDS_MAC,
            "alpha_deg": 0., "temperature_k": TEMPERATURE_K,
            "alpha_label": "zero-pressure-gradient baseline; physical AoA unspecified",
        }],
        "trip": {
            "type": "grit", "position_xc": {"upper": x_c, "lower": x_c},
            "height_margin": 1.1, "max_k_delta99": .8,
            "available_heights_mm": [],
        },
        "analysis": {"source": "profiles", "profiles_file": "profiles.csv"},
    }


def _write_json(path: Path, data: dict) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(data, handle, indent=2, allow_nan=False)
        handle.write("\n")


def generate(output: Path) -> dict:
    """Write three fixed-position configs, shared profiles and geometry metadata."""
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    solution = blasius_solution()
    state = tunnel_state(MACH, REYNOLDS_MAC, REFERENCE_MAC_M, TEMPERATURE_K)
    config = make_config(.1)
    fields = ["condition_id", "eta", "surface", "x_c", "y_m", "u_m_s",
              "nu_m2_s", "ue_m_s", "edge_mach", "laminar_confirmed"]
    count = 0
    with (output / "profiles.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for station in config["wing"]["stations"]:
            for xc in TRIP_POSITIONS:
                x = xc * station["chord_m"]
                scale = math.sqrt(state["nu_inf_m2_s"] * x / state["U_inf_m_s"])
                y = solution["eta"] * scale
                u = solution["u_over_U"] * state["U_inf_m_s"]
                # Validate dimensional profiles before declaring this modeled
                # laminar baseline to the generic profile-input route.
                profile_metrics(y, u, state["nu_inf_m2_s"], state["U_inf_m_s"])
                for surface in ("upper", "lower"):
                    common = {"condition_id": CONDITION_ID, "eta": station["eta"],
                              "surface": surface, "x_c": xc,
                              "nu_m2_s": state["nu_inf_m2_s"],
                              "ue_m_s": state["U_inf_m_s"], "edge_mach": MACH,
                              "laminar_confirmed": "true"}
                    for local_y, velocity in zip(y, u):
                        writer.writerow({**common, "y_m": float(local_y), "u_m_s": float(velocity)})
                        count += 1
    config_names = []
    for xc in TRIP_POSITIONS:
        name = f"config_xc{round(xc*100):02d}.json"
        _write_json(output / name, make_config(xc))
        config_names.append(name)
    metadata = geometry_metadata()
    metadata["blasius_numerics"] = {k: v for k, v in solution.items() if k not in ("eta", "u_over_U")}
    metadata["blasius_numerics"]["equation"] = "f''' + 0.5*f*f''=0; eta=y*sqrt(U/(nu*x)); u/U=f'"
    metadata["blasius_numerics"]["normalization_note"] = (
        "Velocity fractions are divided by the computed f'(10), only after "
        "checking abs(f'(10)-1)<=1e-8. No larger correction is allowed."
    )
    metadata["tunnel_state"] = state
    metadata["temperature_note"] = (
        "Static temperature 288.15 K is an explicit default used to infer "
        "dimensional tunnel speed, density and pressure. With constant-property "
        "Blasius profiles at specified Re_MAC, nu/U=MAC/Re_MAC fixes the sizing "
        "lengths independently of this temperature default."
    )
    metadata["profile_laminar_confirmation"] = (
        "laminar_confirmed=true describes the explicitly assumed undisturbed "
        "Blasius model. It is not experimental confirmation that the wing stays laminar."
    )
    _write_json(output / "geometry_metadata.json", metadata)
    return {
        "output": str(output.resolve()), "config_files": config_names,
        "profile_rows": count, "span_stations": 21,
        "positions_xc": list(TRIP_POSITIONS), "profile_points": len(solution["eta"]),
        "eta_delta99": solution["eta_delta99"],
        "provided_reference_mac_m": REFERENCE_MAC_M,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path,
                        default=Path(__file__).parent / "examples" / "left_wing_flatplate")
    args = parser.parse_args()
    print(json.dumps(generate(args.output), indent=2))


if __name__ == "__main__":
    main()
