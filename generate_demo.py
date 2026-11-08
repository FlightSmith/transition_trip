"""Generate reproducible, SYNTHETIC ONERA-M6-inspired trip-design inputs.

The planform is inspired by published M6 dimensions. All pressures and boundary
layers are analytic/model-generated examples, not experimental ONERA M6 data.
"""
from __future__ import annotations

import argparse
import copy
import csv
import json
from pathlib import Path

import numpy as np

from aero import build_surrogate_profile, edge_state, thwaites, tunnel_state


def synthetic_cp(x_c: np.ndarray, eta: float, alpha_deg: float,
                 surface: str) -> np.ndarray:
    """A mild-recovery pressure curve; no aerodynamic validation is implied."""
    base = -0.22 - 0.48 * (1.0 - 0.3 * eta) * np.exp(-x_c / 0.25) + 0.2 * x_c
    alpha_term = 0.035 * alpha_deg * (1.0 - 0.3 * eta) * np.exp(-x_c / 0.25)
    return base - alpha_term if surface == "upper" else base + alpha_term


def make_config() -> dict:
    stations = [{"eta": eta, "chord_m": 0.8059 * (1.0 - (1.0 - 0.562) * eta)}
                for eta in (0.2, 0.44, 0.65, 0.8, 0.9)]
    conditions = [
        {"id": f"a{alpha}_r{re_million}", "mach": 0.18,
         "reynolds": re_million * 1_000_000.0, "alpha_deg": float(alpha),
         "temperature_k": 288.15}
        for re_million in (3, 6) for alpha in (0, 3, 6)
    ]
    return {
        "schema_version": 1,
        "name": "Synthetic M6-inspired wing: 10% chord grit trip",
        "data_origin": (
            "SYNTHETIC example. Published M6-inspired planform only; Cp is an "
            "analytic mild-recovery curve and is not ONERA test data or CFD. "
            "Profiles are generated from the same Thwaites/quartic surrogate "
            "as the Cp route and are not independent validation. Surface distance "
            "is approximated by x_c times local chord. Initial theta uses an "
            "assumed linear stagnation-to-first-point velocity ramp."
        ),
        "wing": {"semispan_m": 1.1963, "sweep_le_deg": 30.0,
                 "reference_length_m": 0.64607, "stations": stations},
        "conditions": conditions,
        "trip": {
            "type": "grit", "position_xc": {"upper": 0.1, "lower": 0.1},
            "height_margin": 1.1, "max_k_delta99": 0.8,
            "available_heights_mm": [0.05, 0.06, 0.08, 0.10, 0.12, 0.15,
                                     0.18, 0.20, 0.25, 0.30, 0.40, 0.50],
        },
        "analysis": {"source": "cp", "cp_file": "cp.csv"},
    }


def write_json(path: Path, data: dict) -> None:
    path.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n",
                    encoding="utf-8")


def generate(output: Path) -> dict:
    """Write deterministic inputs and return a compact generation summary."""
    output.mkdir(parents=True, exist_ok=True)
    config = make_config()
    x_c = np.unique(np.r_[np.linspace(0.002, 0.35, 180), 0.1])
    trip_index = int(np.flatnonzero(x_c == 0.1)[0])
    cp_fields = ["condition_id", "eta", "surface", "x_c", "s_m", "cp", "theta0_m"]
    profile_fields = ["condition_id", "eta", "surface", "x_c", "y_m",
                      "u_m_s", "nu_m2_s", "ue_m_s", "edge_mach",
                      "laminar_confirmed"]
    cp_count = profile_count = 0
    with (output / "cp.csv").open("w", newline="", encoding="utf-8") as cp_file, \
            (output / "profiles.csv").open("w", newline="", encoding="utf-8") as profile_file:
        cp_writer = csv.DictWriter(cp_file, fieldnames=cp_fields)
        profile_writer = csv.DictWriter(profile_file, fieldnames=profile_fields)
        cp_writer.writeheader()
        profile_writer.writeheader()
        for condition in config["conditions"]:
            state = tunnel_state(condition["mach"], condition["reynolds"],
                                 config["wing"]["reference_length_m"],
                                 condition["temperature_k"])
            for station in config["wing"]["stations"]:
                s_m = x_c * station["chord_m"]
                for surface in ("upper", "lower"):
                    cp = synthetic_cp(x_c, station["eta"], condition["alpha_deg"], surface)
                    edge = edge_state(cp, state, condition["mach"])
                    ue = edge["U_e_m_s"]
                    # Explicit synthetic initial condition: the stagnation-flow
                    # limit theta^2 = 0.075 nu / (dU/ds), with dU/ds = Ue0/s0.
                    theta0 = float(np.sqrt(0.075 * state["nu_inf_m2_s"] * s_m[0] / ue[0]))
                    layer = thwaites(s_m, ue, state["nu_inf_m2_s"], theta0)
                    profile = build_surrogate_profile(
                        float(layer["theta_m"][trip_index]),
                        float(layer["lambda_thwaites"][trip_index]),
                        float(ue[trip_index]), state["nu_inf_m2_s"],
                    )
                    common = {"condition_id": condition["id"], "eta": station["eta"],
                              "surface": surface}
                    for i, x in enumerate(x_c):
                        cp_writer.writerow({**common, "x_c": float(x), "s_m": float(s_m[i]),
                                            "cp": float(cp[i]),
                                            "theta0_m": theta0 if i == 0 else ""})
                        cp_count += 1
                    for y, u, nu in zip(profile["y_m"], profile["u_m_s"], profile["nu_m2_s"]):
                        profile_writer.writerow({
                            **common, "x_c": 0.1, "y_m": float(y), "u_m_s": float(u),
                            "nu_m2_s": float(nu), "ue_m_s": float(ue[trip_index]),
                            "edge_mach": float(edge["M_e"][trip_index]),
                            "laminar_confirmed": "true",
                        })
                        profile_count += 1
    write_json(output / "config.json", config)
    profile_config = copy.deepcopy(config)
    profile_config["name"] = "Synthetic M6-inspired wing: supplied-profile input"
    profile_config["analysis"] = {"source": "profiles", "profiles_file": "profiles.csv"}
    write_json(output / "config_profiles.json", profile_config)
    transonic = copy.deepcopy(config)
    transonic["name"] = "Deliberately invalid transonic Cp fallback demonstration"
    transonic["data_origin"] += (
        " This negative test changes Mach to 0.84 while retaining the synthetic Cp. "
        "It must be rejected by the incompressible boundary-layer fallback."
    )
    for condition in transonic["conditions"]:
        condition["mach"] = 0.84
    write_json(output / "config_transonic_rejected.json", transonic)
    return {"output": str(output.resolve()), "cp_rows": cp_count,
            "profile_rows": profile_count,
            "conditions": len(config["conditions"]),
            "span_stations": len(config["wing"]["stations"]),
            "fixed_x_c": 0.1}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path,
                        default=Path(__file__).parent / "examples" / "m6_synthetic")
    args = parser.parse_args()
    print(json.dumps(generate(args.output), indent=2))


if __name__ == "__main__":
    main()
