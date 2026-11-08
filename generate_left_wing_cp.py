"""Create portable Cp-input examples for the supplied left-wing geometry.

Every generated Cp is ZERO PLACEHOLDER DATA for a flat-plate demonstration,
NOT a pressure measurement or prediction on the swept wing.  Replace cp.csv,
the condition metadata, and the configuration's provenance with actual data
before interpreting the calculation as a wing estimate.  A run of this
generator deliberately replaces its example files in the selected directory.

The placeholder grid starts at x=0 with finite Ue, as in the flat-plate model,
so theta0=0 is appropriate ONLY for this example.  Measured Cp that begins at
x>0 requires a positive, physically justified theta0_m at its first station;
do not copy zero or assume missing upstream boundary-layer history is zero.
Supply surface distance s_m if known instead of the example's s=x_c*chord.
The Cp route is the Thwaites/quartic approximation and is not identical to
the Blasius profile route in generate_left_wing.py.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from generate_left_wing import TRIP_POSITIONS, make_config as make_profile_config


CONDITION_ID = "case_01"
CP_FIELDS = ["condition_id", "eta", "surface", "x_c", "s_m", "cp", "theta0_m"]


def make_config(x_c: float) -> dict:
    """Reuse the user's geometry/references with explicitly placeholder Cp input."""
    config = make_profile_config(x_c)
    config["name"] = f"Left-wing Cp input example: {100*x_c:g}% chord (placeholder Cp=0)"
    config["data_origin"] = (
        "PLACEHOLDER Cp=0 DEMONSTRATION, NOT ACTUAL WING Cp. The supplied "
        "left-wing geometry, MAC=4.629 m, Sref=49 m^2, Mach=0.2 and "
        "Re_MAC=8e6 are retained. All Cp values are zero solely to exercise "
        "the input workflow with a zero-pressure-gradient flat-plate baseline. "
        "Both surfaces are identical by this assumption. The Cp route uses "
        "Thwaites integral thickness and a quartic profile surrogate, not "
        "the numerical Blasius profiles of the separate flat-plate example. "
        "alpha_deg=0 is a placeholder label; physical test angle of attack "
        "is unspecified and must be set when adding real cases. Sweep sets "
        "coordinates only; three-dimensional flow and crossflow are not solved. "
        "Replace cp.csv and this provenance with the actual measurement or "
        "CFD source before using the results as a wing estimate."
    )
    condition = config["conditions"][0]
    condition["id"] = CONDITION_ID
    condition["alpha_deg"] = 0.
    condition["alpha_label"] = "placeholder zero-pressure-gradient case; physical AoA unspecified"
    config["analysis"] = {"source": "cp", "cp_file": "cp.csv"}
    return config


def _write_json(path: Path, value: dict) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(value, handle, indent=2, allow_nan=False)
        handle.write("\n")


def generate(output: Path) -> dict:
    """Write the shared placeholder Cp CSV and 7%, 10%, 12% chord configs."""
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    stations = make_config(.10)["wing"]["stations"]
    # Integer division gives a deterministic 0.002-spaced grid including all
    # three requested positions without almost-duplicate floating-point points.
    x_c = [i / 500 for i in range(176)]
    count = 0
    with (output / "cp.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=CP_FIELDS, lineterminator="\n")
        writer.writeheader()
        for station in stations:
            for surface in ("upper", "lower"):
                for i, x in enumerate(x_c):
                    writer.writerow({
                        "condition_id": CONDITION_ID, "eta": station["eta"],
                        "surface": surface, "x_c": x,
                        "s_m": x * station["chord_m"],
                        "cp": 0., "theta0_m": 0. if i == 0 else "",
                    })
                    count += 1
    config_names = []
    for xc in TRIP_POSITIONS:
        name = f"config_xc{round(xc*100):02d}.json"
        _write_json(output / name, make_config(xc))
        config_names.append(name)
    return {
        "output": str(output.resolve()), "config_files": config_names,
        "cp_file": "cp.csv", "cp_rows": count,
        "span_stations": len(stations), "conditions": 1,
        "positions_xc": list(TRIP_POSITIONS),
        "cp_warning": "ALL Cp VALUES ARE ZERO PLACEHOLDERS; REPLACE WITH ACTUAL WING DATA.",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path,
                        default=Path(__file__).parent / "examples" / "left_wing_cp")
    args = parser.parse_args()
    print(json.dumps(generate(args.output), indent=2))


if __name__ == "__main__":
    main()
