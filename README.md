# Wind-tunnel trip sizing from an undisturbed laminar boundary layer

This Python tool sizes a trip at a specified fraction of local chord on the upper
and lower surfaces of a wing. It evaluates every supplied tunnel condition and
angle of attack, then selects a common exposed height at each spanwise station.
The example uses **10% chord on both surfaces**, as requested.

The calculation is a preliminary engineering estimate. It evaluates a roughness
Reynolds-number criterion; it does not predict natural transition, the distance
from the trip to complete transition, or the drag penalty of a particular trip.
Verify the selected trip in the tunnel using an appropriate transition diagnostic.

## Run the example

Use Python 3.10 or later and install the dependencies listed in `requirements.txt`:

```text
python -m pip install -r requirements.txt
python generate_demo.py
python trip_design.py --config examples/m6_synthetic/config.json --output results/m6_cp
python trip_design.py --config examples/m6_synthetic/config_profiles.json --output results/m6_profiles
```

Open `report.html` in each output directory. The machine-readable results are
`design.csv`, `conditions.csv`, and `summary.json`. If plotting is available,
`overview.png` and the vector version `overview.svg` summarize the design across
the span. Add `--no-plot` to omit those figures. A nonzero exit status means
the requested complete design did not pass all checks or the inputs were invalid;
inspect the report and condition-level reasons.

Input provenance in `summary.json` records SHA-256 hashes and paths relative to
the configuration directory. Inputs outside that directory are identified by
filename and hash, without recording machine-specific absolute paths.

On the prepared Windows workspace, `run_example.ps1` runs both examples using
an available Python installation or the Codex bundled Python runtime. Plotting
dependencies are installed locally in `.runtime`; that directory is excluded
from the portable archive. On another machine, install `requirements.txt` first.

For the supplied synthetic case, both routes select **0.25 mm exposed grit
height at 10% chord on both surfaces**, at all five sampled span stations.
The tested cases cover Mach 0.18, reference Reynolds numbers 3 and 6 million,
and angles of attack 0, 3, and 6 degrees. The maximum selected `k/delta99` is
0.7711, and the minimum selected `Re_k` is 843.5. These are demonstration
results for the artificial input data, not a trip prescription for the real M6.

Run the numerical and input-validation tests from this directory with:

```text
python -m unittest discover -v
```

These check mathematical limiting cases, validation behavior, and complete
synthetic workflows. They are software verification, not experimental validation
of trip effectiveness.

The deliberately invalid high-Mach example should fail:

```text
python trip_design.py --config examples/m6_synthetic/config_transonic_rejected.json --output results/m6_rejected
```

## What the example represents

The planform is inspired by the published ONERA M6 wing: nominal root chord
0.8059 m, taper ratio 0.562, semispan 1.1963 m, leading-edge sweep 30 degrees,
and reference mean aerodynamic chord 0.64607 m. The five evaluated span fractions
are 0.20, 0.44, 0.65, 0.80, and 0.90 of semispan. The rounded original wingtip and
the airfoil's true surface geometry are not modeled.

The six example conditions combine Mach 0.18, reference Reynolds numbers
3 million and 6 million, and angles of attack 0, 3, and 6 degrees. Pressure curves
are smooth, analytic examples with mild adverse pressure gradients. Upper and
lower pressures are symmetric at zero angle of attack, with equal and opposite
synthetic loading increments as incidence increases. They are
**not measured ONERA M6 pressures, a reproduction of an M6 transonic case, or
validated CFD**. Surface distance is approximated by `s_m = x_c * chord_m`.
The initial momentum thickness assumes a linear velocity ramp from a virtual
stagnation point to the first pressure station.

The profile-input example is generated from the same pressure/integral/profile
model as the Cp example. Agreement between those routes is an input-path check,
not independent aerodynamic validation. `laminar_confirmed=true` in those files
means that the synthetic example is defined as laminar; it is not experimental
evidence about a real wing.

## Method and selection policy

At a candidate height `k`, the tool evaluates the undisturbed velocity and
kinematic viscosity at that height:

```text
Re_k(k) = u(k) * k / nu(k)
```

It solves this equation for the requested critical roughness Reynolds number,
applies the configured height margin, and takes the largest required height
across all specified conditions. Available exposed heights, when supplied, are
searched in ascending order for a height meeting every condition's criterion
and thickness limit. Without a hardware list, the continuous common-height
search evaluates the piecewise profile boundaries and exact criterion-crossing
boundaries. This handles supplied profiles whose varying viscosity makes
`Re_k(k)` nonmonotone: merely choosing the largest individual minimum would not
always suffice. An infeasible case is reported
explicitly; it must not be treated as a valid trip recommendation.

The built-in grit criterion uses `Re_k = 600`, a minimum running freestream
Reynolds number `Re_x = U_inf * x / nu_inf = 100000`, and a subsonic edge-Mach
limit of 0.95. These are restricted to preliminary sizing of roughly equiaxed,
sparsely distributed grit. The 0.95 limit is a conservative software applicability
limit, not a physical transition boundary. The NASA grit guidance does not
establish transferable criteria for arbitrary wires, zigzag tapes, adhesive
ridges, or densely packed roughness.

The default example height multiplier of 1.10 and maximum `k/delta99` of 0.80 are
**user-controlled design policies**, not universal transition correlations.
The reported height is the exposed protrusion normal to the surface, including
the effect of adhesive on protrusion, not nominal particle diameter or tape stock
thickness. A stock-size list should therefore describe effective exposed heights.

For a swept wing, sweep locates the trip in planform:

```text
y = eta * semispan
x_global = y * tan(sweep_le) + (x/c) * local_chord
```

For a nonstraight leading edge, each station can instead provide `x_le_m`;
its placement becomes `x_global = x_le_m + (x/c) * local_chord`.

It does not provide a crossflow or attachment-line transition model. The tool
does not apply an unexplained cosine correction to turn strip results into a
three-dimensional boundary layer. Angle of attack identifies each input case;
its effect must already be present in that case's Cp or supplied profiles.
The envelope covers only the discrete supplied cases. Sampling 0, 3, and
6 degrees does not establish feasibility at every intervening angle; add cases
where flow changes rapidly.

## Boundary-layer input choices

**Cp route:** local edge conditions are obtained from compressible isentropic
relations for ideal-gas air. The momentum thickness is then estimated with the
incompressible Thwaites integral using freestream viscosity:

```text
theta(s)^2 = [theta0^2 * Ue(s0)^6 + 0.45 * nu * integral(Ue^5 ds)] / Ue(s)^6
lambda = theta^2 / nu * dUe/ds
```

The fallback requires both freestream Mach and local edge Mach to be at most
0.3, with the local edge limit enforced over the upstream integration history.
It rejects approximate laminar separation at
`lambda <= -0.09` and does not restart an attached layer downstream. A quartic
Pohlhausen profile matched to the Thwaites momentum thickness and shape factor
provides `u(y)` and `delta99`. This is a surrogate profile, particularly in a
pressure gradient, rather than a direct boundary-layer solution. Positive
initial surface distance requires an explicit initial momentum thickness.
Good leading-edge Cp coverage, surface distance, and initial conditions matter.

The implemented Thwaites shape-factor fits are:

```text
H = 2.088 + 0.0731/(lambda + 0.14)         for -0.09 < lambda < 0
H = 2.61 - 3.75*lambda + 5.24*lambda^2   for 0 <= lambda <= 0.1
```

With `z = y/delta`, the reconstruction uses:

```text
u/Ue = 2*z - 2*z^3 + z^4 + (Q/6)*z*(1-z)^3,  0 <= z <= 1
theta/delta      = 37/315 - Q/945 - Q^2/9072
delta_star/delta = 3/10 - Q/120
```

The code solves for `Q` by matching the ratio `H = delta_star/theta`, then
scales `delta` to match the integral momentum thickness. `Q` is a fitted shape
coefficient, **not** the measured local pressure-gradient parameter. The
monotone quartic family is restricted to `-12 <= Q <= 12`; cases outside that
family or the stated Thwaites fit range are rejected. `delta99` is the height
where the reconstructed velocity first reaches `0.99 Ue`, and differs from the
quartic support height `delta`.

**Profile route:** use a measured or computed undisturbed laminar profile at the
specified trip station for every condition, span station, and requested surface.
The tool uses local `u(y)` and `nu(y)` directly, without a Thwaites reconstruction.
Profiles must reach the boundary-layer edge and include the no-slip wall point.
For compressible data, supply local `nu = mu/rho`, not a constant freestream
value. The input velocity is the tangential speed used by the chosen roughness
correlation; the tool cannot infer a 3D velocity vector from Cp. A compressible
profile does not by itself validate a roughness criterion at that Mach number.
The reported momentum thickness and shape factor for a supplied profile are
incompressible integrals truncated at its last sample. They are diagnostics, not
compressible momentum-thickness definitions; the sizing calculation uses the
actual supplied `u(y)` and `nu(y)`.

## Input files and units

All dimensional input uses SI units except hardware heights, which explicitly
use millimeters. Reynolds numbers are dimensionless; angles are degrees.
JSON file references are relative to the configuration file's directory.

Configuration sections:

| Section | Required content |
| --- | --- |
| `schema_version` | `1` |
| `name`, `data_origin` | Optional case label and provenance description; strongly recommended |
| `wing` | `semispan_m`, `sweep_le_deg`, `reference_length_m`, and `stations` |
| `wing.stations` | Objects containing `eta` and streamwise `chord_m`; optional `x_le_m` overrides sweep-based placement |
| `conditions` | Objects containing unique `id`, `mach`, `reynolds`, `alpha_deg`, and optional `temperature_k` |
| `trip` | `type`, `position_xc`; optional `height_margin` (default 1.10), `max_k_delta99` (default 0.80), and `available_heights_mm` |
| `analysis` | `source: "cp"` with `cp_file`, or `source: "profiles"` with `profiles_file` |

`height_margin` must be at least 1 and `max_k_delta99` must be in `(0,1]`.
If no hardware-height list is supplied, the tool reports the continuous required
height after the margin. It also checks whether one uniform height works at
every supplied station and surface; that alternative is reported separately.

`eta` is `y/semispan`, not `y/full_span`. `reynolds` is based on the one declared
`reference_length_m`; it is not a unit Reynolds number. Mach, static temperature,
and Reynolds number determine a self-consistent freestream density/pressure via
the ideal-gas law and Sutherland viscosity. Static temperature defaults to
288.15 K; provide the actual tunnel value for dimensional sizing.

Set the fixed trip positions independently, for example:

```json
"position_xc": {"upper": 0.1, "lower": 0.1}
```

One surface can be omitted to analyze only the other. An optional global
`analysis.laminar_until_xc` rejects trips at or downstream of a known laminar
limit in either input route. In Cp mode, an optional `known_transition_xc`
column supplies a more specific limit for each case/station/surface series.
Populate one or more rows; all nonempty values must be finite and consistent
within that series. The absence of either limit does not establish laminar flow:
Cp mode assumes it, while profile mode requires the explicit confirmation.

Cp CSV columns:

| Column | Meaning |
| --- | --- |
| `condition_id` | Exact condition ID from JSON |
| `eta`, `surface` | Declared span station and `upper` or `lower` |
| `x_c` | Streamwise local chord fraction |
| `s_m` | Optional distance along the modeled surface/streamline from the chosen origin, m |
| `cp` | `(p - p_inf)/(0.5 rho_inf U_inf^2)` |
| `theta0_m` | Momentum thickness at the first row of that pressure series; mandatory and positive for `s_m > 0`; leave other rows blank |
| `known_transition_xc` | Optional known transition location; nonempty entries must be finite and consistent within the series |

Use ordered, increasing `x_c` and `s_m` without duplicates. Cover the upstream
history and the requested trip position; missing data must not be silently
extrapolated. Surface distance and chord fraction are distinct coordinates on
a real wing.

If `s_m` is omitted for the entire series, the code uses `s_m = x_c * chord_m`
and labels that approximation in `conditions.csv`. Provide true surface or
streamline distance where it matters; a mixture of populated and blank `s_m`
cells is rejected. The tool linearly interpolates Cp and surface distance to an
exact requested trip station, using at least three samples through that station.
No pressure or surface-distance extrapolation is performed. A zero initial
distance permits zero initial momentum thickness, but exact stagnation velocity
is not supported by this integrator; use a justified downstream initialization.

Profile CSV columns:

| Column | Meaning |
| --- | --- |
| `condition_id`, `eta`, `surface`, `x_c` | Identify the exact trip station and condition |
| `y_m` | Distance normal to the wall, starting at zero |
| `u_m_s` | Undisturbed tangential velocity at `y_m` |
| `nu_m2_s` | Local kinematic viscosity at `y_m` |
| `ue_m_s` | Profile edge speed, repeated on every row |
| `edge_mach` | Profile edge Mach, repeated on every row |
| `laminar_confirmed` | `true` only when the profile's laminar status is established or explicitly assumed for a synthetic example |

Profiles must be monotone and attached; separated or reversed-flow inputs are
outside the implemented roughness-sizing model. The tool does not silently
extend insufficient wall-normal coverage. Add more angle-of-attack/Mach/Reynolds
cases by adding matching configuration entries and complete CSV data for each.
The profile route requires a profile at the exact requested `x_c` (within
`1e-9`); it does not interpolate profiles between chordwise or spanwise
locations. Samples are interpolated only in the wall-normal direction for
roughness-height sizing.

## Other trip types and calibrated criteria

Use `type: "zigzag"`, `type: "wire"`, or `type: "custom"` only with a calibration appropriate to
the device geometry and tunnel conditions. Add a `calibration` object inside
`trip`, containing `critical_re_k`, `min_re_x`, `max_edge_mach`, and `notes`
describing the source and applicability. This changes the roughness criterion,
not the low-Mach restriction of the Cp-based boundary-layer fallback. No
universal wire or zigzag threshold is supplied.

For a physical model, establish the laminar extent upstream of the trip and
use a boundary-layer method suitable for the sweep, pressure gradients, wall
temperature, and Mach range. Natural transition, crossflow, leading-edge
contamination, shocks, and tunnel turbulence are not predicted by this tool.

## References

- Braslow, Hicks, and Harris,
  [Use of Grit-Type Boundary-Layer-Transition Trips on Wind-Tunnel Models](https://ntrs.nasa.gov/api/citations/19660026829/downloads/19660026829.pdf),
  NASA TN D-3579 (1966). Roughness criterion, running-Reynolds applicability,
  and practical grit installation guidance.
- Braslow and Knox,
  [Simplified Method for Determination of Critical Height of Distributed Roughness Particles for Boundary-Layer Transition at Mach Numbers from 0 to 5](https://ntrs.nasa.gov/citations/19930085292),
  NACA TN-4363 (1958). Basis for evaluating undisturbed properties at roughness
  height; the broader report is not implemented as a universal compressible
  correlation here.
- Thwaites,
  [Approximate Calculation of the Laminar Boundary Layer](https://doi.org/10.1017/S0001925900000184)
  (1949).
- NASA Glenn, [ONERA M6 Wing](https://www.grc.nasa.gov/WWW/wind/valid/m6wing/m6wing.html).
  Published planform and experimental reference context only; no pressure data
  from that experiment are used in this example.
- ONERA/NASA Turbulence Modeling Resource,
  [original M6 geometry](https://turbmodels.larc.nasa.gov/Onerawingnumerics_val/ONERA_M6_Test_Case_TMR.pdf),
  including the 0.8059 m root chord.
