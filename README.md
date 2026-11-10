# Wind-tunnel trip sizing from an undisturbed laminar boundary layer

This Python tool sizes a trip at a specified fraction of local chord on the upper
and lower surfaces of a wing. It evaluates every supplied tunnel condition and
angle of attack, then selects a common exposed height at each spanwise station.
The example uses **10% chord on both surfaces**, as requested.

The calculation is a preliminary engineering estimate. It evaluates a roughness
Reynolds-number criterion. The Cp and supplied-profile routes do not predict
natural transition; the optional XFOIL route uses its own transition prediction
to check that the proposed trip lies in a laminar region. No route predicts the
distance from the trip to complete transition or a particular trip's drag penalty.
Verify the selected trip in the tunnel using an appropriate transition diagnostic.

## Run after downloading from GitHub

No Codex installation is required. On a new Windows computer:

1. Install **Python 3.10 or newer** from the
   [official Windows downloads page](https://www.python.org/downloads/windows/).
   Use a normal Python installation with the `py` launcher or `python` on PATH.
2. Open [the GitHub repository](https://github.com/FlightSmith/transition_trip),
   choose **Code → Download ZIP**, and extract it. Work in the extracted folder
   containing `setup_windows.cmd`, not inside the ZIP.
3. Double-click **`setup_windows.cmd`** and wait for **Setup complete**. This
   creates `.venv` and installs the Python dependencies in the project folder.
4. Download the [official Windows XFOIL 6.99 ZIP](https://web.mit.edu/drela/Public/web/xfoil/XFOIL6.99.zip)
   and extract it so the executable is, for example,
   `C:\tools\XFOIL\xfoil.exe`. Adjust the path below if yours differs.

Open PowerShell in the extracted project folder. First check XFOIL:

```powershell
& "C:\tools\XFOIL\xfoil.exe"
```

At the XFOIL prompt, type **`QUIT`** and press Enter to return to PowerShell.
Then set its location for this PowerShell window and run the NACA 64A204 study:

```powershell
$env:XFOIL_EXE = "C:\tools\XFOIL\xfoil.exe"
.\.venv\Scripts\python.exe run_scaled_naca64a204.py
```

Open `results\left_wing_naca64a204_scale8\report.html` when the run finishes.
The environment-variable setting lasts only for the current PowerShell window;
set it again in a new window, or pass `--xfoil-exe "C:\tools\XFOIL\xfoil.exe"`
on the run command. XFOIL is needed only for solver-based runs, not the included
Cp/profile examples.

For the other workflows, double-click `run_windows.cmd` to choose:

1. M6 synthetic Cp example.
2. Left-wing Blasius example at 7%, 10%, and 12% chord.
3. Left-wing Cp-input demonstration at those three positions.
4. Your own Cp CSV with the supplied left-wing geometry and tunnel condition.
5. Your own JSON configuration, using either Cp or laminar-profile input.
6. A folder of section files and matching Cp or pressure curves.
7. A folder of section files analyzed with an external XFOIL executable.

**macOS/Linux:** run these commands from the extracted project directory:

```sh
python3 setup_environment.py
.venv/bin/python run_trip.py
```

Setup requires an internet connection to install dependencies. Calculations
then run locally. Do not copy `.venv` between computers; run setup on each one.
The setup and menu use project-relative locations and do not require the
original computer's paths. Existing example inputs are included in GitHub;
regenerating them is optional and replaces the selected example files.

## Rerun with a different airfoil

Provide **one normalized two-dimensional airfoil shape**; it is applied uniformly
across the same 1:8 wing with zero twist. XFOIL calculates the new pressures and
boundary layers, so **no existing Cp file is required**. The Mach, Reynolds
numbers, angles, span stations, and trip positions stay as in the scaled study.

For a ready example, run the included NACA 2412 file after the setup above:

```powershell
.\.venv\Scripts\python.exe run_scaled_naca64a204.py --airfoil .\examples\airfoils\naca2412.dat --output .\results\naca2412_scale8
```

For your own file, change just its path and use a separate output folder:

```powershell
.\.venv\Scripts\python.exe run_scaled_naca64a204.py --airfoil "C:\airfoils\my_airfoil.dat" --output .\results\my_airfoil_scale8
```

Open `report.html` in that output folder. Despite the script's historical name,
the calculation uses the supplied shape. Report titles use the file's first
line; add `--airfoil-name "My airfoil"` to override the title. Custom inputs
are copied under the selected output's `inputs/` folder. If `--output` is
omitted, the default is `results/left_wing_scale8_<file-stem>/`; the original
NACA 64A204 result folder is protected from custom-airfoil runs. Reusing a
custom output folder replaces its generated inputs and reports.

The `.dat` format is a **single name line**, followed by at least **10 pairs**
of whitespace-separated `x/c z/c` coordinates. Use a unit chord, with the
contour ordered **trailing edge → upper surface → leading edge → lower surface
→ trailing edge**. Do not include a point-count line or the seven-line pressure
curve header. Do **not** divide the coordinates by 8: the script applies the
model's local chords automatically. The bundled
[`naca2412.dat`](examples/airfoils/naca2412.dat) is a 35-point Selig-format
example from the [UIUC Airfoil Coordinates Database](https://m-selig.ae.illinois.edu/ads/coord/naca2412.dat).

Use `--prepare-only` to inspect the generated geometry without running XFOIL.
To rebuild reports later from saved solver records, select the existing output
and omit `--airfoil`:

```powershell
.\.venv\Scripts\python.exe run_scaled_naca64a204.py --summarize-only --output .\results\my_airfoil_scale8
```

## Section files with Cp or pressure curves

The supplied `synthetic_Cp/` folder can be used directly. It contains two
section shapes and three pressure-coefficient states at each section. Its
ready-to-run settings are in `examples/section_curves/workflow.json`:

```powershell
.\run_windows.cmd --section-config .\examples\section_curves\workflow.json --output .\results\section_curves
```

On macOS/Linux, use `.venv/bin/python run_trip.py` instead of
`.\run_windows.cmd`. The workflow also runs directly as
`python section_workflow.py --config WORKFLOW.json --output OUTPUT_DIR`.
The section workflow JSON is a different input format from the ordinary
`trip_design.py` sizing configuration.

The confirmed example conditions are Mach **0.2**, **Re_MAC = 8 million**,
and **MAC = 4.629 m**, with state 1, 2, and 3 mapped to **0, 1, and 2 degrees**.
The default positions are **7%, 10%, and 12% local chord**, on both surfaces.
Static temperature defaults to **288.15 K** and is editable in the workflow
JSON. State IDs are labels; the angles come from the explicit `states` list.

Choose menu option **6** to enter a folder, tunnel conditions, three state
angles, trip positions, pressure type, and initialization method interactively.
The menu saves a reusable `workflow.json` in the chosen output directory.
For different state IDs/counts or per-state tunnel conditions, edit a copy of
the example workflow JSON instead of using the three-state wizard.

### Raw file format and geometry

- Geometry filenames follow `section_state1_station-1.250`; the two columns
  are dimensional **X, Z in meters**. Pressure filenames follow
  `cp_dist_state2_station-1.250`; their columns are dimensional **X** and **Cp
  or pressure in Pa**, according to `value_kind`. The station suffix is the
  signed global **Y coordinate in meters**. A text extension is optional.
- Numeric data begin on **line 8** by default (`data_start_line` is one-based).
  Each file must contain one two-column curve followed by `*END`. E/D exponent
  notation and `$`/`#` comments are accepted. Header scaling/offset directives
  are not applied: export already-scaled physical numbers. The supplied files
  use identity scaling and zero offsets.
- Keep the ordered **TE → LE → TE contour**. The first branch may be upper or
  lower, but each pressure file must use the same branch convention as its
  section. Geometry determines the upper/lower labels. Do not globally sort X;
  that would mix the two surfaces. Different pressure sampling grids are mapped
  separately on each branch without extrapolating beyond the section geometry.
- `section_state: 1` reuses the state-1 geometry for every pressure state.
  If additional state-specific section files are present, their geometry must
  match; changing geometry requires separate workflows.

The actual supplied sections have `X_LE=2.1875 m`, chord `1.75 m` at
`Y=-1.250 m`, and `X_LE=2.75 m`, chord `1.0 m` at `Y=-5.000 m`. Their vertical
offsets are `0.025 m` and `0.100 m`. These are **different geometry from the
earlier left-wing template**; the importer derives geometry from these files.
Its internal `eta` runs from 0 at the innermost supplied station to 1 at the
outermost supplied station. It does not infer an unprovided full-wing root,
tip, or semispan, and it preserves the signed Y values in the output.

The importer splits the surfaces and computes physical contour arc lengths.
For XFOIL it translates the leading edge to the origin, removes section chord
incidence, and divides both coordinates by the physical LE-to-TE chord length.
The normalized airfoil files are written to `airfoils/` for inspection. The
default `alpha_reference: "section_chord"` treats each specified angle as
relative to that section's chord. With `"global_x"`, the recorded geometric
chord incidence is subtracted before calling XFOIL.

### Pressure conversion and initial thickness

Set `value_kind` to `"cp"` for coefficients. For absolute pressure, use
`"pressure"` or `"absolute_pressure"`; for signed gauge values, use
`"gauge_pressure"`. Pressure conversion requires both reference values:

```text
Cp = (input_pressure - static_pressure_pa) / dynamic_pressure_pa
```

Put `dynamic_pressure_pa` and `static_pressure_pa` in the workflow's global
`pressure_reference` object, or put those keys directly in a particular
`states` entry to override them for that state. Dynamic pressure must be
positive. Samples and static reference must use the **same absolute or gauge
convention**; absolute pressures must be positive. CLI options `--value-kind`,
`--dynamic-pressure-pa`, and `--static-pressure-pa` override the saved settings;
CLI pressure references override both global and per-state values.

These reference pressures only normalize the supplied pressure data. The
dimensional tunnel state is still determined by Mach, Re_MAC, MAC, and static
temperature; pressure references do not replace that state automatically.
Use mutually consistent tunnel measurements and settings.

Cp integration requires an explicit `cp_initialization` choice:

- `{"method": "synthetic_le_zero"}` starts with zero momentum thickness at
  the geometric LE. This is **only the supplied synthetic demonstration
  assumption**, not a valid generic initialization for a real wing. Its Cp=1
  LE values imply finite speed under compressible inversion at Mach 0.2;
  geometric LE and stagnation must not be silently equated.
- `{"method": "stagnation_ramp"}` starts at the first positive surface
  distance with `theta0 = sqrt(0.075 * nu_inf * s0 / Ue0)`. It explicitly
  assumes a linear velocity rise from a stagnation origin to that point.
  It is an approximation requiring a justified origin and upstream flow
  history, particularly at nonzero incidence. The wizard defaults to this
  choice for folders other than the supplied synthetic example.
- `"method": "provided"` additionally requires a positive `theta0_m` in
  the same object, representing the initial momentum thickness at the first
  retained downstream sample. Supply a physically justified value; missing
  upstream history is not zero. A state may override `cp_initialization`.

The true arc distance is retained when the LE point is skipped. Initialization
method, initial position, distance, and thickness are recorded in
`workflow_summary.json`. Neither the conversion nor the geometry import
establishes that a real boundary layer is laminar.

## Optional XFOIL section analysis

Install XFOIL separately; it is **not bundled or installed by the Python
setup**. The [official MIT XFOIL page](https://web.mit.edu/drela/Public/web/xfoil/)
provides a [Windows 6.99 archive](https://web.mit.edu/drela/Public/web/xfoil/XFOIL6.99.zip)
and source distributions. Extract the Windows package and point the tool to
its executable. On macOS/Linux, use a compatible native installation/build;
the Windows executable is not a native Unix program.

```powershell
.\run_windows.cmd --section-config .\examples\section_curves\workflow.json --mode xfoil --xfoil-exe "C:\tools\XFOIL\xfoil.exe" --output .\results\section_xfoil
```

Alternatively choose menu option **7**. The executable can also be selected
with `XFOIL_EXE`, the workflow's `xfoil.executable`, or `xfoil` on PATH, in
that precedence after `--xfoil-exe`. The example uses `ncrit=9`, 200 panels,
200 iterations, and a 60-second timeout per section/state; these are editable
in the workflow's `xfoil` object. Run with `--prepare-only` to inspect normalized
airfoils and the solver plan without executing XFOIL or producing trip designs.

Each section uses its own chord Reynolds number:

```text
Re_local = Re_MAC * section_chord / reference_mac_m
```

XFOIL computes its **own two-dimensional Cp and integral boundary layer from
the section geometry**. It does not prescribe the imported wing Cp to its
viscous solver. The baseline permits natural transition; it does not force
transition at the requested 7%, 10%, or 12% trip position. A usable result
requires convergence, positive attached-flow diagnostics, and a sampling
bracket upstream of the predicted transition. Failures are reported as invalid
cases, without silently substituting the Cp-based model.

XFOIL's DUMP file contains integral quantities, **not a resolved wall-normal
velocity profile**. The tool checks H/Hk format differences, obtains kinematic
shape factor Hk, and reconstructs a quartic profile matched to Hk and momentum
thickness. Transferring the solver thickness and using constant freestream
kinematic viscosity are low-Mach approximations; the dump does not supply
`nu(y)` or a wall-temperature profile. This route remains limited to freestream
and sampled local Mach at most 0.30 and the supported attached quartic shape
family. Laminar status is predicted by XFOIL, not established experimentally.
The [official user primer](https://web.mit.edu/drela/Public/web/xfoil/xfoil_doc.txt)
describes its integral formulation and transition model.

Both section workflows treat the supplied sections independently. They do not
solve swept-wing crossflow, attachment-line effects, a full 3D pressure field,
or unsampled angles/stations. XFOIL and imported-Cp results are different model
paths and should not be presented as independent validation of one another.

### Section workflow outputs

The output folder contains `workflow_summary.json` with geometry, input hashes,
initialization/solver diagnostics, and the effective settings. Curve mode writes
converted `cp.csv`; XFOIL mode writes `boundary_layer.csv`, `xfoil_cp.csv`, and
the individual solver command/output files under `xfoil_runs/`.

Each requested position has a generated sizing configuration and a report
directory such as `xc_0.07/`, `xc_0.1/`, or `xc_0.12/`, containing `report.html`,
`design.csv`, `conditions.csv`, `summary.json`, and optional PNG/SVG figures.
Use `--no-plot` to skip figures. Check the per-position design status and
condition-level reasons; the workflow returns a nonzero exit code if the full
requested envelope is not feasible. Running again into the same output folder
replaces generated tables/reports, so use separate folders to retain runs.

## Original left wing at 1:8 scale with NACA 64A204

`run_scaled_naca64a204.py` prepares and analyzes the original left-wing
planform at **1:8 scale**, using the same **NACA 64A204** section everywhere
with **zero twist**. It evaluates **Mach 0.2**, **Re_MAC = 3, 5, and 8 million**,
and **0, 1, and 2 degrees** incidence, at **7%, 10%, and 12% local chord** on
both surfaces. The 21 span stations and nine tunnel/AoA combinations define
**189 XFOIL operating points**; the three trip positions reuse those solutions.
This thin-airfoil case starts with 280 panels and 200 iterations. Failed solves
are retried with 600 iterations at 280, then 240, then 200 panels; every attempt
is retained. Converged solutions are never replaced to obtain a preferred trip
result. Intermediate boundary-layer march warnings are recorded; the final
coupled convergence and a matching fresh viscous polar point are required.

From the project folder on Windows, after setup and installing XFOIL:

```powershell
.\.venv\Scripts\python.exe run_scaled_naca64a204.py --xfoil-exe "C:\tools\XFOIL\xfoil.exe"
```

On macOS/Linux, use `.venv/bin/python` and a native XFOIL executable. Alternatively
set `XFOIL_EXE` or make `xfoil` available on PATH. Useful rerun options are:

- `--prepare-only`: regenerate model sections, settings, and the solver plan
  without executing XFOIL or producing trip recommendations.
- `--summarize-only`: rebuild the comparison and per-Re reports from saved
  solver records; it does not rerun XFOIL.
- `--no-plot`: skip figures. `--output PATH` selects a separate result folder.

The default comparison is
`results/left_wing_naca64a204_scale8/report.html`, with `comparison.csv` and
`comparison_summary.json` beside it. `re3e6/`, `re5e6/`, and `re8e6/` contain
separate Reynolds-number envelopes over the three angles; `all_conditions/`
contains the combined nine-condition envelope, shared solver evidence, and
per-position reports. A complete report may still describe an infeasible
envelope and return exit code 2; inspect the rejection reasons rather than
treating missing heights as zero.

Every original dimensional coordinate and chord is divided by 8, and reference
area is divided by 64. The model has **MAC 0.578625 m**, **root chord 0.695 m**,
**tip chord 0.125 m**, and **reference area 0.765625 m²**. The specified Reynolds
numbers apply to this **model MAC** and are **not divided by 8**. Local section
Reynolds numbers follow `Re_local = Re_MAC * local_chord / model_MAC`. The
reference area is retained as geometry metadata and does not set trip height.

The 401-point airfoil was generated from the PDAS NACA456 source algorithm at
its native mapping points, with the 64A thickness family and modified meanline
for NACA 64A204. Comparison with the published coordinate table differed by at
most approximately **5.2×10⁻⁷ chord**; this checks geometry, not aerodynamic
accuracy. See the recorded source, method, limitations, and hashes in
[`airfoil_source.json`](examples/left_wing_naca64a204_scale8/airfoil_source.json).
`python generate_naca64a204.py --out DIRECTORY` reproduces the source airfoil
and its verification from the bundled mapping arrays and published table,
without a network connection or Fortran compiler.
The normalized coordinates, generated sections, workflow, and model geometry
are retained in `examples/left_wing_naca64a204_scale8/`. Angles refer to the
nominal NACA chord/global X; the small geometric-LE chord rotation is corrected
when creating the local XFOIL angle, so it is not introduced as spanwise twist.

The default grit criterion requires `Re_x >= 100000`. At Re_MAC = 3 million,
the tip does not meet that requirement at any of these three positions. Thus
the default criterion cannot provide a full-span uniform height covering all
three Reynolds numbers; this is a calibration-applicability limit, not proof
that the tip cannot be tripped. Predicted natural transition can independently
exclude stations, particularly on the upper surface. Retain those invalid
stations in the comparison. The same low-Mach integral/profile reconstruction
limits apply, and independent 2D XFOIL sections do not solve the swept wing's
three-dimensional flow.

## Use your own Cp distribution

The ready-to-copy left-wing templates are in `examples/left_wing_cp/`.
They retain MAC=4.629 m, Mach=0.2, Re_MAC=8 million, signed left-wing coordinates,
and the three requested trip locations. **Their Cp values are all zero and are
only a format demonstration, not actual wing pressures.**

Copy that folder to a separate case folder. Update its `cp.csv`, the condition
metadata, and `data_origin` in the JSON files. Then run a configuration directly:

```powershell
.\run_windows.cmd --config .\my_case\config_xc10.json --output .\results\my_case
```

You can instead supply a CSV on the command line while keeping a configuration
as the geometry/tunnel template:

```powershell
.\run_windows.cmd --config .\examples\left_wing_cp\config_xc10.json --cp-file "C:\test-data\my_cp.csv" --position-xc 0.10 --output .\results\my_cp_run
```

On macOS/Linux, replace `.\run_windows.cmd` with `.venv/bin/python run_trip.py`
and use the appropriate paths. The same options work with `trip_design.py`.
`--cp-file` explicitly selects Cp mode; `--profiles-file` selects the profile
mode. They are mutually exclusive. `--position-xc` changes the position for
every surface listed in the configuration. Omit it to retain separate upper
and lower positions. `--no-plot` skips charts. Input files are never edited by
these overrides, and `summary.json` records the overrides and actual data hashes.

Command-line paths are relative to the current working directory; file paths
inside JSON are relative to that JSON file's directory. Quote paths with spaces.

For real Cp data:

- Match every `condition_id`, `eta`, and `surface` to the JSON configuration.
  Set the actual angle of attack, Mach, reference Reynolds number and static
  temperature; changing labels does not calculate a new pressure distribution.
- The default left-wing template contains 21 stations and one `case_01`.
  Supply all requested stations or edit `wing.stations` to match the stations
  you actually have. There is no spanwise or angle-of-attack interpolation.
- Use `Cp=(p-p_inf)/(0.5*rho_inf*U_inf^2)` for the matching tunnel condition,
  with upstream chordwise coverage through the trip, not just Cp at the trip.
- Provide physical surface distance `s_m` when available. If the profile starts
  downstream, supply a positive, physically justified `theta0_m` in the first
  row of each series. **Do not retain the demo's zero initial thickness for a
  real downstream start, or reset the distance origin to bypass initialization.**
  The exact stagnation point with zero edge velocity is not handled; start at
  a nonstagnant point with a justified initial boundary-layer thickness.
- Both freestream and upstream local Mach must be at most 0.30 in Cp mode.
  Known transition/separation upstream invalidates the laminar trip estimate.

For the supplied left panel, convert global coordinates using:

```text
eta = (Y + 1.333) / (-4.917)
local_chord = 5.56 - 4.56 * eta
X_LE = 9.280 + 3.277 * eta
x_c = (X - X_LE) / local_chord
```

Thus `eta=0` is the specified root and `eta=1` the tip; it is not `abs(Y)/6.25`.
The full CSV schema is documented below. All outputs are estimates subject to
the validity limits of the chosen boundary-layer and trip models.

## Run the original example directly

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

On the original prepared Windows workspace, `run_example.ps1` runs both M6 examples using
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

The separate `generate_left_wing.py` example uses the user-supplied left panel
from `(X_LE,Y,c)=(9.280,-1.333,5.56)` m to `(12.557,-6.250,1.0)` m,
Mach 0.2 and Reynolds number 8 million on the supplied 4.629 m reference MAC.
It generates idealized zero-pressure-gradient Blasius laminar profiles for
7%, 10% and 12% local-chord positions, and runs through the supplied-profile
input path. The upper and lower estimates coincide solely because the same
flat-plate assumption is applied to both. These are not wing pressure data.

```text
python generate_left_wing.py
python trip_design.py --config examples/left_wing_flatplate/config_xc07.json --output results/left_wing_flatplate/xc07
python trip_design.py --config examples/left_wing_flatplate/config_xc10.json --output results/left_wing_flatplate/xc10
python trip_design.py --config examples/left_wing_flatplate/config_xc12.json --output results/left_wing_flatplate/xc12
```

This case uses continuous heights, with a 10% height margin and the 0.8
`k/delta99` cap; no hardware-size inventory is assumed. The supplied dimensions
are treated as physical model dimensions. The exposed trapezoidal panel does
not define the supplied full-wing reference area or MAC; those reference
quantities are retained as supplied. Flat-plate results omit root junction,
wingtip and sweep/crossflow effects. The 7% uniform-height interval is very
narrow and sensitive to profile assumptions; use the spanwise results as a
preliminary estimate, not a validated installation specification.

The original M6 example is described below.

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
An optional station `y_m` sets its signed global spanwise coordinate, allowing
left-side panels and roots offset from the centerline. Without this override,
the spanwise coordinate remains `eta * semispan_m`. In a panel-only analysis,
`semispan_m` is the analyzed root-to-tip span and `eta` runs from that panel's
root to its tip; supply both `y_m` and `x_le_m` to preserve aircraft coordinates.
This station coordinate is distinct from the wall-normal `y_m` in profile CSVs.

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
| `wing.stations` | Objects containing `eta` and streamwise `chord_m`; optional `x_le_m` and signed global `y_m` override placement |
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
temperature, and Mach range. The Cp/profile routes do not predict natural
transition; optional XFOIL provides a model-dependent section prediction.
Crossflow, leading-edge contamination, shocks, and the tunnel's actual
turbulence field are not predicted by this tool.

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
