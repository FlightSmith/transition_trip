"""Summarize a saved 1:8 NACA64A204 wing study without rerunning XFOIL.

Input root must contain all_conditions/config_xc_0.07.json, config_xc_0.1.json,
config_xc_0.12.json, boundary_layer.csv and workflow_summary.json. The nine
conditions are Re_MAC=3e6/5e6/8e6 crossed with alpha=0/1/2 degrees. Derived
per-Re data contain only their own condition IDs; no original BL data changes.
"""

from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import html
import json
import math
import os
from pathlib import Path

import numpy as np

from trip_design import run_design, write_outputs


POSITIONS = (.07, .10, .12)
REYNOLDS = (3_000_000., 5_000_000., 8_000_000.)


def _position_label(position):
    return f"xc_{position:.6f}".rstrip("0").rstrip(".")


def _read_json(path):
    return Path(path).read_text(encoding="utf-8-sig")


def _write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n",
                          encoding="utf-8", newline="\n")


def _write_csv(path, fields, rows):
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _condition_signature(config):
    return sorted((str(c["id"]), float(c["reynolds"]), float(c["alpha_deg"]), float(c["mach"]))
                  for c in config["conditions"])


def _validate_inputs(configs, records):
    first = configs[POSITIONS[0]]
    signature = _condition_signature(first)
    expected = {(re, alpha) for re in REYNOLDS for alpha in (0., 1., 2.)}
    actual = {(re, alpha) for _, re, alpha, _ in signature}
    if len(signature) != 9 or len({s[0] for s in signature}) != 9 or actual != expected:
        raise ValueError("Expected nine unique conditions: Re_MAC=3e6/5e6/8e6 crossed with alpha=0/1/2 deg")
    if any(not math.isclose(mach, .2, abs_tol=1e-12) for _, _, _, mach in signature):
        raise ValueError("This scaled-wing comparison expects Mach 0.20 for every condition")
    stations = sorted(first["wing"]["stations"], key=lambda s: float(s["eta"]))
    if len(stations) != 21 or not np.allclose([s["eta"] for s in stations], np.linspace(0, 1, 21), rtol=0, atol=1e-7):
        raise ValueError("Expected 21 span stations at eta=0,0.05,...,1")
    if not (math.isclose(float(stations[0]["chord_m"]), .695, rel_tol=0, abs_tol=1e-6)
            and math.isclose(float(stations[-1]["chord_m"]), .125, rel_tol=0, abs_tol=1e-6)
            and math.isclose(float(first["wing"]["reference_length_m"]), .578625, rel_tol=0, abs_tol=1e-9)):
        raise ValueError("Expected 1:8 model root/tip chords .695/.125 m and reference MAC .578625 m")
    for position, config in configs.items():
        if _condition_signature(config) != signature or config["wing"] != first["wing"]:
            raise ValueError("All position configurations must use identical geometry and conditions")
        if config["analysis"]["source"] != "integral":
            raise ValueError("This report requires saved integral boundary-layer data")
        if config["trip"]["position_xc"] != {"upper": position, "lower": position}:
            raise ValueError("Configuration positions must match their requested upper/lower x/c")
        if ({k: v for k, v in config["trip"].items() if k != "position_xc"}
                != {k: v for k, v in first["trip"].items() if k != "position_xc"}):
            raise ValueError("All positions must use the same trip calibration and sizing policy")
    known_ids = {entry[0] for entry in signature}
    if any(str(row["condition_id"]) not in known_ids for row in records):
        raise ValueError("Boundary-layer file contains an unknown condition ID")
    missing_ids = known_ids - {str(row["condition_id"]) for row in records}
    if missing_ids:
        raise ValueError(f"Boundary-layer file omits configured conditions: {sorted(missing_ids)}")
    return first


def _sample_record(designs, eta):
    matches = [row for row in designs if abs(float(row["eta"]) - eta) < 1e-7]
    if len(matches) != 1:
        return {"height_mm": None, "status": "not_sampled", "reason": f"eta={eta:g} is not uniquely sampled", "y_span_m": None}
    row = matches[0]
    return dict(row, height_mm=row["height_mm"] if row["status"] == "feasible" else None)


def _comparison_rows(scope, reynolds, position, summary, conditions_count):
    rows = []
    uniform = summary["uniform_height"]
    uniform_height = uniform.get("height_mm") if uniform["status"] == "feasible" and summary["status"] == "feasible" else None
    for surface in ("upper", "lower"):
        designs = sorted((d for d in summary["design"] if d["surface"] == surface), key=lambda d: d["eta"])
        valid = [d for d in designs if d["status"] == "feasible" and d["height_mm"] is not None]
        row = {
            "scope": scope, "reynolds_mac": reynolds, "position_xc": position,
            "surface": surface, "conditions_count": conditions_count,
            "valid_station_count": len(valid), "total_station_count": len(designs),
            "full_surface_coverage": len(valid) == len(designs) and bool(designs),
            "min_valid_height_mm": min((d["height_mm"] for d in valid), default=None),
            "max_valid_height_mm": max((d["height_mm"] for d in valid), default=None),
            "valid_eta_values": ";".join(f"{d['eta']:.8g}" for d in valid),
            "uniform_status": uniform["status"],
            "uniform_all_stations_height_mm": uniform_height,
            "uniform_reason": uniform.get("reason", "") if uniform_height is None else "",
            "detail_report": f"{scope}/{_position_label(position)}/report.html",
        }
        for label, eta in (("root", 0.), ("mid", .5), ("tip", 1.)):
            sampled = _sample_record(designs, eta)
            row.update({f"{label}_eta": eta, f"{label}_y_m": sampled.get("y_span_m"),
                        f"{label}_height_mm": sampled["height_mm"],
                        f"{label}_status": sampled["status"], f"{label}_reason": sampled.get("reason", "")})
        rows.append(row)
    return rows


def _plot_comparison(output, summaries):
    os.environ.setdefault("MPLCONFIGDIR", str((output / ".plot_cache").resolve()))
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return False
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(1, 3, figsize=(15.6, 5.4), sharey=True)
    colors = {3_000_000.: "#176b9b", 5_000_000.: "#c86620", 8_000_000.: "#27805b"}
    maximum = 0.
    for ax, position in zip(axes, POSITIONS):
        counts = []
        for reynolds in REYNOLDS:
            scope = f"re{int(reynolds / 1e6)}e6"
            summary = summaries[(scope, position)]
            pair_counts = []
            for surface, style in (("upper", "-"), ("lower", "--")):
                designs = sorted((d for d in summary["design"] if d["surface"] == surface), key=lambda d: d["eta"])
                x = [d["eta"] for d in designs]
                y = [d["height_mm"] if d["status"] == "feasible" and d["height_mm"] is not None else np.nan for d in designs]
                finite = [v for v in y if np.isfinite(v)]
                maximum = max(maximum, max(finite, default=0.))
                pair_counts.append(len(finite))
                ax.plot(x, y, linestyle=style, marker="o" if surface == "upper" else "s",
                        markersize=3.4, linewidth=1.6, color=colors[reynolds],
                        label=f"Re {reynolds/1e6:g}M, {surface}")
            counts.append(f"{reynolds/1e6:g}M: {pair_counts[0]}/{pair_counts[1]}")
        ax.set(title=f"Trip at {100*position:g}% chord", xlabel="Span fraction eta (root to tip)", xlim=(0, 1))
        ax.set_xticks(np.linspace(0, 1, 6))
        ax.grid(alpha=.18)
        ax.text(.03, .96, "Valid upper/lower stations (of 21)\n" + " | ".join(counts),
                transform=ax.transAxes, va="top", fontsize=8,
                bbox={"facecolor": "white", "edgecolor": "#dae1e7", "alpha": .9, "pad": 5})
    axes[0].set_ylabel("Exposed trip height k (mm)")
    axes[0].set_ylim(0, maximum * 1.23 if maximum > 0 else 1.)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", bbox_to_anchor=(.5, .075), ncol=3, frameon=False, fontsize=9)
    first_summary = summaries[("re3e6", POSITIONS[0])]
    calibration = first_summary["calibration"]
    fig.suptitle("NACA64A204 wing | 1:8 model | Mach 0.20", fontsize=16, fontweight="bold", y=.98)
    fig.text(.5, .918, "At each station: one height across alpha = 0, 1, 2 degrees for the indicated Reynolds number", ha="center", fontsize=10)
    fig.text(.5, .018, "Line gaps are invalid/infeasible stations. Valid portions do not establish full-span coverage.\n"
             f"Re_k = {calibration['critical_re_k']:g}; minimum Re_x = {calibration['min_re_x']:,.0f}; 2D section BL with approximate profile reconstruction.",
             ha="center", fontsize=9, color="#43525c")
    fig.subplots_adjust(top=.84, bottom=.28, left=.07, right=.985, wspace=.1)
    fig.savefig(output / "comparison.png", dpi=190)
    fig.savefig(output / "comparison.svg")
    plt.close(fig)
    return True


def _report_html(output, rows, metadata, chart_created, sizing_summary):
    escape = lambda value: html.escape(str(value))
    def number(value):
        return "—" if value is None else f"{value:.6f}"
    tables = []
    for scope, label in (("re3e6", "Re_MAC=3,000,000; alpha=0,1,2°"),
                         ("re5e6", "Re_MAC=5,000,000; alpha=0,1,2°"),
                         ("re8e6", "Re_MAC=8,000,000; alpha=0,1,2°"),
                         ("all_conditions", "Combined envelope: all nine conditions")):
        body = []
        for row in (r for r in rows if r["scope"] == scope):
            samples = []
            for location in ("root", "mid", "tip"):
                value = row[f"{location}_height_mm"]
                text = number(value) if value is not None else row[f"{location}_status"]
                samples.append(f'<td title="{escape(row[f"{location}_reason"])}">{escape(text)}</td>')
            valid_range = f"{number(row['min_valid_height_mm'])}–{number(row['max_valid_height_mm'])}"
            uniform = number(row["uniform_all_stations_height_mm"]) if row["uniform_all_stations_height_mm"] is not None else row["uniform_status"]
            body.append(f'<tr><td><a href="{escape(row["detail_report"])}">{100*row["position_xc"]:g}%</a></td>'
                        f'<td>{escape(row["surface"])}</td><td>{row["valid_station_count"]}/{row["total_station_count"]}</td>'
                        + ''.join(samples) + f'<td>{valid_range}</td><td>{escape(uniform)}</td></tr>')
        tables.append(f'<h2>{escape(label)}</h2><div class="scroll"><table><thead><tr>'
                      '<th>x/c</th><th>Surface</th><th>Valid stations</th><th>Root k</th><th>Midspan k</th><th>Tip k</th>'
                      '<th>Range on valid stations only</th><th>Uniform k, both surfaces/all stations</th>'
                      '</tr></thead><tbody>' + ''.join(body) + '</tbody></table></div>')
    diag = metadata.get("diagnostics", [])
    converged = sum(d.get("status") == "converged" for d in diag)
    failed = sum(d.get("status") == "failed" for d in diag)
    calibration = sizing_summary["calibration"]
    chart = '<img src="comparison.png" alt="Three chord positions comparing trip height by Reynolds number and surface; lines stop at invalid stations">' if chart_created else '<p>Comparison chart not generated. Install matplotlib and rerun this report with plotting enabled.</p>'
    document = '''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>1:8 swept-wing trip comparison</title><style>
body{font:16px/1.5 system-ui,sans-serif;max-width:1440px;margin:35px auto;padding:0 24px;color:#173447}
h1{font-size:30px;margin-bottom:8px}h2{margin-top:30px;font-size:21px}.note{background:#edf4f7;border-left:4px solid #237ba0;padding:16px 20px}
.warning{background:#fff4e6;border-left-color:#bd7428}.scroll{overflow-x:auto}table{border-collapse:collapse;width:100%;font-size:13px}
th,td{text-align:left;border-bottom:1px solid #d2dde3;padding:9px 8px}th{background:#f3f6f8}img{max-width:100%;height:auto}
a{color:#176b9b}.meta{color:#52616c}code{background:#eef3f6;padding:2px 4px}</style>
<h1>1:8 NACA64A204 swept-wing trip comparison</h1>
<p class="meta">Mach 0.20 · angles 0°, 1°, 2° · Re_MAC 3, 5, 8 million · trip positions 7%, 10%, 12% chord</p>
<p class="note"><strong>Model geometry:</strong> root/tip chord 0.695/0.125 m; reference MAC 0.578625 m. Linear dimensions use 1:8 scaling.
Heights in every table are millimetres. The plotted span covers the 21 supplied stations of the left panel.</p>
<p class="note warning"><strong>Partial coverage is not a full-span design.</strong> Blank heights and line gaps indicate failed applicability or sizing checks.
Ranges use valid stations only. A uniform height is reported only if one height passes every sampled station on both surfaces across the stated condition group.
Taking the largest height from a partial set does not establish a usable uniform trip.</p>
''' + chart + ''.join(tables) + f'''
<h2>Interpretation and data</h2>
<p>The selected trip calibration uses Re_k=u(k)k/nu(k)={calibration['critical_re_k']:g} and requires running Re_x≥{calibration['min_re_x']:,.0f}.
The installed-height multiplier is {sizing_summary['height_margin']:.2f} and the cap is k/delta99≤{sizing_summary['max_k_delta99']:.2f}.
These are preliminary sizing assumptions. XFOIL supplies two-dimensional section boundary-layer integrals; the velocity profile is reconstructed approximately.
Sweep-induced crossflow and three-dimensional transition are not modeled. Solver convergence does not by itself establish an applicable trip design.</p>
<p>Recorded original solver outcomes: {converged} converged, {failed} failed. This comparison reruns sizing from saved integral data; it does not rerun XFOIL.</p>
<p><code>all_conditions/</code> contains the original nine-condition workflow, geometry, solver evidence and shared boundary-layer records.
<code>re3e6/</code>, <code>re5e6/</code> and <code>re8e6/</code> contain filtered copies of those records, three-condition configurations and derived sizing reports.
The supplied Reynolds numbers are referenced to the scaled MAC; they are not divided by 8 again.</p>
<p><a href="comparison.csv">Comparison CSV</a> · <a href="comparison_summary.json">Report metadata</a> ·
<a href="all_conditions/workflow_summary.json">Original workflow summary</a></p></html>'''
    (output / "report.html").write_text(document, encoding="utf-8", newline="\n")


def summarize(output_root, plot=True):
    """Create per-Re and combined-envelope reports from the completed saved study."""
    output = Path(output_root).resolve()
    source = output / "all_conditions"
    configs = {position: json.loads(_read_json(source / f"config_{_position_label(position)}.json")) for position in POSITIONS}
    data_path = source / "boundary_layer.csv"
    for config in configs.values():
        if (source / config["analysis"].get("integral_file", "")).resolve() != data_path.resolve():
            raise ValueError("Every source configuration must reference all_conditions/boundary_layer.csv")
    with data_path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        fields = reader.fieldnames
        records = list(reader)
    if not fields or not records or "condition_id" not in fields:
        raise ValueError("Saved boundary-layer CSV is missing its header or data")
    first = _validate_inputs(configs, records)
    workflow_path = source / "workflow_summary.json"
    workflow = json.loads(_read_json(workflow_path))
    summaries, comparison, group_outcomes = {}, [], []
    for reynolds in (*REYNOLDS, "all"):
        scope = f"re{int(reynolds / 1e6)}e6" if reynolds != "all" else "all_conditions"
        directory = output / scope
        directory.mkdir(exist_ok=True)
        selected_conditions = [c for c in first["conditions"] if reynolds == "all" or float(c["reynolds"]) == reynolds]
        ids = {str(c["id"]) for c in selected_conditions}
        if reynolds != "all":
            filtered = [row for row in records if str(row["condition_id"]) in ids]
            _write_csv(directory / "boundary_layer.csv", fields, filtered)
        for position in POSITIONS:
            label = _position_label(position)
            config_path = directory / f"config_{label}.json"
            if reynolds != "all":
                cfg = copy.deepcopy(configs[position])
                cfg["conditions"] = selected_conditions
                cfg["analysis"] = dict(cfg["analysis"], source="integral", integral_file="boundary_layer.csv")
                cfg["name"] = f"1:8 NACA64A204 wing; Re_MAC={reynolds/1e6:g}M; {100*position:g}% chord"
                cfg["data_origin"] = (str(cfg.get("data_origin", "")) + " Filtered from the saved nine-condition boundary-layer dataset; no solver rerun.")
                _write_json(config_path, cfg)
            cfg, summary, condition_rows = run_design(config_path)
            write_outputs(directory / label, cfg, summary, condition_rows, plot=plot)
            summaries[(scope, position)] = summary
            comparison.extend(_comparison_rows(scope, reynolds, position, summary, len(selected_conditions)))
            group_outcomes.append({"scope": scope, "position_xc": position, "status": summary["status"],
                                   "uniform_height": summary["uniform_height"], "condition_count": len(selected_conditions)})
    csv_fields = list(dict.fromkeys(key for row in comparison for key in row))
    _write_csv(output / "comparison.csv", csv_fields, comparison)
    chart_created = _plot_comparison(output, summaries) if plot else False
    sources = [data_path, workflow_path, *(source / f"config_{_position_label(x)}.json" for x in POSITIONS)]
    audit = {
        "description": "Derived sizing summaries from saved BL integrals; no XFOIL execution",
        "scale": "1:8", "reference_mac_m": .578625, "root_chord_m": .695, "tip_chord_m": .125,
        "chart_created": chart_created, "comparison_rows": len(comparison), "groups": group_outcomes,
        "uniform_definition": "One common height passing all supplied stations on both surfaces and all selected conditions",
        "range_definition": "Minimum/maximum over valid stations only; not a full-panel design if coverage is partial",
        "inputs": [{"path": path.relative_to(output).as_posix(), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()} for path in sources],
    }
    _write_json(output / "comparison_summary.json", audit)
    _report_html(output, comparison, workflow, chart_created, summaries[("re3e6", POSITIONS[0])])
    return {"output": str(output), "report": str(output / "report.html"),
            "comparison_csv": str(output / "comparison.csv"), "comparison_rows": len(comparison),
            "chart_created": chart_created, "groups": group_outcomes}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output_root", type=Path, help="Study folder containing all_conditions/")
    parser.add_argument("--no-plot", action="store_true")
    args = parser.parse_args(argv)
    result = summarize(args.output_root, plot=not args.no_plot)
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
