#!/usr/bin/env python3
"""Run deterministic verification cases and write aggregate/machine artifacts."""

from __future__ import annotations

import csv
import hashlib
import html
import json
import math
from pathlib import Path
from typing import Iterable

import numpy as np

from wood_hill import (
    VOIGT_LABELS,
    finite_difference_tangent,
    load_material,
    run_component_path,
    run_uniaxial_path,
)
from wood_hill.conventions import (
    rotation_about_global_t,
    strain_material_to_global,
    stress_material_to_global,
)
from wood_hill.extensions import extension_roadmap


ROOT = Path(__file__).resolve().parents[1]
PARAMETERS = ROOT / "parameters" / "demo_wood.json"
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"
RESULTS.mkdir(exist_ok=True)
FIGURES.mkdir(exist_ok=True)


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str] | None = None) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty table: {path}")
    names = list(rows[0]) if fieldnames is None else fieldnames
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=names)
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 16), b""):
            digest.update(block)
    return digest.hexdigest()


def path_records(case: str, direction: str, points, component: int, segment: str = "monotonic"):
    records = []
    for point in points:
        result = point.result
        state = point.state
        records.append(
            {
                "case": case,
                "segment": segment,
                "direction": direction,
                "step": point.step,
                "controlled_strain": point.target_strain,
                "conjugate_stress_MPa": state.stress[component],
                "equivalent_stress_MPa": material.hill.equivalent_stress(state.stress),
                "yield_radius_MPa": material.radius(state.kappa),
                "kappa": state.kappa,
                "plastic_multiplier": result.plastic_multiplier,
                "yield_residual_MPa": result.yield_residual,
                "newton_residual_scaled": result.newton_residual,
                "plastic_dissipation_MPa": state.plastic_dissipation,
                "elastic_energy_MPa": material.elastic_energy(state),
                "lateral_stress_norm_MPa": point.lateral_stress_norm,
                "elastic_step": int(result.elastic),
                "return_iterations": result.iterations,
            }
        )
    return records


def summarize_path(case: str, direction: str, points, component: int) -> dict[str, object]:
    plastic = [point for point in points if not point.result.elastic]
    onset = plastic[0].state.stress[component] if plastic else math.nan
    return {
        "case": case,
        "direction": direction,
        "steps": len(points),
        "final_controlled_strain": points[-1].target_strain,
        "final_conjugate_stress_MPa": points[-1].state.stress[component],
        "yield_onset_stress_MPa": onset,
        "final_kappa": points[-1].state.kappa,
        "cumulative_plastic_dissipation_MPa": points[-1].state.plastic_dissipation,
        "minimum_plastic_dissipation_increment_MPa": min(
            point.result.plastic_dissipation_increment for point in points
        ),
        "minimum_elastic_energy_MPa": min(material.elastic_energy(point.state) for point in points),
        "max_abs_plastic_yield_residual_MPa": max(
            (abs(point.result.yield_residual) for point in plastic), default=0.0
        ),
        "max_scaled_newton_residual": max(
            (point.result.newton_residual for point in plastic), default=0.0
        ),
        "max_lateral_stress_norm_MPa": max(point.lateral_stress_norm for point in points),
        "elastic_steps": sum(point.result.elastic for point in points),
        "plastic_steps": sum(not point.result.elastic for point in points),
        "max_return_iterations": max(point.result.iterations for point in points),
    }


def svg_plot(
    path: Path,
    series: list[tuple[str, np.ndarray, np.ndarray, str]],
    title: str,
    x_label: str,
    y_label: str,
    log_x: bool = False,
    log_y: bool = False,
) -> None:
    width, height = 820, 540
    left, right, top, bottom = 92, 28, 58, 78
    plot_w, plot_h = width - left - right, height - top - bottom
    all_x = np.concatenate([entry[1] for entry in series]).astype(float)
    all_y = np.concatenate([entry[2] for entry in series]).astype(float)
    if log_x:
        all_x = np.log10(all_x)
    if log_y:
        all_y = np.log10(all_y)
    x_min, x_max = float(np.min(all_x)), float(np.max(all_x))
    y_min, y_max = float(np.min(all_y)), float(np.max(all_y))
    x_pad = 0.04 * max(x_max - x_min, 1e-12)
    y_pad = 0.08 * max(y_max - y_min, 1e-12)
    x_min, x_max = x_min - x_pad, x_max + x_pad
    y_min, y_max = y_min - y_pad, y_max + y_pad

    def xy(x: float, y: float) -> tuple[float, float]:
        if log_x:
            x = math.log10(x)
        if log_y:
            y = math.log10(y)
        px = left + (x - x_min) / (x_max - x_min) * plot_w
        py = top + (y_max - y) / (y_max - y_min) * plot_h
        return px, py

    lines = [
        '<svg xmlns="http://www.w3.org/2000/svg" width="820" height="540" viewBox="0 0 820 540">',
        '<rect width="100%" height="100%" fill="white"/>',
        f'<text x="410" y="30" text-anchor="middle" font-family="sans-serif" font-size="19">{html.escape(title)}</text>',
        f'<rect x="{left}" y="{top}" width="{plot_w}" height="{plot_h}" fill="none" stroke="#333"/>',
    ]
    for i in range(6):
        fraction = i / 5
        x_value = x_min + fraction * (x_max - x_min)
        y_value = y_min + fraction * (y_max - y_min)
        px = left + fraction * plot_w
        py = top + (1 - fraction) * plot_h
        lines.extend(
            [
                f'<line x1="{px:.2f}" y1="{top}" x2="{px:.2f}" y2="{top + plot_h}" stroke="#e6e6e6"/>',
                f'<line x1="{left}" y1="{py:.2f}" x2="{left + plot_w}" y2="{py:.2f}" stroke="#e6e6e6"/>',
                f'<text x="{px:.2f}" y="{top + plot_h + 23}" text-anchor="middle" font-family="sans-serif" font-size="12">{(10**x_value if log_x else x_value):.3g}</text>',
                f'<text x="{left - 10}" y="{py + 4:.2f}" text-anchor="end" font-family="sans-serif" font-size="12">{(10**y_value if log_y else y_value):.3g}</text>',
            ]
        )
    for label, xs, ys, color in series:
        points = " ".join(f"{px:.2f},{py:.2f}" for px, py in (xy(float(x), float(y)) for x, y in zip(xs, ys)))
        lines.append(
            f'<polyline points="{points}" fill="none" stroke="{color}" stroke-width="2.3" stroke-linejoin="round"/>'
        )
    lines.extend(
        [
            f'<text x="{left + plot_w / 2}" y="{height - 24}" text-anchor="middle" font-family="sans-serif" font-size="15">{html.escape(x_label)}</text>',
            f'<text x="23" y="{top + plot_h / 2}" text-anchor="middle" font-family="sans-serif" font-size="15" transform="rotate(-90 23 {top + plot_h / 2})">{html.escape(y_label)}</text>',
        ]
    )
    legend_x, legend_y = left + 14, top + 20
    for i, (label, _, _, color) in enumerate(series):
        y = legend_y + 22 * i
        lines.append(f'<line x1="{legend_x}" y1="{y}" x2="{legend_x + 25}" y2="{y}" stroke="{color}" stroke-width="3"/>')
        lines.append(f'<text x="{legend_x + 34}" y="{y + 4}" font-family="sans-serif" font-size="13">{html.escape(label)}</text>')
    lines.append("</svg>")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


material = load_material(PARAMETERS)
colors = {"L": "#2166ac", "R": "#b2182b", "T": "#4d9221", "RT": "#762a83", "TL": "#e08214", "LR": "#008080"}

# L/R/T material-axis uniaxial-stress paths.
uniaxial_targets = {
    "L": np.linspace(0.0, -0.010, 161)[1:],
    "R": np.linspace(0.0, -0.025, 161)[1:],
    "T": np.linspace(0.0, -0.030, 161)[1:],
}
all_records: list[dict[str, object]] = []
summaries: list[dict[str, object]] = []
uniaxial_points = {}
for axis, direction in enumerate(("L", "R", "T")):
    points = run_uniaxial_path(material, axis, uniaxial_targets[direction])
    uniaxial_points[direction] = points
    all_records.extend(path_records(f"uniaxial_{direction}", direction, points, axis))
    summaries.append(summarize_path(f"uniaxial_{direction}", direction, points, axis))

# Pure engineering shear paths.
shear_targets = {
    "RT": np.linspace(0.0, 0.080, 161)[1:],
    "TL": np.linspace(0.0, 0.025, 161)[1:],
    "LR": np.linspace(0.0, 0.022, 161)[1:],
}
shear_points = {}
for component, direction in zip((3, 4, 5), ("RT", "TL", "LR")):
    points = run_component_path(material, component, shear_targets[direction])
    shear_points[direction] = points
    all_records.extend(path_records(f"pure_shear_{direction}", direction, points, component))
    summaries.append(summarize_path(f"pure_shear_{direction}", direction, points, component))

# R-axis load, fully elastic unload, then reload beyond the prior maximum.
load_targets = np.linspace(0.0, -0.025, 101)[1:]
unload_targets = np.linspace(-0.025, -0.013, 61)[1:]
reload_targets = np.linspace(-0.013, -0.030, 86)[1:]
cycle_targets = np.concatenate((load_targets, unload_targets, reload_targets))
cycle_points = run_uniaxial_path(material, 1, cycle_targets)
for point in cycle_points:
    segment = "load" if point.step <= len(load_targets) else (
        "unload" if point.step <= len(load_targets) + len(unload_targets) else "reload"
    )
    all_records.extend(path_records("R_load_unload_reload", "R", [point], 1, segment=segment))
summaries.append(summarize_path("R_load_unload_reload", "R", cycle_points, 1))

write_csv(RESULTS / "loading_paths.csv", all_records)
write_csv(RESULTS / "path_summary.csv", summaries)

# Elastic and Hill surface checks.
checks = material.elasticity.checks()
calibration = material.hill.calibration_errors()
material_checks = [
    {"check": key, "value": value, "acceptance": "diagnostic"}
    for key, value in checks.items()
]
material_checks.extend(
    {"check": f"hill_directional_calibration_error_{key}_MPa", "value": value, "acceptance": "abs <= 1e-10"}
    for key, value in calibration.items()
)
material_checks.extend(
    [
        {
            "check": "hill_matrix_min_eigenvalue",
            "value": float(np.linalg.eigvalsh(material.hill.matrix)[0]),
            "acceptance": ">= -1e-10 (positive semidefinite)",
        },
        {
            "check": "hydrostatic_equivalent_stress_MPa",
            "value": material.hill.equivalent_stress(np.ones(6) * np.array([12, 12, 12, 0, 0, 0])),
            "acceptance": "abs <= 1e-10 (pressure insensitive)",
        },
    ]
)
write_csv(RESULTS / "material_checks.csv", material_checks)

# Analytic algorithmic tangent versus independent central finite differences.
tangent_cases = []
elastic_increment = np.array([1e-5, -2e-5, 1e-5, 2e-5, 0.0, 0.0])
plastic_increment = np.array([0.0, -0.015, 0.0, 0.003, 0.002, 0.0])
preloaded = run_uniaxial_path(material, 1, np.linspace(0, -0.018, 81)[1:])[-1].state
continued_increment = np.array([-0.0005, -0.003, 0.0002, 0.001, -0.0004, 0.0003])
for name, old_state, increment in (
    ("elastic_origin", material.initial_state(), elastic_increment),
    ("plastic_origin", material.initial_state(), plastic_increment),
    ("plastic_preloaded", preloaded, continued_increment),
):
    result = material.integrate(old_state, increment)
    finite_difference = finite_difference_tangent(material, old_state, increment, step=1e-8)
    difference = result.tangent - finite_difference
    tangent_cases.append(
        {
            "case": name,
            "branch": "elastic" if result.elastic else "plastic",
            "fd_step": 1e-8,
            "relative_frobenius_error": float(np.linalg.norm(difference) / np.linalg.norm(finite_difference)),
            "max_absolute_error_MPa": float(np.max(np.abs(difference))),
            "analytic_symmetry_error_MPa": float(np.max(np.abs(result.tangent - result.tangent.T))),
            "analytic_min_eigenvalue_MPa": float(np.linalg.eigvalsh(result.tangent)[0]),
        }
    )

# Rotation covariance and global finite-difference tangent at 37 degrees.
rotation = rotation_about_global_t(np.deg2rad(37.0))
material_increment = np.array([-0.008, -0.012, 0.002, 0.040, 0.005, -0.003])
global_increment = strain_material_to_global(material_increment, rotation)
direct = material.integrate(material.initial_state(), material_increment)
rotated, global_stress, global_tangent = material.integrate_global(
    material.initial_state(), global_increment, rotation
)
expected_global_stress = stress_material_to_global(direct.state.stress, rotation)
global_fd = np.empty((6, 6))
for i in range(6):
    perturb = np.zeros(6)
    perturb[i] = 1e-8
    plus = material.integrate_global(material.initial_state(), global_increment + perturb, rotation)[1]
    minus = material.integrate_global(material.initial_state(), global_increment - perturb, rotation)[1]
    global_fd[:, i] = (plus - minus) / (2e-8)
rotation_metrics = {
    "angle_degrees": 37.0,
    "material_state_stress_covariance_error_MPa": float(np.max(np.abs(rotated.state.stress - direct.state.stress))),
    "global_stress_transform_error_MPa": float(np.max(np.abs(global_stress - expected_global_stress))),
    "work_conjugacy_error_MPa": float(
        abs(global_stress @ global_increment - direct.state.stress @ material_increment)
    ),
    "global_tangent_relative_fd_error": float(np.linalg.norm(global_tangent - global_fd) / np.linalg.norm(global_fd)),
}
tangent_cases.append(
    {
        "case": "rotated_37_degrees_global",
        "branch": "elastic" if direct.elastic else "plastic",
        "fd_step": 1e-8,
        "relative_frobenius_error": rotation_metrics["global_tangent_relative_fd_error"],
        "max_absolute_error_MPa": float(np.max(np.abs(global_tangent - global_fd))),
        "analytic_symmetry_error_MPa": float(np.max(np.abs(global_tangent - global_tangent.T))),
        "analytic_min_eigenvalue_MPa": float(np.linalg.eigvalsh(global_tangent)[0]),
    }
)
write_csv(RESULTS / "tangent_checks.csv", tangent_cases)

# Non-proportional step-size convergence.
def mixed_path_endpoint(number_of_steps: int):
    state = material.initial_state()
    previous = np.zeros(6)
    max_residual = 0.0
    for t in np.linspace(0.0, 1.0, number_of_steps + 1)[1:]:
        current = np.array(
            [-0.008 * t, -0.018 * t, 0.003 * t, 0.050 * t * t,
             0.006 * np.sin(np.pi * t), -0.004 * t * t]
        )
        result = material.integrate(state, current - previous)
        state = result.state
        previous = current
        if not result.elastic:
            max_residual = max(max_residual, abs(result.yield_residual))
    return state, max_residual


reference_state, _ = mixed_path_endpoint(2560)
convergence_rows = []
step_counts = (10, 20, 40, 80, 160, 320, 640)
for count in step_counts:
    state, max_residual = mixed_path_endpoint(count)
    convergence_rows.append(
        {
            "steps": count,
            "step_size": 1.0 / count,
            "stress_error_MPa": float(np.linalg.norm(state.stress - reference_state.stress)),
            "kappa_error": abs(state.kappa - reference_state.kappa),
            "max_abs_yield_residual_MPa": max_residual,
            "final_plastic_dissipation_MPa": state.plastic_dissipation,
        }
    )
fit_rows = convergence_rows[1:-1]
convergence_order = float(
    np.polyfit(
        np.log([row["step_size"] for row in fit_rows]),
        np.log([row["stress_error_MPa"] for row in fit_rows]),
        1,
    )[0]
)
write_csv(RESULTS / "step_convergence.csv", convergence_rows)

# Dedicated load/unload/reload checks.
load_end = cycle_points[len(load_targets) - 1]
unload_end = cycle_points[len(load_targets) + len(unload_targets) - 1]
cycle_checks = [
    {"check": "kappa_at_load_end", "value": load_end.state.kappa, "acceptance": "diagnostic"},
    {"check": "kappa_at_unload_end", "value": unload_end.state.kappa, "acceptance": "unchanged from load end"},
    {
        "check": "unload_kappa_change",
        "value": unload_end.state.kappa - load_end.state.kappa,
        "acceptance": "abs <= 1e-12",
    },
    {
        "check": "elastic_unload_steps",
        "value": sum(point.result.elastic for point in cycle_points[len(load_targets):len(load_targets) + len(unload_targets)]),
        "acceptance": f"{len(unload_targets)} expected",
    },
    {
        "check": "minimum_kappa_increment",
        "value": float(min(np.diff([0.0] + [point.state.kappa for point in cycle_points]))),
        "acceptance": ">= -1e-14",
    },
]
write_csv(RESULTS / "cycle_checks.csv", cycle_checks)

# Figures are original line drawings made from the solver results.
strengths = {"L": 40.0, "R": 6.0, "T": 6.0, "RT": 3.0, "TL": 4.0, "LR": 4.5}
svg_plot(
    FIGURES / "uniaxial_material_axes.svg",
    [
        (
            direction,
            np.array([point.target_strain for point in uniaxial_points[direction]]),
            np.array([point.state.stress[i] / strengths[direction] for point in uniaxial_points[direction]]),
            colors[direction],
        )
        for i, direction in enumerate(("L", "R", "T"))
    ],
    "M0 material-axis uniaxial-stress paths",
    "Axial strain",
    "Stress / initial directional yield stress",
)
svg_plot(
    FIGURES / "pure_shear_paths.svg",
    [
        (
            direction,
            np.array([point.target_strain for point in shear_points[direction]]),
            np.array([point.state.stress[i] / strengths[direction] for point in shear_points[direction]]),
            colors[direction],
        )
        for i, direction in zip((3, 4, 5), ("RT", "TL", "LR"))
    ],
    "M0 pure engineering-shear paths",
    "Engineering shear strain",
    "Shear stress / initial shear yield stress",
)
svg_plot(
    FIGURES / "load_unload_reload.svg",
    [
        (
            "R load-unload-reload",
            np.array([point.target_strain for point in cycle_points]),
            np.array([point.state.stress[1] for point in cycle_points]),
            colors["R"],
        )
    ],
    "R-axis load, elastic unload, and reload",
    "R axial strain",
    "R stress (MPa)",
)
svg_plot(
    FIGURES / "step_convergence.svg",
    [
        (
            "Stress endpoint error",
            np.array([row["step_size"] for row in convergence_rows]),
            np.array([row["stress_error_MPa"] for row in convergence_rows]),
            "#2166ac",
        )
    ],
    f"Non-proportional path convergence (observed order {convergence_order:.2f})",
    "Normalized step size",
    "Endpoint stress error (MPa)",
    log_x=True,
    log_y=True,
)

max_plastic_residual = max(row["max_abs_plastic_yield_residual_MPa"] for row in summaries)
max_newton_residual = max(row["max_scaled_newton_residual"] for row in summaries)
min_dissipation = min(row["minimum_plastic_dissipation_increment_MPa"] for row in summaries)
max_lateral = max(row["max_lateral_stress_norm_MPa"] for row in summaries)
summary = {
    "project_status": "validated-prototype",
    "implemented_scope": "M0 material point only",
    "not_implemented": [contract.stage for contract in extension_roadmap()] + [
        "finite_element_model",
        "Marc_or_HYPELA2_subroutine",
        "contact",
        "damage",
        "densification",
        "engineering_design_validation",
    ],
    "parameter_status": "custom_algorithm_demonstration_not_calibrated",
    "voigt_convention": {
        "order": list(VOIGT_LABELS),
        "stress_shear": "physical",
        "strain_shear": "engineering",
        "work_pair": "sigma_vector dot epsilon_vector",
    },
    "solver": {
        "return_mapping": "fully implicit backward Euler with safeguarded Newton",
        "tangent": "Jacobian-consistent algorithmic tangent, independently central-difference verified",
        "max_abs_plastic_yield_residual_MPa": max_plastic_residual,
        "max_scaled_newton_residual": max_newton_residual,
        "max_return_iterations": max(row["max_return_iterations"] for row in summaries),
        "max_fd_tangent_relative_error": max(row["relative_frobenius_error"] for row in tangent_cases),
    },
    "elasticity": checks,
    "rotation": rotation_metrics,
    "paths": {
        "maximum_lateral_stress_norm_MPa": max_lateral,
        "minimum_plastic_dissipation_increment_MPa": min_dissipation,
        "cycle_unload_kappa_change": unload_end.state.kappa - load_end.state.kappa,
        "nonproportional_step_convergence_order": convergence_order,
        "step_10_to_640_stress_error_reduction": convergence_rows[0]["stress_error_MPa"] / convergence_rows[-1]["stress_error_MPa"],
    },
    "claim_boundary": [
        "Algorithm verification is not material calibration.",
        "Hill M0 is pressure-insensitive and symmetric in tension/compression.",
        "No Xu embedment FE, Marc/HYPELA2, contact, densification, damage, or safety claim is reproduced.",
    ],
}
(RESULTS / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

roadmap_rows = []
for contract in extension_roadmap():
    roadmap_rows.append(
        {
            "stage": contract.stage,
            "status": contract.status,
            "proposed_state_variables": "; ".join(contract.proposed_state_variables) or "none",
            "required_verification_gates": "; ".join(contract.required_verification_gates),
        }
    )
write_csv(RESULTS / "extension_roadmap.csv", roadmap_rows)

artifact_paths = sorted(
    [path for path in RESULTS.glob("*") if path.name != "artifact_manifest.csv"]
    + list(FIGURES.glob("*.svg"))
)
manifest_rows = [
    {
        "path": str(path.relative_to(ROOT)),
        "bytes": path.stat().st_size,
        "sha256": sha256(path),
    }
    for path in artifact_paths
]
write_csv(RESULTS / "artifact_manifest.csv", manifest_rows)

print(
    json.dumps(
        {
            "status": summary["project_status"],
            "max_yield_residual_MPa": max_plastic_residual,
            "max_tangent_relative_error": summary["solver"]["max_fd_tangent_relative_error"],
            "convergence_order": convergence_order,
        },
        indent=2,
    )
)

