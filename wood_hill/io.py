"""Strict JSON parameter loading for the demonstration material."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from .elasticity import OrthotropicElasticity
from .hill import HillYield
from .model import HillMaterial


def load_parameter_document(path: str | Path) -> dict[str, Any]:
    with Path(path).open("r", encoding="utf-8") as stream:
        document = json.load(stream)
    if document.get("parameter_status") != "custom_algorithm_demonstration_not_calibrated":
        raise ValueError("parameter document must declare its non-calibrated demonstration status")
    return document


def material_from_document(document: dict[str, Any]) -> HillMaterial:
    elasticity = OrthotropicElasticity(**document["orthotropic_elasticity_MPa"])
    hill = HillYield(**document["hill_initial_yield_MPa"])
    controls = document.get("solver_controls", {})
    return HillMaterial(
        elasticity=elasticity,
        hill=hill,
        hardening_modulus=float(document["linear_isotropic_hardening_MPa"]),
        relative_tolerance=float(controls.get("relative_tolerance", 1e-11)),
        max_iterations=int(controls.get("max_iterations", 40)),
    )


def load_material(path: str | Path) -> HillMaterial:
    return material_from_document(load_parameter_document(path))

