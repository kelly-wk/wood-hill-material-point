"""Material-point loading drivers used by tests and verification scripts."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .model import HillMaterial, MaterialState, UpdateResult


@dataclass(frozen=True)
class PathPoint:
    step: int
    target_strain: float
    state: MaterialState
    result: UpdateResult
    lateral_stress_norm: float


def uniaxial_stress_increment(
    material: HillMaterial,
    state: MaterialState,
    axis: int,
    axial_increment: float,
    lateral_tolerance: float = 1e-10,
    max_iterations: int = 20,
) -> UpdateResult:
    """Apply an axial strain increment while solving for zero lateral stress.

    The driver is valid for material-axis L/R/T tests (``axis`` 0, 1, or 2).
    Each lateral Newton iteration re-integrates from the unchanged old state.
    """
    if axis not in (0, 1, 2):
        raise ValueError("uniaxial axis must be 0 (L), 1 (R), or 2 (T)")
    lateral = [i for i in (0, 1, 2) if i != axis]
    dstrain = np.zeros(6)
    dstrain[axis] = axial_increment
    # Elastic Poisson response provides a close first estimate.
    target_stress_increment = np.zeros(6)
    target_stress_increment[axis] = material.elasticity.stiffness[axis, axis] * axial_increment
    elastic_guess = material.elasticity.compliance @ target_stress_increment
    dstrain[lateral] = elastic_guess[lateral]

    result = material.integrate(state, dstrain)
    stress_scale = max(material.initial_radius, abs(result.state.stress[axis]), 1.0)
    for _ in range(max_iterations):
        residual = result.state.stress[lateral]
        if np.linalg.norm(residual, ord=np.inf) <= lateral_tolerance * stress_scale:
            return result
        jacobian = result.tangent[np.ix_(lateral, lateral)]
        correction = np.linalg.solve(jacobian, -residual)
        accepted = False
        old_norm = float(np.linalg.norm(residual, ord=np.inf))
        alpha = 1.0
        for _ in range(16):
            candidate = dstrain.copy()
            candidate[lateral] += alpha * correction
            candidate_result = material.integrate(state, candidate)
            candidate_norm = float(
                np.linalg.norm(candidate_result.state.stress[lateral], ord=np.inf)
            )
            if candidate_norm < old_norm:
                dstrain = candidate
                result = candidate_result
                accepted = True
                break
            alpha *= 0.5
        if not accepted:
            raise RuntimeError("uniaxial lateral-stress Newton line search failed")
    raise RuntimeError("uniaxial lateral-stress solve did not converge")


def run_uniaxial_path(
    material: HillMaterial,
    axis: int,
    strain_targets: np.ndarray,
    state: MaterialState | None = None,
) -> list[PathPoint]:
    targets = np.asarray(strain_targets, dtype=float)
    if targets.ndim != 1 or not np.all(np.isfinite(targets)):
        raise ValueError("strain targets must be a finite one-dimensional array")
    current = material.initial_state() if state is None else state
    previous_target = 0.0
    points: list[PathPoint] = []
    for step, target in enumerate(targets, start=1):
        result = uniaxial_stress_increment(
            material, current, axis, float(target - previous_target)
        )
        current = result.state
        lateral = [i for i in (0, 1, 2) if i != axis]
        points.append(
            PathPoint(
                step=step,
                target_strain=float(target),
                state=current,
                result=result,
                lateral_stress_norm=float(np.linalg.norm(current.stress[lateral], ord=np.inf)),
            )
        )
        previous_target = float(target)
    return points


def run_component_path(
    material: HillMaterial,
    component: int,
    strain_targets: np.ndarray,
    state: MaterialState | None = None,
) -> list[PathPoint]:
    if component not in range(6):
        raise ValueError("component must be a Voigt index 0..5")
    targets = np.asarray(strain_targets, dtype=float)
    current = material.initial_state() if state is None else state
    previous_target = 0.0
    points: list[PathPoint] = []
    for step, target in enumerate(targets, start=1):
        dstrain = np.zeros(6)
        dstrain[component] = target - previous_target
        result = material.integrate(current, dstrain)
        current = result.state
        points.append(
            PathPoint(
                step=step,
                target_strain=float(target),
                state=current,
                result=result,
                lateral_stress_norm=0.0,
            )
        )
        previous_target = float(target)
    return points

