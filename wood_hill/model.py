"""Implicit material-point integration for the validated M0 prototype."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .conventions import (
    strain_global_to_material,
    strain_material_to_global,
    stress_material_to_global,
    transformation_matrices,
    validate_rotation,
)
from .elasticity import OrthotropicElasticity
from .hill import HillYield


def _zero6() -> np.ndarray:
    return np.zeros(6, dtype=float)


@dataclass(frozen=True)
class MaterialState:
    """State stored in material axes; strain shears are engineering shears."""

    stress: np.ndarray = field(default_factory=_zero6)
    plastic_strain: np.ndarray = field(default_factory=_zero6)
    total_strain: np.ndarray = field(default_factory=_zero6)
    kappa: float = 0.0
    plastic_dissipation: float = 0.0

    def __post_init__(self) -> None:
        for name in ("stress", "plastic_strain", "total_strain"):
            value = np.asarray(getattr(self, name), dtype=float)
            if value.shape != (6,) or not np.all(np.isfinite(value)):
                raise ValueError(f"{name} must be a finite vector with shape (6,)")
            value = value.copy()
            value.setflags(write=False)
            object.__setattr__(self, name, value)
        if not np.isfinite(self.kappa) or self.kappa < 0:
            raise ValueError("kappa must be non-negative")
        if not np.isfinite(self.plastic_dissipation) or self.plastic_dissipation < -1e-12:
            raise ValueError("plastic dissipation must be non-negative")


@dataclass(frozen=True)
class UpdateResult:
    state: MaterialState
    tangent: np.ndarray
    elastic: bool
    iterations: int
    plastic_multiplier: float
    yield_residual: float
    newton_residual: float
    plastic_dissipation_increment: float


class HillMaterial:
    """M0: orthotropic elasticity + associated Hill plasticity + linear hardening."""

    def __init__(
        self,
        elasticity: OrthotropicElasticity,
        hill: HillYield,
        hardening_modulus: float,
        relative_tolerance: float = 1e-11,
        max_iterations: int = 40,
    ) -> None:
        if not np.isfinite(hardening_modulus) or hardening_modulus < 0:
            raise ValueError("hardening modulus must be non-negative")
        if relative_tolerance <= 0 or max_iterations < 2:
            raise ValueError("invalid nonlinear solver controls")
        self.elasticity = elasticity
        self.hill = hill
        self.hardening_modulus = float(hardening_modulus)
        self.relative_tolerance = float(relative_tolerance)
        self.max_iterations = int(max_iterations)

    @property
    def initial_radius(self) -> float:
        return float(self.hill.reference_stress)

    def radius(self, kappa: float) -> float:
        return self.initial_radius + self.hardening_modulus * kappa

    def yield_function(self, stress: np.ndarray, kappa: float) -> float:
        return self.hill.equivalent_stress(stress) - self.radius(kappa)

    def initial_state(self) -> MaterialState:
        return MaterialState()

    def _residual_and_jacobian(
        self, stress: np.ndarray, delta_lambda: float, trial: np.ndarray, kappa_old: float
    ) -> tuple[np.ndarray, np.ndarray]:
        c = self.elasticity.stiffness
        n = self.hill.gradient(stress)
        dn = self.hill.hessian(stress)
        residual = np.empty(7, dtype=float)
        residual[:6] = stress - trial + delta_lambda * (c @ n)
        residual[6] = self.hill.equivalent_stress(stress) - self.radius(
            kappa_old + delta_lambda
        )
        jacobian = np.empty((7, 7), dtype=float)
        jacobian[:6, :6] = np.eye(6) + delta_lambda * (c @ dn)
        jacobian[:6, 6] = c @ n
        jacobian[6, :6] = n
        jacobian[6, 6] = -self.hardening_modulus
        return residual, jacobian

    def _scaled_residual_norm(self, residual: np.ndarray, trial: np.ndarray) -> float:
        scale = max(self.initial_radius, float(np.linalg.norm(trial)), 1.0)
        return float(np.linalg.norm(residual / scale, ord=np.inf))

    def integrate(self, state: MaterialState, strain_increment: np.ndarray) -> UpdateResult:
        dstrain = np.asarray(strain_increment, dtype=float)
        if dstrain.shape != (6,) or not np.all(np.isfinite(dstrain)):
            raise ValueError("strain increment must be a finite vector with shape (6,)")
        c = self.elasticity.stiffness
        trial = state.stress + c @ dstrain
        trial_f = self.yield_function(trial, state.kappa)
        yield_scale = max(self.initial_radius, self.radius(state.kappa), 1.0)
        tolerance = self.relative_tolerance * yield_scale

        if trial_f <= tolerance:
            new_state = MaterialState(
                stress=trial,
                plastic_strain=state.plastic_strain,
                total_strain=state.total_strain + dstrain,
                kappa=state.kappa,
                plastic_dissipation=state.plastic_dissipation,
            )
            return UpdateResult(
                state=new_state,
                tangent=c.copy(),
                elastic=True,
                iterations=0,
                plastic_multiplier=0.0,
                yield_residual=float(trial_f),
                newton_residual=0.0,
                plastic_dissipation_increment=0.0,
            )

        n_trial = self.hill.gradient(trial)
        denominator = float(n_trial @ c @ n_trial + self.hardening_modulus)
        if denominator <= 0:
            raise FloatingPointError("non-positive return-mapping denominator")
        delta_lambda = max(0.0, trial_f / denominator)
        stress = trial - delta_lambda * (c @ n_trial)

        converged = False
        final_norm = float("inf")
        iteration = 0
        for iteration in range(1, self.max_iterations + 1):
            residual, jacobian = self._residual_and_jacobian(
                stress, delta_lambda, trial, state.kappa
            )
            final_norm = self._scaled_residual_norm(residual, trial)
            if final_norm <= self.relative_tolerance:
                converged = True
                break
            step = np.linalg.solve(jacobian, -residual)
            merit = final_norm
            accepted = False
            step_length = 1.0
            for _ in range(24):
                candidate_lambda = delta_lambda + step_length * step[6]
                candidate_stress = stress + step_length * step[:6]
                if candidate_lambda < 0 or self.hill.equivalent_stress(candidate_stress) <= 1e-14:
                    step_length *= 0.5
                    continue
                candidate_residual, _ = self._residual_and_jacobian(
                    candidate_stress, candidate_lambda, trial, state.kappa
                )
                candidate_merit = self._scaled_residual_norm(candidate_residual, trial)
                if candidate_merit < merit * (1.0 - 1e-4 * step_length) or candidate_merit < self.relative_tolerance:
                    stress = candidate_stress
                    delta_lambda = candidate_lambda
                    accepted = True
                    break
                step_length *= 0.5
            if not accepted:
                raise RuntimeError("return mapping line search failed")

        if not converged:
            residual, _ = self._residual_and_jacobian(stress, delta_lambda, trial, state.kappa)
            final_norm = self._scaled_residual_norm(residual, trial)
            if final_norm > self.relative_tolerance:
                raise RuntimeError(
                    f"return mapping did not converge after {self.max_iterations} iterations"
                )

        residual, jacobian = self._residual_and_jacobian(
            stress, delta_lambda, trial, state.kappa
        )
        rhs = np.vstack((c, np.zeros((1, 6))))
        sensitivity = np.linalg.solve(jacobian, rhs)
        consistent_tangent = sensitivity[:6, :]
        # Roundoff can leave an O(1e-12) skew part for an associative model.
        consistent_tangent = 0.5 * (consistent_tangent + consistent_tangent.T)

        n = self.hill.gradient(stress)
        plastic_increment = delta_lambda * n
        dissipation_increment = float(stress @ plastic_increment)
        dissipation_tolerance = 1e-11 * max(1.0, abs(dissipation_increment))
        if dissipation_increment < -dissipation_tolerance:
            raise FloatingPointError("negative incremental plastic dissipation")
        dissipation_increment = max(0.0, dissipation_increment)
        new_state = MaterialState(
            stress=stress,
            plastic_strain=state.plastic_strain + plastic_increment,
            total_strain=state.total_strain + dstrain,
            kappa=state.kappa + delta_lambda,
            plastic_dissipation=state.plastic_dissipation + dissipation_increment,
        )
        yield_residual = self.yield_function(stress, new_state.kappa)
        return UpdateResult(
            state=new_state,
            tangent=consistent_tangent,
            elastic=False,
            iterations=iteration,
            plastic_multiplier=float(delta_lambda),
            yield_residual=float(yield_residual),
            newton_residual=float(self._scaled_residual_norm(residual, trial)),
            plastic_dissipation_increment=dissipation_increment,
        )

    def integrate_global(
        self,
        state_material: MaterialState,
        strain_increment_global: np.ndarray,
        rotation: np.ndarray,
    ) -> tuple[UpdateResult, np.ndarray, np.ndarray]:
        """Update in material axes and return global stress and tangent."""
        q = validate_rotation(rotation)
        dstrain_material = strain_global_to_material(strain_increment_global, q)
        result = self.integrate(state_material, dstrain_material)
        stress_global = stress_material_to_global(result.state.stress, q)
        t_eps, t_sig = transformation_matrices(q)
        tangent_global = t_sig @ result.tangent @ t_eps
        return result, stress_global, tangent_global

    def elastic_energy(self, state: MaterialState) -> float:
        elastic_strain = state.total_strain - state.plastic_strain
        value = 0.5 * float(elastic_strain @ self.elasticity.stiffness @ elastic_strain)
        return max(0.0, value)


def finite_difference_tangent(
    material: HillMaterial,
    state: MaterialState,
    strain_increment: np.ndarray,
    step: float = 1e-8,
) -> np.ndarray:
    """Central-difference derivative of the complete update mapping.

    This intentionally re-integrates from the same old state for every column.
    It is an independent verification oracle for the returned algorithmic tangent.
    """
    if step <= 0:
        raise ValueError("finite-difference step must be positive")
    base = np.asarray(strain_increment, dtype=float)
    tangent = np.empty((6, 6), dtype=float)
    for i in range(6):
        perturbation = np.zeros(6)
        perturbation[i] = step
        plus = material.integrate(state, base + perturbation).state.stress
        minus = material.integrate(state, base - perturbation).state.stress
        tangent[:, i] = (plus - minus) / (2.0 * step)
    return tangent

