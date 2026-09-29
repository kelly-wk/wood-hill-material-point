"""Hill (1948) quadratic equivalent stress in material axes."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class HillYield:
    """Quadratic, pressure-insensitive, tension/compression-symmetric surface.

    Directional values are initial yield stresses in the material axes:
    ``X_L, Y_R, Z_T, S_RT, S_TL, S_LR``. The matrix is normalized to
    ``reference_stress`` so ``q(sigma_directional_yield)`` equals that
    reference stress for every supplied direction.
    """

    X_L: float
    Y_R: float
    Z_T: float
    S_RT: float
    S_TL: float
    S_LR: float
    reference_stress: float | None = None

    def __post_init__(self) -> None:
        strengths = np.array(
            [self.X_L, self.Y_R, self.Z_T, self.S_RT, self.S_TL, self.S_LR],
            dtype=float,
        )
        if not np.all(np.isfinite(strengths)) or np.any(strengths <= 0):
            raise ValueError("all Hill directional yield stresses must be positive")
        ref = self.X_L if self.reference_stress is None else self.reference_stress
        if not np.isfinite(ref) or ref <= 0:
            raise ValueError("reference stress must be positive")
        object.__setattr__(self, "reference_stress", float(ref))
        eigenvalues = np.linalg.eigvalsh(self.matrix)
        tolerance = 1e-11 * max(1.0, np.max(np.abs(eigenvalues)))
        if eigenvalues[0] < -tolerance:
            raise ValueError(
                "directional yield stresses produce a non-convex Hill quadratic"
            )

    @property
    def coefficients(self) -> dict[str, float]:
        x, y, z = self.X_L, self.Y_R, self.Z_T
        scale = self.reference_stress**2
        return {
            "F": 0.5 * scale * (1.0 / y**2 + 1.0 / z**2 - 1.0 / x**2),
            "G": 0.5 * scale * (1.0 / z**2 + 1.0 / x**2 - 1.0 / y**2),
            "H": 0.5 * scale * (1.0 / x**2 + 1.0 / y**2 - 1.0 / z**2),
            "L": 0.5 * scale / self.S_RT**2,
            "M": 0.5 * scale / self.S_TL**2,
            "N": 0.5 * scale / self.S_LR**2,
        }

    @property
    def matrix(self) -> np.ndarray:
        c = self.coefficients
        f, g, h, l, m, n = (c[k] for k in ("F", "G", "H", "L", "M", "N"))
        a = np.zeros((6, 6), dtype=float)
        a[:3, :3] = np.array(
            [[g + h, -h, -g], [-h, f + h, -f], [-g, -f, f + g]]
        )
        a[3, 3] = 2.0 * l
        a[4, 4] = 2.0 * m
        a[5, 5] = 2.0 * n
        return a

    def equivalent_stress(self, stress: np.ndarray) -> float:
        s = np.asarray(stress, dtype=float)
        value = float(s @ self.matrix @ s)
        tolerance = 1e-12 * max(1.0, float(np.linalg.norm(s)) ** 2)
        if value < -tolerance:
            raise FloatingPointError("negative Hill quadratic value")
        return float(np.sqrt(max(0.0, value)))

    def gradient(self, stress: np.ndarray) -> np.ndarray:
        s = np.asarray(stress, dtype=float)
        q = self.equivalent_stress(s)
        if q <= np.finfo(float).eps:
            raise FloatingPointError("Hill gradient is undefined at zero equivalent stress")
        return self.matrix @ s / q

    def hessian(self, stress: np.ndarray) -> np.ndarray:
        s = np.asarray(stress, dtype=float)
        q = self.equivalent_stress(s)
        if q <= np.finfo(float).eps:
            raise FloatingPointError("Hill Hessian is undefined at zero equivalent stress")
        n = self.matrix @ s / q
        return self.matrix / q - np.outer(n, n) / q

    def calibration_errors(self) -> dict[str, float]:
        strengths = [self.X_L, self.Y_R, self.Z_T, self.S_RT, self.S_TL, self.S_LR]
        errors: dict[str, float] = {}
        for i, (label, value) in enumerate(zip(("L", "R", "T", "RT", "TL", "LR"), strengths)):
            stress = np.zeros(6)
            stress[i] = value
            errors[label] = self.equivalent_stress(stress) - self.reference_stress
        return errors

