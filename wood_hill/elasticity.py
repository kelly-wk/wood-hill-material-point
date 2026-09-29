"""Small-strain orthotropic elasticity in material L/R/T axes."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class OrthotropicElasticity:
    """Nine independent engineering constants.

    Poisson ratio ``nu_ij`` means minus strain in ``j`` caused by uniaxial
    stress in ``i``. Reciprocal ratios are computed from
    ``nu_ij / E_i = nu_ji / E_j`` rather than accepted independently.
    Moduli and stresses use the same units (MPa in the demonstration).
    """

    E_L: float
    E_R: float
    E_T: float
    nu_LR: float
    nu_LT: float
    nu_RT: float
    G_RT: float
    G_TL: float
    G_LR: float

    def __post_init__(self) -> None:
        values = np.array(
            [self.E_L, self.E_R, self.E_T, self.G_RT, self.G_TL, self.G_LR],
            dtype=float,
        )
        if not np.all(np.isfinite(values)) or np.any(values <= 0):
            raise ValueError("all Young and shear moduli must be positive finite values")
        poisson = np.array([self.nu_LR, self.nu_LT, self.nu_RT], dtype=float)
        if not np.all(np.isfinite(poisson)):
            raise ValueError("Poisson ratios must be finite")
        # Building the matrices performs the decisive thermodynamic check.
        _ = self.stiffness

    @property
    def reciprocal_poisson(self) -> dict[str, float]:
        return {
            "nu_RL": self.nu_LR * self.E_R / self.E_L,
            "nu_TL": self.nu_LT * self.E_T / self.E_L,
            "nu_TR": self.nu_RT * self.E_T / self.E_R,
        }

    @property
    def compliance(self) -> np.ndarray:
        s = np.zeros((6, 6), dtype=float)
        s[0, 0] = 1.0 / self.E_L
        s[1, 1] = 1.0 / self.E_R
        s[2, 2] = 1.0 / self.E_T
        s[0, 1] = s[1, 0] = -self.nu_LR / self.E_L
        s[0, 2] = s[2, 0] = -self.nu_LT / self.E_L
        s[1, 2] = s[2, 1] = -self.nu_RT / self.E_R
        s[3, 3] = 1.0 / self.G_RT
        s[4, 4] = 1.0 / self.G_TL
        s[5, 5] = 1.0 / self.G_LR
        return s

    @property
    def stiffness(self) -> np.ndarray:
        s = self.compliance
        if not np.allclose(s, s.T, rtol=0.0, atol=1e-13):
            raise ValueError("orthotropic compliance violates reciprocity")
        eigenvalues = np.linalg.eigvalsh(s)
        if eigenvalues[0] <= 0:
            raise ValueError(
                "orthotropic compliance is not positive definite; check Poisson ratios"
            )
        c = np.linalg.inv(s)
        if not np.allclose(c, c.T, rtol=1e-12, atol=1e-10):
            raise ValueError("computed stiffness is not symmetric")
        return c

    def checks(self) -> dict[str, float]:
        s = self.compliance
        c = self.stiffness
        rec = self.reciprocal_poisson
        reciprocity_errors = [
            abs(self.nu_LR / self.E_L - rec["nu_RL"] / self.E_R),
            abs(self.nu_LT / self.E_L - rec["nu_TL"] / self.E_T),
            abs(self.nu_RT / self.E_R - rec["nu_TR"] / self.E_T),
        ]
        return {
            "compliance_min_eigenvalue": float(np.linalg.eigvalsh(s)[0]),
            "stiffness_min_eigenvalue": float(np.linalg.eigvalsh(c)[0]),
            "stiffness_condition_number": float(np.linalg.cond(c)),
            "stiffness_symmetry_error": float(np.max(np.abs(c - c.T))),
            "max_reciprocity_error": float(max(reciprocity_errors)),
        }

