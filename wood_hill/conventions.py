"""Work-conjugate Voigt conventions and material-axis rotations.

Stress order is [LL, RR, TT, RT, TL, LR] with physical shear stress.
Strain order is [LL, RR, TT, RT, TL, LR] with engineering shear strain.
Consequently ``stress @ strain`` is the tensor double contraction.
"""

from __future__ import annotations

import numpy as np


VOIGT_LABELS = ("LL", "RR", "TT", "RT", "TL", "LR")


def stress_to_tensor(stress: np.ndarray) -> np.ndarray:
    s = np.asarray(stress, dtype=float)
    if s.shape != (6,):
        raise ValueError("stress must have shape (6,)")
    return np.array(
        [[s[0], s[5], s[4]], [s[5], s[1], s[3]], [s[4], s[3], s[2]]],
        dtype=float,
    )


def stress_from_tensor(tensor: np.ndarray) -> np.ndarray:
    s = np.asarray(tensor, dtype=float)
    if s.shape != (3, 3):
        raise ValueError("stress tensor must have shape (3, 3)")
    if not np.allclose(s, s.T, rtol=0.0, atol=1e-12):
        raise ValueError("stress tensor must be symmetric")
    return np.array([s[0, 0], s[1, 1], s[2, 2], s[1, 2], s[2, 0], s[0, 1]])


def strain_to_tensor(strain: np.ndarray) -> np.ndarray:
    e = np.asarray(strain, dtype=float)
    if e.shape != (6,):
        raise ValueError("strain must have shape (6,)")
    return np.array(
        [
            [e[0], 0.5 * e[5], 0.5 * e[4]],
            [0.5 * e[5], e[1], 0.5 * e[3]],
            [0.5 * e[4], 0.5 * e[3], e[2]],
        ],
        dtype=float,
    )


def strain_from_tensor(tensor: np.ndarray) -> np.ndarray:
    e = np.asarray(tensor, dtype=float)
    if e.shape != (3, 3):
        raise ValueError("strain tensor must have shape (3, 3)")
    if not np.allclose(e, e.T, rtol=0.0, atol=1e-12):
        raise ValueError("strain tensor must be symmetric")
    return np.array([e[0, 0], e[1, 1], e[2, 2], 2 * e[1, 2], 2 * e[2, 0], 2 * e[0, 1]])


def validate_rotation(rotation: np.ndarray, tolerance: float = 1e-12) -> np.ndarray:
    q = np.asarray(rotation, dtype=float)
    if q.shape != (3, 3):
        raise ValueError("rotation must have shape (3, 3)")
    if not np.allclose(q.T @ q, np.eye(3), rtol=0.0, atol=tolerance):
        raise ValueError("rotation columns must form an orthonormal material basis")
    if np.linalg.det(q) < 1 - tolerance:
        raise ValueError("rotation must be proper (determinant +1)")
    return q


def strain_global_to_material(strain_global: np.ndarray, rotation: np.ndarray) -> np.ndarray:
    """Rotate global engineering strain into material L/R/T axes.

    ``rotation`` has the material basis vectors as columns in global coordinates.
    """
    q = validate_rotation(rotation)
    e_global = strain_to_tensor(strain_global)
    return strain_from_tensor(q.T @ e_global @ q)


def strain_material_to_global(strain_material: np.ndarray, rotation: np.ndarray) -> np.ndarray:
    q = validate_rotation(rotation)
    e_material = strain_to_tensor(strain_material)
    return strain_from_tensor(q @ e_material @ q.T)


def stress_global_to_material(stress_global: np.ndarray, rotation: np.ndarray) -> np.ndarray:
    q = validate_rotation(rotation)
    s_global = stress_to_tensor(stress_global)
    return stress_from_tensor(q.T @ s_global @ q)


def stress_material_to_global(stress_material: np.ndarray, rotation: np.ndarray) -> np.ndarray:
    q = validate_rotation(rotation)
    s_material = stress_to_tensor(stress_material)
    return stress_from_tensor(q @ s_material @ q.T)


def transformation_matrices(rotation: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return ``T_eps`` and ``T_sig`` for material/global transformations.

    ``eps_m = T_eps @ eps_g`` and ``sig_g = T_sig @ sig_m``. Therefore an
    algorithmic material tangent transforms as ``D_g = T_sig @ D_m @ T_eps``.
    """
    q = validate_rotation(rotation)
    t_eps = np.column_stack(
        [strain_global_to_material(np.eye(6)[i], q) for i in range(6)]
    )
    t_sig = np.column_stack(
        [stress_material_to_global(np.eye(6)[i], q) for i in range(6)]
    )
    return t_eps, t_sig


def rotation_about_global_t(angle_radians: float) -> np.ndarray:
    """Rotate the material L/R plane about the global third axis."""
    c = float(np.cos(angle_radians))
    s = float(np.sin(angle_radians))
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])

