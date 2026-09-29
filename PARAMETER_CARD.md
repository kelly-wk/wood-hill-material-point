# Parameter card: illustrative M0 material

Status: `custom_algorithm_demonstration_not_calibrated`

These numbers exist only to exercise the solver. They do not represent a species, grade, moisture state, temperature, strain rate, product, standard, Xu et al. specimen, or safe design input.

## Orthotropic elasticity

| Parameter | Value | Unit |
|---|---:|---|
| `E_L` | 12000 | MPa |
| `E_R` | 900 | MPa |
| `E_T` | 600 | MPa |
| `nu_LR` | 0.35 | – |
| `nu_LT` | 0.40 | – |
| `nu_RT` | 0.30 | – |
| `G_RT` | 55 | MPa |
| `G_TL` | 650 | MPa |
| `G_LR` | 750 | MPa |

Reverse Poisson ratios are derived, not input: `nu_RL = 0.02625`, `nu_TL = 0.02000`, and `nu_TR = 0.20000`.

## Hill initial surface

| Direction | Initial yield | Unit |
|---|---:|---|
| L uniaxial `X_L` | 40.0 | MPa |
| R uniaxial `Y_R` | 6.0 | MPa |
| T uniaxial `Z_T` | 6.0 | MPa |
| RT shear `S_RT` | 3.0 | MPa |
| TL shear `S_TL` | 4.0 | MPa |
| LR shear `S_LR` | 4.5 | MPa |
| reference radius `sigma_0` | 40.0 | MPa |
| isotropic hardening `H_iso` | 2500.0 | MPa |

The surface treats positive and negative stresses identically. The word “yield” here identifies the M0 numerical surface, not measured onset, damage, or a design strength.

## Numerical controls

| Control | Value |
|---|---:|
| Relative residual tolerance | `1e-11` |
| Maximum local Newton iterations | `40` |
| Tangent verification perturbation | `1e-8` strain |

The JSON document is authoritative for execution. This card exists for human review.

