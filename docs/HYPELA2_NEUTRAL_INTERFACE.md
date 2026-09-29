# Solver-neutral material-subroutine handoff

This note is pseudocode and a state contract. It is **not** a Marc/HYPELA2 subroutine, does not use a proprietary call signature, and was not compiled or tested in Marc.

## Host assumptions to resolve before any port

The host adapter must independently confirm:

- stress/strain component order and engineering-versus-tensor shear convention;
- whether the host or material routine owns global/material rotation;
- small-strain versus finite-strain flags and objective-rate responsibilities;
- stress and state initialization on restart;
- tangent layout, symmetry expectation, and units;
- whether the increment is total, mechanical, or excludes thermal/hygroscopic parts;
- cutback/error signalling and thread safety.

None of these is inferred from proprietary solver documentation.

## M0 state vector

If the host stores stress separately, a minimal neutral state vector is:

| Zero-based slots | Meaning | Convention |
|---|---|---|
| `0:6` | plastic strain | `[epsp_LL, epsp_RR, epsp_TT, gammap_RT, gammap_TL, gammap_LR]`; engineering plastic shear |
| `6` | `kappa` | non-negative accumulated equivalent plastic multiplier |
| `7` | cumulative plastic dissipation | MPa |

Total strain is retained by the Python reference implementation for energy diagnostics; a host that already provides total strain need not duplicate it. M1/M2/M3 variables are absent because those stages are not implemented. A future schema change must be versioned rather than silently reinterpreting M0 slots.

## Neutral update pseudocode

```text
inputs:
    stress_n in material axes
    state_n = {plastic_strain_n, kappa_n, plastic_dissipation_n}
    engineering_strain_increment in material axes
    parameter card and tolerances

trial = stress_n + C * engineering_strain_increment
f_trial = q_Hill(trial) - (sigma0 + H_iso * kappa_n)

if f_trial <= tolerance:
    stress_np1 = trial
    state_np1 = state_n
    D_alg = C
else:
    unknown x = [stress_np1(6), delta_lambda]
    initialize delta_lambda from a normal-projection estimate
    repeat safeguarded Newton:
        n = d q_Hill / d stress at stress_np1
        dn = d² q_Hill / d stress² at stress_np1
        r_sigma = stress_np1 - trial + delta_lambda * C * n
        r_f = q_Hill(stress_np1)
              - sigma0 - H_iso * (kappa_n + delta_lambda)
        solve exact 7x7 residual Jacobian for correction
        backtrack until residual decreases and delta_lambda >= 0
    plastic_strain_np1 = plastic_strain_n + delta_lambda * n
    kappa_np1 = kappa_n + delta_lambda
    plastic_dissipation_np1 = plastic_dissipation_n
        + stress_np1 dot (delta_lambda * n)
    solve the converged 7x7 Jacobian against [C; 0]
    D_alg = top six rows of that sensitivity

outputs:
    stress_np1, state_np1, D_alg, local status/diagnostics
```

## Required port-level tests

A port is not validated until it reproduces, in the host convention:

1. six elastic basis increments and the complete stiffness matrix;
2. L/R/T uniaxial-stress paths with lateral stresses driven to zero;
3. RT/TL/LR pure shear;
4. an arbitrary rotated material frame with work invariance;
5. load–unload–reload history and non-decreasing `kappa`;
6. yield residual and non-negative dissipation bounds;
7. central finite-difference tangent at elastic, virgin-plastic, and preloaded-plastic states;
8. one-element global Newton convergence and step cutback behavior.

Passing Python material-point tests alone cannot establish any Marc/HYPELA2 or FE claim.
