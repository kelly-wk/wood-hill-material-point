# Data card: generated material-point verification records

## Summary

This project contains no empirical dataset. Every CSV row is generated deterministically by the local M0 solver from the invented parameter card. No external source artifact is transformed into a training or benchmark dataset.

Project status is **validated-prototype**. Generated records demonstrate internal algorithm behavior; they are not material calibration and are not an engineering safety model.

## Record unit and files

The main unit in [`results/loading_paths.csv`](results/loading_paths.csv) is one accepted material-point increment. Fields include case/segment, controlled strain, its conjugate stress, Hill equivalent stress, yield radius, hardening history, plastic multiplier, yield residual, Newton residual, accumulated plastic dissipation, stored elastic energy, lateral-stress residual, branch flag, and iteration count.

Aggregate files:

- `path_summary.csv`: one row per L/R/T, shear, or cyclic path;
- `material_checks.csv`: stiffness, reciprocity, Hill convexity/calibration, and hydrostatic checks;
- `tangent_checks.csv`: analytic-versus-finite-difference tangent metrics;
- `step_convergence.csv`: non-proportional step-refinement endpoints;
- `cycle_checks.csv`: unloading and monotonic-history assertions;
- `extension_roadmap.csv`: M1/M2/M3 contracts explicitly marked `not_implemented`;
- `summary.json`: acceptance thresholds and claim boundary.

## Generation

Run `PYTHONPATH=. python3 scripts/run_verification.py`. No randomness is used. Floating-point results were produced with Python 3.9.6 and NumPy 2.0.2 on the recorded environment. Tiny final-digit changes may occur with a different BLAS/LAPACK implementation; tests use tolerances well above roundoff.

## Voigt and units

- order: `LL, RR, TT, RT, TL, LR`;
- normal/shear stresses and all moduli: MPa;
- normal and engineering shear strains: dimensionless;
- plastic work density: MPa, numerically equal to MJ/m³;
- compression: negative in the generated axial paths.

## Parameter origin

All values in [`parameters/demo_wood.json`](parameters/demo_wood.json) are custom algorithm-demonstration values. They were not copied from Xu et al., a standard, a species database, or a product sheet. They are selected to create a positive-definite orthotropic stiffness and convex, strongly anisotropic Hill surface.

## Appropriate use

- unit and regression testing for return mapping;
- studying work-conjugate Voigt conversions;
- verifying material-axis rotations;
- teaching grouped verification of yield residual, dissipation, tangent, and step convergence;
- serving as a neutral starting point for an independently validated FE port.

## Inappropriate use

- estimating any real wood property or connection capacity;
- FE validation, contact, screw, dowel, embedment, insertion, or withdrawal claims;
- failure, damage, densification, fracture, or localization prediction;
- structural design, code compliance, or safety decisions;
- claiming M1/M2/M3 or Marc/HYPELA2 is implemented.

## Privacy, copyright, and licensing

There are no people, locations, experimental specimens, or exact source rows in the generated records. Identity-bearing draft text, source figures, tables, and proprietary manual/code content are excluded. Generated records and original project code/documentation are MIT-licensed; cited publications retain their own rights.
