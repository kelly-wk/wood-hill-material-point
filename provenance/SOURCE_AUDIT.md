# Public Source and Provenance Audit

Audit date: 2026-09-29 (Asia/Tokyo)

## Scope decision

This repository implements a deliberately narrow M0 reference model:

- small-strain orthotropic elasticity;
- Hill-type quadratic anisotropic plasticity;
- associated flow and linear isotropic hardening;
- backward-Euler return mapping;
- a Jacobian-consistent algorithmic tangent.

Tension/compression asymmetry, densification history, damage, component finite elements, contact, proprietary solver integration, and physical calibration remain outside the implemented scope.

## Numerical provenance

The executable parameter card is labeled `custom_algorithm_demonstration_not_calibrated`. No numerical value in it is asserted to come from Hill, Xu et al., Oudjene–Khelifa, a design standard, a species database, or a product. Generated CSV, JSON, and SVG artifacts are deterministic outputs of [`scripts/run_verification.py`](../scripts/run_verification.py).

## Primary method sources

- Hill's original quadratic anisotropic yield theory: <https://doi.org/10.1098/rspa.1948.0045>
- Simo and Taylor on consistent tangents for return-mapping algorithms: <https://doi.org/10.1016/0045-7825(85)90070-2>
- Xu et al. on a timber dowel embedment FE context using anisotropic elastoplasticity: <https://doi.org/10.1007/s10086-012-1296-0>
- Oudjene and Khelifa on anisotropic wood plasticity coupled to densification: <https://doi.org/10.1016/j.conbuildmat.2009.06.034>

These sources motivate the theory and longer-term roadmap. This project does not reproduce their figures, tables, experiments, parameters, finite-element models, or validation claims.

## License and disclosure boundary

New code, documentation, generated records, and original SVG figures are MIT-licensed. Cited publications and proprietary solver documentation retain their own rights and are not redistributed. The repository contains no raw experiments, personal identifiers, local machine paths, credentials, commercial solver files, or calibrated design values.
