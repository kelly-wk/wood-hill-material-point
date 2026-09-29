"""Explicitly unimplemented extension contracts.

These records prevent roadmap items from being mistaken for solver features.
They do not alter M0 stress integration.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class ExtensionContract:
    stage: str
    status: str
    proposed_state_variables: tuple[str, ...]
    required_verification_gates: tuple[str, ...]


def extension_roadmap() -> tuple[ExtensionContract, ...]:
    return (
        ExtensionContract(
            stage="M1",
            status="not_implemented",
            proposed_state_variables=(),
            required_verification_gates=(
                "independent L/R/T tension and compression calibration",
                "convexity check for tension-compression-asymmetric surface",
                "cyclic reversal and finite-difference tangent tests",
            ),
        ),
        ExtensionContract(
            stage="M2",
            status="not_implemented",
            proposed_state_variables=("densification_history",),
            required_verification_gates=(
                "transverse compression-unloading experiments",
                "non-negative dissipation for evolving stiffness",
                "state transfer across insertion and withdrawal analyses",
            ),
        ),
        ExtensionContract(
            stage="M3",
            status="not_implemented",
            proposed_state_variables=("damage_L", "damage_R", "damage_T", "damage_shear"),
            required_verification_gates=(
                "objective regularization or fracture-energy calibration",
                "mesh sensitivity study",
                "independent crack-location or full-field validation",
            ),
        ),
    )

