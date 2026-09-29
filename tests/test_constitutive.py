from __future__ import annotations

import json
import math
import unittest
from pathlib import Path

import numpy as np

from wood_hill import (
    HillMaterial,
    HillYield,
    OrthotropicElasticity,
    finite_difference_tangent,
    load_material,
    run_component_path,
    run_uniaxial_path,
)
from wood_hill.conventions import (
    rotation_about_global_t,
    strain_from_tensor,
    strain_global_to_material,
    strain_material_to_global,
    strain_to_tensor,
    stress_from_tensor,
    stress_material_to_global,
    stress_to_tensor,
)
from wood_hill.extensions import extension_roadmap


ROOT = Path(__file__).resolve().parents[1]


class ConstitutiveTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.material = load_material(ROOT / "parameters" / "demo_wood.json")

    def test_voigt_tensor_round_trip_and_work_conjugacy(self) -> None:
        stress = np.array([3.0, -2.0, 1.0, 0.7, -0.4, 0.2])
        strain = np.array([0.01, -0.02, 0.03, 0.04, -0.05, 0.06])
        np.testing.assert_allclose(stress_from_tensor(stress_to_tensor(stress)), stress)
        np.testing.assert_allclose(strain_from_tensor(strain_to_tensor(strain)), strain)
        tensor_work = float(np.sum(stress_to_tensor(stress) * strain_to_tensor(strain)))
        self.assertAlmostEqual(float(stress @ strain), tensor_work, places=14)

    def test_rotation_preserves_work(self) -> None:
        rotation = rotation_about_global_t(math.radians(37.0))
        stress_m = np.array([3.0, -2.0, 1.0, 0.7, -0.4, 0.2])
        strain_m = np.array([0.01, -0.02, 0.03, 0.04, -0.05, 0.06])
        stress_g = stress_material_to_global(stress_m, rotation)
        strain_g = strain_material_to_global(strain_m, rotation)
        self.assertAlmostEqual(float(stress_m @ strain_m), float(stress_g @ strain_g), places=13)
        np.testing.assert_allclose(strain_global_to_material(strain_g, rotation), strain_m, atol=1e-14)

    def test_elastic_matrix_is_symmetric_positive_definite_and_reciprocal(self) -> None:
        checks = self.material.elasticity.checks()
        self.assertGreater(checks["compliance_min_eigenvalue"], 0.0)
        self.assertGreater(checks["stiffness_min_eigenvalue"], 0.0)
        self.assertLess(checks["stiffness_symmetry_error"], 1e-10)
        self.assertLess(checks["max_reciprocity_error"], 1e-15)
        unit_L_stress = np.array([1.0, 0, 0, 0, 0, 0])
        strain = self.material.elasticity.compliance @ unit_L_stress
        self.assertAlmostEqual(-strain[1] / strain[0], 0.35, places=13)
        self.assertAlmostEqual(-strain[2] / strain[0], 0.40, places=13)

    def test_non_positive_definite_elasticity_is_rejected(self) -> None:
        with self.assertRaises(ValueError):
            OrthotropicElasticity(1000, 1000, 1000, 0.9, 0.9, 0.9, 100, 100, 100)

    def test_hill_directional_calibration_and_hydrostatic_null(self) -> None:
        errors = self.material.hill.calibration_errors()
        self.assertLess(max(abs(value) for value in errors.values()), 1e-10)
        hydrostatic = np.array([12.0, 12.0, 12.0, 0.0, 0.0, 0.0])
        self.assertLess(self.material.hill.equivalent_stress(hydrostatic), 1e-10)
        self.assertGreaterEqual(np.linalg.eigvalsh(self.material.hill.matrix)[0], -1e-10)

    def test_non_convex_hill_calibration_is_rejected(self) -> None:
        with self.assertRaises(ValueError):
            HillYield(40.0, 6.0, 5.0, 3.0, 4.0, 4.5)

    def test_elastic_update_is_exact(self) -> None:
        state = self.material.initial_state()
        increment = np.array([1e-5, -2e-5, 1e-5, 2e-5, 0.0, 0.0])
        result = self.material.integrate(state, increment)
        self.assertTrue(result.elastic)
        np.testing.assert_allclose(
            result.state.stress,
            self.material.elasticity.stiffness @ increment,
            rtol=1e-14,
            atol=1e-14,
        )
        np.testing.assert_allclose(result.tangent, self.material.elasticity.stiffness)
        self.assertEqual(result.state.kappa, 0.0)

    def test_plastic_return_mapping_consistency_and_dissipation(self) -> None:
        increment = np.array([0.0, -0.015, 0.0, 0.003, 0.002, 0.0])
        result = self.material.integrate(self.material.initial_state(), increment)
        self.assertFalse(result.elastic)
        self.assertGreater(result.plastic_multiplier, 0.0)
        self.assertLess(abs(result.yield_residual), 1e-8)
        self.assertLess(result.newton_residual, 1e-10)
        self.assertGreaterEqual(result.plastic_dissipation_increment, 0.0)
        self.assertGreaterEqual(self.material.elastic_energy(result.state), 0.0)
        self.assertLessEqual(result.iterations, 6)

    def test_algorithmic_tangent_matches_finite_difference_elastic_and_plastic(self) -> None:
        cases = (
            np.array([1e-5, -2e-5, 1e-5, 2e-5, 0.0, 0.0]),
            np.array([0.0, -0.015, 0.0, 0.003, 0.002, 0.0]),
        )
        for increment in cases:
            with self.subTest(increment=increment.tolist()):
                result = self.material.integrate(self.material.initial_state(), increment)
                fd = finite_difference_tangent(
                    self.material, self.material.initial_state(), increment, step=1e-8
                )
                relative = np.linalg.norm(result.tangent - fd) / np.linalg.norm(fd)
                self.assertLess(relative, 1e-7)
                self.assertLess(np.max(np.abs(result.tangent - result.tangent.T)), 1e-9)

    def test_material_axis_uniaxial_paths(self) -> None:
        targets = (-0.010, -0.025, -0.030)
        strengths = (40.0, 6.0, 6.0)
        for axis, (target, strength) in enumerate(zip(targets, strengths)):
            with self.subTest(axis=axis):
                points = run_uniaxial_path(
                    self.material, axis, np.linspace(0.0, target, 121)[1:]
                )
                plastic = [point for point in points if not point.result.elastic]
                self.assertTrue(plastic)
                onset = abs(plastic[0].state.stress[axis])
                # The reported onset is the first plastic point on a finite grid.
                self.assertLess(abs(onset - strength), 0.25)
                lateral = [i for i in (0, 1, 2) if i != axis]
                self.assertLess(
                    max(np.max(np.abs(point.state.stress[lateral])) for point in points),
                    1e-8,
                )
                self.assertLess(
                    max(abs(point.result.yield_residual) for point in plastic), 1e-8
                )

    def test_pure_shear_paths(self) -> None:
        targets = (0.080, 0.025, 0.022)
        strengths = (3.0, 4.0, 4.5)
        for component, target, strength in zip((3, 4, 5), targets, strengths):
            with self.subTest(component=component):
                points = run_component_path(
                    self.material, component, np.linspace(0.0, target, 121)[1:]
                )
                plastic = [point for point in points if not point.result.elastic]
                self.assertTrue(plastic)
                self.assertLess(abs(plastic[0].state.stress[component] - strength), 0.05)
                self.assertGreaterEqual(min(point.result.plastic_dissipation_increment for point in points), 0.0)

    def test_load_unload_reload_preserves_history_on_elastic_unload(self) -> None:
        load = np.linspace(0.0, -0.025, 81)[1:]
        unload = np.linspace(-0.025, -0.013, 49)[1:]
        reload = np.linspace(-0.013, -0.030, 69)[1:]
        points = run_uniaxial_path(self.material, 1, np.concatenate((load, unload, reload)))
        load_end = points[len(load) - 1].state
        unload_points = points[len(load):len(load) + len(unload)]
        unload_end = unload_points[-1].state
        self.assertTrue(all(point.result.elastic for point in unload_points))
        self.assertAlmostEqual(load_end.kappa, unload_end.kappa, places=13)
        kappas = np.array([0.0] + [point.state.kappa for point in points])
        self.assertGreaterEqual(np.min(np.diff(kappas)), -1e-14)
        self.assertGreater(points[-1].state.kappa, load_end.kappa)

    def test_rotated_material_axis_covariance_and_global_tangent(self) -> None:
        rotation = rotation_about_global_t(math.radians(37.0))
        increment_m = np.array([-0.008, -0.012, 0.002, 0.040, 0.005, -0.003])
        increment_g = strain_material_to_global(increment_m, rotation)
        direct = self.material.integrate(self.material.initial_state(), increment_m)
        rotated, stress_g, tangent_g = self.material.integrate_global(
            self.material.initial_state(), increment_g, rotation
        )
        np.testing.assert_allclose(rotated.state.stress, direct.state.stress, atol=1e-11)
        np.testing.assert_allclose(
            stress_g, stress_material_to_global(direct.state.stress, rotation), atol=1e-11
        )
        global_fd = np.empty((6, 6))
        for i in range(6):
            perturb = np.zeros(6)
            perturb[i] = 1e-8
            plus = self.material.integrate_global(
                self.material.initial_state(), increment_g + perturb, rotation
            )[1]
            minus = self.material.integrate_global(
                self.material.initial_state(), increment_g - perturb, rotation
            )[1]
            global_fd[:, i] = (plus - minus) / 2e-8
        self.assertLess(np.linalg.norm(tangent_g - global_fd) / np.linalg.norm(global_fd), 1e-7)

    def test_nonproportional_step_refinement_reduces_endpoint_error(self) -> None:
        def endpoint(n: int):
            state = self.material.initial_state()
            previous = np.zeros(6)
            for t in np.linspace(0.0, 1.0, n + 1)[1:]:
                current = np.array(
                    [-0.008 * t, -0.018 * t, 0.003 * t, 0.050 * t * t,
                     0.006 * np.sin(np.pi * t), -0.004 * t * t]
                )
                state = self.material.integrate(state, current - previous).state
                previous = current
            return state

        reference = endpoint(320)
        coarse = endpoint(20)
        medium = endpoint(40)
        coarse_error = np.linalg.norm(coarse.stress - reference.stress)
        medium_error = np.linalg.norm(medium.stress - reference.stress)
        self.assertLess(medium_error, 0.6 * coarse_error)

    def test_extension_contracts_are_explicitly_unimplemented(self) -> None:
        contracts = extension_roadmap()
        self.assertEqual([entry.stage for entry in contracts], ["M1", "M2", "M3"])
        self.assertTrue(all(entry.status == "not_implemented" for entry in contracts))
        self.assertTrue(all(entry.required_verification_gates for entry in contracts))


if __name__ == "__main__":
    unittest.main()
