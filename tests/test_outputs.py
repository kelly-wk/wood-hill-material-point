from __future__ import annotations

import csv
import hashlib
import json
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 16), b""):
            digest.update(block)
    return digest.hexdigest()


class OutputTests(unittest.TestCase):
    def test_required_artifacts_exist_and_are_nonempty(self) -> None:
        required = [
            "README.md",
            "REPORT.md",
            "DATA_CARD.md",
            "LICENSE",
            "Makefile",
            "requirements.txt",
            "parameters/demo_wood.json",
            "provenance/SOURCE_AUDIT.md",
            "docs/HYPELA2_NEUTRAL_INTERFACE.md",
            "results/summary.json",
            "results/material_checks.csv",
            "results/path_summary.csv",
            "results/loading_paths.csv",
            "results/tangent_checks.csv",
            "results/step_convergence.csv",
            "results/cycle_checks.csv",
            "results/extension_roadmap.csv",
            "results/artifact_manifest.csv",
            "figures/uniaxial_material_axes.svg",
            "figures/pure_shear_paths.svg",
            "figures/load_unload_reload.svg",
            "figures/step_convergence.svg",
        ]
        missing = [path for path in required if not (ROOT / path).is_file()]
        self.assertFalse(missing, f"missing artifacts: {missing}")
        self.assertTrue(all((ROOT / path).stat().st_size > 0 for path in required))

    def test_summary_meets_numerical_acceptance_bounds(self) -> None:
        summary = json.loads((ROOT / "results" / "summary.json").read_text())
        self.assertEqual(summary["project_status"], "validated-prototype")
        self.assertEqual(summary["implemented_scope"], "M0 material point only")
        self.assertEqual(summary["parameter_status"], "custom_algorithm_demonstration_not_calibrated")
        self.assertLess(summary["solver"]["max_abs_plastic_yield_residual_MPa"], 1e-8)
        self.assertLess(summary["solver"]["max_scaled_newton_residual"], 1e-10)
        self.assertLess(summary["solver"]["max_fd_tangent_relative_error"], 1e-7)
        self.assertLess(summary["rotation"]["work_conjugacy_error_MPa"], 1e-12)
        self.assertLess(summary["rotation"]["global_tangent_relative_fd_error"], 1e-7)
        self.assertGreater(summary["paths"]["nonproportional_step_convergence_order"], 0.8)
        self.assertLess(summary["paths"]["nonproportional_step_convergence_order"], 1.3)
        self.assertGreater(summary["paths"]["step_10_to_640_stress_error_reduction"], 50)
        self.assertEqual(summary["paths"]["minimum_plastic_dissipation_increment_MPa"], 0.0)
        for stage in ("M1", "M2", "M3", "Marc_or_HYPELA2_subroutine"):
            self.assertIn(stage, summary["not_implemented"])

    def test_machine_tables_cover_required_cases(self) -> None:
        with (ROOT / "results" / "path_summary.csv").open(newline="") as stream:
            rows = list(csv.DictReader(stream))
        cases = {row["case"] for row in rows}
        expected = {
            "uniaxial_L", "uniaxial_R", "uniaxial_T", "pure_shear_RT",
            "pure_shear_TL", "pure_shear_LR", "R_load_unload_reload",
        }
        self.assertEqual(cases, expected)
        self.assertTrue(all(float(row["minimum_plastic_dissipation_increment_MPa"]) >= 0 for row in rows))
        self.assertTrue(all(float(row["max_abs_plastic_yield_residual_MPa"]) < 1e-8 for row in rows))

    def test_artifact_manifest_hashes(self) -> None:
        with (ROOT / "results" / "artifact_manifest.csv").open(newline="") as stream:
            rows = list(csv.DictReader(stream))
        self.assertTrue(rows)
        for row in rows:
            path = ROOT / row["path"]
            self.assertTrue(path.is_file(), row["path"])
            self.assertEqual(int(row["bytes"]), path.stat().st_size)
            self.assertEqual(row["sha256"], sha256(path))

    def test_parameter_card_declares_custom_noncalibrated_values(self) -> None:
        document = json.loads((ROOT / "parameters" / "demo_wood.json").read_text())
        self.assertEqual(
            document["parameter_status"],
            "custom_algorithm_demonstration_not_calibrated",
        )
        self.assertIn("invented", document["disclaimer"].lower())

    def test_documentation_preserves_claim_boundary(self) -> None:
        text = "\n".join(
            (ROOT / path).read_text(encoding="utf-8")
            for path in ("README.md", "REPORT.md", "DATA_CARD.md", "provenance/SOURCE_AUDIT.md")
        )
        required_phrases = (
            "validated-prototype",
            "not material calibration",
            "not implemented",
            "not an engineering safety model",
            "tension/compression symmetric",
        )
        for phrase in required_phrases:
            self.assertIn(phrase.lower(), text.lower())

    def test_no_private_source_or_identity_locator_is_packaged(self) -> None:
        forbidden_suffixes = {".pdf", ".doc", ".docx", ".zip"}
        self.assertFalse(
            [path for path in ROOT.rglob("*") if path.suffix.lower() in forbidden_suffixes]
        )
        for path in ROOT.rglob("*"):
            if (not path.is_file() or path.suffix in {".svg", ".pyc"}
                    or "__pycache__" in path.parts):
                continue
            text = path.read_text(encoding="utf-8", errors="ignore")
            self.assertNotIn("drive.google.com/" + "file/d/", text)
            self.assertNotIn("/" + "Users/", text)


if __name__ == "__main__":
    unittest.main()
