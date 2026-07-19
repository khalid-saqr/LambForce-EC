from pathlib import Path
import yaml

from lambforce_ec.protocol import (
    REQUIRED_PARAMETER_COLUMNS,
    validate_parameter_registry,
    validate_traceability_matrix,
)

ROOT = Path(__file__).resolve().parents[1]


def test_traceability_matrix_is_complete():
    matrix = yaml.safe_load((ROOT / "protocol/readme_traceability.yaml").read_text())
    counts = validate_traceability_matrix(matrix)
    assert counts["implemented"] >= 30
    assert counts["blocked_on_archive"] == 2
    assert {row["requirement_id"] for row in matrix["requirements"]} == {
        f"R{i:02d}" for i in range(1, 35)
    }


def test_parameter_registry_has_full_readme_schema():
    assert validate_parameter_registry(ROOT / "registry/parameter_registry.csv") >= 25
    header = (ROOT / "registry/parameter_registry.csv").read_text().splitlines()[0].split(",")
    assert REQUIRED_PARAMETER_COLUMNS <= set(header)


def test_decision_gate_registry():
    gates = yaml.safe_load((ROOT / "protocol/decision_gates.yaml").read_text())
    assert gates["protocol_version"] == "2.0.0"
    assert len(gates["gates"]) == 9
    assert gates["gates"]["cross_artery_generality"]["minimum_reference_arteries"] == 4
    assert len(gates["allowed_classifications"]) == 5
