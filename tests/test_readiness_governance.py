from pathlib import Path
import copy
import csv

import pytest

from lambforce_ec.exceptions import ValidationError
from lambforce_ec.io import config_from_mapping
from lambforce_ec.protocol import (
    load_yaml,
    validate_parameter_registry,
    validate_reference_arteries,
    validate_source_registry,
    validate_traceability_matrix,
)
from lambforce_ec.published_source import (
    validate_published_inputs,
    validate_published_source_binding,
)
from lambforce_ec.validation import sha256_file

ROOT = Path(__file__).resolve().parents[1]


def test_parameter_registry_semantics(tmp_path):
    registry = ROOT / "registry/parameter_registry.csv"
    assert validate_parameter_registry(registry) == 36
    rows = list(csv.DictReader(registry.open(encoding="utf-8", newline="")))
    rows[-2]["independent_or_derived"] = "indepent"
    bad = tmp_path / "bad.csv"
    with bad.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    with pytest.raises(ValidationError):
        validate_parameter_registry(bad)


def test_reference_artery_registry_matches_frozen_readme_values():
    registry = load_yaml(ROOT / "configs/reference_arteries.yaml")
    assert validate_reference_arteries(registry) == 6
    bad = copy.deepcopy(registry)
    bad["arteries"][0]["radius_m"] = 0.014
    with pytest.raises(ValidationError):
        validate_reference_arteries(bad)


def test_published_input_registry_semantics():
    registry = load_yaml(ROOT / "registry/published_v2_hydrodynamics.yaml")
    assert validate_published_inputs(registry) == 6
    bad = copy.deepcopy(registry)
    bad["arteries"][0]["pressure_gradient_scale_pa_per_m"] = 9001.0
    with pytest.raises(ValidationError):
        validate_published_inputs(bad)
    source_registry = load_yaml(ROOT / "registry/source_registry.yaml")
    validate_published_source_binding(
        registry, sha256_file(ROOT / "registry/published_v2_hydrodynamics.yaml"), source_registry
    )


def test_source_registry_semantics():
    registry = load_yaml(ROOT / "registry/source_registry.yaml")
    assert validate_source_registry(registry) == 2
    bad = copy.deepcopy(registry)
    bad["sources"][1]["claim_bearing"] = True
    with pytest.raises(ValidationError):
        validate_source_registry(bad)


def test_traceability_resolves_readme_code_and_tests():
    matrix = load_yaml(ROOT / "protocol/readme_traceability.yaml")
    counts = validate_traceability_matrix(matrix, repository_root=ROOT)
    assert counts == {"implemented": 32, "blocked_on_source_reproduction": 2}
    typo = copy.deepcopy(matrix)
    typo["requirements"][8]["status"] = "implemed"
    with pytest.raises(ValidationError):
        validate_traceability_matrix(typo, repository_root=ROOT)
    stale = copy.deepcopy(matrix)
    stale["readme_git_blob_sha"] = "0" * 40
    with pytest.raises(ValidationError):
        validate_traceability_matrix(stale, repository_root=ROOT)


def test_unknown_or_prohibited_configuration_fields_are_rejected():
    with pytest.raises(ValidationError):
        config_from_mapping({"material": {"lamb_force_transfer_efficiency": 0.5}})
    with pytest.raises(ValidationError):
        config_from_mapping({"outcome_selected_prestress": 1.0})
