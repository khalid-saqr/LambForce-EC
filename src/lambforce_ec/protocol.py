from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping
import csv
import yaml

from .exceptions import ProvenanceError, ValidationError
from .models import ArteryRecord

REQUIRED_PARAMETER_COLUMNS = {
    "parameter_id",
    "symbol",
    "description",
    "si_unit",
    "frozen_value",
    "frozen_lower",
    "frozen_upper",
    "frozen_set",
    "source_doi",
    "source_url",
    "source_cell_type",
    "source_vascular_bed",
    "measurement_method",
    "source_strength",
    "independent_or_derived",
    "correlation_group",
    "primary_or_secondary",
    "claim_bearing",
    "transformation_rule",
    "notes",
}


def load_yaml(path: str | Path) -> dict[str, Any]:
    with Path(path).open("r", encoding="utf-8") as stream:
        value = yaml.safe_load(stream)
    if not isinstance(value, dict):
        raise ValidationError(f"{path} must contain a YAML mapping.")
    return value


def validate_parameter_registry(path: str | Path) -> int:
    with Path(path).open("r", encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream)
        missing = REQUIRED_PARAMETER_COLUMNS - set(reader.fieldnames or [])
        if missing:
            raise ValidationError(f"Parameter registry missing columns: {', '.join(sorted(missing))}")
        rows = list(reader)
    if not rows:
        raise ValidationError("Parameter registry is empty.")
    ids = [row["parameter_id"] for row in rows]
    if len(ids) != len(set(ids)):
        raise ValidationError("Parameter registry contains duplicate parameter_id values.")
    return len(rows)


def source_records(registry: Mapping[str, Any]) -> list[dict[str, Any]]:
    records = registry.get("sources")
    if not isinstance(records, list):
        raise ValidationError("Source registry must contain a sources list.")
    return records


def assert_claim_bearing_source(record: ArteryRecord, registry: Mapping[str, Any] | None) -> None:
    if registry is None:
        raise ProvenanceError("Claim-bearing execution requires the versioned source registry.")
    matches = [
        item
        for item in source_records(registry)
        if item.get("source_identifier") == record.source_identifier
        and item.get("source_version") == record.source_version
        and record.artery_id in item.get("artery_ids", [])
    ]
    if len(matches) != 1:
        raise ProvenanceError("Source identifier/version/artery does not resolve uniquely.")
    item = matches[0]
    if item.get("archive_status") != "verified":
        raise ProvenanceError("Hydrodynamic archive is not yet marked verified.")
    expected = item.get("archive_sha256")
    if not isinstance(expected, str) or expected != record.source_checksum:
        raise ProvenanceError("Artery source checksum does not match the verified archive registry.")


def validate_traceability_matrix(matrix: Mapping[str, Any]) -> dict[str, int]:
    requirements = matrix.get("requirements")
    if not isinstance(requirements, list) or not requirements:
        raise ValidationError("Traceability matrix must contain requirements.")
    seen: set[str] = set()
    counts = {"implemented": 0, "blocked_on_archive": 0}
    required_keys = {
        "requirement_id",
        "readme_section",
        "requirement",
        "implementation",
        "tests",
        "outputs_or_gate",
        "status",
    }
    for item in requirements:
        if not isinstance(item, dict) or required_keys - set(item):
            raise ValidationError("Traceability requirement is incomplete.")
        rid = str(item["requirement_id"])
        if rid in seen:
            raise ValidationError(f"Duplicate traceability requirement_id: {rid}")
        seen.add(rid)
        status = str(item["status"])
        if status not in counts:
            raise ValidationError(f"Unsupported traceability status: {status}")
        counts[status] += 1
    return counts
