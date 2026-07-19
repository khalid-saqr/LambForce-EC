from __future__ import annotations

import ast
import csv
import hashlib
import os
import re
import sysconfig
from pathlib import Path
from typing import Any, Mapping

import numpy as np
import yaml

from .exceptions import ValidationError

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
REFERENCE_ARTERIES = {
    "aortic_root": (0.0150, 22.03),
    "thoracic_aorta": (0.0120, 17.62),
    "femoral": (0.0040, 5.87),
    "carotid": (0.0035, 5.14),
    "iliac": (0.0045, 6.61),
    "brachial": (0.0020, 2.94),
}
REFERENCE_MODES = {"historical_v2", "verified"}
_SHA256 = re.compile(r"^[0-9a-f]{64}$")
_GIT_SHA = re.compile(r"^[0-9a-f]{40}$")


def data_roots() -> list[Path]:
    candidates = [
        Path(__file__).resolve().parents[2],
        Path(sysconfig.get_path("data")) / "share/lambforce_ec",
    ]
    if os.environ.get("LAMBFORCE_EC_DATA_ROOT"):
        candidates.insert(0, Path(os.environ["LAMBFORCE_EC_DATA_ROOT"]))
    result: list[Path] = []
    for path in candidates:
        resolved = path.resolve()
        if resolved not in result:
            result.append(resolved)
    return result


def resolve_data_path(relative_path: str | Path, explicit: str | Path | None = None) -> Path:
    if explicit is not None:
        path = Path(explicit)
        if path.is_file():
            return path
        raise ValidationError(f"Required file does not exist: {path}")
    relative = Path(relative_path)
    for root in data_roots():
        candidate = root / relative
        if candidate.is_file():
            return candidate
    raise ValidationError(f"Could not resolve {relative} from the registered data roots.")


def load_yaml(path: str | Path) -> dict[str, Any]:
    value = yaml.safe_load(Path(path).read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValidationError(f"{path} must contain a YAML mapping.")
    return value


def load_default_yaml(relative_path: str) -> dict[str, Any]:
    return load_yaml(resolve_data_path(relative_path))


def _number(value: str, label: str) -> float | None:
    if not value.strip():
        return None
    try:
        result = float(value)
    except ValueError as exc:
        raise ValidationError(f"{label} must be numeric.") from exc
    if not np.isfinite(result):
        raise ValidationError(f"{label} must be finite.")
    return result


def validate_parameter_registry(path: str | Path) -> int:
    with Path(path).open(encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream)
        missing = REQUIRED_PARAMETER_COLUMNS - set(reader.fieldnames or [])
        if missing:
            raise ValidationError(
                f"Parameter registry missing columns: {', '.join(sorted(missing))}"
            )
        rows = list(reader)
    ids = [row["parameter_id"].strip() for row in rows]
    if not rows or any(not value for value in ids) or len(ids) != len(set(ids)):
        raise ValidationError(
            "Parameter registry must be non-empty with unique parameter_id values."
        )
    grades = {"A", "B", "C", "D", "A-C", "B-C", "C-D", "ground_truth"}
    for row in rows:
        pid = row["parameter_id"].strip()
        for key in ("symbol", "description", "si_unit", "measurement_method", "correlation_group"):
            if not row[key].strip():
                raise ValidationError(f"{pid}.{key} must be non-empty.")
        if row["source_strength"].strip() not in grades:
            raise ValidationError(f"{pid}.source_strength is invalid.")
        independence = row["independent_or_derived"].strip()
        if independence not in {"independent", "derived"}:
            raise ValidationError(f"{pid}.independent_or_derived is invalid.")
        if row["primary_or_secondary"].strip() not in {"primary", "secondary", "feasibility"}:
            raise ValidationError(f"{pid}.primary_or_secondary is invalid.")
        if row["claim_bearing"].strip().lower() not in {"true", "false"}:
            raise ValidationError(f"{pid}.claim_bearing must be Boolean.")
        for key in ("frozen_value", "frozen_lower", "frozen_upper"):
            _number(row[key], f"{pid}.{key}")
        for value in filter(None, row["frozen_set"].split(";")):
            _number(value, f"{pid}.frozen_set")
        if independence == "derived" and not row["transformation_rule"].strip():
            raise ValidationError(f"{pid} is derived but has no transformation_rule.")
        if independence == "independent" and row["source_strength"].strip() != "D":
            if not (row["source_doi"].strip() or row["source_url"].strip()):
                raise ValidationError(f"{pid} has no source DOI or URL.")
    return len(rows)


def validate_reference_arteries(registry: Mapping[str, Any]) -> int:
    source, arteries = registry.get("source"), registry.get("arteries")
    if not isinstance(source, dict) or source.get("doi") != "10.1038/s41598-026-47474-x":
        raise ValidationError("Reference arteries must cite the immutable paper DOI.")
    if not isinstance(arteries, list) or len(arteries) != 6:
        raise ValidationError("Reference registry must contain exactly six arteries.")
    if {item.get("artery_id") for item in arteries if isinstance(item, dict)} != set(
        REFERENCE_ARTERIES
    ):
        raise ValidationError("Reference artery IDs differ from the frozen set.")
    for item in arteries:
        artery_id = str(item["artery_id"])
        radius, alpha = REFERENCE_ARTERIES[artery_id]
        if not np.isclose(float(item.get("radius_m")), radius, rtol=0, atol=1e-15):
            raise ValidationError(f"{artery_id}.radius_m differs from README.")
        if not np.isclose(float(item.get("womersley_alpha")), alpha, rtol=0, atol=1e-12):
            raise ValidationError(f"{artery_id}.womersley_alpha differs from README.")
        amplitudes = item.get("pressure_harmonic_amplitudes")
        if not isinstance(amplitudes, list) or len(amplitudes) != 6:
            raise ValidationError(f"{artery_id} requires six finite harmonic amplitudes.")
        if any(not np.isfinite(float(value)) for value in amplitudes):
            raise ValidationError(f"{artery_id} harmonic amplitudes must be finite.")
    return 6


def source_records(registry: Mapping[str, Any]) -> list[dict[str, Any]]:
    value = registry.get("sources")
    if not isinstance(value, list):
        raise ValidationError("Source registry must contain a sources list.")
    return value


def _complete_payload_map(value: Any, arteries: set[str]) -> bool:
    if not isinstance(value, dict) or set(value) != REFERENCE_MODES:
        return False
    for mode in REFERENCE_MODES:
        row = value[mode]
        if not isinstance(row, dict) or set(row) != arteries:
            return False
        if any(not _SHA256.match(str(digest)) for digest in row.values()):
            return False
    return True


def validate_source_registry(registry: Mapping[str, Any]) -> int:
    if not isinstance(registry.get("registry_version"), str):
        raise ValidationError("Source registry requires registry_version.")
    required = {
        "source_identifier",
        "source_version",
        "artery_ids",
        "source_status",
        "repository_commit_sha",
        "published_notebook_path",
        "published_notebook_blob_sha",
        "published_input_registry_sha256",
        "reproduction_commit_sha",
        "record_payload_sha256_by_mode_and_artery",
        "qualification_status",
        "role",
        "claim_bearing",
    }
    statuses = {"frozen_published_source", "verified_reproduction", "synthetic_non_claim_bearing"}
    qualifications = {
        "awaiting_six_artery_reproduction",
        "awaiting_historical_oracle",
        "verified",
        "not_applicable",
    }
    seen: set[tuple[str, str]] = set()
    for item in source_records(registry):
        if not isinstance(item, dict) or required - set(item):
            raise ValidationError("Source-registry record is incomplete.")
        key = (str(item["source_identifier"]), str(item["source_version"]))
        if key in seen:
            raise ValidationError("Duplicate source identifier/version.")
        seen.add(key)
        arteries = item["artery_ids"]
        if not isinstance(arteries, list) or not arteries or len(arteries) != len(set(arteries)):
            raise ValidationError("artery_ids must be a non-empty unique list.")
        if (
            item["source_status"] not in statuses
            or item["qualification_status"] not in qualifications
        ):
            raise ValidationError("Invalid published-source status.")
        if not isinstance(item["claim_bearing"], bool):
            raise ValidationError("claim_bearing must be Boolean.")
        role = item["role"]
        if role not in {"published_hydrodynamic_ground_truth", "software_validation_only"}:
            raise ValidationError("Invalid source role.")
        if role == "software_validation_only":
            if item["claim_bearing"] or item["source_status"] != "synthetic_non_claim_bearing":
                raise ValidationError("Software-validation source cannot become claim-bearing.")
            continue
        if set(arteries) != set(REFERENCE_ARTERIES) or item["claim_bearing"] is not True:
            raise ValidationError("Published source must cover the frozen six arteries.")
        for field in ("repository_commit_sha", "published_notebook_blob_sha"):
            if not isinstance(item[field], str) or not _GIT_SHA.match(item[field]):
                raise ValidationError(f"Published source requires {field}.")
        if (
            not isinstance(item["published_notebook_path"], str)
            or not item["published_notebook_path"].strip()
        ):
            raise ValidationError("Published source requires published_notebook_path.")
        if not isinstance(item["published_input_registry_sha256"], str) or not _SHA256.match(
            item["published_input_registry_sha256"]
        ):
            raise ValidationError("Published source requires published_input_registry_sha256.")
        if item["source_version"] != item["published_notebook_blob_sha"]:
            raise ValidationError("source_version must equal the frozen notebook blob SHA.")
        if item["source_status"] == "verified_reproduction":
            if item["qualification_status"] != "verified":
                raise ValidationError("Verified reproduction requires verified qualification.")
            if not isinstance(item["reproduction_commit_sha"], str) or not _GIT_SHA.match(
                item["reproduction_commit_sha"]
            ):
                raise ValidationError("Verified reproduction requires reproduction_commit_sha.")
            if not _complete_payload_map(
                item["record_payload_sha256_by_mode_and_artery"], set(arteries)
            ):
                raise ValidationError(
                    "Verified reproduction requires all twelve record payload hashes."
                )
    return len(source_records(registry))


def git_blob_sha(path: str | Path) -> str:
    data = Path(path).read_bytes()
    return hashlib.sha1(f"blob {len(data)}\0".encode() + data).hexdigest()


def _path(token: str, root: Path) -> Path | None:
    token = token.strip().split("::", 1)[0]
    if "/" in token:
        return root / token.split()[0]
    aliases = {
        name: f"src/lambforce_ec/{name}.py"
        for name in (
            "workflow",
            "loads",
            "models",
            "analysis",
            "constitutive",
            "harmonics",
            "structural",
            "validation",
            "protocol",
            "governance",
            "io",
            "published_source",
            "published_contract",
            "published_solver",
        )
    }
    key = token.split(".", 1)[0]
    return root / aliases[key] if key in aliases else None


def validate_traceability_matrix(
    matrix: Mapping[str, Any], repository_root: str | Path | None = None
) -> dict[str, int]:
    rows = matrix.get("requirements")
    required = {
        "requirement_id",
        "readme_section",
        "requirement",
        "implementation",
        "tests",
        "outputs_or_gate",
        "status",
    }
    if not isinstance(rows, list) or any(
        not isinstance(row, dict) or required - set(row) for row in rows
    ):
        raise ValidationError("Traceability matrix is incomplete.")
    ids = [str(row["requirement_id"]) for row in rows]
    if set(ids) != {f"R{i:02d}" for i in range(1, 35)} or len(ids) != len(set(ids)):
        raise ValidationError("Traceability matrix must contain exactly R01 through R34.")
    counts = {"implemented": 0, "blocked_on_source_reproduction": 0}
    for row in rows:
        status = str(row["status"])
        if status not in counts:
            raise ValidationError(f"Unsupported traceability status: {status}")
        counts[status] += 1
    if counts != {"implemented": 32, "blocked_on_source_reproduction": 2}:
        raise ValidationError(
            "Traceability counts must remain 32 implemented and 2 source-blocked."
        )
    root = Path(repository_root) if repository_root is not None else None
    readme = root / "README.md" if root else resolve_data_path("README.md")
    if git_blob_sha(readme) != matrix.get("readme_git_blob_sha"):
        raise ValidationError("Traceability README blob SHA does not match README.")
    if root:
        for row in rows:
            paths = [
                path
                for path in (_path(token, root) for token in str(row["implementation"]).split(";"))
                if path
            ]
            if not paths or any(not path.is_file() for path in paths):
                raise ValidationError(
                    f"{row['requirement_id']} references a missing implementation path."
                )
            for reference in str(row["tests"]).split(";"):
                parts = reference.strip().split("::", 1)
                test_path = root / parts[0]
                if not test_path.is_file():
                    raise ValidationError(
                        f"{row['requirement_id']} references a missing test path."
                    )
                if len(parts) == 2:
                    tree = ast.parse(test_path.read_text(encoding="utf-8"))
                    symbols = {
                        node.name
                        for node in tree.body
                        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef))
                    }
                    if parts[1] not in symbols:
                        raise ValidationError(
                            f"{row['requirement_id']} references a missing test symbol."
                        )
    return counts
