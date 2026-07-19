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

from .exceptions import ProvenanceError, ValidationError
from .models import ArteryRecord
from .provenance import verify_record_payload

REQUIRED_PARAMETER_COLUMNS = {
    "parameter_id", "symbol", "description", "si_unit", "frozen_value", "frozen_lower",
    "frozen_upper", "frozen_set", "source_doi", "source_url", "source_cell_type",
    "source_vascular_bed", "measurement_method", "source_strength",
    "independent_or_derived", "correlation_group", "primary_or_secondary", "claim_bearing",
    "transformation_rule", "notes",
}
REFERENCE_ARTERIES = {
    "aortic_root": (0.0150, 22.03), "thoracic_aorta": (0.0120, 17.62),
    "femoral": (0.0040, 5.87), "carotid": (0.0035, 5.14),
    "iliac": (0.0045, 6.61), "brachial": (0.0020, 2.94),
}
_SHA256 = re.compile(r"^[0-9a-f]{64}$")
_GIT_SHA = re.compile(r"^[0-9a-f]{40}$")


def data_roots() -> list[Path]:
    candidates = [Path(__file__).resolve().parents[2], Path(sysconfig.get_path("data")) / "share/lambforce_ec"]
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
            raise ValidationError(f"Parameter registry missing columns: {', '.join(sorted(missing))}")
        rows = list(reader)
    ids = [row["parameter_id"].strip() for row in rows]
    if not rows or any(not value for value in ids) or len(ids) != len(set(ids)):
        raise ValidationError("Parameter registry must be non-empty with unique parameter_id values.")
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
    if {item.get("artery_id") for item in arteries if isinstance(item, dict)} != set(REFERENCE_ARTERIES):
        raise ValidationError("Reference artery IDs differ from the frozen set.")
    for item in arteries:
        artery_id = str(item["artery_id"])
        radius, alpha = REFERENCE_ARTERIES[artery_id]
        if not np.isclose(float(item.get("radius_m")), radius, rtol=0, atol=1e-15):
            raise ValidationError(f"{artery_id}.radius_m differs from README.")
        if not np.isclose(float(item.get("womersley_alpha")), alpha, rtol=0, atol=1e-12):
            raise ValidationError(f"{artery_id}.womersley_alpha differs from README.")
        amplitudes = item.get("pressure_harmonic_amplitudes")
        if not isinstance(amplitudes, list) or len(amplitudes) != 6 or any(not np.isfinite(float(v)) for v in amplitudes):
            raise ValidationError(f"{artery_id} requires six finite harmonic amplitudes.")
    return 6


def source_records(registry: Mapping[str, Any]) -> list[dict[str, Any]]:
    value = registry.get("sources")
    if not isinstance(value, list):
        raise ValidationError("Source registry must contain a sources list.")
    return value


def validate_source_registry(registry: Mapping[str, Any]) -> int:
    if not isinstance(registry.get("registry_version"), str):
        raise ValidationError("Source registry requires registry_version.")
    seen: set[tuple[str, str]] = set()
    statuses = {"awaiting_archive", "verified", "synthetic_non_claim_bearing"}
    required = {"source_identifier", "source_version", "artery_ids", "archive_status", "archive_sha256",
                "source_member_sha256_by_artery", "conversion_manifest_sha256_by_artery",
                "converter_commit_sha", "role", "claim_bearing"}
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
        if item["archive_status"] not in statuses or not isinstance(item["claim_bearing"], bool):
            raise ValidationError("Invalid source status or claim_bearing value.")
        role = item["role"]
        if role not in {"immutable_hydrodynamic_ground_truth", "software_validation_only"}:
            raise ValidationError("Invalid source role.")
        if item["claim_bearing"] and (role != "immutable_hydrodynamic_ground_truth" or set(arteries) != set(REFERENCE_ARTERIES)):
            raise ValidationError("Claim-bearing source must be the immutable six-artery archive.")
        if role == "software_validation_only" and item["claim_bearing"]:
            raise ValidationError("Software-validation source cannot be claim-bearing.")
        if item["archive_status"] == "verified":
            if not isinstance(item["archive_sha256"], str) or not _SHA256.match(item["archive_sha256"]):
                raise ValidationError("Verified source requires archive_sha256.")
            if not isinstance(item["converter_commit_sha"], str) or not _GIT_SHA.match(item["converter_commit_sha"]):
                raise ValidationError("Verified source requires converter_commit_sha.")
            for field in ("source_member_sha256_by_artery", "conversion_manifest_sha256_by_artery"):
                values = item[field]
                if not isinstance(values, dict) or set(values) != set(arteries) or any(not _SHA256.match(str(v)) for v in values.values()):
                    raise ValidationError(f"Verified source requires a complete {field} map.")
    return len(source_records(registry))


def _harmonics(value: Any, label: str) -> None:
    if not isinstance(value, dict):
        raise ProvenanceError(f"archive_harmonics.{label} must be a mapping.")
    magnitude, phase = value.get("magnitude"), value.get("phase_rad")
    if not isinstance(magnitude, list) or not isinstance(phase, list) or not magnitude or len(magnitude) != len(phase):
        raise ProvenanceError(f"archive_harmonics.{label} arrays are invalid.")
    if any(not np.isfinite(float(v)) for v in magnitude + phase):
        raise ProvenanceError(f"archive_harmonics.{label} contains non-finite values.")


def validate_claim_bearing_hydrodynamic_contract(record: ArteryRecord) -> None:
    if record.artery_id not in REFERENCE_ARTERIES:
        raise ProvenanceError("Claim-bearing execution is restricted to the six registered arteries.")
    if record.lamb_density_isotropic_n_m3 is None:
        raise ProvenanceError("Claim-bearing record requires the archived isotropic field.")
    contract = record.metadata.get("hydrodynamic_contract")
    required = {"rho_kg_m3", "nu_zz_m2_s", "reference_area_m2", "fluid_integration_depth_m",
                "womersley_alpha", "harmonic_rms_tolerance", "source_waveform_identifier",
                "source_archive_member", "radial_collocation_order", "archive_harmonics"}
    if not isinstance(contract, dict) or required - set(contract):
        raise ProvenanceError("Claim-bearing hydrodynamic contract is incomplete.")
    for key in ("rho_kg_m3", "nu_zz_m2_s", "reference_area_m2", "fluid_integration_depth_m", "womersley_alpha", "harmonic_rms_tolerance"):
        if not np.isfinite(float(contract[key])) or float(contract[key]) <= 0:
            raise ProvenanceError(f"hydrodynamic_contract.{key} must be positive.")
    if any(not isinstance(contract[key], str) or not contract[key].strip() for key in ("source_waveform_identifier", "source_archive_member")):
        raise ProvenanceError("Hydrodynamic source identifiers must be non-empty.")
    radial = np.asarray(record.radial_coordinate_m, dtype=float)
    if int(contract["radial_collocation_order"]) != radial.size:
        raise ProvenanceError("radial_collocation_order differs from stored grid.")
    tolerance = max(1e-12, 1e-10 * record.radius_m)
    if not np.isclose(radial[-1], record.radius_m, rtol=0, atol=tolerance):
        raise ProvenanceError("Radial grid must terminate at the artery wall.")
    if not np.isclose(radial[0], record.radius_m - float(contract["fluid_integration_depth_m"]), rtol=0, atol=tolerance):
        raise ProvenanceError("Radial grid does not match fluid integration depth.")
    archive = contract["archive_harmonics"]
    labels = {"signed_lamb_load", "isotropic_lamb_load", "wall_shear_stress"}
    if not isinstance(archive, dict) or set(archive) != labels:
        raise ProvenanceError("archive_harmonics must contain signed, isotropic, and WSS records.")
    for label in labels:
        _harmonics(archive[label], label)


def assert_claim_bearing_source(record: ArteryRecord, registry: Mapping[str, Any] | None) -> None:
    verify_record_payload(record)
    validate_claim_bearing_hydrodynamic_contract(record)
    if registry is None:
        raise ProvenanceError("Claim-bearing execution requires the source registry.")
    validate_source_registry(registry)
    matches = [item for item in source_records(registry) if item.get("source_identifier") == record.source_identifier and item.get("source_version") == record.source_version and record.artery_id in item.get("artery_ids", [])]
    if len(matches) != 1:
        raise ProvenanceError("Source identifier/version/artery does not resolve uniquely.")
    item = matches[0]
    checks = (
        (item.get("claim_bearing") is True, "Source is not claim-bearing."),
        (item.get("role") == "immutable_hydrodynamic_ground_truth", "Source role is not immutable ground truth."),
        (item.get("archive_status") == "verified", "Archive is not verified."),
        (item.get("archive_sha256") == record.source_archive_sha256, "Archive checksum mismatch."),
        (item["source_member_sha256_by_artery"].get(record.artery_id) == record.source_member_sha256, "Member checksum mismatch."),
        (item["conversion_manifest_sha256_by_artery"].get(record.artery_id) == record.conversion_manifest_sha256, "Manifest checksum mismatch."),
        (item.get("converter_commit_sha") == record.converter_commit_sha, "Converter commit mismatch."),
    )
    for passed, message in checks:
        if not passed:
            raise ProvenanceError(message)


def git_blob_sha(path: str | Path) -> str:
    data = Path(path).read_bytes()
    return hashlib.sha1(f"blob {len(data)}\0".encode() + data).hexdigest()


def _path(token: str, root: Path) -> Path | None:
    token = token.strip().split("::", 1)[0]
    if "/" in token:
        return root / token.split()[0]
    aliases = {name: f"src/lambforce_ec/{name}.py" for name in ("workflow", "loads", "models", "analysis", "constitutive", "harmonics", "structural", "validation", "protocol", "io")}
    return root / aliases[token.split(".", 1)[0]] if token.split(".", 1)[0] in aliases else None


def validate_traceability_matrix(matrix: Mapping[str, Any], repository_root: str | Path | None = None) -> dict[str, int]:
    rows = matrix.get("requirements")
    required = {"requirement_id", "readme_section", "requirement", "implementation", "tests", "outputs_or_gate", "status"}
    if not isinstance(rows, list) or any(not isinstance(row, dict) or required - set(row) for row in rows):
        raise ValidationError("Traceability matrix is incomplete.")
    ids = [str(row["requirement_id"]) for row in rows]
    if set(ids) != {f"R{i:02d}" for i in range(1, 35)} or len(ids) != len(set(ids)):
        raise ValidationError("Traceability matrix must contain exactly R01 through R34.")
    counts = {"implemented": 0, "blocked_on_archive": 0}
    for row in rows:
        status = str(row["status"])
        if status not in counts:
            raise ValidationError(f"Unsupported traceability status: {status}")
        counts[status] += 1
    if counts != {"implemented": 32, "blocked_on_archive": 2}:
        raise ValidationError("Traceability counts must remain 32 implemented and 2 blocked.")
    root = Path(repository_root) if repository_root is not None else None
    readme = root / "README.md" if root else resolve_data_path("README.md")
    if git_blob_sha(readme) != matrix.get("readme_git_blob_sha"):
        raise ValidationError("Traceability README blob SHA does not match README.")
    if root:
        for row in rows:
            paths = [path for path in (_path(token, root) for token in str(row["implementation"]).split(";")) if path]
            if not paths or any(not path.is_file() for path in paths):
                raise ValidationError(f"{row['requirement_id']} references a missing implementation path.")
            for reference in str(row["tests"]).split(";"):
                parts = reference.strip().split("::", 1)
                test_path = root / parts[0]
                if not test_path.is_file():
                    raise ValidationError(f"{row['requirement_id']} references a missing test path.")
                if len(parts) == 2:
                    tree = ast.parse(test_path.read_text(encoding="utf-8"))
                    symbols = {node.name for node in tree.body if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef))}
                    if parts[1] not in symbols:
                        raise ValidationError(f"{row['requirement_id']} references a missing test symbol.")
    return counts
