from __future__ import annotations

from typing import Any, Mapping
import re

import numpy as np

from .exceptions import ProvenanceError
from .governance import (
    REFERENCE_ARTERIES,
    REFERENCE_MODES,
    REQUIRED_PARAMETER_COLUMNS,
    data_roots,
    git_blob_sha,
    load_default_yaml,
    load_yaml,
    resolve_data_path,
    source_records,
    validate_parameter_registry,
    validate_reference_arteries,
    validate_source_registry,
    validate_traceability_matrix,
)
from .models import ArteryRecord
from .provenance import verify_record_payload

_SHA256 = re.compile(r"^[0-9a-f]{64}$")
_GIT_SHA = re.compile(r"^[0-9a-f]{40}$")


def _harmonics(value: Any, label: str) -> None:
    if not isinstance(value, dict):
        raise ProvenanceError(f"generated_harmonics.{label} must be a mapping.")
    magnitude, phase = value.get("magnitude"), value.get("phase_rad")
    if not isinstance(magnitude, list) or not isinstance(phase, list) or not magnitude:
        raise ProvenanceError(f"generated_harmonics.{label} arrays are invalid.")
    if len(magnitude) != len(phase) or any(not np.isfinite(float(v)) for v in magnitude + phase):
        raise ProvenanceError(f"generated_harmonics.{label} contains invalid values.")


def validate_claim_bearing_hydrodynamic_contract(record: ArteryRecord) -> None:
    if record.artery_id not in REFERENCE_ARTERIES:
        raise ProvenanceError(
            "Claim-bearing execution is restricted to the six registered arteries."
        )
    if record.lamb_density_isotropic_n_m3 is None:
        raise ProvenanceError("Claim-bearing record requires the published-source isotropic field.")
    published = record.metadata.get("published_source")
    required_published = {
        "repository",
        "repository_commit_sha",
        "published_notebook_path",
        "published_notebook_blob_sha",
        "published_input_registry_sha256",
        "reproduction_commit_sha",
        "reproduction_configuration_sha256",
        "reproduction_mode",
    }
    if not isinstance(published, dict) or required_published - set(published):
        raise ProvenanceError("Published-source provenance is incomplete.")
    if published["reproduction_mode"] not in REFERENCE_MODES:
        raise ProvenanceError("Unknown published-source reproduction mode.")
    if any(
        not _GIT_SHA.match(str(published[key]))
        for key in (
            "repository_commit_sha",
            "published_notebook_blob_sha",
            "reproduction_commit_sha",
        )
    ):
        raise ProvenanceError("Published-source Git identities are invalid.")
    if any(
        not _SHA256.match(str(published[key]))
        for key in ("published_input_registry_sha256", "reproduction_configuration_sha256")
    ):
        raise ProvenanceError("Published-source SHA-256 identities are invalid.")
    contract = record.metadata.get("hydrodynamic_contract")
    required = {
        "rho_kg_m3",
        "nu_zz_m2_s",
        "reference_area_m2",
        "fluid_integration_depth_m",
        "womersley_alpha",
        "published_womersley_alpha",
        "harmonic_rms_tolerance",
        "source_waveform_identifier",
        "source_notebook_path",
        "source_repository_commit_sha",
        "source_notebook_blob_sha",
        "radial_collocation_order",
        "reproduction_mode",
        "generated_harmonics",
    }
    if not isinstance(contract, dict) or required - set(contract):
        raise ProvenanceError("Claim-bearing hydrodynamic contract is incomplete.")
    for key in (
        "rho_kg_m3",
        "nu_zz_m2_s",
        "reference_area_m2",
        "fluid_integration_depth_m",
        "womersley_alpha",
        "published_womersley_alpha",
        "harmonic_rms_tolerance",
    ):
        if not np.isfinite(float(contract[key])) or float(contract[key]) <= 0:
            raise ProvenanceError(f"hydrodynamic_contract.{key} must be positive.")
    expected_alpha = REFERENCE_ARTERIES[record.artery_id][1]
    if not np.isclose(
        float(contract["published_womersley_alpha"]), expected_alpha, rtol=0.0, atol=1e-12
    ):
        raise ProvenanceError("Published Womersley alpha differs from the frozen Table 1 value.")
    if contract["reproduction_mode"] != published["reproduction_mode"]:
        raise ProvenanceError("Hydrodynamic and provenance reproduction modes disagree.")
    if contract["source_repository_commit_sha"] != published["repository_commit_sha"]:
        raise ProvenanceError("Hydrodynamic source commit disagrees with provenance.")
    if contract["source_notebook_blob_sha"] != published["published_notebook_blob_sha"]:
        raise ProvenanceError("Hydrodynamic notebook blob disagrees with provenance.")
    radial = np.asarray(record.radial_coordinate_m, dtype=float)
    if int(contract["radial_collocation_order"]) != radial.size:
        raise ProvenanceError("radial_collocation_order differs from stored grid.")
    tolerance = max(1e-12, 1e-10 * record.radius_m)
    if not np.isclose(radial[-1], record.radius_m, rtol=0, atol=tolerance):
        raise ProvenanceError("Radial grid must terminate at the artery wall.")
    if not np.isclose(
        radial[0],
        record.radius_m - float(contract["fluid_integration_depth_m"]),
        rtol=0,
        atol=tolerance,
    ):
        raise ProvenanceError("Radial grid does not match the published control-volume depth.")
    harmonics = contract["generated_harmonics"]
    labels = {"signed_lamb_load", "isotropic_lamb_load", "wall_shear_stress"}
    if not isinstance(harmonics, dict) or set(harmonics) != labels:
        raise ProvenanceError(
            "generated_harmonics must contain signed, isotropic, and WSS records."
        )
    for label in labels:
        _harmonics(harmonics[label], label)


def assert_claim_bearing_source(record: ArteryRecord, registry: Mapping[str, Any] | None) -> None:
    verify_record_payload(record)
    validate_claim_bearing_hydrodynamic_contract(record)
    if registry is None:
        raise ProvenanceError("Claim-bearing execution requires the source registry.")
    validate_source_registry(registry)
    matches = [
        item
        for item in source_records(registry)
        if item.get("source_identifier") == record.source_identifier
        and item.get("source_version") == record.source_version
        and record.artery_id in item.get("artery_ids", [])
    ]
    if len(matches) != 1:
        raise ProvenanceError(
            "Published source identifier/version/artery does not resolve uniquely."
        )
    item = matches[0]
    published = record.metadata["published_source"]
    mode = published["reproduction_mode"]
    payload_map = item.get("record_payload_sha256_by_mode_and_artery", {})
    expected_payload = (
        payload_map.get(mode, {}).get(record.artery_id) if isinstance(payload_map, dict) else None
    )
    checks = (
        (item.get("claim_bearing") is True, "Source is not claim-bearing."),
        (
            item.get("role") == "published_hydrodynamic_ground_truth",
            "Source role is not published ground truth.",
        ),
        (
            item.get("source_status") == "verified_reproduction",
            "Six-artery reproduction is not verified.",
        ),
        (
            item.get("qualification_status") == "verified",
            "Historical-source qualification is not verified.",
        ),
        (
            item.get("repository_commit_sha") == published["repository_commit_sha"],
            "Source repository commit mismatch.",
        ),
        (
            item.get("published_notebook_blob_sha") == published["published_notebook_blob_sha"],
            "Notebook blob mismatch.",
        ),
        (
            item.get("published_input_registry_sha256")
            == published["published_input_registry_sha256"],
            "Published-input registry mismatch.",
        ),
        (
            item.get("reproduction_commit_sha") == published["reproduction_commit_sha"],
            "Reproduction commit mismatch.",
        ),
        (expected_payload == record.record_payload_sha256, "Registered record payload mismatch."),
    )
    for passed, message in checks:
        if not passed:
            raise ProvenanceError(message)


__all__ = [
    "REFERENCE_ARTERIES",
    "REFERENCE_MODES",
    "REQUIRED_PARAMETER_COLUMNS",
    "assert_claim_bearing_source",
    "data_roots",
    "git_blob_sha",
    "load_default_yaml",
    "load_yaml",
    "resolve_data_path",
    "source_records",
    "validate_claim_bearing_hydrodynamic_contract",
    "validate_parameter_registry",
    "validate_reference_arteries",
    "validate_source_registry",
    "validate_traceability_matrix",
]
