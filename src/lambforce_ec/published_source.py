from __future__ import annotations

import json
from pathlib import Path
import re
from typing import Any, Mapping

import numpy as np
import yaml

from .exceptions import ProvenanceError, ValidationError
from .io import load_artery_npz, save_artery_npz
from .loads import integrate_radial_density
from .models import ArteryRecord
from .provenance import compute_record_payload_sha256
from .published_contract import (
    REFERENCE_ARTERY_IDS,
    REFERENCE_MODES,
    PublishedArteryCase,
    ReproductionProfile,
    artery_case,
    reproduction_profile,
    validate_published_inputs,
    validate_published_source_binding,
)
from .published_solver import PublishedWomersleySolver, classical_womersley, compute_fields
from .validation import canonical_json_bytes, sha256_bytes, sha256_file

_GIT_SHA = re.compile(r"^[0-9a-f]{40}$")
_SHA256 = re.compile(r"^[0-9a-f]{64}$")


def _harmonics(values: np.ndarray) -> dict[str, list[float]]:
    coefficients = np.fft.rfft(np.asarray(values, float)) / values.size
    return {
        "magnitude": np.abs(coefficients).tolist(),
        "phase_rad": np.angle(coefficients).tolist(),
    }


def _relative_rms(value: np.ndarray, reference: np.ndarray) -> float:
    numerator = float(np.sqrt(np.mean((np.asarray(value) - np.asarray(reference)) ** 2)))
    denominator = max(float(np.sqrt(np.mean(np.asarray(reference) ** 2))), 1e-30)
    return numerator / denominator


def reproduce_artery(
    case: PublishedArteryCase,
    mapping: Mapping[str, Any],
    input_registry_sha256: str,
    reproduction_commit_sha: str,
    mode: str,
    profile: ReproductionProfile,
) -> tuple[ArteryRecord, dict[str, Any]]:
    if not _SHA256.match(input_registry_sha256) or not _GIT_SHA.match(reproduction_commit_sha):
        raise ValidationError("Reproduction identities must be complete SHA values.")
    source = mapping["source"]
    anisotropic = compute_fields(case, mapping, mode, profile, isotropic=False)
    isotropic = compute_fields(case, mapping, mode, profile, isotropic=True)
    source_identity = {
        "repository": source["repository"],
        "repository_commit_sha": source["repository_commit_sha"],
        "published_notebook_path": source["published_notebook_path"],
        "published_notebook_blob_sha": source["published_notebook_blob_sha"],
        "published_input_registry_sha256": input_registry_sha256,
    }
    source_snapshot_sha256 = sha256_bytes(canonical_json_bytes(source_identity))
    configuration = {
        "mode": mode,
        "radial_order": profile.radial_order,
        "time_points": profile.time_points,
        "quadrature_nodes": profile.quadrature_nodes,
        "anisotropy": mapping["anisotropy"],
        "fixed_order_interpolation_and_reconstruction": True,
    }
    configuration_sha256 = sha256_bytes(canonical_json_bytes(configuration))
    source_identity_sha256 = sha256_bytes(
        canonical_json_bytes({**source_identity, "kind": "published_source_snapshot"})
    )
    total = integrate_radial_density(
        anisotropic["radial_coordinate_m"], anisotropic["force_density_n_m3"]
    )
    isotropic_load = integrate_radial_density(
        isotropic["radial_coordinate_m"], isotropic["force_density_n_m3"]
    )
    contract = {
        "rho_kg_m3": float(mapping["fluid"]["density_kg_m3"]),
        "nu_zz_m2_s": float(mapping["fluid"]["kinematic_viscosity_m2_s"]),
        "reference_area_m2": float(mapping["control_volume"]["reference_area_m2"]),
        "fluid_integration_depth_m": float(anisotropic["fluid_integration_depth_m"]),
        "womersley_alpha": float(anisotropic["alpha"]),
        "published_womersley_alpha": case.published_womersley_alpha,
        "radial_collocation_order": profile.quadrature_nodes,
        "harmonic_rms_tolerance": float(mapping["numerical"]["harmonic_rms_tolerance"]),
        "source_waveform_identifier": (
            f"{source['repository']}:{source['published_notebook_path']}:{case.artery_id}"
        ),
        "source_notebook_path": source["published_notebook_path"],
        "source_repository_commit_sha": source["repository_commit_sha"],
        "source_notebook_blob_sha": source["published_notebook_blob_sha"],
        "reproduction_mode": mode,
        "generated_harmonics": {
            "signed_lamb_load": _harmonics(total),
            "isotropic_lamb_load": _harmonics(isotropic_load),
            "wall_shear_stress": _harmonics(anisotropic["wall_shear_stress_pa"]),
        },
    }
    record = ArteryRecord(
        artery_id=case.artery_id,
        artery_name=case.artery_name,
        radius_m=case.radius_m,
        omega0_rad_s=float(anisotropic["omega0_rad_s"]),
        radial_coordinate_m=np.asarray(anisotropic["radial_coordinate_m"]),
        time_s=np.asarray(anisotropic["time_s"]),
        lamb_density_signed_n_m3=np.asarray(anisotropic["force_density_n_m3"]),
        wall_shear_stress_pa=np.asarray(anisotropic["wall_shear_stress_pa"]),
        lamb_density_isotropic_n_m3=np.asarray(isotropic["force_density_n_m3"]),
        source_identifier=f"github:{source['repository']}#{source['published_notebook_path']}",
        source_version=str(source["published_notebook_blob_sha"]),
        source_checksum=source_snapshot_sha256,
        source_member_sha256=source_identity_sha256,
        conversion_manifest_sha256=configuration_sha256,
        converter_commit_sha=reproduction_commit_sha,
        metadata={
            "published_source": {
                **source_identity,
                "reproduction_commit_sha": reproduction_commit_sha,
                "reproduction_configuration_sha256": configuration_sha256,
                "reproduction_mode": mode,
            },
            "hydrodynamic_contract": contract,
        },
    )
    record.validate()
    record.record_payload_sha256 = compute_record_payload_sha256(record)
    metrics = {
        "max_backward_residual": anisotropic["max_backward_residual"],
        "computed_womersley_alpha": float(anisotropic["alpha"]),
        "published_womersley_alpha": case.published_womersley_alpha,
        "womersley_alpha_relative_difference": abs(
            float(anisotropic["alpha"]) - case.published_womersley_alpha
        )
        / case.published_womersley_alpha,
        "differentiation_polynomial_error": anisotropic["differentiation_polynomial_error"],
        "isotropic_classical_linf_error": isotropic["isotropic_classical_linf_error"],
        "signed_lamb_peak_abs_pa": float(np.max(np.abs(total))),
        "isotropic_lamb_peak_abs_pa": float(np.max(np.abs(isotropic_load))),
        "anisotropy_increment_peak_abs_pa": float(np.max(np.abs(total - isotropic_load))),
        "wss_peak_abs_pa": float(np.max(np.abs(anisotropic["wall_shear_stress_pa"]))),
    }
    return record, metrics


def reproduce_all_six(
    published_inputs_path: str | Path,
    output_directory: str | Path,
    reproduction_commit_sha: str,
    profile_name: str = "publication",
) -> dict[str, Any]:
    input_path = Path(published_inputs_path)
    mapping = yaml.safe_load(input_path.read_text(encoding="utf-8"))
    if not isinstance(mapping, dict):
        raise ValidationError("Published-input registry must contain a mapping.")
    validate_published_inputs(mapping)
    input_sha = sha256_file(input_path)
    profile = reproduction_profile(mapping, profile_name)
    output = Path(output_directory)
    rows: dict[str, Any] = {}
    records: dict[tuple[str, str], ArteryRecord] = {}
    for item in mapping["arteries"]:
        case = artery_case(item)
        rows[case.artery_id] = {"published_inputs": dict(item), "modes": {}}
        for mode in REFERENCE_MODES:
            record, metrics = reproduce_artery(
                case, mapping, input_sha, reproduction_commit_sha, mode, profile
            )
            target = output / "records" / mode / f"{case.artery_id}.npz"
            save_artery_npz(record, target)
            loaded = load_artery_npz(target)
            records[(case.artery_id, mode)] = loaded
            rows[case.artery_id]["modes"][mode] = {
                "record_path": str(target.relative_to(output)),
                "record_payload_sha256": loaded.record_payload_sha256,
                **metrics,
            }
        historical = records[(case.artery_id, "historical_v2")]
        verified = records[(case.artery_id, "verified")]
        h_total = integrate_radial_density(
            historical.radial_coordinate_m, historical.lamb_density_signed_n_m3
        )
        v_total = integrate_radial_density(
            verified.radial_coordinate_m, verified.lamb_density_signed_n_m3
        )
        h_iso = integrate_radial_density(
            historical.radial_coordinate_m, historical.lamb_density_isotropic_n_m3
        )
        v_iso = integrate_radial_density(
            verified.radial_coordinate_m, verified.lamb_density_isotropic_n_m3
        )
        rows[case.artery_id]["historical_vs_verified"] = {
            "signed_lamb_relative_rms": _relative_rms(v_total, h_total),
            "isotropic_lamb_relative_rms": _relative_rms(v_iso, h_iso),
            "wss_relative_rms": _relative_rms(
                verified.wall_shear_stress_pa, historical.wall_shear_stress_pa
            ),
        }
    report = {
        "status": "REPRODUCED_AWAITING_HISTORICAL_V2_ORACLE",
        "claim_bearing": False,
        "source": mapping["source"],
        "published_input_registry_sha256": input_sha,
        "reproduction_commit_sha": reproduction_commit_sha,
        "profile": profile_name,
        "numerical_profile": profile.__dict__,
        "modes": list(REFERENCE_MODES),
        "arteries": rows,
        "historical_reference": {
            "status": "pending_cold_execution",
            "required_oracle": "output-stripped cold execution of frozen picoNewton_v2.ipynb",
            "notebook_blob_sha": mapping["source"]["published_notebook_blob_sha"],
        },
    }
    output.mkdir(parents=True, exist_ok=True)
    (output / "reproduction_report.json").write_text(
        json.dumps(report, indent=2, sort_keys=True), encoding="utf-8"
    )
    return report


def verify_reproduction_directory(
    reproduction_directory: str | Path,
    published_inputs_path: str | Path,
    source_registry: Mapping[str, Any],
) -> dict[str, Any]:
    root = Path(reproduction_directory)
    report_path = root / "reproduction_report.json"
    if not report_path.is_file():
        raise ProvenanceError("Reproduction directory is missing reproduction_report.json.")
    report = json.loads(report_path.read_text(encoding="utf-8"))
    input_path = Path(published_inputs_path)
    input_sha = sha256_file(input_path)
    mapping = yaml.safe_load(input_path.read_text(encoding="utf-8"))
    validate_published_source_binding(mapping, input_sha, source_registry)
    if (
        report.get("published_input_registry_sha256") != input_sha
        or report.get("source") != mapping["source"]
    ):
        raise ProvenanceError("Reproduction report source identity is stale.")
    count = 0
    max_residual = 0.0
    for artery_id in REFERENCE_ARTERY_IDS:
        for mode in REFERENCE_MODES:
            entry = report["arteries"][artery_id]["modes"][mode]
            record = load_artery_npz(root / entry["record_path"])
            if record.record_payload_sha256 != entry["record_payload_sha256"]:
                raise ProvenanceError("Record checksum differs from reproduction report.")
            published = record.metadata.get("published_source", {})
            expected = {
                "repository_commit_sha": mapping["source"]["repository_commit_sha"],
                "published_notebook_blob_sha": mapping["source"]["published_notebook_blob_sha"],
                "published_input_registry_sha256": input_sha,
                "reproduction_mode": mode,
            }
            if any(published.get(key) != value for key, value in expected.items()):
                raise ProvenanceError(f"{artery_id}/{mode} published-source metadata mismatch.")
            max_residual = max(max_residual, float(entry["max_backward_residual"]))
            count += 1
    status = "PASS_WORKFLOW_READY_HISTORICAL_ORACLE_PENDING"
    return {
        "status": status,
        "claim_bearing": False,
        "records_verified": count,
        "arteries_verified": 6,
        "modes_verified": list(REFERENCE_MODES),
        "max_backward_residual": max_residual,
        "published_input_registry_sha256": input_sha,
        "historical_reference_status": report.get("historical_reference", {}).get("status"),
    }


def save_reproduction_verification(report: Mapping[str, Any], path: str | Path) -> Path:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(dict(report), indent=2, sort_keys=True), encoding="utf-8")
    return target


__all__ = [
    "PublishedWomersleySolver",
    "classical_womersley",
    "reproduce_all_six",
    "save_reproduction_verification",
    "validate_published_inputs",
    "validate_published_source_binding",
    "verify_reproduction_directory",
]
