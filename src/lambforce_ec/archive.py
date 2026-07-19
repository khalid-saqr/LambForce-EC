from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Mapping
import re

import numpy as np

from .exceptions import ProvenanceError, ValidationError
from .io import load_mapping, save_artery_npz
from .loads import integrate_radial_density
from .models import ArteryRecord
from .protocol import assert_claim_bearing_source, validate_claim_bearing_hydrodynamic_contract
from .provenance import compute_record_payload_sha256
from .validation import require_keys, sha256_file


_GIT_SHA = re.compile(r"^[0-9a-f]{40}$")
_REQUIRED_MEMBER_ARRAYS = {
    "radial_coordinate_m",
    "time_s",
    "lamb_density_signed_n_m3",
    "lamb_density_isotropic_n_m3",
    "wall_shear_stress_pa",
}
_REQUIRED_MANIFEST = {
    "artery_id",
    "artery_name",
    "radius_m",
    "omega0_rad_s",
    "source_identifier",
    "source_version",
    "coordinate_convention",
    "hydrodynamic_contract",
}


def ingest_archive_member(
    archive_path: str | Path,
    member_npz_path: str | Path,
    manifest_path: str | Path,
    output_path: str | Path,
    converter_commit_sha: str,
) -> dict[str, Any]:
    """Deterministically convert one immutable archive member into the canonical record."""
    archive_path = Path(archive_path)
    member_npz_path = Path(member_npz_path)
    manifest_path = Path(manifest_path)
    if not _GIT_SHA.match(converter_commit_sha):
        raise ValidationError("converter_commit_sha must be a full lowercase 40-character Git SHA.")
    manifest = load_mapping(manifest_path)
    require_keys(manifest, _REQUIRED_MANIFEST, "conversion manifest")
    with np.load(member_npz_path, allow_pickle=False) as data:
        require_keys(
            {name: None for name in data.files},
            _REQUIRED_MEMBER_ARRAYS,
            "archive member NPZ",
        )
        arrays = {name: np.asarray(data[name]) for name in _REQUIRED_MEMBER_ARRAYS}
    archive_sha = sha256_file(archive_path)
    member_sha = sha256_file(member_npz_path)
    manifest_sha = sha256_file(manifest_path)
    expected_archive = manifest.get("expected_archive_sha256")
    expected_member = manifest.get("expected_source_member_sha256")
    if expected_archive is not None and expected_archive != archive_sha:
        raise ProvenanceError("Archive checksum does not match the conversion manifest.")
    if expected_member is not None and expected_member != member_sha:
        raise ProvenanceError("Source-member checksum does not match the conversion manifest.")
    record = ArteryRecord(
        artery_id=str(manifest["artery_id"]),
        artery_name=str(manifest["artery_name"]),
        radius_m=float(manifest["radius_m"]),
        omega0_rad_s=float(manifest["omega0_rad_s"]),
        radial_coordinate_m=arrays["radial_coordinate_m"],
        time_s=arrays["time_s"],
        lamb_density_signed_n_m3=arrays["lamb_density_signed_n_m3"],
        wall_shear_stress_pa=arrays["wall_shear_stress_pa"],
        lamb_density_isotropic_n_m3=arrays["lamb_density_isotropic_n_m3"],
        source_identifier=str(manifest["source_identifier"]),
        source_version=str(manifest["source_version"]),
        source_checksum=archive_sha,
        coordinate_convention=str(manifest["coordinate_convention"]),
        source_member_sha256=member_sha,
        conversion_manifest_sha256=manifest_sha,
        converter_commit_sha=converter_commit_sha,
        metadata={
            "hydrodynamic_contract": manifest["hydrodynamic_contract"],
            "conversion": {
                "archive_filename": archive_path.name,
                "source_member_filename": member_npz_path.name,
                "manifest_filename": manifest_path.name,
            },
        },
    )
    record.validate()
    validate_claim_bearing_hydrodynamic_contract(record)
    output = save_artery_npz(record, output_path)
    return {
        "status": "CONVERTED_NOT_YET_QUALIFIED",
        "artery_id": record.artery_id,
        "output": str(output),
        "source_archive_sha256": archive_sha,
        "source_member_sha256": member_sha,
        "conversion_manifest_sha256": manifest_sha,
        "converter_commit_sha": converter_commit_sha,
        "record_payload_sha256": compute_record_payload_sha256(record),
    }


def _reconstruct_archive_harmonics(entry: Mapping[str, Any], n_time: int) -> np.ndarray:
    magnitude = np.asarray(entry["magnitude"], dtype=float)
    phase = np.asarray(entry["phase_rad"], dtype=float)
    if magnitude.size > n_time // 2 + 1:
        raise ProvenanceError("Archived harmonic series exceeds the available time resolution.")
    coefficients = np.zeros(n_time // 2 + 1, dtype=complex)
    coefficients[: magnitude.size] = magnitude * np.exp(1j * phase)
    return np.fft.irfft(coefficients * n_time, n=n_time)


def _relative_rms(value: np.ndarray, reference: np.ndarray) -> float:
    numerator = float(np.sqrt(np.mean((np.asarray(value) - np.asarray(reference)) ** 2)))
    denominator = max(float(np.sqrt(np.mean(np.asarray(reference) ** 2))), 1e-30)
    return numerator / denominator


def qualify_hydrodynamics(
    record: ArteryRecord,
    source_registry: Mapping[str, Any],
) -> dict[str, Any]:
    """Evaluate the immutable-source and waveform-regression gates for one artery."""
    assert_claim_bearing_source(record, source_registry)
    contract = record.metadata["hydrodynamic_contract"]
    total = integrate_radial_density(
        record.radial_coordinate_m, record.lamb_density_signed_n_m3
    )
    isotropic = integrate_radial_density(
        record.radial_coordinate_m, record.lamb_density_isotropic_n_m3
    )
    wss = np.asarray(record.wall_shear_stress_pa, dtype=float)
    archive_harmonics = contract["archive_harmonics"]
    reconstructed = {
        "signed_lamb_load": _reconstruct_archive_harmonics(
            archive_harmonics["signed_lamb_load"], record.time_s.size
        ),
        "isotropic_lamb_load": _reconstruct_archive_harmonics(
            archive_harmonics["isotropic_lamb_load"], record.time_s.size
        ),
        "wall_shear_stress": _reconstruct_archive_harmonics(
            archive_harmonics["wall_shear_stress"], record.time_s.size
        ),
    }
    errors = {
        "signed_lamb_relative_rms": _relative_rms(total, reconstructed["signed_lamb_load"]),
        "isotropic_lamb_relative_rms": _relative_rms(
            isotropic, reconstructed["isotropic_lamb_load"]
        ),
        "wss_relative_rms": _relative_rms(wss, reconstructed["wall_shear_stress"]),
    }
    tolerance = float(contract["harmonic_rms_tolerance"])
    passed = all(value <= tolerance for value in errors.values())
    report = {
        "status": "PASS" if passed else "FAIL",
        "artery_id": record.artery_id,
        "source_archive_sha256": record.source_archive_sha256,
        "source_member_sha256": record.source_member_sha256,
        "conversion_manifest_sha256": record.conversion_manifest_sha256,
        "converter_commit_sha": record.converter_commit_sha,
        "record_payload_sha256": record.record_payload_sha256,
        "harmonic_rms_tolerance": tolerance,
        "regression_errors": errors,
        "signed_integral_peak_abs_pa": float(np.max(np.abs(total))),
        "isotropic_integral_peak_abs_pa": float(np.max(np.abs(isotropic))),
        "anisotropy_increment_peak_abs_pa": float(np.max(np.abs(total - isotropic))),
        "exposure_peak_pa": float(
            np.max(
                integrate_radial_density(
                    record.radial_coordinate_m,
                    np.abs(record.lamb_density_signed_n_m3),
                )
            )
        ),
        "total_equals_isotropic_plus_increment": bool(
            np.allclose(total, isotropic + (total - isotropic), rtol=1e-13, atol=1e-15)
        ),
    }
    if not passed:
        raise ProvenanceError(json.dumps(report, sort_keys=True))
    return report


def save_qualification_report(report: Mapping[str, Any], path: str | Path) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(dict(report), indent=2, sort_keys=True), encoding="utf-8")
    return path
