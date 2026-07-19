from __future__ import annotations

from typing import Any

import numpy as np

from .exceptions import ProvenanceError
from .models import ArteryRecord
from .validation import checksum_arrays


def record_payload_arrays(record: ArteryRecord) -> dict[str, np.ndarray]:
    arrays = {
        "radial_coordinate_m": np.asarray(record.radial_coordinate_m),
        "time_s": np.asarray(record.time_s),
        "lamb_density_signed_n_m3": np.asarray(record.lamb_density_signed_n_m3),
        "wall_shear_stress_pa": np.asarray(record.wall_shear_stress_pa),
    }
    if record.lamb_density_isotropic_n_m3 is not None:
        arrays["lamb_density_isotropic_n_m3"] = np.asarray(record.lamb_density_isotropic_n_m3)
    return arrays


def record_payload_metadata(record: ArteryRecord) -> dict[str, Any]:
    """Immutable metadata included in the converted-record checksum."""
    return {
        "artery_id": record.artery_id,
        "artery_name": record.artery_name,
        "radius_m": float(record.radius_m),
        "omega0_rad_s": float(record.omega0_rad_s),
        "source_identifier": record.source_identifier,
        "source_version": record.source_version,
        "source_archive_sha256": record.source_archive_sha256,
        "source_member_sha256": record.source_member_sha256,
        "conversion_manifest_sha256": record.conversion_manifest_sha256,
        "converter_commit_sha": record.converter_commit_sha,
        "coordinate_convention": record.coordinate_convention,
        "metadata": record.metadata,
    }


def compute_record_payload_sha256(record: ArteryRecord) -> str:
    return checksum_arrays(record_payload_arrays(record), record_payload_metadata(record))


def verify_record_payload(record: ArteryRecord) -> str:
    if record.record_payload_sha256 is None:
        raise ProvenanceError("Artery record is missing record_payload_sha256.")
    actual = compute_record_payload_sha256(record)
    if actual != record.record_payload_sha256:
        raise ProvenanceError(
            "Artery record payload checksum mismatch; arrays or immutable metadata were altered."
        )
    return actual
