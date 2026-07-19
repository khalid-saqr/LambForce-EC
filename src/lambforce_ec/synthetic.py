from __future__ import annotations

import numpy as np

from .models import ArteryRecord
from .validation import canonical_json_bytes, checksum_arrays, sha256_bytes


def make_synthetic_artery(
    artery_id: str = "synthetic_validation_artery",
    radius_m: float = 4.0e-3,
    frequency_hz: float = 1.2,
    n_radial: int = 48,
    n_time: int = 128,
) -> ArteryRecord:
    """Generate a non-claim-bearing signed field for software verification."""
    omega0 = 2 * np.pi * frequency_hz
    period = 1 / frequency_hz
    time = np.arange(n_time) * period / n_time
    depth = 10e-6
    radial = np.linspace(radius_m - depth, radius_m, n_radial)
    y = (radius_m - radial) / depth
    envelope = np.exp(-4 * y)[:, None]
    phase = omega0 * time[None, :]
    total = 5.0e3 * envelope * (
        np.sin(phase) + 0.35 * np.sin(2 * phase + 0.4) + 0.15 * np.sin(4 * phase - 0.2)
    )
    isotropic = 0.30 * 5.0e3 * envelope * np.sin(phase - 0.1)
    wss = 2.0 * np.sin(phase.ravel() - 0.25) + 0.4 * np.sin(2 * phase.ravel() + 0.5)
    metadata = {
        "generator": "lambforce_ec.synthetic.make_synthetic_artery",
        "purpose": "software verification only; not a physiological result",
        "frequency_hz": frequency_hz,
        "archive_status": "synthetic_non_claim_bearing",
    }
    arrays = {
        "radial_coordinate_m": radial,
        "time_s": time,
        "lamb_density_signed_n_m3": total,
        "wall_shear_stress_pa": wss,
        "lamb_density_isotropic_n_m3": isotropic,
    }
    archive_checksum = checksum_arrays(arrays, {**metadata, "kind": "synthetic_archive"})
    member_checksum = checksum_arrays(arrays, {**metadata, "kind": "synthetic_member"})
    manifest_checksum = sha256_bytes(canonical_json_bytes({**metadata, "kind": "manifest"}))
    record = ArteryRecord(
        artery_id=artery_id,
        artery_name="Synthetic validation artery",
        radius_m=radius_m,
        omega0_rad_s=omega0,
        radial_coordinate_m=radial,
        time_s=time,
        lamb_density_signed_n_m3=total,
        wall_shear_stress_pa=wss,
        lamb_density_isotropic_n_m3=isotropic,
        source_identifier="synthetic://lambforce-ec/software-validation",
        source_version="2.1.0-readiness",
        source_checksum=archive_checksum,
        source_member_sha256=member_checksum,
        conversion_manifest_sha256=manifest_checksum,
        converter_commit_sha="0" * 40,
        metadata=metadata,
    )
    record.validate()
    return record
