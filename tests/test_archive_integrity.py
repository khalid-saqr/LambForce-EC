import numpy as np
import pytest

from lambforce_ec.archive import qualify_hydrodynamics
from lambforce_ec.exceptions import ProvenanceError
from lambforce_ec.io import load_artery_npz, save_artery_npz
from lambforce_ec.loads import integrate_radial_density
from lambforce_ec.models import ArteryRecord
from lambforce_ec.protocol import assert_claim_bearing_source
from lambforce_ec.provenance import compute_record_payload_sha256
from lambforce_ec.synthetic import make_synthetic_artery


def _harmonics(values: np.ndarray) -> dict[str, list[float]]:
    coefficients = np.fft.rfft(np.asarray(values, dtype=float)) / values.size
    return {
        "magnitude": np.abs(coefficients).tolist(),
        "phase_rad": np.angle(coefficients).tolist(),
    }


def _claim_record() -> tuple[ArteryRecord, dict]:
    base = make_synthetic_artery(
        artery_id="femoral",
        radius_m=4.0e-3,
        n_radial=16,
        n_time=32,
    )
    total = integrate_radial_density(
        base.radial_coordinate_m, base.lamb_density_signed_n_m3
    )
    isotropic = integrate_radial_density(
        base.radial_coordinate_m, base.lamb_density_isotropic_n_m3
    )
    contract = {
        "rho_kg_m3": 1060.0,
        "nu_zz_m2_s": 3.5e-6,
        "reference_area_m2": 1.0e-10,
        "fluid_integration_depth_m": float(
            base.radial_coordinate_m[-1] - base.radial_coordinate_m[0]
        ),
        "womersley_alpha": 5.87,
        "radial_collocation_order": base.radial_coordinate_m.size,
        "harmonic_rms_tolerance": 1.0e-12,
        "source_waveform_identifier": "test:femoral-waveform",
        "source_archive_member": "femoral.npz",
        "archive_harmonics": {
            "signed_lamb_load": _harmonics(total),
            "isotropic_lamb_load": _harmonics(isotropic),
            "wall_shear_stress": _harmonics(base.wall_shear_stress_pa),
        },
    }
    record = ArteryRecord(
        artery_id="femoral",
        artery_name="Femoral",
        radius_m=base.radius_m,
        omega0_rad_s=base.omega0_rad_s,
        radial_coordinate_m=base.radial_coordinate_m.copy(),
        time_s=base.time_s.copy(),
        lamb_density_signed_n_m3=base.lamb_density_signed_n_m3.copy(),
        wall_shear_stress_pa=base.wall_shear_stress_pa.copy(),
        lamb_density_isotropic_n_m3=base.lamb_density_isotropic_n_m3.copy(),
        source_identifier="doi:10.1038/s41598-026-47474-x#immutable-hydrodynamic-archive",
        source_version="test-archive-v1",
        source_checksum="a" * 64,
        source_member_sha256="b" * 64,
        conversion_manifest_sha256="c" * 64,
        converter_commit_sha="d" * 40,
        metadata={"hydrodynamic_contract": contract},
    )
    record.record_payload_sha256 = compute_record_payload_sha256(record)
    registry = {
        "registry_version": "test",
        "sources": [
            {
                "source_identifier": record.source_identifier,
                "source_version": record.source_version,
                "artery_ids": [
                    "aortic_root",
                    "thoracic_aorta",
                    "femoral",
                    "carotid",
                    "iliac",
                    "brachial",
                ],
                "archive_status": "verified",
                "archive_sha256": record.source_archive_sha256,
                "source_member_sha256_by_artery": {
                    artery: (record.source_member_sha256 if artery == "femoral" else "e" * 64)
                    for artery in (
                        "aortic_root",
                        "thoracic_aorta",
                        "femoral",
                        "carotid",
                        "iliac",
                        "brachial",
                    )
                },
                "conversion_manifest_sha256_by_artery": {
                    artery: (
                        record.conversion_manifest_sha256 if artery == "femoral" else "f" * 64
                    )
                    for artery in (
                        "aortic_root",
                        "thoracic_aorta",
                        "femoral",
                        "carotid",
                        "iliac",
                        "brachial",
                    )
                },
                "converter_commit_sha": record.converter_commit_sha,
                "role": "immutable_hydrodynamic_ground_truth",
                "claim_bearing": True,
            }
        ],
    }
    return record, registry


def test_claim_bearing_requires_full_archive_binding():
    record, registry = _claim_record()
    assert_claim_bearing_source(record, registry)
    report = qualify_hydrodynamics(record, registry)
    assert report["status"] == "PASS"
    assert max(report["regression_errors"].values()) < 1e-12


def test_tampered_record_payload_is_rejected(tmp_path):
    record, _ = _claim_record()
    path = save_artery_npz(record, tmp_path / "record.npz")
    with np.load(path, allow_pickle=False) as data:
        payload = {name: np.array(data[name], copy=True) for name in data.files}
    payload["lamb_density_signed_n_m3"][0, 0] += 1.0
    tampered = tmp_path / "tampered.npz"
    np.savez_compressed(tampered, **payload)
    with pytest.raises(ProvenanceError):
        load_artery_npz(tampered)


def test_missing_isotropic_field_is_never_replaced_by_zero():
    record, registry = _claim_record()
    record.lamb_density_isotropic_n_m3 = None
    record.record_payload_sha256 = compute_record_payload_sha256(record)
    with pytest.raises(ProvenanceError):
        assert_claim_bearing_source(record, registry)


def test_metadata_archive_checksum_alone_cannot_authorize_tampered_arrays():
    record, registry = _claim_record()
    record.lamb_density_signed_n_m3[0, 0] += 1.0
    with pytest.raises(ProvenanceError):
        assert_claim_bearing_source(record, registry)
