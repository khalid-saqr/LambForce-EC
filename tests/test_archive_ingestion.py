import numpy as np
import yaml

from lambforce_ec.archive import ingest_archive_member
from lambforce_ec.io import load_artery_npz
from lambforce_ec.loads import integrate_radial_density
from lambforce_ec.synthetic import make_synthetic_artery


def _harmonics(values):
    coeff = np.fft.rfft(np.asarray(values, dtype=float)) / len(values)
    return {
        "magnitude": np.abs(coeff).tolist(),
        "phase_rad": np.angle(coeff).tolist(),
    }


def test_archive_ingestion_is_deterministic_and_checksums_every_layer(tmp_path):
    base = make_synthetic_artery(
        artery_id="femoral", radius_m=4.0e-3, n_radial=16, n_time=32
    )
    archive = tmp_path / "immutable-archive.bin"
    archive.write_bytes(b"immutable archive fixture")
    member = tmp_path / "femoral-member.npz"
    np.savez_compressed(
        member,
        radial_coordinate_m=base.radial_coordinate_m,
        time_s=base.time_s,
        lamb_density_signed_n_m3=base.lamb_density_signed_n_m3,
        lamb_density_isotropic_n_m3=base.lamb_density_isotropic_n_m3,
        wall_shear_stress_pa=base.wall_shear_stress_pa,
    )
    total = integrate_radial_density(
        base.radial_coordinate_m, base.lamb_density_signed_n_m3
    )
    isotropic = integrate_radial_density(
        base.radial_coordinate_m, base.lamb_density_isotropic_n_m3
    )
    manifest = {
        "artery_id": "femoral",
        "artery_name": "Femoral",
        "radius_m": base.radius_m,
        "omega0_rad_s": base.omega0_rad_s,
        "source_identifier": "doi:10.1038/s41598-026-47474-x#immutable-hydrodynamic-archive",
        "source_version": "fixture-v1",
        "coordinate_convention": "outward_normal_positive",
        "hydrodynamic_contract": {
            "rho_kg_m3": 1060.0,
            "nu_zz_m2_s": 3.5e-6,
            "reference_area_m2": 1.0e-10,
            "fluid_integration_depth_m": float(
                base.radial_coordinate_m[-1] - base.radial_coordinate_m[0]
            ),
            "womersley_alpha": 5.87,
            "radial_collocation_order": 16,
            "harmonic_rms_tolerance": 1.0e-12,
            "source_waveform_identifier": "fixture:femoral",
            "source_archive_member": member.name,
            "archive_harmonics": {
                "signed_lamb_load": _harmonics(total),
                "isotropic_lamb_load": _harmonics(isotropic),
                "wall_shear_stress": _harmonics(base.wall_shear_stress_pa),
            },
        },
    }
    manifest_path = tmp_path / "manifest.yaml"
    manifest_path.write_text(yaml.safe_dump(manifest), encoding="utf-8")

    first = ingest_archive_member(
        archive, member, manifest_path, tmp_path / "first.npz", "1" * 40
    )
    second = ingest_archive_member(
        archive, member, manifest_path, tmp_path / "second.npz", "1" * 40
    )
    assert first["record_payload_sha256"] == second["record_payload_sha256"]
    assert first["source_archive_sha256"] == second["source_archive_sha256"]
    assert first["source_member_sha256"] == second["source_member_sha256"]
    loaded = load_artery_npz(tmp_path / "first.npz")
    assert loaded.record_payload_sha256 == first["record_payload_sha256"]
    assert loaded.source_member_sha256 == first["source_member_sha256"]
