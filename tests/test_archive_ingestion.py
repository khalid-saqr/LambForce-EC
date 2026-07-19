from pathlib import Path

from lambforce_ec.published_source import reproduce_all_six, verify_reproduction_directory
from lambforce_ec.protocol import load_yaml

ROOT = Path(__file__).resolve().parents[1]


def test_all_six_published_source_reproduction_is_deterministic(tmp_path):
    inputs = ROOT / "registry/published_v2_hydrodynamics.yaml"
    first = reproduce_all_six(inputs, tmp_path / "first", "1" * 40, "verification")
    second = reproduce_all_six(inputs, tmp_path / "second", "1" * 40, "verification")
    assert first["status"] == "REPRODUCED_AWAITING_HISTORICAL_V2_ORACLE"
    assert set(first["arteries"]) == {
        "aortic_root", "thoracic_aorta", "femoral", "carotid", "iliac", "brachial"
    }
    hashes_first = {
        (artery, mode): row["record_payload_sha256"]
        for artery, data in first["arteries"].items()
        for mode, row in data["modes"].items()
    }
    hashes_second = {
        (artery, mode): row["record_payload_sha256"]
        for artery, data in second["arteries"].items()
        for mode, row in data["modes"].items()
    }
    assert len(hashes_first) == 12
    assert hashes_first == hashes_second
    registry = load_yaml(ROOT / "registry/source_registry.yaml")
    verification = verify_reproduction_directory(tmp_path / "first", inputs, registry)
    assert verification["records_verified"] == 12
    assert verification["status"] == "PASS_WORKFLOW_READY_HISTORICAL_ORACLE_PENDING"
    assert verification["claim_bearing"] is False
