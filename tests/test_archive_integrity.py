from pathlib import Path
import copy

import numpy as np
import pytest
import yaml

from lambforce_ec.exceptions import ProvenanceError
from lambforce_ec.io import load_artery_npz
from lambforce_ec.protocol import assert_claim_bearing_source, load_yaml
from lambforce_ec.published_source import reproduce_all_six

ROOT = Path(__file__).resolve().parents[1]


def _reproduced_record(tmp_path):
    inputs = ROOT / "registry/published_v2_hydrodynamics.yaml"
    report = reproduce_all_six(inputs, tmp_path / "run", "2" * 40, "verification")
    path = tmp_path / "run" / report["arteries"]["femoral"]["modes"]["verified"]["record_path"]
    return load_artery_npz(path), report


def _verified_registry(record, report):
    source = load_yaml(ROOT / "registry/source_registry.yaml")
    registry = copy.deepcopy(source)
    item = registry["sources"][0]
    item["source_status"] = "verified_reproduction"
    item["qualification_status"] = "verified"
    item["reproduction_commit_sha"] = record.metadata["published_source"]["reproduction_commit_sha"]
    payload = {
        mode: {
            artery: report["arteries"][artery]["modes"][mode]["record_payload_sha256"]
            for artery in report["arteries"]
        }
        for mode in ("historical_v2", "verified")
    }
    item["record_payload_sha256_by_mode_and_artery"] = payload
    return registry


def test_claim_bearing_requires_verified_published_source(tmp_path):
    record, report = _reproduced_record(tmp_path)
    frozen = load_yaml(ROOT / "registry/source_registry.yaml")
    with pytest.raises(ProvenanceError):
        assert_claim_bearing_source(record, frozen)
    assert_claim_bearing_source(record, _verified_registry(record, report))


def test_tampered_record_payload_is_rejected(tmp_path):
    record, report = _reproduced_record(tmp_path)
    source_path = tmp_path / "run" / report["arteries"]["femoral"]["modes"]["verified"]["record_path"]
    with np.load(source_path, allow_pickle=False) as data:
        payload = {name: np.array(data[name], copy=True) for name in data.files}
    payload["lamb_density_signed_n_m3"][0, 0] += 1.0
    tampered = tmp_path / "tampered.npz"
    np.savez_compressed(tampered, **payload)
    with pytest.raises(ProvenanceError):
        load_artery_npz(tampered)


def test_source_identity_is_not_an_archive_placeholder():
    registry = yaml.safe_load((ROOT / "registry/source_registry.yaml").read_text())
    item = registry["sources"][0]
    assert item["source_identifier"] == "github:khalid-saqr/picoNewton#picoNewton_v2.ipynb"
    assert item["published_notebook_blob_sha"] == "9d61c237cda75df338ce0383038f7765c886f503"
    assert "archive_status" not in item
    assert "archive_sha256" not in item
