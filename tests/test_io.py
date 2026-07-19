import json
import numpy as np

from lambforce_ec.io import load_artery_npz, save_artery_npz, save_result
from lambforce_ec.models import RunConfig
from lambforce_ec.synthetic import make_synthetic_artery
from lambforce_ec.workflow import run_case


def test_npz_roundtrip_and_manifest(tmp_path):
    record = make_synthetic_artery(n_time=32)
    input_path = save_artery_npz(record, tmp_path / "input.npz")
    loaded = load_artery_npz(input_path)
    assert loaded.artery_id == record.artery_id
    assert np.allclose(loaded.lamb_density_signed_n_m3, record.lamb_density_signed_n_m3)
    result = run_case(loaded, config=RunConfig(nx=8, nz=8))
    output = save_result(result, tmp_path / "result")
    manifest = json.loads((output / "manifest.json").read_text())
    assert len(manifest["fields_sha256"]) == 64
    assert manifest["metadata"]["configuration_checksum"]
