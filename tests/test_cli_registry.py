import yaml
import pytest

from lambforce_ec.cli import main
from lambforce_ec.exceptions import ValidationError
from lambforce_ec.io import config_from_mapping
from lambforce_ec.registry import create_solver, solver_registry


def test_solver_registry_and_invalid_solver():
    registry = solver_registry()
    assert set(registry) == {"periodic_spectral_2d", "bounded_fd_2d", "lumped_0d"}
    assert create_solver("periodic_spectral_2d").solver_id == "periodic_spectral_2d"
    with pytest.raises(ValidationError):
        create_solver("missing")


def test_cli_end_to_end(tmp_path):
    input_path = tmp_path / "input.npz"
    assert main(["synthetic", str(input_path)]) == 0
    assert main(["validate", str(input_path)]) == 0
    config_path = tmp_path / "config.yaml"
    config = {
        "run": {"nx": 8, "nz": 8, "load_distribution": "localized_bound"}
    }
    config_path.write_text(yaml.safe_dump(config), encoding="utf-8")
    geometry, material, run = config_from_mapping(config)
    assert geometry.area_m2 > 0 and material.parameter_set_id and run.nx == 8
    output = tmp_path / "result"
    assert main([
        "run", "--input", str(input_path), "--config", str(config_path), "--output", str(output)
    ]) == 0
    assert (output / "manifest.json").exists()
