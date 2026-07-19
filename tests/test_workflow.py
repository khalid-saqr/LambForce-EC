import numpy as np

from lambforce_ec.models import RunConfig
from lambforce_ec.synthetic import make_synthetic_artery
from lambforce_ec.workflow import run_case, run_registered_comparison


def test_production_workflow_outputs_and_comparison():
    record = make_synthetic_artery(n_radial=24, n_time=64)
    config = RunConfig(nx=16, nz=16, load_distribution="localized_bound")
    result = run_case(record, config=config)
    assert result.arrays["membrane_tension_max_principal_n_m"].shape == (64, 16, 16)
    assert result.metadata["source_checksum"] == record.source_checksum
    assert result.metadata["normal_residual_relative_l2_max"] < 1e-10
    assert result.metadata["normal_resultant_conservation_relative_max"] < 1e-10
    assert np.min(result.arrays["strain_energy_density_j_m3"]) >= -1e-20
    comparison = run_registered_comparison(record, config=config)
    delta = comparison["membrane_tension_increment_lamb_n_m"]
    assert np.linalg.norm(delta) > 0
    assert np.isfinite(comparison["incremental_lamb_wss_ratio"])


def test_generic_nonreference_artery_is_accepted():
    record = make_synthetic_artery(artery_id="new_external_artery_2031", n_time=32)
    result = run_case(record, config=RunConfig(nx=8, nz=8))
    assert result.metadata["artery_id"] == "new_external_artery_2031"
