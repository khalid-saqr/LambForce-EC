from dataclasses import replace
import numpy as np

from lambforce_ec.models import RunConfig
from lambforce_ec.synthetic import make_synthetic_artery
from lambforce_ec.workflow import run_case, run_registered_comparison


def test_registered_comparison_has_direct_increments():
    record = make_synthetic_artery(n_time=32)
    comparison = run_registered_comparison(record, config=RunConfig(nx=8, nz=8))
    for key in (
        "membrane_tension_increment_lamb_n_m",
        "membrane_tension_increment_isotropic_n_m",
        "membrane_tension_increment_anisotropy_n_m",
    ):
        assert comparison[key].shape == (32, 8, 8)
        assert np.all(np.isfinite(comparison[key]))


def test_wss_present_and_slip_limit():
    record = make_synthetic_artery(n_time=32)
    bonded = RunConfig(nx=8, nz=8, load_case="wss_only")
    slip = replace(bonded, membrane_cortex_coupling="tangential_slip_limit")
    bonded_result = run_case(record, config=bonded)
    slip_result = run_case(record, config=slip)
    assert np.linalg.norm(bonded_result.arrays["tension_xz_wss_n_m"]) > 0
    assert np.all(slip_result.arrays["tension_xz_wss_n_m"] == 0)
    assert np.linalg.norm(slip_result.arrays["displacement_tangential_m"]) > 0


def test_required_output_fields_and_summaries():
    record = make_synthetic_artery(n_time=32)
    result = run_case(record, config=RunConfig(nx=8, nz=8))
    required = {
        "membrane_tension_max_principal_n_m",
        "strain_max_principal",
        "displacement_normal_cell_m",
        "curvature_change_m_inv",
        "glycocalyx_strain_normal",
        "glycocalyx_reaction_stress_pa",
        "strain_energy_density_j_m3",
        "tension_loading_rate_n_m_s",
    }
    assert required <= set(result.arrays)
    summary = result.summaries["membrane_tension_max_principal_n_m"]
    for key in (
        "peak_signed", "peak_absolute", "rms", "cycle_mean", "percentile_95_absolute",
        "peak_to_peak", "time_of_peak_s", "index_of_peak", "spatial_concentration",
        "signed_asymmetry", "harmonic_magnitude", "harmonic_phase_rad", "harmonic_gain",
    ):
        assert key in summary
    assert "work_per_cycle_j" in result.metadata
