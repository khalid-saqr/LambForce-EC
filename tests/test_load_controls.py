import numpy as np

from lambforce_ec.loads import extract_scalar_loads, select_load_case, spatial_kernel
from lambforce_ec.models import Geometry, LOAD_CASES, RunConfig
from lambforce_ec.synthetic import make_synthetic_artery
from lambforce_ec.workflow import required_run_matrix


def test_signed_and_exposure_integrals_are_distinct():
    record = make_synthetic_artery(n_time=32)
    radial_shape = np.linspace(-1.0, 1.0, record.radial_coordinate_m.size)[:, None]
    temporal = np.sin(record.omega0_rad_s * record.time_s)[None, :]
    record.lamb_density_signed_n_m3 = 8.0e3 * radial_shape * temporal + 1.0e3 * temporal
    loads = extract_scalar_loads(record)
    assert np.any(loads["lamb_signed_pa"] < 0)
    assert np.all(loads["lamb_exposure_pa"] >= np.abs(loads["lamb_signed_pa"]) - 1e-14)
    assert np.max(loads["lamb_exposure_pa"] - np.abs(loads["lamb_signed_pa"])) > 1e-4


def test_exposure_is_non_wss_diagnostic():
    loads = extract_scalar_loads(make_synthetic_artery(n_time=32))
    normal, tangential = select_load_case(loads, "exposure_diagnostic")
    assert np.all(normal >= 0)
    assert np.all(tangential == 0)


def test_total_equals_isotropic_plus_increment():
    loads = extract_scalar_loads(make_synthetic_artery(n_time=32))
    assert np.allclose(
        loads["lamb_signed_pa"],
        loads["lamb_isotropic_pa"] + loads["lamb_anisotropy_increment_pa"],
        rtol=1e-13,
        atol=1e-15,
    )


def test_required_load_cases_are_complete():
    rows = required_run_matrix()
    assert {row["load_case"] for row in rows} == LOAD_CASES
    assert len(rows) == 12


def test_sign_controls_and_zero_normal_control():
    loads = extract_scalar_loads(make_synthetic_artery(n_time=32))
    inward, tau_in = select_load_case(loads, "inward_only")
    outward, tau_out = select_load_case(loads, "outward_only")
    zero, tau_zero = select_load_case(loads, "zero_normal_load")
    assert np.all(inward <= 0) and np.all(outward >= 0)
    assert np.allclose(inward + outward, loads["lamb_signed_pa"])
    assert np.all(zero == 0)
    assert np.allclose(tau_in, loads["wss_pa"])
    assert np.allclose(tau_out, loads["wss_pa"])
    assert np.allclose(tau_zero, loads["wss_pa"])


def test_all_load_distributions_conserve_resultant():
    geometry = Geometry()
    thickness = np.linspace(0.2e-6, 0.8e-6, 12 * 10).reshape(12, 10)
    for distribution in ("uniform_apical", "localized_bound", "glycocalyx_resolved"):
        config = RunConfig(nx=12, nz=10, load_distribution=distribution)
        kernel = spatial_kernel(
            geometry, config, thickness if distribution == "glycocalyx_resolved" else None
        )
        dx = geometry.length_x_m / config.nx
        dz = geometry.length_z_m / config.nz
        assert abs(np.sum(kernel) * dx * dz / geometry.area_m2 - 1) < 1e-12
