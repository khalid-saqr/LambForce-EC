import numpy as np

from lambforce_ec.loads import integrate_radial_density, spatial_kernel
from lambforce_ec.models import Geometry, RunConfig


def test_signed_radial_integration_and_exposure_scale():
    r = np.linspace(0.0, 10e-6, 101)
    density = np.full((r.size, 7), 2.5e3)
    q = integrate_radial_density(r, density)
    assert np.allclose(q, 2.5e-2, rtol=1e-13)


def test_all_spatial_kernels_conserve_area_resultant():
    geometry = Geometry()
    for distribution in ("uniform_apical", "localized_bound"):
        config = RunConfig(nx=24, nz=20, load_distribution=distribution)
        kernel = spatial_kernel(geometry, config)
        dx = geometry.length_x_m / config.nx
        dz = geometry.length_z_m / config.nz
        assert abs(np.sum(kernel) * dx * dz / geometry.area_m2 - 1) < 1e-12
    thickness = np.linspace(0.2e-6, 0.8e-6, 24 * 20).reshape(24, 20)
    config = RunConfig(nx=24, nz=20, load_distribution="glycocalyx_resolved")
    kernel = spatial_kernel(geometry, config, thickness)
    dx = geometry.length_x_m / config.nx
    dz = geometry.length_z_m / config.nz
    assert abs(np.sum(kernel) * dx * dz / geometry.area_m2 - 1) < 1e-12
