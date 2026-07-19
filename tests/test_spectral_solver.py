import numpy as np

from lambforce_ec.constitutive import plate_foundation_properties
from lambforce_ec.models import Geometry, MaterialState
from lambforce_ec.solvers import SpectralPlateSolver


def test_uniform_static_matches_foundation_limit_and_conserves_force():
    geometry = Geometry()
    material = MaterialState()
    q = np.full((32, 32), 0.02)
    solution = SpectralPlateSolver().solve(q, geometry, material)
    _, kf, kg = plate_foundation_properties(geometry, material, 0.0)
    assert np.allclose(solution.displacement_cell_m, 0.02 / kf, rtol=1e-12, atol=1e-30)
    assert np.allclose(solution.displacement_apical_top_m, 0.02 / kf + 0.02 / kg, rtol=1e-12)
    assert solution.residual_relative_l2 < 1e-12
    assert abs(solution.reaction_resultant_n / solution.applied_resultant_n - 1) < 1e-12
    assert solution.work_measure_j > 0


def test_single_fourier_mode_matches_analytic_operator():
    geometry = Geometry()
    material = MaterialState()
    nx = nz = 32
    x = np.arange(nx) * geometry.length_x_m / nx
    z = np.arange(nz) * geometry.length_z_m / nz
    shape = np.cos(2 * np.pi * x[:, None] / geometry.length_x_m) * np.cos(
        4 * np.pi * z[None, :] / geometry.length_z_m
    )
    q = 0.03 * shape
    solution = SpectralPlateSolver().solve(q, geometry, material)
    d, kf, _ = plate_foundation_properties(geometry, material, 0.0)
    k2 = (2 * np.pi / geometry.length_x_m) ** 2 + (4 * np.pi / geometry.length_z_m) ** 2
    expected = q / (d * k2**2 + kf)
    assert np.linalg.norm(solution.displacement_cell_m.real - expected) / np.linalg.norm(expected) < 1e-11
