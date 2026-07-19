import numpy as np

from lambforce_ec.models import Geometry, MaterialState, RunConfig
from lambforce_ec.solvers import BoundedFDPlateSolver, SpectralPlateSolver
from lambforce_ec.synthetic import make_synthetic_artery
from lambforce_ec.workflow import run_case


def test_zero_load_and_positive_energy_limits():
    geometry = Geometry()
    material = MaterialState().for_rheology("elastic")
    zero = np.zeros((8, 8))
    spectral = SpectralPlateSolver().solve(zero, geometry, material)
    bounded = BoundedFDPlateSolver().solve(zero, geometry, material)
    assert np.allclose(spectral.displacement_cell_m, 0)
    assert np.allclose(bounded.displacement_cell_m, 0)
    q = np.full((8, 8), 0.02)
    assert SpectralPlateSolver().solve(q, geometry, material).work_measure_j > 0
    assert BoundedFDPlateSolver().solve(q, geometry, material).work_measure_j > 0


def test_elastic_primary_has_zero_dissipation():
    record = make_synthetic_artery(n_time=32)
    result = run_case(record, config=RunConfig(nx=8, nz=8, rheology_mode="elastic"))
    assert all(value == 0 for value in result.metadata["average_dissipated_power_w_by_harmonic"])
