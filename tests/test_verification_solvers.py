import numpy as np

from lambforce_ec.models import Geometry, MaterialState
from lambforce_ec.solvers import BoundedFDPlateSolver, LumpedFoundationSolver


def test_lumped_and_bounded_solver_residuals():
    geometry = Geometry()
    material = MaterialState()
    lumped = LumpedFoundationSolver().solve(0.02, geometry, material)
    assert abs(lumped.reaction_resultant_n / lumped.applied_resultant_n - 1) < 1e-12
    q = np.full((10, 12), 0.02)
    bounded = BoundedFDPlateSolver().solve(q, geometry, material)
    assert bounded.residual_relative_l2 < 1e-10
