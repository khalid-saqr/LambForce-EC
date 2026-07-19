"""Analytical lumped foundation benchmark.

This candidate is retained only as a zero-dimensional verification baseline.
It cannot supply spatial membrane fields and therefore cannot be the primary
claim-bearing solver.
"""
from __future__ import annotations

from .common import Geometry, MaterialState, complex_properties, tangential_response


def solve_uniform_normal(q_pa: complex, geometry: Geometry, material: MaterialState, omega_rad_s: float = 0.0) -> dict[str, complex]:
    _, kf = complex_properties(geometry, material, omega_rad_s)
    w = q_pa / kf
    force = q_pa * geometry.area_m2
    reaction = kf * w * geometry.area_m2
    return {
        'normal_displacement_m': w,
        'applied_force_n': force,
        'reaction_force_n': reaction,
        'foundation_stiffness_n_m3': kf,
    }


def solve_tangential(tau_w_pa: complex, geometry: Geometry, material: MaterialState, omega_rad_s: float = 0.0) -> dict[str, complex]:
    return tangential_response(tau_w_pa, geometry, material, omega_rad_s)
