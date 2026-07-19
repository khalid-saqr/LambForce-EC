from __future__ import annotations

from typing import Any, Mapping

import numpy as np
import scipy.sparse as sp
from scipy.special import jv

from .exceptions import ValidationError
from .published_contract import PublishedArteryCase, REFERENCE_MODES, ReproductionProfile


class PublishedWomersleySolver:
    """Independent Chebyshev implementation of the frozen harmonic system."""

    def __init__(self, radial_order: int, mode: str):
        if radial_order < 30 or mode not in REFERENCE_MODES:
            raise ValidationError("Invalid radial order or reproduction mode.")
        self.n = radial_order + 1
        self.mode = mode
        k = np.arange(self.n)
        x = np.cos(np.pi * k / radial_order)
        c = np.ones(self.n)
        c[0] = c[-1] = 2
        c *= (-1.0) ** k
        if mode == "verified":
            grid = np.tile(x, (self.n, 1)).T
            difference = grid - grid.T
            derivative = np.outer(c, 1 / c) / (difference + np.eye(self.n))
        else:
            grid = np.tile(x, (self.n, 1))
            difference = grid - grid.T + np.eye(self.n)
            derivative = np.outer(c, 1 / c) / difference
        derivative -= np.diag(np.sum(derivative, axis=1))
        self.r = (1 - x) / 2
        self.d = -2 * derivative
        safe_r = self.r.copy()
        safe_r[0] = 1 if mode == "verified" else 1e-12
        d1 = sp.csr_matrix(self.d)
        d2 = sp.csr_matrix(self.d @ self.d)
        self.l0 = d2 + sp.diags(1 / safe_r) @ d1
        self.l1 = self.l0 - sp.diags(1 / safe_r**2)

    def derivative_polynomial_error(self) -> float:
        return float(np.max(np.abs(self.d @ self.r**2 - 2 * self.r)))

    def solve(
        self, alpha: float, harmonic: int, forcing: float, beta: float, gamma: float, delta: float
    ):
        identity = sp.eye(self.n, format="csr")
        azz = ((1j * harmonic * alpha**2) * identity - self.l0).tolil()
        azt = (-beta * self.l1).tolil()
        atz = (-gamma * self.l0).tolil()
        att = ((1j * harmonic * alpha**2) * identity - delta * self.l1).tolil()
        bz = forcing * np.ones(self.n, complex)
        bt = np.zeros(self.n, complex)
        azz[0, :], azt[0, :], bz[0] = self.d[0], 0, 0
        atz[0, :], att[0, :], att[0, 0], bt[0] = 0, 0, 1, 0
        azz[-1, :], azz[-1, -1], azt[-1, :], bz[-1] = 0, 1, 0, 0
        atz[-1, :], att[-1, :], att[-1, -1], bt[-1] = 0, 0, 1, 0
        matrix = sp.vstack([sp.hstack([azz, azt]), sp.hstack([atz, att])]).toarray()
        rhs = np.concatenate([bz, bt])
        solution = np.linalg.solve(matrix, rhs)
        residual = matrix @ solution - rhs
        backward = np.linalg.norm(residual, np.inf) / max(
            np.linalg.norm(matrix, np.inf) * np.linalg.norm(solution, np.inf)
            + np.linalg.norm(rhs, np.inf),
            1e-30,
        )
        return solution[: self.n], solution[self.n :], float(backward)

    def vorticity(self, axial: np.ndarray, azimuthal: np.ndarray):
        omega_theta = -(self.d @ axial)
        derivative = self.d @ (self.r * azimuthal)
        omega_z = np.empty_like(azimuthal)
        omega_z[1:] = derivative[1:] / self.r[1:]
        omega_z[0] = 2 * (self.d @ azimuthal)[0]
        return omega_z, omega_theta


def classical_womersley(radial: np.ndarray, alpha: float) -> np.ndarray:
    kappa = alpha * np.sqrt(-1j)
    return 1 / (1j * alpha**2) * (1 - jv(0, kappa * radial) / jv(0, kappa))


def _interpolate(nodes: np.ndarray, field: np.ndarray, points: np.ndarray) -> np.ndarray:
    weights = (-1.0) ** np.arange(nodes.size)
    weights[0] *= 0.5
    weights[-1] *= 0.5
    output = np.empty((points.size, field.shape[1]), complex)
    for point_index, point in enumerate(points):
        difference = point - nodes
        exact = np.flatnonzero(difference == 0)
        if exact.size:
            output[point_index] = field[int(exact[0])]
            continue
        numerator = np.zeros(field.shape[1], complex)
        denominator = 0.0
        for node_index in range(nodes.size):
            factor = weights[node_index] / difference[node_index]
            numerator += factor * field[node_index]
            denominator += factor
        output[point_index] = numerator / denominator
    return output


def _reconstruct(coefficients: np.ndarray, basis: np.ndarray) -> np.ndarray:
    values = np.zeros((coefficients.shape[0], basis.shape[1]), complex)
    for index in range(coefficients.shape[1]):
        values += coefficients[:, index, None] * basis[index, None, :]
    return np.real(values)


def _reconstruct_scalar(coefficients: np.ndarray, basis: np.ndarray) -> np.ndarray:
    values = np.zeros(basis.shape[1], complex)
    for index in range(coefficients.size):
        values += coefficients[index] * basis[index]
    return np.real(values)


def compute_fields(
    case: PublishedArteryCase,
    mapping: Mapping[str, Any],
    mode: str,
    profile: ReproductionProfile,
    isotropic: bool,
) -> dict[str, Any]:
    fluid, control, anisotropy = mapping["fluid"], mapping["control_volume"], mapping["anisotropy"]
    density = float(fluid["density_kg_m3"])
    viscosity = float(fluid["kinematic_viscosity_m2_s"])
    frequency = float(fluid["fundamental_frequency_hz"])
    omega0 = 2 * np.pi * frequency
    dynamic_viscosity = density * viscosity
    alpha = case.radius_m * np.sqrt(omega0 / viscosity)
    velocity_scale = case.pressure_gradient_scale_pa_per_m * case.radius_m**2 / dynamic_viscosity
    beta = 0 if isotropic else float(anisotropy["beta"])
    gamma = 0 if isotropic else float(anisotropy["gamma"])
    delta = 1 if isotropic else float(anisotropy["delta"])
    solver = PublishedWomersleySolver(profile.radial_order, mode)
    fields = {name: [] for name in ("uz", "ut", "oz", "ot")}
    residuals = []
    for harmonic, coefficient in enumerate(case.harmonic_coefficients, 1):
        uz, ut, residual = solver.solve(alpha, harmonic, coefficient, beta, gamma, delta)
        oz, ot = solver.vorticity(uz, ut)
        for name, value in zip(fields, (uz, ut, oz, ot)):
            fields[name].append(value)
        residuals.append(residual)
    harmonic_fields = {name: np.stack(values, axis=1) for name, values in fields.items()}
    cycle = np.arange(profile.time_points) / profile.time_points
    basis = np.exp(1j * 2 * np.pi * np.outer(np.arange(1, 7), cycle))
    depth = float(control["reference_volume_m3"]) / float(control["reference_area_m2"])
    query = np.linspace(1 - depth / case.radius_m, 1, profile.quadrature_nodes)
    near = {name: _interpolate(solver.r, values, query) for name, values in harmonic_fields.items()}
    real = {name: _reconstruct(values, basis) for name, values in near.items()}
    if mode == "historical_v2":
        lamb = _reconstruct(near["ut"] * near["oz"] - near["uz"] * near["ot"], basis)
    else:
        lamb = real["ut"] * real["oz"] - real["uz"] * real["ot"]
    force_density = density * velocity_scale**2 * lamb / case.radius_m
    d_uz = solver.d @ harmonic_fields["uz"]
    d_ut = solver.d @ harmonic_fields["ut"]
    shear_h = (
        (d_uz[-1] + beta * (d_ut[-1] - harmonic_fields["ut"][-1]))
        * dynamic_viscosity
        * velocity_scale
        / case.radius_m
    )
    return {
        "radial_coordinate_m": query * case.radius_m,
        "time_s": cycle / frequency,
        "force_density_n_m3": force_density,
        "wall_shear_stress_pa": _reconstruct_scalar(shear_h, basis),
        "alpha": float(alpha),
        "omega0_rad_s": float(omega0),
        "max_backward_residual": max(residuals),
        "differentiation_polynomial_error": solver.derivative_polynomial_error(),
        "isotropic_classical_linf_error": float(
            np.max(np.abs(harmonic_fields["uz"][:, 0] - classical_womersley(solver.r, alpha)))
        )
        if isotropic
        else None,
        "fluid_integration_depth_m": depth,
    }
