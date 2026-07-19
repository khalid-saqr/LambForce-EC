"""Bounded finite-difference plate verification prototype.

The operator uses the square of a Dirichlet Laplacian, corresponding to a
Navier-type bounded plate. This is a boundary-condition verification model,
not the primary periodic monolayer representation.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve

from .common import Geometry, MaterialState, complex_properties


@dataclass
class FDSolution:
    x_m: np.ndarray
    z_m: np.ndarray
    load_pa: np.ndarray
    displacement_m: np.ndarray
    total_reaction_pa: np.ndarray
    residual_relative_l2: float
    applied_resultant_n: complex
    reaction_resultant_n: complex
    work_measure_j: float
    matrix_nnz: int


def _dirichlet_laplacian(nx: int, nz: int, dx: float, dz: float) -> sparse.csr_matrix:
    tx = sparse.diags(
        [np.ones(nx - 1), -2.0 * np.ones(nx), np.ones(nx - 1)],
        offsets=[-1, 0, 1],
        format='csr',
    ) / dx**2
    tz = sparse.diags(
        [np.ones(nz - 1), -2.0 * np.ones(nz), np.ones(nz - 1)],
        offsets=[-1, 0, 1],
        format='csr',
    ) / dz**2
    return sparse.kron(sparse.eye(nz, format='csr'), tx, format='csr') + sparse.kron(
        tz, sparse.eye(nx, format='csr'), format='csr'
    )


def solve(load_pa: np.ndarray, geometry: Geometry, material: MaterialState, omega_rad_s: float = 0.0) -> FDSolution:
    q = np.asarray(load_pa, dtype=complex)
    if q.ndim != 2 or min(q.shape) < 4:
        raise ValueError('load_pa must be a 2D interior-grid array.')
    nx, nz = q.shape
    dx = geometry.length_x_m / (nx + 1)
    dz = geometry.length_z_m / (nz + 1)
    x = (np.arange(nx) + 1) * dx
    z = (np.arange(nz) + 1) * dz
    d, kf = complex_properties(geometry, material, omega_rad_s)

    lap = _dirichlet_laplacian(nx, nz, dx, dz)
    bih = lap @ lap
    identity = sparse.eye(nx * nz, dtype=complex, format='csr')
    a = (d * bih + kf * identity).tocsr()
    qvec = q.reshape(nx * nz, order='F')
    wvec = spsolve(a, qvec)
    reaction = a @ wvec
    residual = reaction - qvec
    residual_norm = np.linalg.norm(residual) / max(np.linalg.norm(qvec), 1e-30)
    darea = dx * dz
    f_applied = np.sum(qvec) * darea
    f_reaction = np.sum(reaction) * darea
    work = 0.5 * float(np.real(np.vdot(wvec, reaction) * darea))

    return FDSolution(
        x_m=x,
        z_m=z,
        load_pa=q,
        displacement_m=wvec.reshape((nx, nz), order='F'),
        total_reaction_pa=reaction.reshape((nx, nz), order='F'),
        residual_relative_l2=float(residual_norm),
        applied_resultant_n=complex(f_applied),
        reaction_resultant_n=complex(f_reaction),
        work_measure_j=work,
        matrix_nnz=int(a.nnz),
    )


def uniform_load(q_pa: complex, nx: int, nz: int) -> np.ndarray:
    return np.full((nx, nz), q_pa, dtype=complex)


def cosine_load(q_mean_pa: complex, q_amplitude_pa: complex, geometry: Geometry, nx: int, nz: int, mode_x: int = 1, mode_z: int = 1) -> np.ndarray:
    x = (np.arange(nx) + 1) * geometry.length_x_m / (nx + 1)
    z = (np.arange(nz) + 1) * geometry.length_z_m / (nz + 1)
    shape = np.cos(2.0 * np.pi * mode_x * x[:, None] / geometry.length_x_m) * np.cos(
        2.0 * np.pi * mode_z * z[None, :] / geometry.length_z_m
    )
    return q_mean_pa + q_amplitude_pa * shape
