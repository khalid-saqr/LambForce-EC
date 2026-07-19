"""Periodic Fourier-spectral plate-on-foundation prototype.

Equation:
    D nabla^4 w + k_f w = q(x,z)

The periodic domain represents a repeating endothelial monolayer unit cell.
The method supports real static and complex harmonic amplitudes.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np

from .common import Geometry, MaterialState, complex_properties


@dataclass
class SpectralSolution:
    x_m: np.ndarray
    z_m: np.ndarray
    load_pa: np.ndarray
    displacement_m: np.ndarray
    foundation_reaction_pa: np.ndarray
    bending_reaction_pa: np.ndarray
    total_reaction_pa: np.ndarray
    residual_relative_l2: float
    applied_resultant_n: complex
    reaction_resultant_n: complex
    work_measure_j: float
    average_dissipated_power_w: float


def solve(load_pa: np.ndarray, geometry: Geometry, material: MaterialState, omega_rad_s: float = 0.0) -> SpectralSolution:
    q = np.asarray(load_pa, dtype=complex)
    if q.ndim != 2 or min(q.shape) < 4:
        raise ValueError('load_pa must be a 2D array with at least 4 points per axis.')
    nx, nz = q.shape
    dx = geometry.length_x_m / nx
    dz = geometry.length_z_m / nz
    x = np.arange(nx) * dx
    z = np.arange(nz) * dz

    d, kf = complex_properties(geometry, material, omega_rad_s)
    kx = 2.0 * np.pi * np.fft.fftfreq(nx, d=dx)
    kz = 2.0 * np.pi * np.fft.fftfreq(nz, d=dz)
    k2 = kx[:, None] ** 2 + kz[None, :] ** 2
    denom = d * k2**2 + kf
    if np.any(np.abs(denom) == 0):
        raise ZeroDivisionError('Singular spectral operator.')

    qhat = np.fft.fft2(q)
    what = qhat / denom
    w = np.fft.ifft2(what)
    rf_hat = kf * what
    rb_hat = d * k2**2 * what
    rf = np.fft.ifft2(rf_hat)
    rb = np.fft.ifft2(rb_hat)
    rt = rf + rb
    residual = rt - q
    residual_norm = np.linalg.norm(residual.ravel()) / max(np.linalg.norm(q.ravel()), 1e-30)

    darea = dx * dz
    f_applied = np.sum(q) * darea
    f_reaction = np.sum(rt) * darea
    work_measure = 0.5 * float(np.real(np.vdot(w, q) * darea))
    dissipated = 0.0
    if omega_rad_s > 0.0:
        dissipated = 0.5 * omega_rad_s * float(np.imag(np.vdot(w, q) * darea))

    return SpectralSolution(
        x_m=x,
        z_m=z,
        load_pa=q,
        displacement_m=w,
        foundation_reaction_pa=rf,
        bending_reaction_pa=rb,
        total_reaction_pa=rt,
        residual_relative_l2=float(residual_norm),
        applied_resultant_n=complex(f_applied),
        reaction_resultant_n=complex(f_reaction),
        work_measure_j=work_measure,
        average_dissipated_power_w=dissipated,
    )


def uniform_load(q_pa: complex, nx: int, nz: int) -> np.ndarray:
    return np.full((nx, nz), q_pa, dtype=complex)


def cosine_load(q_mean_pa: complex, q_amplitude_pa: complex, geometry: Geometry, nx: int, nz: int, mode_x: int = 1, mode_z: int = 1) -> np.ndarray:
    x = np.arange(nx) * geometry.length_x_m / nx
    z = np.arange(nz) * geometry.length_z_m / nz
    shape = np.cos(2.0 * np.pi * mode_x * x[:, None] / geometry.length_x_m) * np.cos(
        2.0 * np.pi * mode_z * z[None, :] / geometry.length_z_m
    )
    return q_mean_pa + q_amplitude_pa * shape


def projected_cosine_amplitude(field: np.ndarray, geometry: Geometry, mode_x: int = 1, mode_z: int = 1) -> complex:
    nx, nz = field.shape
    x = np.arange(nx) * geometry.length_x_m / nx
    z = np.arange(nz) * geometry.length_z_m / nz
    shape = np.cos(2.0 * np.pi * mode_x * x[:, None] / geometry.length_x_m) * np.cos(
        2.0 * np.pi * mode_z * z[None, :] / geometry.length_z_m
    )
    return complex(np.vdot(shape, field) / np.vdot(shape, shape))
