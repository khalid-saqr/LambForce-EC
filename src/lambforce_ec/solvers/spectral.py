from __future__ import annotations

from dataclasses import dataclass
import numpy as np

from ..constitutive import plate_foundation_properties, sls_complex_modulus
from ..models import Geometry, MaterialState


@dataclass
class SpectralSolution:
    x_m: np.ndarray
    z_m: np.ndarray
    load_pa: np.ndarray
    displacement_cell_m: np.ndarray
    displacement_apical_top_m: np.ndarray
    glycocalyx_strain: np.ndarray
    glycocalyx_reaction_pa: np.ndarray
    foundation_reaction_pa: np.ndarray
    bending_reaction_pa: np.ndarray
    total_reaction_pa: np.ndarray
    curvature_x_m_inv: np.ndarray
    curvature_z_m_inv: np.ndarray
    curvature_xz_m_inv: np.ndarray
    strain_x: np.ndarray
    strain_z: np.ndarray
    shear_strain_xz: np.ndarray
    tension_x_n_m: np.ndarray
    tension_z_n_m: np.ndarray
    tension_xz_n_m: np.ndarray
    residual_relative_l2: float
    applied_resultant_n: complex
    reaction_resultant_n: complex
    work_measure_j: float
    average_dissipated_power_w: float


class SpectralPlateSolver:
    solver_id = "periodic_spectral_2d"
    solver_version = "1.0.0"

    def solve(
        self,
        load_pa: np.ndarray,
        geometry: Geometry,
        material: MaterialState,
        omega_rad_s: float = 0.0,
    ) -> SpectralSolution:
        q = np.asarray(load_pa, dtype=complex)
        if q.ndim != 2 or min(q.shape) < 4:
            raise ValueError("load_pa must be a 2D array with at least four points per axis.")
        geometry.validate()
        material.validate()
        nx, nz = q.shape
        dx = geometry.length_x_m / nx
        dz = geometry.length_z_m / nz
        x = np.arange(nx) * dx
        z = np.arange(nz) * dz

        d, kf, kg = plate_foundation_properties(geometry, material, omega_rad_s)
        ec = sls_complex_modulus(material.cortex, omega_rad_s)
        eg = sls_complex_modulus(material.glycocalyx, omega_rad_s)
        nu = material.poisson_ratio

        kx = 2 * np.pi * np.fft.fftfreq(nx, d=dx)
        kz = 2 * np.pi * np.fft.fftfreq(nz, d=dz)
        kx2 = kx[:, None]
        kz2 = kz[None, :]
        wave_number_sq = kx2**2 + kz2**2
        denominator = d * wave_number_sq**2 + kf
        if np.any(np.abs(denominator) < 1e-30):
            raise ZeroDivisionError("Singular spectral operator.")

        qhat = np.fft.fft2(q)
        what = qhat / denominator
        w_cell = np.fft.ifft2(what)
        foundation = np.fft.ifft2(kf * what)
        bending = np.fft.ifft2(d * wave_number_sq**2 * what)
        reaction = foundation + bending
        residual = reaction - q
        residual_relative_l2 = np.linalg.norm(residual.ravel()) / max(
            np.linalg.norm(q.ravel()), 1e-30
        )

        wxx = np.fft.ifft2(-(kx2**2) * what)
        wzz = np.fft.ifft2(-(kz2**2) * what)
        wxz = np.fft.ifft2(-(kx2 * kz2) * what)
        curvature_x = -wxx
        curvature_z = -wzz
        curvature_xz = -2 * wxz

        surface_z = geometry.cortex_thickness_m / 2
        strain_x = surface_z * curvature_x
        strain_z = surface_z * curvature_z
        shear_strain_xz = surface_z * curvature_xz
        plane_stress_factor = ec / (1 - nu**2)
        stress_x = plane_stress_factor * (strain_x + nu * strain_z)
        stress_z = plane_stress_factor * (strain_z + nu * strain_x)
        stress_xz = ec / (2 * (1 + nu)) * shear_strain_xz
        tension_x = stress_x * geometry.cortex_thickness_m
        tension_z = stress_z * geometry.cortex_thickness_m
        tension_xz = stress_xz * geometry.cortex_thickness_m

        glycocalyx_strain = q / eg
        glycocalyx_displacement = q / kg
        w_top = w_cell + glycocalyx_displacement

        darea = dx * dz
        applied = np.sum(q) * darea
        reacted = np.sum(reaction) * darea
        work = 0.5 * float(np.real(np.vdot(w_top, q) * darea))
        dissipated = 0.0
        if omega_rad_s > 0:
            dissipated = 0.5 * omega_rad_s * float(np.imag(np.vdot(w_top, q) * darea))

        return SpectralSolution(
            x_m=x,
            z_m=z,
            load_pa=q,
            displacement_cell_m=w_cell,
            displacement_apical_top_m=w_top,
            glycocalyx_strain=glycocalyx_strain,
            glycocalyx_reaction_pa=q,
            foundation_reaction_pa=foundation,
            bending_reaction_pa=bending,
            total_reaction_pa=reaction,
            curvature_x_m_inv=curvature_x,
            curvature_z_m_inv=curvature_z,
            curvature_xz_m_inv=curvature_xz,
            strain_x=strain_x,
            strain_z=strain_z,
            shear_strain_xz=shear_strain_xz,
            tension_x_n_m=tension_x,
            tension_z_n_m=tension_z,
            tension_xz_n_m=tension_xz,
            residual_relative_l2=float(residual_relative_l2),
            applied_resultant_n=complex(applied),
            reaction_resultant_n=complex(reacted),
            work_measure_j=work,
            average_dissipated_power_w=max(0.0, dissipated),
        )
