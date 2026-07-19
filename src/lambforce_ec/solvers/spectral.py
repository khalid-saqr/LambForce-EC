from __future__ import annotations

from dataclasses import dataclass
import numpy as np
from scipy.sparse.linalg import LinearOperator, gmres

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
    iterative_solver_info: int = 0


class SpectralPlateSolver:
    solver_id = "periodic_spectral_2d"
    solver_version = "2.0.0-phase0"

    def solve(
        self,
        load_pa: np.ndarray,
        geometry: Geometry,
        material: MaterialState,
        omega_rad_s: float = 0.0,
        foundation_stiffness_map_n_m3: np.ndarray | None = None,
        numerical_tolerance_relative: float = 1e-10,
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

        d, kf_scalar, kg = plate_foundation_properties(geometry, material, omega_rad_s)
        ec = sls_complex_modulus(material.cortex, omega_rad_s)
        eg = sls_complex_modulus(material.glycocalyx, omega_rad_s)
        nu = material.poisson_ratio

        kx = 2 * np.pi * np.fft.fftfreq(nx, d=dx)
        kz = 2 * np.pi * np.fft.fftfreq(nz, d=dz)
        kx2 = kx[:, None]
        kz2 = kz[None, :]
        wave_number_sq = kx2**2 + kz2**2
        biharmonic_symbol = wave_number_sq**2

        iterative_info = 0
        if foundation_stiffness_map_n_m3 is None:
            kmap = np.full(q.shape, kf_scalar, dtype=complex)
            denominator = d * biharmonic_symbol + kf_scalar
            if np.any(np.abs(denominator) < 1e-30):
                raise ZeroDivisionError("Singular spectral operator.")
            qhat = np.fft.fft2(q)
            what = qhat / denominator
            w_cell = np.fft.ifft2(what)
        else:
            kmap = np.asarray(foundation_stiffness_map_n_m3, dtype=complex)
            if kmap.shape != q.shape or np.any(np.abs(kmap) < 1e-30):
                raise ValueError("foundation_stiffness_map_n_m3 has an invalid shape or value.")
            mean_k = complex(np.mean(kmap))
            denominator = d * biharmonic_symbol + mean_k

            def apply(vector: np.ndarray) -> np.ndarray:
                field = vector.reshape(q.shape)
                bending = np.fft.ifft2(d * biharmonic_symbol * np.fft.fft2(field))
                return (bending + kmap * field).reshape(-1)

            def precondition(vector: np.ndarray) -> np.ndarray:
                field = vector.reshape(q.shape)
                return np.fft.ifft2(np.fft.fft2(field) / denominator).reshape(-1)

            operator = LinearOperator((q.size, q.size), matvec=apply, dtype=complex)
            preconditioner = LinearOperator((q.size, q.size), matvec=precondition, dtype=complex)
            try:
                solution, iterative_info = gmres(
                    operator,
                    q.reshape(-1),
                    M=preconditioner,
                    rtol=numerical_tolerance_relative,
                    atol=0.0,
                    restart=min(100, q.size),
                    maxiter=500,
                )
            except TypeError:  # SciPy 1.10/1.11 compatibility: ``tol`` preceded ``rtol``.
                solution, iterative_info = gmres(
                    operator,
                    q.reshape(-1),
                    M=preconditioner,
                    tol=numerical_tolerance_relative,
                    atol=0.0,
                    restart=min(100, q.size),
                    maxiter=500,
                )
            if iterative_info != 0:
                raise RuntimeError(f"Variable-foundation GMRES did not converge: info={iterative_info}")
            w_cell = solution.reshape(q.shape)
            what = np.fft.fft2(w_cell)

        foundation = kmap * w_cell
        bending = np.fft.ifft2(d * biharmonic_symbol * what)
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
        is_viscoelastic = any(
            component.ratio_e0_einf > 1
            for component in (material.cortex, material.cytosol, material.glycocalyx, material.nucleus)
        )
        if omega_rad_s > 0 and is_viscoelastic:
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
            iterative_solver_info=int(iterative_info),
        )
