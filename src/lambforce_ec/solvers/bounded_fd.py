from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve

from ..constitutive import plate_foundation_properties, sls_complex_modulus
from ..models import Geometry, MaterialState


@dataclass
class BoundedFDSolution:
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
    matrix_nnz: int
    boundary_condition: str


class BoundedFDPlateSolver:
    solver_id = "bounded_fd_2d"
    solver_version = "2.1.0-readiness"

    @staticmethod
    def _d2(n: int, h: float) -> sparse.csr_matrix:
        return sparse.diags(
            [np.ones(n - 1), -2 * np.ones(n), np.ones(n - 1)],
            [-1, 0, 1],
            format="csr",
        ) / h**2

    @staticmethod
    def _d4(n: int, h: float, boundary_condition: str) -> sparse.csr_matrix:
        """Fourth derivative with ghost-node elimination.

        `compliant_edge` is the simply-supported bound: w=0 and w''=0.
        `clamped_edge` is the clamped bound: w=0 and w'=0.
        """
        matrix = sparse.diags(
            [
                np.ones(n - 2),
                -4 * np.ones(n - 1),
                6 * np.ones(n),
                -4 * np.ones(n - 1),
                np.ones(n - 2),
            ],
            [-2, -1, 0, 1, 2],
            format="lil",
        )
        if boundary_condition == "compliant_edge":
            matrix[0, 0] = 5.0
            matrix[-1, -1] = 5.0
        elif boundary_condition == "clamped_edge":
            # Quartic ghost-node elimination using w=0 and w'=0 at each edge.
            # The first ghost value is 6*w1 - 2*w2 + w3/3, which yields a
            # fourth-derivative boundary row exact for clamped quartic fields.
            matrix[0, 0] = 12.0
            matrix[0, 1] = -6.0
            matrix[0, 2] = 4.0 / 3.0
            matrix[-1, -1] = 12.0
            matrix[-1, -2] = -6.0
            matrix[-1, -3] = 4.0 / 3.0
        else:
            raise ValueError(f"Unsupported bounded boundary condition: {boundary_condition}")
        return matrix.tocsr() / h**4

    @staticmethod
    def _curvatures(
        w: np.ndarray, dx: float, dz: float
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        padded = np.pad(w, ((1, 1), (1, 1)), mode="constant")
        wxx = (
            padded[2:, 1:-1]
            - 2 * padded[1:-1, 1:-1]
            + padded[:-2, 1:-1]
        ) / dx**2
        wzz = (
            padded[1:-1, 2:]
            - 2 * padded[1:-1, 1:-1]
            + padded[1:-1, :-2]
        ) / dz**2
        wxz = (
            padded[2:, 2:]
            - padded[2:, :-2]
            - padded[:-2, 2:]
            + padded[:-2, :-2]
        ) / (4 * dx * dz)
        return -wxx, -wzz, -2 * wxz

    def solve(
        self,
        load_pa: np.ndarray,
        geometry: Geometry,
        material: MaterialState,
        omega_rad_s: float = 0.0,
        foundation_stiffness_map_n_m3: np.ndarray | None = None,
        boundary_condition: str = "compliant_edge",
        numerical_tolerance_relative: float = 1e-10,
    ) -> BoundedFDSolution:
        q = np.asarray(load_pa, dtype=complex)
        if q.ndim != 2 or min(q.shape) < 6:
            raise ValueError("load_pa must be a 2D interior-grid array with at least six points.")
        if numerical_tolerance_relative <= 0:
            raise ValueError("numerical_tolerance_relative must be positive.")
        geometry.validate()
        material.validate()
        nx, nz = q.shape
        dx = geometry.length_x_m / (nx + 1)
        dz = geometry.length_z_m / (nz + 1)
        x = (np.arange(nx) + 1) * dx
        z = (np.arange(nz) + 1) * dz
        d, kf, kg = plate_foundation_properties(geometry, material, omega_rad_s)
        ec = sls_complex_modulus(material.cortex, omega_rad_s)
        eg = sls_complex_modulus(material.glycocalyx, omega_rad_s)
        nu = material.poisson_ratio

        d2x = self._d2(nx, dx)
        d2z = self._d2(nz, dz)
        d4x = self._d4(nx, dx, boundary_condition)
        d4z = self._d4(nz, dz, boundary_condition)
        ix = sparse.eye(nx, format="csr")
        iz = sparse.eye(nz, format="csr")
        biharmonic = (
            sparse.kron(iz, d4x, format="csr")
            + 2 * sparse.kron(d2z, d2x, format="csr")
            + sparse.kron(d4z, ix, format="csr")
        )
        if foundation_stiffness_map_n_m3 is None:
            kmap = np.full(q.shape, kf, dtype=complex)
        else:
            kmap = np.asarray(foundation_stiffness_map_n_m3, dtype=complex)
            if kmap.shape != q.shape or np.any(np.abs(kmap) < 1e-30):
                raise ValueError("foundation_stiffness_map_n_m3 has an invalid shape or value.")
        operator = (
            d * biharmonic
            + sparse.diags(kmap.reshape(-1, order="F"), format="csr")
        ).tocsr()
        qvec = q.reshape(nx * nz, order="F")
        wvec = spsolve(operator, qvec)
        bending_vec = d * (biharmonic @ wvec)
        foundation_vec = kmap.reshape(-1, order="F") * wvec
        reaction = bending_vec + foundation_vec
        residual = reaction - qvec
        residual_norm = np.linalg.norm(residual) / max(np.linalg.norm(qvec), 1e-30)
        if not np.isfinite(residual_norm) or residual_norm > numerical_tolerance_relative:
            raise RuntimeError(
                "Bounded finite-difference solve failed the configured residual tolerance: "
                f"{residual_norm:.3e} > {numerical_tolerance_relative:.3e}."
            )
        w = wvec.reshape((nx, nz), order="F")
        curvature_x, curvature_z, curvature_xz = self._curvatures(w, dx, dz)
        surface_z = geometry.cortex_thickness_m / 2
        strain_x = surface_z * curvature_x
        strain_z = surface_z * curvature_z
        shear_strain_xz = surface_z * curvature_xz
        plane_stress_factor = ec / (1 - nu**2)
        tension_x = (
            plane_stress_factor
            * (strain_x + nu * strain_z)
            * geometry.cortex_thickness_m
        )
        tension_z = (
            plane_stress_factor
            * (strain_z + nu * strain_x)
            * geometry.cortex_thickness_m
        )
        tension_xz = (
            ec
            / (2 * (1 + nu))
            * shear_strain_xz
            * geometry.cortex_thickness_m
        )
        glycocalyx_strain = q / eg
        w_top = w + q / kg
        darea = dx * dz
        applied = np.sum(qvec) * darea
        reacted = np.sum(reaction) * darea
        work = 0.5 * float(np.real(np.vdot(w_top, q) * darea))
        dissipated = 0.0
        is_viscoelastic = any(
            component.ratio_e0_einf > 1
            for component in (
                material.cortex,
                material.cytosol,
                material.glycocalyx,
                material.nucleus,
            )
        )
        if omega_rad_s > 0 and is_viscoelastic:
            dissipated = 0.5 * omega_rad_s * float(
                np.imag(np.vdot(w_top, q) * darea)
            )
        return BoundedFDSolution(
            x_m=x,
            z_m=z,
            load_pa=q,
            displacement_cell_m=w,
            displacement_apical_top_m=w_top,
            glycocalyx_strain=glycocalyx_strain,
            glycocalyx_reaction_pa=q,
            foundation_reaction_pa=foundation_vec.reshape((nx, nz), order="F"),
            bending_reaction_pa=bending_vec.reshape((nx, nz), order="F"),
            total_reaction_pa=reaction.reshape((nx, nz), order="F"),
            curvature_x_m_inv=curvature_x,
            curvature_z_m_inv=curvature_z,
            curvature_xz_m_inv=curvature_xz,
            strain_x=strain_x,
            strain_z=strain_z,
            shear_strain_xz=shear_strain_xz,
            tension_x_n_m=tension_x,
            tension_z_n_m=tension_z,
            tension_xz_n_m=tension_xz,
            residual_relative_l2=float(residual_norm),
            applied_resultant_n=complex(applied),
            reaction_resultant_n=complex(reacted),
            work_measure_j=work,
            average_dissipated_power_w=max(0.0, dissipated),
            matrix_nnz=int(operator.nnz),
            boundary_condition=boundary_condition,
        )
