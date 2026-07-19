from __future__ import annotations

from dataclasses import dataclass
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve

from ..constitutive import plate_foundation_properties
from ..models import Geometry, MaterialState


@dataclass
class BoundedFDSolution:
    x_m: np.ndarray
    z_m: np.ndarray
    load_pa: np.ndarray
    displacement_cell_m: np.ndarray
    total_reaction_pa: np.ndarray
    residual_relative_l2: float
    applied_resultant_n: complex
    reaction_resultant_n: complex
    work_measure_j: float
    matrix_nnz: int


class BoundedFDPlateSolver:
    solver_id = "bounded_fd_2d"
    solver_version = "1.0.0"

    @staticmethod
    def _dirichlet_laplacian(nx: int, nz: int, dx: float, dz: float) -> sparse.csr_matrix:
        tx = sparse.diags(
            [np.ones(nx - 1), -2 * np.ones(nx), np.ones(nx - 1)], [-1, 0, 1], format="csr"
        ) / dx**2
        tz = sparse.diags(
            [np.ones(nz - 1), -2 * np.ones(nz), np.ones(nz - 1)], [-1, 0, 1], format="csr"
        ) / dz**2
        return sparse.kron(sparse.eye(nz, format="csr"), tx, format="csr") + sparse.kron(
            tz, sparse.eye(nx, format="csr"), format="csr"
        )

    def solve(
        self,
        load_pa: np.ndarray,
        geometry: Geometry,
        material: MaterialState,
        omega_rad_s: float = 0.0,
    ) -> BoundedFDSolution:
        q = np.asarray(load_pa, dtype=complex)
        if q.ndim != 2 or min(q.shape) < 4:
            raise ValueError("load_pa must be a 2D interior-grid array.")
        nx, nz = q.shape
        dx = geometry.length_x_m / (nx + 1)
        dz = geometry.length_z_m / (nz + 1)
        x = (np.arange(nx) + 1) * dx
        z = (np.arange(nz) + 1) * dz
        d, kf, _ = plate_foundation_properties(geometry, material, omega_rad_s)
        lap = self._dirichlet_laplacian(nx, nz, dx, dz)
        operator = (d * (lap @ lap) + kf * sparse.eye(nx * nz, dtype=complex, format="csr")).tocsr()
        qvec = q.reshape(nx * nz, order="F")
        wvec = spsolve(operator, qvec)
        reaction = operator @ wvec
        residual = reaction - qvec
        residual_norm = np.linalg.norm(residual) / max(np.linalg.norm(qvec), 1e-30)
        darea = dx * dz
        applied = np.sum(qvec) * darea
        reacted = np.sum(reaction) * darea
        work = 0.5 * float(np.real(np.vdot(wvec, qvec) * darea))
        return BoundedFDSolution(
            x_m=x,
            z_m=z,
            load_pa=q,
            displacement_cell_m=wvec.reshape((nx, nz), order="F"),
            total_reaction_pa=reaction.reshape((nx, nz), order="F"),
            residual_relative_l2=float(residual_norm),
            applied_resultant_n=complex(applied),
            reaction_resultant_n=complex(reacted),
            work_measure_j=work,
            matrix_nnz=int(operator.nnz),
        )
