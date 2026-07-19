from __future__ import annotations

from dataclasses import dataclass
import numpy as np

from ..constitutive import plate_foundation_properties
from ..models import Geometry, MaterialState


@dataclass
class LumpedSolution:
    load_pa: complex
    displacement_cell_m: complex
    displacement_apical_top_m: complex
    foundation_reaction_pa: complex
    applied_resultant_n: complex
    reaction_resultant_n: complex
    work_measure_j: float


class LumpedFoundationSolver:
    solver_id = "lumped_0d"
    solver_version = "1.0.0"

    def solve(
        self,
        load_pa: complex,
        geometry: Geometry,
        material: MaterialState,
        omega_rad_s: float = 0.0,
    ) -> LumpedSolution:
        _, kf, kg = plate_foundation_properties(geometry, material, omega_rad_s)
        q = complex(load_pa)
        w_cell = q / kf
        w_top = w_cell + q / kg
        force = q * geometry.area_m2
        work = 0.5 * float(np.real(np.conjugate(w_top) * force))
        return LumpedSolution(
            load_pa=q,
            displacement_cell_m=w_cell,
            displacement_apical_top_m=w_top,
            foundation_reaction_pa=kf * w_cell,
            applied_resultant_n=force,
            reaction_resultant_n=kf * w_cell * geometry.area_m2,
            work_measure_j=work,
        )
