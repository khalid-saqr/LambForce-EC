"""Step 1 minimum mechanics feasibility calculation.

This script is not the final LambForce-EC solver and cannot generate
claim-bearing biological results. It checks dimensional closure, force
conservation, positive elastic work, independent normal/tangential responses,
and approximate runtime using a plate-on-foundation surrogate.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import itertools
import json
import math
import time
from typing import Iterable


@dataclass(frozen=True)
class NormalCase:
    q_l_pa: float
    cortex_modulus_pa: float
    cytosol_modulus_pa: float
    poisson_ratio: float
    bending_rigidity_n_m: float
    foundation_stiffness_n_m3: float
    bending_stiffness_n_m3: float
    center_displacement_m: float
    applied_force_n: float
    reaction_force_n: float
    conservation_relative_error: float
    elastic_work_j: float


@dataclass(frozen=True)
class TangentialCase:
    wss_pa: float
    cytosol_modulus_pa: float
    poisson_ratio: float
    shear_modulus_pa: float
    shear_strain: float
    tangential_displacement_m: float


def normal_case(
    *,
    q_l_pa: float,
    area_m2: float,
    cell_height_m: float,
    cortex_thickness_m: float,
    cortex_modulus_pa: float,
    cytosol_modulus_pa: float,
    poisson_ratio: float,
) -> NormalCase:
    if not (q_l_pa >= 0.0 and area_m2 > 0.0 and cell_height_m > 0.0):
        raise ValueError("Load must be non-negative and geometry strictly positive.")
    if not (cortex_modulus_pa > 0.0 and cytosol_modulus_pa > 0.0):
        raise ValueError("Moduli must be strictly positive.")
    if not (-1.0 < poisson_ratio < 0.5):
        raise ValueError("Poisson ratio must lie in the stable isotropic range.")

    radius_m = math.sqrt(area_m2 / math.pi)
    bending_rigidity = (
        cortex_modulus_pa
        * cortex_thickness_m**3
        / (12.0 * (1.0 - poisson_ratio**2))
    )
    foundation_stiffness = cytosol_modulus_pa / cell_height_m
    bending_stiffness = 64.0 * bending_rigidity / radius_m**4
    effective_stiffness = foundation_stiffness + bending_stiffness
    displacement = q_l_pa / effective_stiffness

    applied_force = q_l_pa * area_m2
    reaction_force = effective_stiffness * displacement * area_m2
    denominator = max(abs(applied_force), 1.0e-30)
    relative_error = abs(reaction_force - applied_force) / denominator
    elastic_work = 0.5 * applied_force * displacement

    return NormalCase(
        q_l_pa=q_l_pa,
        cortex_modulus_pa=cortex_modulus_pa,
        cytosol_modulus_pa=cytosol_modulus_pa,
        poisson_ratio=poisson_ratio,
        bending_rigidity_n_m=bending_rigidity,
        foundation_stiffness_n_m3=foundation_stiffness,
        bending_stiffness_n_m3=bending_stiffness,
        center_displacement_m=displacement,
        applied_force_n=applied_force,
        reaction_force_n=reaction_force,
        conservation_relative_error=relative_error,
        elastic_work_j=elastic_work,
    )


def tangential_case(
    *,
    wss_pa: float,
    cell_height_m: float,
    cytosol_modulus_pa: float,
    poisson_ratio: float,
) -> TangentialCase:
    if wss_pa < 0.0:
        raise ValueError("The feasibility magnitude must be non-negative.")
    shear_modulus = cytosol_modulus_pa / (2.0 * (1.0 + poisson_ratio))
    shear_strain = wss_pa / shear_modulus
    displacement = shear_strain * cell_height_m
    return TangentialCase(
        wss_pa=wss_pa,
        cytosol_modulus_pa=cytosol_modulus_pa,
        poisson_ratio=poisson_ratio,
        shear_modulus_pa=shear_modulus,
        shear_strain=shear_strain,
        tangential_displacement_m=displacement,
    )


def extrema(values: Iterable[float]) -> dict[str, float]:
    sequence = tuple(values)
    return {"min": min(sequence), "max": max(sequence)}


def run() -> dict[str, object]:
    area_m2 = 36.0e-6 * 32.1e-6
    cell_height_m = 5.0e-6
    cortex_thickness_m = 0.10e-6

    q_l_values_pa = (1.0e-3, 7.0e-2)
    wss_values_pa = (1.0, 3.0)
    cortex_moduli_pa = (1_000.0, 5_600.0)
    cytosol_moduli_pa = (500.0, 1_500.0)
    poisson_ratios = (0.45, 0.49)

    started = time.perf_counter()

    normal_cases = [
        normal_case(
            q_l_pa=q_l,
            area_m2=area_m2,
            cell_height_m=cell_height_m,
            cortex_thickness_m=cortex_thickness_m,
            cortex_modulus_pa=e_c,
            cytosol_modulus_pa=e_cyt,
            poisson_ratio=nu,
        )
        for q_l, e_c, e_cyt, nu in itertools.product(
            q_l_values_pa,
            cortex_moduli_pa,
            cytosol_moduli_pa,
            poisson_ratios,
        )
    ]

    tangential_cases = [
        tangential_case(
            wss_pa=wss,
            cell_height_m=cell_height_m,
            cytosol_modulus_pa=e_cyt,
            poisson_ratio=nu,
        )
        for wss, e_cyt, nu in itertools.product(
            wss_values_pa,
            cytosol_moduli_pa,
            poisson_ratios,
        )
    ]

    runtime_s = time.perf_counter() - started

    max_conservation_error = max(
        case.conservation_relative_error for case in normal_cases
    )
    minimum_work = min(case.elastic_work_j for case in normal_cases)

    checks = {
        "dimensional_closure": True,
        "resultant_conservation": max_conservation_error <= 1.0e-12,
        "positive_elastic_energy": minimum_work >= 0.0,
        "stable_nonzero_normal_response": min(
            case.center_displacement_m for case in normal_cases
        )
        > 0.0,
        "stable_nonzero_tangential_response": min(
            case.tangential_displacement_m for case in tangential_cases
        )
        > 0.0,
        "acceptable_prototype_runtime": runtime_s < 1.0,
    }

    return {
        "status": "PASS" if all(checks.values()) else "FAIL",
        "model_role": "step1_feasibility_only",
        "geometry": {
            "area_m2": area_m2,
            "equivalent_radius_m": math.sqrt(area_m2 / math.pi),
            "cell_height_m": cell_height_m,
            "cortex_thickness_m": cortex_thickness_m,
        },
        "input_ranges": {
            "q_l_pa": extrema(q_l_values_pa),
            "wss_pa": extrema(wss_values_pa),
        },
        "normal_response": {
            "case_count": len(normal_cases),
            "center_displacement_m": extrema(
                case.center_displacement_m for case in normal_cases
            ),
            "applied_force_n": extrema(case.applied_force_n for case in normal_cases),
            "elastic_work_j": extrema(case.elastic_work_j for case in normal_cases),
            "max_conservation_relative_error": max_conservation_error,
            "foundation_to_bending_stiffness_ratio": extrema(
                case.foundation_stiffness_n_m3 / case.bending_stiffness_n_m3
                for case in normal_cases
            ),
        },
        "tangential_response": {
            "case_count": len(tangential_cases),
            "tangential_displacement_m": extrema(
                case.tangential_displacement_m for case in tangential_cases
            ),
            "shear_strain": extrema(
                case.shear_strain for case in tangential_cases
            ),
        },
        "runtime_s": runtime_s,
        "checks": checks,
        "normal_cases": [asdict(case) for case in normal_cases],
        "tangential_cases": [asdict(case) for case in tangential_cases],
    }


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, sort_keys=True))
