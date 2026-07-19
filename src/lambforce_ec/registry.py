from __future__ import annotations

from .exceptions import ValidationError
from .models import (
    HARMONIC_CONTROLS,
    LATERAL_SUPPORTS,
    LOAD_CASES,
    LOAD_DISTRIBUTIONS,
    MEMBRANE_CORTEX_COUPLINGS,
    NUCLEAR_REPRESENTATIONS,
    RHEOLOGY_MODES,
)
from .solvers import BoundedFDPlateSolver, LumpedFoundationSolver, SpectralPlateSolver


def solver_registry() -> dict[str, type]:
    return {
        SpectralPlateSolver.solver_id: SpectralPlateSolver,
        BoundedFDPlateSolver.solver_id: BoundedFDPlateSolver,
        LumpedFoundationSolver.solver_id: LumpedFoundationSolver,
    }


def protocol_surface_registry() -> dict[str, tuple[str, ...]]:
    return {
        "load_case": tuple(sorted(LOAD_CASES)),
        "harmonic_control": tuple(sorted(HARMONIC_CONTROLS)),
        "load_distribution": tuple(sorted(LOAD_DISTRIBUTIONS)),
        "lateral_support": tuple(sorted(LATERAL_SUPPORTS)),
        "membrane_cortex_coupling": tuple(sorted(MEMBRANE_CORTEX_COUPLINGS)),
        "nuclear_representation": tuple(sorted(NUCLEAR_REPRESENTATIONS)),
        "rheology_mode": tuple(sorted(RHEOLOGY_MODES)),
        "prestress_state": ("zero",),
    }


def create_solver(solver_id: str):
    try:
        return solver_registry()[solver_id]()
    except KeyError as exc:
        raise ValidationError(f"Unknown solver_id: {solver_id}") from exc
