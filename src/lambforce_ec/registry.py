from __future__ import annotations

from .exceptions import ValidationError
from .solvers import BoundedFDPlateSolver, LumpedFoundationSolver, SpectralPlateSolver


def solver_registry() -> dict[str, type]:
    """Return registered solver implementations without binding studies to classes."""
    return {
        SpectralPlateSolver.solver_id: SpectralPlateSolver,
        BoundedFDPlateSolver.solver_id: BoundedFDPlateSolver,
        LumpedFoundationSolver.solver_id: LumpedFoundationSolver,
    }


def create_solver(solver_id: str):
    try:
        return solver_registry()[solver_id]()
    except KeyError as exc:
        raise ValidationError(f"Unknown solver_id: {solver_id}") from exc
