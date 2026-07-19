"""LambForce-EC production mechanics package."""

from .models import ArteryRecord, Geometry, MaterialState, RunConfig, SLSMaterial
from .workflow import SimulationResult, required_run_matrix, run_case, run_registered_comparison

__version__ = "0.4.0"

__all__ = [
    "ArteryRecord",
    "Geometry",
    "MaterialState",
    "RunConfig",
    "SLSMaterial",
    "SimulationResult",
    "run_case",
    "run_registered_comparison",
    "required_run_matrix",
]
