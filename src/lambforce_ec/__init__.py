"""LambForce-EC production mechanics and published-source reproduction package."""

from .models import ArteryRecord, Geometry, MaterialState, RunConfig, SLSMaterial
from .published_source import (
    reproduce_all_six,
    reproduce_artery,
    validate_published_inputs,
    verify_reproduction_directory,
)
from .workflow import SimulationResult, required_run_matrix, run_case, run_registered_comparison

__version__ = "0.6.0"

__all__ = [
    "ArteryRecord",
    "Geometry",
    "MaterialState",
    "RunConfig",
    "SLSMaterial",
    "SimulationResult",
    "reproduce_all_six",
    "reproduce_artery",
    "validate_published_inputs",
    "verify_reproduction_directory",
    "run_case",
    "run_registered_comparison",
    "required_run_matrix",
]
