from __future__ import annotations

from typing import Protocol, runtime_checkable
import numpy as np

from .models import Geometry, MaterialState


@runtime_checkable
class NormalMechanicsSolver(Protocol):
    """Replaceable normal-mechanics interface used by study workflows."""

    solver_id: str
    solver_version: str

    def solve(
        self,
        load_pa: np.ndarray,
        geometry: Geometry,
        material: MaterialState,
        omega_rad_s: float = 0.0,
    ) -> object: ...


@runtime_checkable
class ArteryInputAdapter(Protocol):
    """Interface for immutable hydrodynamic archives or exported field records."""

    def load(self, source: str) -> object: ...


@runtime_checkable
class SpatialLoadModel(Protocol):
    """Interface for force-conserving apical load distributions."""

    model_id: str

    def kernel(self, geometry: Geometry, nx: int, nz: int) -> np.ndarray: ...
