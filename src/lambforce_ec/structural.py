from __future__ import annotations

from itertools import product
import numpy as np

from .constitutive import sls_complex_modulus
from .exceptions import ValidationError
from .models import (
    GLYCOCALYX_STATES,
    Geometry,
    MaterialState,
    RunConfig,
    LATERAL_SUPPORTS,
    MEMBRANE_CORTEX_COUPLINGS,
    NUCLEAR_REPRESENTATIONS,
    LOAD_DISTRIBUTIONS,
)


def validate_correlated_glycocalyx(geometry: Geometry, material: MaterialState) -> None:
    geometry.validate()
    material.validate()
    if geometry.glycocalyx_state_id != material.glycocalyx_state_id:
        raise ValidationError("Geometry and material glycocalyx_state_id values must match.")
    expected_h, expected_e = GLYCOCALYX_STATES[geometry.glycocalyx_state_id]
    if not np.isclose(geometry.glycocalyx_thickness_m, expected_h, rtol=0, atol=1e-15):
        raise ValidationError(
            "glycocalyx_thickness_m violates the frozen correlated glycocalyx state."
        )
    if not np.isclose(material.glycocalyx.einf_pa, expected_e, rtol=0, atol=1e-12):
        raise ValidationError("glycocalyx modulus violates the frozen correlated glycocalyx state.")


def grid_coordinates(
    geometry: Geometry, config: RunConfig
) -> tuple[np.ndarray, np.ndarray, float, float]:
    geometry.validate()
    config.validate()
    if config.lateral_support == "periodic_monolayer":
        dx = geometry.length_x_m / config.nx
        dz = geometry.length_z_m / config.nz
        x = np.arange(config.nx) * dx
        z = np.arange(config.nz) * dz
    else:
        dx = geometry.length_x_m / (config.nx + 1)
        dz = geometry.length_z_m / (config.nz + 1)
        x = (np.arange(config.nx) + 1) * dx
        z = (np.arange(config.nz) + 1) * dz
    return x, z, dx, dz


def nuclear_mask(geometry: Geometry, config: RunConfig) -> np.ndarray:
    x, z, _, _ = grid_coordinates(geometry, config)
    x0 = 0.5 * geometry.length_x_m
    z0 = 0.5 * geometry.length_z_m
    ax = 0.5 * geometry.nucleus_axis_x_m
    az = 0.5 * geometry.nucleus_axis_z_m
    return ((x[:, None] - x0) / ax) ** 2 + ((z[None, :] - z0) / az) ** 2 <= 1.0


def foundation_stiffness_map(
    geometry: Geometry,
    material: MaterialState,
    config: RunConfig,
    omega_rad_s: float,
) -> np.ndarray:
    active = material.for_rheology(config.rheology_mode)
    base = sls_complex_modulus(active.cytosol, omega_rad_s) / geometry.cell_height_m
    values = np.full((config.nx, config.nz), base, dtype=complex)
    if config.nuclear_representation == "stiff_nuclear_region":
        nucleus = sls_complex_modulus(active.nucleus, omega_rad_s) / geometry.nucleus_height_m
        values[nuclear_mask(geometry, config)] = nucleus
    return values


def structural_ensemble() -> list[dict[str, str]]:
    """Return the complete frozen structural Cartesian product.

    Phase 0 chooses the full product rather than an outcome-selected subset. The lateral support
    determines the compatible solver at configuration generation time.
    """
    rows: list[dict[str, str]] = []
    for load, support, coupling, nucleus in product(
        sorted(LOAD_DISTRIBUTIONS),
        sorted(LATERAL_SUPPORTS),
        sorted(MEMBRANE_CORTEX_COUPLINGS),
        sorted(NUCLEAR_REPRESENTATIONS),
    ):
        rows.append(
            {
                "load_distribution": load,
                "lateral_support": support,
                "membrane_cortex_coupling": coupling,
                "nuclear_representation": nucleus,
                "prestress_state": "zero",
                "solver_id": (
                    "periodic_spectral_2d" if support == "periodic_monolayer" else "bounded_fd_2d"
                ),
            }
        )
    return rows
