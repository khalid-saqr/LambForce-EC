from __future__ import annotations

import numpy as np
from scipy.integrate import trapezoid

from .exceptions import ValidationError
from .models import ArteryRecord, Geometry, RunConfig
from .validation import assert_conserved


def integrate_radial_density(radial_coordinate_m: np.ndarray, density_n_m3: np.ndarray) -> np.ndarray:
    r = np.asarray(radial_coordinate_m, dtype=float)
    f = np.asarray(density_n_m3, dtype=float)
    if r.ndim != 1 or f.ndim != 2 or f.shape[0] != r.size:
        raise ValidationError("density must have shape (n_radial, n_time).")
    return trapezoid(f, x=r, axis=0)


def extract_scalar_loads(record: ArteryRecord) -> dict[str, np.ndarray]:
    record.validate()
    total = integrate_radial_density(record.radial_coordinate_m, record.lamb_density_signed_n_m3)
    exposure = integrate_radial_density(
        record.radial_coordinate_m, np.abs(record.lamb_density_signed_n_m3)
    )
    if record.lamb_density_isotropic_n_m3 is None:
        isotropic = np.zeros_like(total)
    else:
        isotropic = integrate_radial_density(
            record.radial_coordinate_m, record.lamb_density_isotropic_n_m3
        )
    return {
        "lamb_signed_pa": total,
        "lamb_isotropic_pa": isotropic,
        "lamb_anisotropy_increment_pa": total - isotropic,
        "lamb_exposure_pa": exposure,
        "wss_pa": np.asarray(record.wall_shear_stress_pa, dtype=float),
    }


def select_load_case(loads: dict[str, np.ndarray], load_case: str) -> tuple[np.ndarray, np.ndarray]:
    zero = np.zeros_like(loads["wss_pa"])
    if load_case == "unloaded":
        return zero, zero
    if load_case == "wss_only":
        return zero, loads["wss_pa"]
    if load_case == "lamb_signed_only":
        return loads["lamb_signed_pa"], zero
    if load_case == "wss_plus_lamb_signed":
        return loads["lamb_signed_pa"], loads["wss_pa"]
    if load_case == "isotropic_lamb":
        return loads["lamb_isotropic_pa"], loads["wss_pa"]
    if load_case == "anisotropy_increment":
        return loads["lamb_anisotropy_increment_pa"], loads["wss_pa"]
    if load_case == "exposure_diagnostic":
        return loads["lamb_exposure_pa"], zero
    raise ValidationError(f"Unknown load_case: {load_case}")


def periodic_grid(geometry: Geometry, nx: int, nz: int) -> tuple[np.ndarray, np.ndarray, float, float]:
    geometry.validate()
    x = np.arange(nx) * geometry.length_x_m / nx
    z = np.arange(nz) * geometry.length_z_m / nz
    return x, z, geometry.length_x_m / nx, geometry.length_z_m / nz


def spatial_kernel(
    geometry: Geometry,
    config: RunConfig,
    glycocalyx_thickness_field_m: np.ndarray | None = None,
) -> np.ndarray:
    config.validate()
    x, z, dx, dz = periodic_grid(geometry, config.nx, config.nz)
    if config.load_distribution == "uniform_apical":
        kernel = np.ones((config.nx, config.nz), dtype=float)
    elif config.load_distribution == "localized_bound":
        # Periodic minimum-image Gaussian; widths are preregistered geometry fractions.
        x0, z0 = 0.5 * geometry.length_x_m, 0.5 * geometry.length_z_m
        dxp = np.minimum(np.abs(x[:, None] - x0), geometry.length_x_m - np.abs(x[:, None] - x0))
        dzp = np.minimum(np.abs(z[None, :] - z0), geometry.length_z_m - np.abs(z[None, :] - z0))
        sx = config.localized_sigma_x_fraction * geometry.length_x_m
        sz = config.localized_sigma_z_fraction * geometry.length_z_m
        kernel = np.exp(-0.5 * ((dxp / sx) ** 2 + (dzp / sz) ** 2))
    elif config.load_distribution == "glycocalyx_resolved":
        if glycocalyx_thickness_field_m is None:
            raise ValidationError(
                "glycocalyx_resolved requires a positive glycocalyx_thickness_field_m array."
            )
        h = np.asarray(glycocalyx_thickness_field_m, dtype=float)
        if h.shape != (config.nx, config.nz) or not np.all(np.isfinite(h)) or np.any(h <= 0):
            raise ValidationError("glycocalyx thickness field has an invalid shape or value.")
        # Common-displacement spring limit: local transmitted traction scales as E/h.
        kernel = 1.0 / h
    else:
        raise ValidationError(f"Unknown load_distribution: {config.load_distribution}")

    integral = np.sum(kernel) * dx * dz
    if integral <= 0 or not np.isfinite(integral):
        raise ValidationError("spatial kernel has non-positive integral.")
    kernel = kernel * geometry.area_m2 / integral
    assert_conserved(
        applied=geometry.area_m2,
        reaction=np.sum(kernel) * dx * dz,
        tolerance_relative=config.conservation_tolerance_relative,
        label="spatial kernel",
    )
    return kernel


def distribute_scalar_waveform(q_pa: np.ndarray, kernel: np.ndarray) -> np.ndarray:
    q = np.asarray(q_pa, dtype=float)
    k = np.asarray(kernel, dtype=float)
    if q.ndim != 1 or k.ndim != 2:
        raise ValidationError("q_pa must be one-dimensional and kernel two-dimensional.")
    return q[:, None, None] * k[None, :, :]
