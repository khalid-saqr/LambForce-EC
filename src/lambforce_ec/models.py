from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any
import re
import numpy as np

from .exceptions import ValidationError

_SHA256 = re.compile(r"^[0-9a-f]{64}$")


@dataclass(frozen=True)
class Geometry:
    length_x_m: float = 36.0e-6
    length_z_m: float = 32.1e-6
    cell_height_m: float = 5.0e-6
    cortex_thickness_m: float = 0.10e-6
    glycocalyx_thickness_m: float = 0.50e-6

    @property
    def area_m2(self) -> float:
        return self.length_x_m * self.length_z_m

    def validate(self) -> None:
        values = {
            "length_x_m": self.length_x_m,
            "length_z_m": self.length_z_m,
            "cell_height_m": self.cell_height_m,
            "cortex_thickness_m": self.cortex_thickness_m,
            "glycocalyx_thickness_m": self.glycocalyx_thickness_m,
        }
        for name, value in values.items():
            if not np.isfinite(value) or value <= 0:
                raise ValidationError(f"{name} must be finite and positive.")
        if self.cortex_thickness_m >= self.cell_height_m:
            raise ValidationError("cortex_thickness_m must be smaller than cell_height_m.")


@dataclass(frozen=True)
class SLSMaterial:
    einf_pa: float
    ratio_e0_einf: float = 2.0
    relaxation_time_s: float = 1.0

    def validate(self, name: str) -> None:
        if not np.isfinite(self.einf_pa) or self.einf_pa <= 0:
            raise ValidationError(f"{name}.einf_pa must be finite and positive.")
        if not np.isfinite(self.ratio_e0_einf) or self.ratio_e0_einf < 1:
            raise ValidationError(f"{name}.ratio_e0_einf must be at least one.")
        if not np.isfinite(self.relaxation_time_s) or self.relaxation_time_s <= 0:
            raise ValidationError(f"{name}.relaxation_time_s must be finite and positive.")


@dataclass(frozen=True)
class MaterialState:
    cortex: SLSMaterial = field(default_factory=lambda: SLSMaterial(1_000.0, 2.0, 0.05))
    cytosol: SLSMaterial = field(default_factory=lambda: SLSMaterial(500.0, 2.0, 2.0))
    glycocalyx: SLSMaterial = field(default_factory=lambda: SLSMaterial(390.0, 2.0, 0.10))
    nucleus: SLSMaterial = field(default_factory=lambda: SLSMaterial(5_000.0, 2.0, 0.25))
    poisson_ratio: float = 0.45
    parameter_set_id: str = "reference_reduced_model_v1"

    def validate(self) -> None:
        self.cortex.validate("cortex")
        self.cytosol.validate("cytosol")
        self.glycocalyx.validate("glycocalyx")
        self.nucleus.validate("nucleus")
        if not (-1.0 < self.poisson_ratio < 0.5):
            raise ValidationError("poisson_ratio is outside the stable isotropic range.")
        if not self.parameter_set_id.strip():
            raise ValidationError("parameter_set_id must be non-empty.")


@dataclass
class ArteryRecord:
    artery_id: str
    artery_name: str
    radius_m: float
    omega0_rad_s: float
    radial_coordinate_m: np.ndarray
    time_s: np.ndarray
    lamb_density_signed_n_m3: np.ndarray
    wall_shear_stress_pa: np.ndarray
    lamb_density_isotropic_n_m3: np.ndarray | None
    source_identifier: str
    source_version: str
    source_checksum: str
    coordinate_convention: str = "outward_normal_positive"
    metadata: dict[str, Any] = field(default_factory=dict)

    def validate(self) -> None:
        if not self.artery_id.strip() or not self.artery_name.strip():
            raise ValidationError("artery_id and artery_name must be non-empty.")
        if not np.isfinite(self.radius_m) or self.radius_m <= 0:
            raise ValidationError("radius_m must be finite and positive.")
        if not np.isfinite(self.omega0_rad_s) or self.omega0_rad_s <= 0:
            raise ValidationError("omega0_rad_s must be finite and positive.")
        if not _SHA256.match(self.source_checksum):
            raise ValidationError("source_checksum must be a lowercase SHA-256 hexadecimal digest.")
        r = np.asarray(self.radial_coordinate_m, dtype=float)
        t = np.asarray(self.time_s, dtype=float)
        f = np.asarray(self.lamb_density_signed_n_m3, dtype=float)
        tau = np.asarray(self.wall_shear_stress_pa, dtype=float)
        if r.ndim != 1 or t.ndim != 1 or r.size < 2 or t.size < 4:
            raise ValidationError("radial_coordinate_m and time_s must be one-dimensional arrays.")
        if not np.all(np.isfinite(r)) or not np.all(np.isfinite(t)):
            raise ValidationError("coordinates must be finite.")
        if not np.all(np.diff(r) > 0):
            raise ValidationError("radial_coordinate_m must be strictly increasing.")
        if not np.all(np.diff(t) > 0):
            raise ValidationError("time_s must be strictly increasing and endpoint-exclusive.")
        if r[-1] > self.radius_m * (1 + 1e-10):
            raise ValidationError("radial coordinates extend beyond radius_m.")
        if f.shape != (r.size, t.size):
            raise ValidationError("lamb_density_signed_n_m3 must have shape (n_radial, n_time).")
        if tau.shape != (t.size,):
            raise ValidationError("wall_shear_stress_pa must have shape (n_time,).")
        if not np.all(np.isfinite(f)) or not np.all(np.isfinite(tau)):
            raise ValidationError("hydrodynamic arrays must be finite.")
        if self.lamb_density_isotropic_n_m3 is not None:
            iso = np.asarray(self.lamb_density_isotropic_n_m3, dtype=float)
            if iso.shape != f.shape or not np.all(np.isfinite(iso)):
                raise ValidationError("isotropic Lamb density must match the signed-field shape.")
        dt = np.diff(t)
        if np.max(np.abs(dt - dt[0])) > max(1e-12, 1e-8 * abs(dt[0])):
            raise ValidationError("time_s must be uniformly sampled for harmonic reconstruction.")
        expected_period = 2 * np.pi / self.omega0_rad_s
        sampled_period = dt[0] * t.size
        if abs(sampled_period - expected_period) / expected_period > 5e-3:
            raise ValidationError("time grid period is inconsistent with omega0_rad_s.")


@dataclass(frozen=True)
class RunConfig:
    nx: int = 32
    nz: int = 32
    solver_id: str = "periodic_spectral_2d"
    structural_model_id: str = "reduced_plate_foundation_series_glycocalyx_v1"
    load_distribution: str = "uniform_apical"
    load_case: str = "wss_plus_lamb_signed"
    harmonic_control: str = "full_waveform"
    localized_sigma_x_fraction: float = 0.15
    localized_sigma_z_fraction: float = 0.15
    conservation_tolerance_relative: float = 1e-12
    protocol_version: str = "1.0.0-step1"
    parameter_freeze_version: str = "2.0.0"

    def validate(self) -> None:
        if self.nx < 4 or self.nz < 4:
            raise ValidationError("nx and nz must each be at least four.")
        if self.solver_id not in {"periodic_spectral_2d", "bounded_fd_2d", "lumped_0d"}:
            raise ValidationError(f"Unsupported solver_id: {self.solver_id}")
        if self.load_distribution not in {"uniform_apical", "localized_bound", "glycocalyx_resolved"}:
            raise ValidationError(f"Unsupported load_distribution: {self.load_distribution}")
        if self.load_case not in {
            "unloaded", "wss_only", "lamb_signed_only", "wss_plus_lamb_signed",
            "isotropic_lamb", "anisotropy_increment", "exposure_diagnostic",
        }:
            raise ValidationError(f"Unsupported load_case: {self.load_case}")
        if self.harmonic_control not in {"fundamental_only", "harmonics_le_2", "full_waveform"}:
            raise ValidationError(f"Unsupported harmonic_control: {self.harmonic_control}")
        if not (0 < self.localized_sigma_x_fraction <= 0.5):
            raise ValidationError("localized_sigma_x_fraction must be in (0, 0.5].")
        if not (0 < self.localized_sigma_z_fraction <= 0.5):
            raise ValidationError("localized_sigma_z_fraction must be in (0, 0.5].")
        if self.conservation_tolerance_relative <= 0:
            raise ValidationError("conservation_tolerance_relative must be positive.")
