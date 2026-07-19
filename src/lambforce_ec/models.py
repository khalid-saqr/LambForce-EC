from __future__ import annotations

from dataclasses import dataclass, field, replace
from typing import Any
import re

import numpy as np

from .exceptions import ValidationError


_SHA256 = re.compile(r"^[0-9a-f]{64}$")
_GIT_SHA = re.compile(r"^[0-9a-f]{40}$")

GLYCOCALYX_STATES: dict[str, tuple[float, float]] = {
    "thin_stiff": (0.11e-6, 1000.0),
    "reference": (0.50e-6, 390.0),
    "thick_soft": (1.00e-6, 25.0),
}

LOAD_CASES = {
    "unloaded",
    "wss_only",
    "lamb_signed_only",
    "wss_plus_lamb_signed",
    "isotropic_lamb",
    "anisotropy_increment",
    "exposure_diagnostic",
    "inward_only",
    "outward_only",
    "zero_normal_load",
}
LATERAL_SUPPORTS = {"periodic_monolayer", "compliant_edge", "clamped_edge"}
MEMBRANE_CORTEX_COUPLINGS = {"perfectly_bonded", "tangential_slip_limit"}
NUCLEAR_REPRESENTATIONS = {"homogeneous_cell_body", "stiff_nuclear_region"}
RHEOLOGY_MODES = {"elastic", "sls"}
HARMONIC_CONTROLS = {"fundamental_only", "harmonics_le_2", "full_waveform"}
LOAD_DISTRIBUTIONS = {"uniform_apical", "localized_bound", "glycocalyx_resolved"}


@dataclass(frozen=True)
class Geometry:
    length_x_m: float = 36.0e-6
    length_z_m: float = 32.1e-6
    cell_height_m: float = 5.0e-6
    cortex_thickness_m: float = 0.10e-6
    glycocalyx_thickness_m: float = 0.50e-6
    glycocalyx_state_id: str = "reference"
    nucleus_axis_x_m: float = 8.0e-6
    nucleus_axis_z_m: float = 6.0e-6
    nucleus_height_m: float = 2.50e-6

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
            "nucleus_axis_x_m": self.nucleus_axis_x_m,
            "nucleus_axis_z_m": self.nucleus_axis_z_m,
            "nucleus_height_m": self.nucleus_height_m,
        }
        for name, value in values.items():
            if not np.isfinite(value) or value <= 0:
                raise ValidationError(f"{name} must be finite and positive.")
        if self.cortex_thickness_m >= self.cell_height_m:
            raise ValidationError("cortex_thickness_m must be smaller than cell_height_m.")
        if self.nucleus_height_m > self.cell_height_m:
            raise ValidationError("nucleus_height_m cannot exceed cell_height_m.")
        if self.nucleus_axis_x_m >= self.length_x_m or self.nucleus_axis_z_m >= self.length_z_m:
            raise ValidationError("nuclear in-plane axes must fit inside the endothelial footprint.")
        if self.glycocalyx_state_id not in GLYCOCALYX_STATES:
            raise ValidationError(f"Unknown glycocalyx_state_id: {self.glycocalyx_state_id}")


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

    def elastic_limit(self) -> "SLSMaterial":
        return replace(self, ratio_e0_einf=1.0)


@dataclass(frozen=True)
class MaterialState:
    cortex: SLSMaterial = field(default_factory=lambda: SLSMaterial(1_000.0, 2.0, 0.05))
    cytosol: SLSMaterial = field(default_factory=lambda: SLSMaterial(500.0, 2.0, 2.0))
    glycocalyx: SLSMaterial = field(default_factory=lambda: SLSMaterial(390.0, 2.0, 0.10))
    nucleus: SLSMaterial = field(default_factory=lambda: SLSMaterial(5_000.0, 2.0, 0.25))
    poisson_ratio: float = 0.45
    parameter_set_id: str = "reference_reduced_model_v2"
    glycocalyx_state_id: str = "reference"

    def validate(self) -> None:
        self.cortex.validate("cortex")
        self.cytosol.validate("cytosol")
        self.glycocalyx.validate("glycocalyx")
        self.nucleus.validate("nucleus")
        if not (-1.0 < self.poisson_ratio < 0.5):
            raise ValidationError("poisson_ratio is outside the stable isotropic range.")
        if not self.parameter_set_id.strip():
            raise ValidationError("parameter_set_id must be non-empty.")
        if self.glycocalyx_state_id not in GLYCOCALYX_STATES:
            raise ValidationError(f"Unknown glycocalyx_state_id: {self.glycocalyx_state_id}")

    def for_rheology(self, mode: str) -> "MaterialState":
        if mode not in RHEOLOGY_MODES:
            raise ValidationError(f"Unsupported rheology_mode: {mode}")
        if mode == "sls":
            return self
        return replace(
            self,
            cortex=self.cortex.elastic_limit(),
            cytosol=self.cytosol.elastic_limit(),
            glycocalyx=self.glycocalyx.elastic_limit(),
            nucleus=self.nucleus.elastic_limit(),
        )


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
    source_member_sha256: str | None = None
    conversion_manifest_sha256: str | None = None
    converter_commit_sha: str | None = None
    record_payload_sha256: str | None = None
    metadata: dict[str, Any] = field(default_factory=dict)

    @property
    def source_archive_sha256(self) -> str:
        """Explicit alias retained for backward-compatible callers."""
        return self.source_checksum

    def validate(self) -> None:
        if not self.artery_id.strip() or not self.artery_name.strip():
            raise ValidationError("artery_id and artery_name must be non-empty.")
        if not np.isfinite(self.radius_m) or self.radius_m <= 0:
            raise ValidationError("radius_m must be finite and positive.")
        if not np.isfinite(self.omega0_rad_s) or self.omega0_rad_s <= 0:
            raise ValidationError("omega0_rad_s must be finite and positive.")
        if not _SHA256.match(self.source_checksum):
            raise ValidationError("source_checksum must be a lowercase SHA-256 hexadecimal digest.")
        for name, value in (
            ("source_member_sha256", self.source_member_sha256),
            ("conversion_manifest_sha256", self.conversion_manifest_sha256),
            ("record_payload_sha256", self.record_payload_sha256),
        ):
            if value is not None and not _SHA256.match(value):
                raise ValidationError(f"{name} must be a lowercase SHA-256 hexadecimal digest.")
        if self.converter_commit_sha is not None and not _GIT_SHA.match(self.converter_commit_sha):
            raise ValidationError("converter_commit_sha must be a full lowercase 40-character Git SHA.")
        if self.coordinate_convention != "outward_normal_positive":
            raise ValidationError("coordinate_convention must be outward_normal_positive.")
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
    structural_model_id: str = "auto"
    load_distribution: str = "uniform_apical"
    lateral_support: str = "periodic_monolayer"
    membrane_cortex_coupling: str = "perfectly_bonded"
    nuclear_representation: str = "homogeneous_cell_body"
    prestress_state: str = "zero"
    rheology_mode: str = "elastic"
    load_case: str = "wss_plus_lamb_signed"
    harmonic_control: str = "full_waveform"
    localized_sigma_x_fraction: float = 0.15
    localized_sigma_z_fraction: float = 0.15
    conservation_tolerance_relative: float = 1e-12
    numerical_tolerance_relative: float = 1e-10
    protocol_version: str = "2.0.0"
    parameter_freeze_version: str = "2.0.0"
    claim_bearing: bool = False

    @property
    def resolved_structural_model_id(self) -> str:
        if self.structural_model_id != "auto":
            return self.structural_model_id
        return ":".join(
            (
                self.load_distribution,
                self.lateral_support,
                self.membrane_cortex_coupling,
                self.nuclear_representation,
                self.prestress_state,
            )
        )

    def validate(self) -> None:
        if self.nx < 6 or self.nz < 6:
            raise ValidationError("nx and nz must each be at least six.")
        if self.solver_id not in {"periodic_spectral_2d", "bounded_fd_2d", "lumped_0d"}:
            raise ValidationError(f"Unsupported solver_id: {self.solver_id}")
        if self.load_distribution not in LOAD_DISTRIBUTIONS:
            raise ValidationError(f"Unsupported load_distribution: {self.load_distribution}")
        if self.lateral_support not in LATERAL_SUPPORTS:
            raise ValidationError(f"Unsupported lateral_support: {self.lateral_support}")
        if self.membrane_cortex_coupling not in MEMBRANE_CORTEX_COUPLINGS:
            raise ValidationError(
                f"Unsupported membrane_cortex_coupling: {self.membrane_cortex_coupling}"
            )
        if self.nuclear_representation not in NUCLEAR_REPRESENTATIONS:
            raise ValidationError(f"Unsupported nuclear_representation: {self.nuclear_representation}")
        if self.prestress_state != "zero":
            raise ValidationError(
                "Only zero prestress is frozen. Nonzero prestress requires a new sourced freeze version."
            )
        if self.rheology_mode not in RHEOLOGY_MODES:
            raise ValidationError(f"Unsupported rheology_mode: {self.rheology_mode}")
        if self.load_case not in LOAD_CASES:
            raise ValidationError(f"Unsupported load_case: {self.load_case}")
        if self.harmonic_control not in HARMONIC_CONTROLS:
            raise ValidationError(f"Unsupported harmonic_control: {self.harmonic_control}")
        if self.lateral_support == "periodic_monolayer" and self.solver_id == "bounded_fd_2d":
            raise ValidationError("periodic_monolayer requires periodic_spectral_2d.")
        if self.lateral_support != "periodic_monolayer" and self.solver_id == "periodic_spectral_2d":
            raise ValidationError("bounded lateral supports require bounded_fd_2d.")
        if self.solver_id == "lumped_0d" and self.claim_bearing:
            raise ValidationError("lumped_0d is an analytical baseline and cannot be claim-bearing.")
        if not (0 < self.localized_sigma_x_fraction <= 0.5):
            raise ValidationError("localized_sigma_x_fraction must be in (0, 0.5].")
        if not (0 < self.localized_sigma_z_fraction <= 0.5):
            raise ValidationError("localized_sigma_z_fraction must be in (0, 0.5].")
        if self.conservation_tolerance_relative <= 0 or self.numerical_tolerance_relative <= 0:
            raise ValidationError("numerical tolerances must be positive.")
        if self.protocol_version != "2.0.0" or self.parameter_freeze_version != "2.0.0":
            raise ValidationError("Protocol and parameter freeze versions must be 2.0.0.")
