"""Common definitions for Step 2 solver prototypes.

These prototypes are numerical-method selection artifacts. They are not the
claim-bearing LambForce-EC biological solver.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
import numpy as np


@dataclass(frozen=True)
class Geometry:
    length_x_m: float = 36.0e-6
    length_z_m: float = 32.1e-6
    cell_height_m: float = 5.0e-6
    cortex_thickness_m: float = 0.10e-6

    @property
    def area_m2(self) -> float:
        return self.length_x_m * self.length_z_m


@dataclass(frozen=True)
class MaterialState:
    cortex_einf_pa: float = 1_000.0
    cortex_ratio_e0_einf: float = 2.0
    cortex_tau_s: float = 0.05
    cytosol_einf_pa: float = 500.0
    cytosol_ratio_e0_einf: float = 2.0
    cytosol_tau_s: float = 2.0
    poisson_ratio: float = 0.45

    def validate(self) -> None:
        if self.cortex_einf_pa <= 0 or self.cytosol_einf_pa <= 0:
            raise ValueError('Elastic moduli must be positive.')
        if self.cortex_ratio_e0_einf < 1 or self.cytosol_ratio_e0_einf < 1:
            raise ValueError('E0/Einf ratios must be at least one.')
        if self.cortex_tau_s <= 0 or self.cytosol_tau_s <= 0:
            raise ValueError('Relaxation times must be positive.')
        if not (-1.0 < self.poisson_ratio < 0.5):
            raise ValueError('Poisson ratio outside stable isotropic range.')


def sls_modulus(einf_pa: float, ratio_e0_einf: float, tau_s: float, omega_rad_s: float) -> complex:
    """Standard-linear-solid complex modulus under exp(i omega t)."""
    e0 = einf_pa * ratio_e0_einf
    if omega_rad_s == 0.0:
        return complex(einf_pa)
    z = 1j * omega_rad_s * tau_s
    return einf_pa + (e0 - einf_pa) * z / (1.0 + z)


def complex_properties(geometry: Geometry, material: MaterialState, omega_rad_s: float) -> tuple[complex, complex]:
    """Return plate bending rigidity D and foundation stiffness k_f."""
    material.validate()
    ec = sls_modulus(
        material.cortex_einf_pa,
        material.cortex_ratio_e0_einf,
        material.cortex_tau_s,
        omega_rad_s,
    )
    ecyt = sls_modulus(
        material.cytosol_einf_pa,
        material.cytosol_ratio_e0_einf,
        material.cytosol_tau_s,
        omega_rad_s,
    )
    nu = material.poisson_ratio
    d = ec * geometry.cortex_thickness_m**3 / (12.0 * (1.0 - nu**2))
    kf = ecyt / geometry.cell_height_m
    return d, kf


def tangential_response(tau_w_pa: complex, geometry: Geometry, material: MaterialState, omega_rad_s: float) -> dict[str, complex]:
    ecyt = sls_modulus(
        material.cytosol_einf_pa,
        material.cytosol_ratio_e0_einf,
        material.cytosol_tau_s,
        omega_rad_s,
    )
    shear_modulus = ecyt / (2.0 * (1.0 + material.poisson_ratio))
    shear_strain = tau_w_pa / shear_modulus
    return {
        'shear_modulus_pa': shear_modulus,
        'shear_strain': shear_strain,
        'tangential_displacement_m': shear_strain * geometry.cell_height_m,
    }


def relative_error(value: complex | float, reference: complex | float, floor: float = 1e-30) -> float:
    return float(abs(value - reference) / max(abs(reference), floor))


def complex_to_pair(value: complex) -> dict[str, float]:
    return {'real': float(np.real(value)), 'imag': float(np.imag(value))}
