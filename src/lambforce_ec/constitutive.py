from __future__ import annotations

import numpy as np
from .models import Geometry, MaterialState, SLSMaterial


def sls_complex_modulus(material: SLSMaterial, omega_rad_s: float) -> complex:
    """Standard-linear-solid modulus for the exp(i*omega*t) convention."""
    material.validate("material")
    if omega_rad_s < 0 or not np.isfinite(omega_rad_s):
        raise ValueError("omega_rad_s must be finite and non-negative.")
    e0 = material.einf_pa * material.ratio_e0_einf
    if omega_rad_s == 0:
        return complex(material.einf_pa)
    z = 1j * omega_rad_s * material.relaxation_time_s
    return material.einf_pa + (e0 - material.einf_pa) * z / (1 + z)


def plate_foundation_properties(
    geometry: Geometry, material: MaterialState, omega_rad_s: float
) -> tuple[complex, complex, complex]:
    geometry.validate()
    material.validate()
    ec = sls_complex_modulus(material.cortex, omega_rad_s)
    ecyt = sls_complex_modulus(material.cytosol, omega_rad_s)
    eg = sls_complex_modulus(material.glycocalyx, omega_rad_s)
    nu = material.poisson_ratio
    bending_rigidity_n_m = ec * geometry.cortex_thickness_m**3 / (12 * (1 - nu**2))
    foundation_stiffness_n_m3 = ecyt / geometry.cell_height_m
    glycocalyx_stiffness_n_m3 = eg / geometry.glycocalyx_thickness_m
    return bending_rigidity_n_m, foundation_stiffness_n_m3, glycocalyx_stiffness_n_m3


def shear_modulus(young_modulus_pa: complex, poisson_ratio: float) -> complex:
    return young_modulus_pa / (2 * (1 + poisson_ratio))
