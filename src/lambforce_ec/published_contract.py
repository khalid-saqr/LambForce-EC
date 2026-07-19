from __future__ import annotations

from dataclasses import dataclass
import re
from typing import Any, Mapping

import numpy as np

from .exceptions import ValidationError

REFERENCE_MODES = ("historical_v2", "verified")
REFERENCE_ARTERY_IDS = (
    "aortic_root", "thoracic_aorta", "femoral", "carotid", "iliac", "brachial",
)
_GIT_SHA = re.compile(r"^[0-9a-f]{40}$")
_SHA256 = re.compile(r"^[0-9a-f]{64}$")

EXPECTED_SOURCE = {
    "repository": "khalid-saqr/picoNewton",
    "repository_commit_sha": "4c3c36db0578373cc4e48d9d8c7e8a85944ed1cb",
    "published_notebook_path": "picoNewton_v2.ipynb",
    "published_notebook_blob_sha": "9d61c237cda75df338ce0383038f7765c886f503",
    "paper_doi": "10.1038/s41598-026-47474-x",
}
EXPECTED_FLUID = {
    "density_kg_m3": 1060.0,
    "kinematic_viscosity_m2_s": 3.5e-6,
    "fundamental_frequency_hz": 1.2,
}
EXPECTED_CONTROL_VOLUME = {"reference_area_m2": 1.0e-10, "reference_volume_m3": 1.0e-15}
EXPECTED_ANISOTROPY = {
    "beta": 0.1, "gamma": 0.1, "delta": 1.0,
    "isotropic_beta": 0.0, "isotropic_gamma": 0.0, "isotropic_delta": 1.0,
}
EXPECTED_NUMERICAL = {
    "radial_order": 150, "time_points": 2048, "quadrature_nodes": 256,
    "harmonics_retained": 6, "harmonic_rms_tolerance": 1.0e-3,
}
EXPECTED_ARTERIES = {
    "aortic_root": ("Aortic Root", 0.015, 22.03, 9000.0, (1, .82, .54, .33, .24, .17)),
    "thoracic_aorta": ("Thoracic Aorta", .012, 17.62, 7000.0, (1, .76, .45, .28, .20, .12)),
    "femoral": ("Femoral", .004, 5.87, 6000.0, (1, .58, .10, -.17, .05, .04)),
    "carotid": ("Carotid", .0035, 5.14, 6500.0, (1, .63, .31, .15, .10, .06)),
    "iliac": ("Iliac", .0045, 6.61, 5500.0, (1, .51, .12, -.11, .05, .03)),
    "brachial": ("Brachial", .002, 2.94, 4000.0, (1, .49, .16, -.05, .02, .01)),
}


@dataclass(frozen=True)
class PublishedArteryCase:
    artery_id: str
    artery_name: str
    radius_m: float
    published_womersley_alpha: float
    pressure_gradient_scale_pa_per_m: float
    harmonic_coefficients: tuple[float, ...]


@dataclass(frozen=True)
class ReproductionProfile:
    radial_order: int
    time_points: int
    quadrature_nodes: int


def _require_exact(actual: Mapping[str, Any], expected: Mapping[str, Any], label: str) -> None:
    for key, frozen in expected.items():
        if key not in actual:
            raise ValidationError(f"{label}.{key} is missing.")
        value = actual[key]
        if isinstance(frozen, float):
            if not np.isclose(float(value), frozen, rtol=0, atol=max(1e-18, abs(frozen) * 1e-14)):
                raise ValidationError(f"{label}.{key} differs from the frozen published value.")
        elif value != frozen:
            raise ValidationError(f"{label}.{key} differs from the frozen published value.")


def validate_published_inputs(mapping: Mapping[str, Any]) -> int:
    required = {"registry_version", "source", "fluid", "control_volume", "anisotropy", "numerical", "arteries"}
    if required - set(mapping):
        raise ValidationError("Published-input registry is incomplete.")
    for section, expected in (
        ("source", EXPECTED_SOURCE), ("fluid", EXPECTED_FLUID),
        ("control_volume", EXPECTED_CONTROL_VOLUME), ("anisotropy", EXPECTED_ANISOTROPY),
        ("numerical", EXPECTED_NUMERICAL),
    ):
        if not isinstance(mapping[section], dict):
            raise ValidationError(f"Published-input {section} must be a mapping.")
        _require_exact(mapping[section], expected, section)
    numerical = mapping["numerical"]
    if numerical.get("verified_nonlinear_ordering") != "reconstruct_real_fields_then_multiply":
        raise ValidationError("Verified nonlinear ordering differs from the frozen formulation.")
    if numerical.get("historical_nonlinear_ordering") != "multiply_corresponding_harmonic_fields_then_reconstruct":
        raise ValidationError("Historical nonlinear ordering differs from the regression path.")
    arteries = mapping["arteries"]
    if not isinstance(arteries, list) or len(arteries) != 6:
        raise ValidationError("Published-input registry must contain exactly six arteries.")
    if {str(row.get("artery_id")) for row in arteries if isinstance(row, dict)} != set(EXPECTED_ARTERIES):
        raise ValidationError("Published artery IDs differ from the frozen six.")
    for row in arteries:
        artery_id = str(row["artery_id"])
        name, radius, alpha, pressure, harmonics = EXPECTED_ARTERIES[artery_id]
        if row.get("artery_name") != name:
            raise ValidationError(f"{artery_id}.artery_name differs from the frozen source.")
        numeric = (float(row["radius_m"]), float(row["published_womersley_alpha"]), float(row["pressure_gradient_scale_pa_per_m"]))
        if not np.allclose(numeric, (radius, alpha, pressure), rtol=0, atol=1e-12):
            raise ValidationError(f"{artery_id} scalar inputs differ from the frozen source.")
        coefficients = row.get("harmonic_coefficients")
        if not isinstance(coefficients, list) or not np.allclose(coefficients, harmonics, rtol=0, atol=1e-15):
            raise ValidationError(f"{artery_id}.harmonic_coefficients differ from the frozen source.")
    return 6


def artery_case(row: Mapping[str, Any]) -> PublishedArteryCase:
    return PublishedArteryCase(
        str(row["artery_id"]), str(row["artery_name"]), float(row["radius_m"]),
        float(row["published_womersley_alpha"]), float(row["pressure_gradient_scale_pa_per_m"]),
        tuple(float(value) for value in row["harmonic_coefficients"]),
    )


def reproduction_profile(mapping: Mapping[str, Any], name: str) -> ReproductionProfile:
    if name == "verification":
        return ReproductionProfile(32, 128, 48)
    if name == "publication":
        numerical = mapping["numerical"]
        return ReproductionProfile(
            int(numerical["radial_order"]), int(numerical["time_points"]),
            int(numerical["quadrature_nodes"]),
        )
    raise ValidationError(f"Unsupported reproduction profile: {name}")


def validate_published_source_binding(
    mapping: Mapping[str, Any], input_sha256: str, source_registry: Mapping[str, Any]
) -> None:
    validate_published_inputs(mapping)
    if not _SHA256.match(input_sha256):
        raise ValidationError("Published-input registry checksum is invalid.")
    identifier = f"github:{mapping['source']['repository']}#{mapping['source']['published_notebook_path']}"
    sources = source_registry.get("sources") if isinstance(source_registry, Mapping) else None
    matches = [row for row in (sources or []) if isinstance(row, dict) and row.get("source_identifier") == identifier]
    if len(matches) != 1:
        raise ValidationError("Published source does not resolve uniquely in source_registry.yaml.")
    row = matches[0]
    checks = {
        "source_version": mapping["source"]["published_notebook_blob_sha"],
        "repository_commit_sha": mapping["source"]["repository_commit_sha"],
        "published_notebook_path": mapping["source"]["published_notebook_path"],
        "published_notebook_blob_sha": mapping["source"]["published_notebook_blob_sha"],
        "published_input_registry_sha256": input_sha256,
        "role": "published_hydrodynamic_ground_truth", "claim_bearing": True,
    }
    for key, expected in checks.items():
        if row.get(key) != expected:
            raise ValidationError(f"source_registry.yaml {key} differs from the frozen source.")
