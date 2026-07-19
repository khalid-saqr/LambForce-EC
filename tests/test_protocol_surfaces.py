import csv
from pathlib import Path
import pytest

from lambforce_ec.exceptions import ValidationError
from lambforce_ec.models import Geometry, MaterialState, RunConfig
from lambforce_ec.registry import protocol_surface_registry
from lambforce_ec.structural import validate_correlated_glycocalyx

ROOT = Path(__file__).resolve().parents[1]


def test_protocol_surface_registry():
    registry = protocol_surface_registry()
    assert len(registry["load_case"]) == 10
    assert len(registry["load_distribution"]) == 3
    assert len(registry["lateral_support"]) == 3
    assert len(registry["membrane_cortex_coupling"]) == 2
    assert len(registry["nuclear_representation"]) == 2
    assert registry["prestress_state"] == ("zero",)


def test_protocol_versions_are_locked():
    RunConfig().validate()
    with pytest.raises(ValidationError):
        RunConfig(protocol_version="1.0.0-step1").validate()
    with pytest.raises(ValidationError):
        RunConfig(parameter_freeze_version="2.0.1").validate()


def test_nonzero_prestress_is_rejected():
    with pytest.raises(ValidationError):
        RunConfig(prestress_state="fitted").validate()


def test_geometry_contains_frozen_nuclear_dimensions():
    geometry = Geometry()
    geometry.validate()
    assert geometry.nucleus_axis_x_m == 8e-6
    assert geometry.nucleus_axis_z_m == 6e-6
    assert geometry.nucleus_height_m == 2.5e-6


def test_glycocalyx_correlation_is_enforced():
    validate_correlated_glycocalyx(Geometry(), MaterialState())
    with pytest.raises(ValidationError):
        validate_correlated_glycocalyx(
            Geometry(glycocalyx_thickness_m=1.0e-6), MaterialState()
        )


def test_elastic_mode_removes_relaxation_dependence():
    material = MaterialState().for_rheology("elastic")
    assert material.cortex.ratio_e0_einf == 1
    assert material.cytosol.ratio_e0_einf == 1
    assert material.glycocalyx.ratio_e0_einf == 1
    assert material.nucleus.ratio_e0_einf == 1


def test_parameter_registry_derivation_rules():
    with (ROOT / "registry/parameter_registry.csv").open(newline="", encoding="utf-8") as stream:
        rows = {row["parameter_id"]: row for row in csv.DictReader(stream)}
    assert rows["shear_modulus"]["independent_or_derived"] == "derived"
    assert rows["bulk_modulus"]["independent_or_derived"] == "derived"
    assert rows["shear_modulus"]["transformation_rule"] == "E/[2(1+nu)]"
