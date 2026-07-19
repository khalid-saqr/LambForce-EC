from dataclasses import replace
import numpy as np

from lambforce_ec.models import Geometry, MaterialState, RunConfig
from lambforce_ec.structural import foundation_stiffness_map, structural_ensemble
from lambforce_ec.synthetic import make_synthetic_artery
from lambforce_ec.workflow import run_case


def test_structural_ensemble_is_full_product():
    ensemble = structural_ensemble()
    assert len(ensemble) == 36
    assert len({tuple(sorted(row.items())) for row in ensemble}) == 36
    assert {row["lateral_support"] for row in ensemble} == {
        "periodic_monolayer", "compliant_edge", "clamped_edge"
    }


def test_stiff_nuclear_region_changes_foundation_map():
    geometry = Geometry()
    material = MaterialState()
    config = RunConfig(nx=16, nz=16, nuclear_representation="stiff_nuclear_region")
    kmap = foundation_stiffness_map(geometry, material, config, 0.0)
    assert np.max(kmap.real) > np.min(kmap.real)
    assert np.count_nonzero(kmap.real == np.max(kmap.real)) > 0


def test_all_lateral_supports_execute():
    record = make_synthetic_artery(n_radial=12, n_time=16)
    periodic = RunConfig(nx=8, nz=8, solver_id="periodic_spectral_2d", lateral_support="periodic_monolayer")
    compliant = replace(periodic, solver_id="bounded_fd_2d", lateral_support="compliant_edge")
    clamped = replace(periodic, solver_id="bounded_fd_2d", lateral_support="clamped_edge")
    results = [run_case(record, config=cfg) for cfg in (periodic, compliant, clamped)]
    for result in results:
        assert result.arrays["displacement_normal_cell_m"].shape == (16, 8, 8)
        assert result.metadata["normal_residual_relative_l2_max"] < 1e-8
    assert not np.allclose(
        results[1].arrays["displacement_normal_cell_m"],
        results[2].arrays["displacement_normal_cell_m"],
        rtol=1e-8,
        atol=1e-30,
    )
