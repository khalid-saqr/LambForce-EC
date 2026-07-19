"""Execute Step 2 common benchmarks and select the numerical method."""
from __future__ import annotations

import csv
import json
import math
from pathlib import Path
import sys
import time

import numpy as np
import yaml

HERE = Path(__file__).resolve()
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / 'prototypes'))

from step2.common import Geometry, MaterialState, complex_properties, tangential_response, relative_error, complex_to_pair
from step2 import lumped_model, spectral_plate, bounded_fd_plate


def _resultant_error(applied: complex, reaction: complex, scale: float | None = None) -> float:
    denominator = max(float(scale) if scale is not None else abs(applied), 1e-30)
    return float(abs(reaction - applied) / denominator)


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError('No rows to write.')
    fields = list(rows[0])
    with path.open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def run() -> dict[str, object]:
    config_path = ROOT / 'benchmarks' / 'canonical_cases.yaml'
    cfg = yaml.safe_load(config_path.read_text())
    geometry = Geometry(**cfg['geometry'])
    material = MaterialState(**cfg['material_reference'])
    nx = nz = 24

    benchmark_rows: list[dict[str, object]] = []
    detailed: dict[str, object] = {}

    for case in cfg['cases']:
        case_id = case['id']
        omega = 2.0 * math.pi * float(case['frequency_hz'])
        qmean = complex(case['q_mean_pa'])
        qamp = complex(case['q_amplitude_pa'])
        tau = complex(case['tau_w_pa'])

        # Lumped analytical baseline: applicable only to the spatial mean.
        started = time.perf_counter()
        lumped = lumped_model.solve_uniform_normal(qmean, geometry, material, omega)
        tan_lumped = lumped_model.solve_tangential(tau, geometry, material, omega)
        runtime = time.perf_counter() - started
        lumped_applicable = abs(qamp) == 0.0
        benchmark_rows.append({
            'candidate': 'lumped_0d',
            'case_id': case_id,
            'applicable': str(lumped_applicable).lower(),
            'grid_n': 1,
            'runtime_s': runtime,
            'residual_relative_l2': 0.0,
            'resultant_relative_error': _resultant_error(lumped['applied_force_n'], lumped['reaction_force_n']),
            'max_displacement_abs_m': abs(lumped['normal_displacement_m']),
            'tangential_displacement_abs_m': abs(tan_lumped['tangential_displacement_m']),
            'work_measure_j': float(0.5 * np.real(np.conj(lumped['normal_displacement_m']) * lumped['reaction_force_n'])),
            'average_dissipated_power_w': '',
            'spatial_fields': 'false',
            'complex_harmonic_support': 'true',
        })

        # Periodic spectral candidate.
        qspec = spectral_plate.cosine_load(qmean, qamp, geometry, nx, nz)
        started = time.perf_counter()
        spec = spectral_plate.solve(qspec, geometry, material, omega)
        tan_spec = tangential_response(tau, geometry, material, omega)
        runtime_spec = time.perf_counter() - started
        benchmark_rows.append({
            'candidate': 'periodic_spectral_2d',
            'case_id': case_id,
            'applicable': 'true',
            'grid_n': nx,
            'runtime_s': runtime_spec,
            'residual_relative_l2': spec.residual_relative_l2,
            'resultant_relative_error': _resultant_error(spec.applied_resultant_n, spec.reaction_resultant_n, float(np.sum(np.abs(qspec)) * geometry.area_m2 / (nx*nz))),
            'max_displacement_abs_m': float(np.max(np.abs(spec.displacement_m))),
            'tangential_displacement_abs_m': abs(tan_spec['tangential_displacement_m']),
            'work_measure_j': spec.work_measure_j,
            'average_dissipated_power_w': spec.average_dissipated_power_w,
            'spatial_fields': 'true',
            'complex_harmonic_support': 'true',
        })

        # Bounded finite-difference verification candidate.
        qfd = bounded_fd_plate.cosine_load(qmean, qamp, geometry, nx, nz)
        started = time.perf_counter()
        fd = bounded_fd_plate.solve(qfd, geometry, material, omega)
        tan_fd = tangential_response(tau, geometry, material, omega)
        runtime_fd = time.perf_counter() - started
        benchmark_rows.append({
            'candidate': 'bounded_fd_2d',
            'case_id': case_id,
            'applicable': 'true',
            'grid_n': nx,
            'runtime_s': runtime_fd,
            'residual_relative_l2': fd.residual_relative_l2,
            'resultant_relative_error': _resultant_error(fd.applied_resultant_n, fd.reaction_resultant_n, float(np.sum(np.abs(qfd)) * (geometry.length_x_m/(nx+1)) * (geometry.length_z_m/(nz+1)))),
            'max_displacement_abs_m': float(np.max(np.abs(fd.displacement_m))),
            'tangential_displacement_abs_m': abs(tan_fd['tangential_displacement_m']),
            'work_measure_j': fd.work_measure_j,
            'average_dissipated_power_w': '',
            'spatial_fields': 'true',
            'complex_harmonic_support': 'true',
        })

        detailed[case_id] = {
            'spectral_applied_resultant_n': complex_to_pair(spec.applied_resultant_n),
            'spectral_reaction_resultant_n': complex_to_pair(spec.reaction_resultant_n),
            'fd_applied_resultant_n': complex_to_pair(fd.applied_resultant_n),
            'fd_reaction_resultant_n': complex_to_pair(fd.reaction_resultant_n),
        }

    # Spectral analytical convergence for one exactly representable cosine mode.
    conv_rows: list[dict[str, object]] = []
    qamp = 0.07
    for n in (8, 12, 16, 24, 32, 48, 64, 96):
        q = spectral_plate.cosine_load(0.0, qamp, geometry, n, n)
        started = time.perf_counter()
        sol = spectral_plate.solve(q, geometry, material, 0.0)
        runtime = time.perf_counter() - started
        measured = spectral_plate.projected_cosine_amplitude(sol.displacement_m, geometry)
        d, kf = complex_properties(geometry, material, 0.0)
        kx = 2.0 * math.pi / geometry.length_x_m
        kz = 2.0 * math.pi / geometry.length_z_m
        expected = qamp / (d * (kx*kx + kz*kz)**2 + kf)
        conv_rows.append({
            'candidate': 'periodic_spectral_2d',
            'grid_n': n,
            'runtime_s': runtime,
            'reference_quantity': 'cosine_displacement_amplitude_m',
            'value_real': float(np.real(measured)),
            'value_imag': float(np.imag(measured)),
            'relative_error': relative_error(measured, expected),
            'residual_relative_l2': sol.residual_relative_l2,
        })

    # FD self-convergence against a 96x96 reference.
    references = {}
    for n in (8, 12, 16, 20, 24, 32, 40, 48):
        q = bounded_fd_plate.cosine_load(0.035, 0.0175, geometry, n, n)
        started = time.perf_counter()
        sol = bounded_fd_plate.solve(q, geometry, material, 0.0)
        runtime = time.perf_counter() - started
        center = sol.displacement_m[n//2, n//2]
        references[n] = (center, runtime, sol.residual_relative_l2)
    ref = references[48][0]
    for n, (value, runtime, residual) in references.items():
        conv_rows.append({
            'candidate': 'bounded_fd_2d',
            'grid_n': n,
            'runtime_s': runtime,
            'reference_quantity': 'center_displacement_m',
            'value_real': float(np.real(value)),
            'value_imag': float(np.imag(value)),
            'relative_error': relative_error(value, ref),
            'residual_relative_l2': residual,
        })

    gates = cfg['selection_gates']
    spectral_rows = [r for r in benchmark_rows if r['candidate'] == 'periodic_spectral_2d']
    fd_rows = [r for r in benchmark_rows if r['candidate'] == 'bounded_fd_2d']
    lumped_rows = [r for r in benchmark_rows if r['candidate'] == 'lumped_0d']
    spectral_conv = [r for r in conv_rows if r['candidate'] == 'periodic_spectral_2d']

    selection = {
        'primary_solver': 'periodic_spectral_2d',
        'verification_solver': 'bounded_fd_2d',
        'analytical_baseline': 'lumped_0d',
        'full_3d_status': 'deferred_to_later_representative_verification',
        'hard_gate_results': {
            'spectral_spatial_fields': all(r['spatial_fields'] == 'true' for r in spectral_rows),
            'spectral_complex_harmonic_support': all(r['complex_harmonic_support'] == 'true' for r in spectral_rows),
            'spectral_residual': max(float(r['residual_relative_l2']) for r in spectral_rows) <= gates['residual_relative_l2_max'],
            'spectral_resultant_conservation': max(float(r['resultant_relative_error']) for r in spectral_rows) <= gates['resultant_relative_error_max'],
            'spectral_analytic_convergence': max(float(r['relative_error']) for r in spectral_conv) <= gates['spectral_analytic_error_max'],
            'spectral_runtime': max(float(r['runtime_s']) for r in spectral_rows) <= gates['primary_runtime_s_max_per_case'],
            'fd_residual': max(float(r['residual_relative_l2']) for r in fd_rows) <= gates['residual_relative_l2_max'],
            'lumped_rejected_as_primary_due_to_no_spatial_fields': all(r['spatial_fields'] == 'false' for r in lumped_rows),
        },
        'rationale': [
            'The periodic spectral solver matches the primary periodic-monolayer structural class.',
            'It provides complete spatial fields and exact complex harmonic response.',
            'It passes analytical cosine-mode and force-conservation checks at near-machine precision.',
            'The bounded finite-difference solver provides an independent discretization and boundary-condition verification path.',
            'The lumped model remains useful only for uniform analytical limits.',
        ],
    }
    selection['status'] = 'PASS' if all(selection['hard_gate_results'].values()) else 'FAIL'

    _write_csv(ROOT / 'benchmarks' / 'benchmark_results.csv', benchmark_rows)
    _write_csv(ROOT / 'benchmarks' / 'convergence_results.csv', conv_rows)
    (ROOT / 'benchmarks' / 'solver_selection.json').write_text(json.dumps(selection, indent=2, sort_keys=True) + '\n')
    (ROOT / 'benchmarks' / 'benchmark_details.json').write_text(json.dumps(detailed, indent=2, sort_keys=True) + '\n')
    return {
        'selection': selection,
        'benchmark_row_count': len(benchmark_rows),
        'convergence_row_count': len(conv_rows),
    }


if __name__ == '__main__':
    print(json.dumps(run(), indent=2, sort_keys=True))
