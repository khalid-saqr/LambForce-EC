from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
from typing import Any
import numpy as np

from .analysis import field_summary, incremental_field, norm_ratio, principal_max
from .constitutive import shear_modulus, sls_complex_modulus
from .harmonics import apply_harmonic_control, forward_real_series, reconstruct_coefficients
from .loads import extract_scalar_loads, select_load_case, spatial_kernel
from .models import ArteryRecord, Geometry, MaterialState, RunConfig
from .solvers import BoundedFDPlateSolver, LumpedFoundationSolver, SpectralPlateSolver
from .validation import assert_conserved, checksum_arrays


@dataclass
class SimulationResult:
    arrays: dict[str, np.ndarray]
    metadata: dict[str, Any]
    summaries: dict[str, dict[str, Any]]


def _solve_periodic_harmonics(
    q_time_pa: np.ndarray,
    tau_time_pa: np.ndarray,
    time_s: np.ndarray,
    kernel: np.ndarray,
    geometry: Geometry,
    material: MaterialState,
    config: RunConfig,
) -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    q_series = apply_harmonic_control(forward_real_series(q_time_pa, time_s), config.harmonic_control)
    tau_series = apply_harmonic_control(forward_real_series(tau_time_pa, time_s), config.harmonic_control)
    solver = SpectralPlateSolver()
    nh = q_series.coefficients.shape[0]
    shape = (nh, config.nx, config.nz)
    names = [
        "load_normal_pa", "displacement_normal_cell_m", "displacement_normal_apical_top_m",
        "glycocalyx_strain_normal", "glycocalyx_reaction_stress_pa", "curvature_x_m_inv",
        "curvature_z_m_inv", "curvature_xz_m_inv", "strain_x", "strain_z",
        "shear_strain_xz_bending", "tension_x_n_m", "tension_z_n_m", "tension_xz_bending_n_m",
    ]
    coeffs = {name: np.zeros(shape, dtype=complex) for name in names}
    residuals: list[float] = []
    conservation_errors: list[float] = []
    work_by_harmonic: list[float] = []
    dissipated_power_by_harmonic: list[float] = []

    for h in range(nh):
        qh = q_series.coefficients[h] * kernel
        solution = solver.solve(qh, geometry, material, float(q_series.omega_rad_s[h]))
        assert_conserved(
            solution.applied_resultant_n,
            solution.reaction_resultant_n,
            config.conservation_tolerance_relative,
            f"normal harmonic {h}",
        )
        values = {
            "load_normal_pa": solution.load_pa,
            "displacement_normal_cell_m": solution.displacement_cell_m,
            "displacement_normal_apical_top_m": solution.displacement_apical_top_m,
            "glycocalyx_strain_normal": solution.glycocalyx_strain,
            "glycocalyx_reaction_stress_pa": solution.glycocalyx_reaction_pa,
            "curvature_x_m_inv": solution.curvature_x_m_inv,
            "curvature_z_m_inv": solution.curvature_z_m_inv,
            "curvature_xz_m_inv": solution.curvature_xz_m_inv,
            "strain_x": solution.strain_x,
            "strain_z": solution.strain_z,
            "shear_strain_xz_bending": solution.shear_strain_xz,
            "tension_x_n_m": solution.tension_x_n_m,
            "tension_z_n_m": solution.tension_z_n_m,
            "tension_xz_bending_n_m": solution.tension_xz_n_m,
        }
        for name, value in values.items():
            coeffs[name][h] = value
        residuals.append(solution.residual_relative_l2)
        denominator = max(abs(solution.applied_resultant_n), 1e-30)
        conservation_errors.append(float(abs(solution.reaction_resultant_n - solution.applied_resultant_n) / denominator))
        is_nyquist = q_series.n_time % 2 == 0 and h == nh - 1
        single_sided_weight = 1.0 if h == 0 or is_nyquist else 4.0
        work_by_harmonic.append(single_sided_weight * solution.work_measure_j)
        dissipated_power_by_harmonic.append(
            single_sided_weight * solution.average_dissipated_power_w
        )

    arrays = {name: reconstruct_coefficients(value, q_series.n_time) for name, value in coeffs.items()}

    # Homogeneous tangential series in the registered reduced shear limit.
    tangential_displacement_coeff = np.zeros(nh, dtype=complex)
    cytosol_shear_strain_coeff = np.zeros(nh, dtype=complex)
    glycocalyx_shear_strain_coeff = np.zeros(nh, dtype=complex)
    wss_tension_resultant_coeff = np.zeros(nh, dtype=complex)
    for h in range(nh):
        omega = float(tau_series.omega_rad_s[h])
        ecyt = sls_complex_modulus(material.cytosol, omega)
        eg = sls_complex_modulus(material.glycocalyx, omega)
        gcyt = shear_modulus(ecyt, material.poisson_ratio)
        gg = shear_modulus(eg, material.poisson_ratio)
        tau_h = tau_series.coefficients[h]
        cytosol_shear_strain_coeff[h] = tau_h / gcyt
        glycocalyx_shear_strain_coeff[h] = tau_h / gg
        tangential_displacement_coeff[h] = (
            tau_h * geometry.cell_height_m / gcyt
            + tau_h * geometry.glycocalyx_thickness_m / gg
        )
        # Direct traction resultant across the apical cortex thickness; no fitted gain.
        wss_tension_resultant_coeff[h] = tau_h * geometry.cortex_thickness_m

    arrays["wall_shear_stress_pa"] = reconstruct_coefficients(tau_series.coefficients, tau_series.n_time)
    arrays["displacement_tangential_m"] = reconstruct_coefficients(
        tangential_displacement_coeff, tau_series.n_time
    )
    arrays["cytosol_shear_strain"] = reconstruct_coefficients(
        cytosol_shear_strain_coeff, tau_series.n_time
    )
    arrays["glycocalyx_shear_strain"] = reconstruct_coefficients(
        glycocalyx_shear_strain_coeff, tau_series.n_time
    )
    wss_tension = reconstruct_coefficients(wss_tension_resultant_coeff, tau_series.n_time)
    arrays["tension_xz_wss_n_m"] = wss_tension[:, None, None] * np.ones(
        (1, config.nx, config.nz)
    )
    arrays["tension_xz_total_n_m"] = (
        arrays["tension_xz_bending_n_m"] + arrays["tension_xz_wss_n_m"]
    )
    arrays["membrane_tension_max_principal_n_m"] = principal_max(
        arrays["tension_x_n_m"], arrays["tension_z_n_m"], arrays["tension_xz_total_n_m"]
    )
    arrays["strain_max_principal"] = principal_max(
        arrays["strain_x"],
        arrays["strain_z"],
        arrays["shear_strain_xz_bending"] / 2,
    )

    # Positive reduced-model elastic-energy density using relaxed moduli.
    ec = material.cortex.einf_pa
    ecyt = material.cytosol.einf_pa
    eg = material.glycocalyx.einf_pa
    nu = material.poisson_ratio
    epsx = arrays["strain_x"]
    epsz = arrays["strain_z"]
    gamma = arrays["shear_strain_xz_bending"]
    cortex_energy = 0.5 * ec / (1 - nu**2) * (epsx**2 + epsz**2 + 2 * nu * epsx * epsz)
    cortex_energy += 0.5 * ec / (2 * (1 + nu)) * gamma**2
    foundation_energy = 0.5 * ecyt * (
        arrays["displacement_normal_cell_m"] / geometry.cell_height_m
    ) ** 2
    glyco_energy = 0.5 * eg * arrays["glycocalyx_strain_normal"] ** 2
    arrays["strain_energy_density_j_m3"] = cortex_energy + foundation_energy + glyco_energy

    dt = float(time_s[1] - time_s[0])
    tension = arrays["membrane_tension_max_principal_n_m"]
    arrays["tension_loading_rate_n_m_s"] = np.gradient(tension, dt, axis=0, edge_order=2)
    metadata = {
        "solver_id": solver.solver_id,
        "solver_version": solver.solver_version,
        "harmonic_omega_rad_s": q_series.omega_rad_s.tolist(),
        "normal_residual_relative_l2_max": float(max(residuals, default=0.0)),
        "normal_resultant_conservation_relative_max": float(max(conservation_errors, default=0.0)),
        "work_measure_j_by_harmonic": work_by_harmonic,
        "average_dissipated_power_w_by_harmonic": dissipated_power_by_harmonic,
        "work_per_cycle_j": float(sum(dissipated_power_by_harmonic) * (time_s.size * dt)),
    }
    return arrays, metadata


def _run_reduced_verification(
    q_time_pa: np.ndarray,
    geometry: Geometry,
    material: MaterialState,
    config: RunConfig,
) -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    q0 = complex(np.mean(q_time_pa))
    if config.solver_id == "lumped_0d":
        solution = LumpedFoundationSolver().solve(q0, geometry, material)
        arrays = {
            "displacement_normal_cell_m": np.asarray([solution.displacement_cell_m.real]),
            "displacement_normal_apical_top_m": np.asarray([solution.displacement_apical_top_m.real]),
        }
        metadata = {
            "solver_id": "lumped_0d",
            "solver_version": "1.0.0",
            "applied_resultant_n": solution.applied_resultant_n.real,
            "reaction_resultant_n": solution.reaction_resultant_n.real,
        }
        return arrays, metadata
    kernel = spatial_kernel(geometry, config)
    solution = BoundedFDPlateSolver().solve(q0 * kernel, geometry, material)
    arrays = {"displacement_normal_cell_m": solution.displacement_cell_m.real[None, :, :]}
    metadata = {
        "solver_id": "bounded_fd_2d",
        "solver_version": "1.0.0",
        "residual_relative_l2": solution.residual_relative_l2,
        "applied_resultant_n": solution.applied_resultant_n.real,
        "reaction_resultant_n": solution.reaction_resultant_n.real,
        "matrix_nnz": solution.matrix_nnz,
    }
    return arrays, metadata


def run_case(
    record: ArteryRecord,
    geometry: Geometry | None = None,
    material: MaterialState | None = None,
    config: RunConfig | None = None,
    glycocalyx_thickness_field_m: np.ndarray | None = None,
) -> SimulationResult:
    geometry = geometry or Geometry()
    material = material or MaterialState()
    config = config or RunConfig()
    record.validate()
    geometry.validate()
    material.validate()
    config.validate()

    scalar_loads = extract_scalar_loads(record)
    q_time, tau_time = select_load_case(scalar_loads, config.load_case)
    if config.solver_id == "periodic_spectral_2d":
        kernel = spatial_kernel(geometry, config, glycocalyx_thickness_field_m)
        arrays, solver_metadata = _solve_periodic_harmonics(
            q_time, tau_time, record.time_s, kernel, geometry, material, config
        )
        arrays["spatial_kernel"] = kernel
    else:
        arrays, solver_metadata = _run_reduced_verification(q_time, geometry, material, config)

    arrays["time_s"] = np.asarray(record.time_s, dtype=float)
    arrays["lamb_scalar_selected_pa"] = np.asarray(q_time, dtype=float)
    arrays["wss_scalar_selected_pa"] = np.asarray(tau_time, dtype=float)
    metadata = {
        "artery_id": record.artery_id,
        "artery_name": record.artery_name,
        "protocol_version": config.protocol_version,
        "parameter_freeze_version": config.parameter_freeze_version,
        "parameter_set_id": material.parameter_set_id,
        "structural_model_id": config.structural_model_id,
        "load_distribution": config.load_distribution,
        "load_case": config.load_case,
        "harmonic_control": config.harmonic_control,
        "mesh_or_basis_resolution": {"nx": config.nx, "nz": config.nz},
        "coordinate_system": record.coordinate_convention,
        "source_identifier": record.source_identifier,
        "source_version": record.source_version,
        "source_checksum": record.source_checksum,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        **solver_metadata,
    }
    metadata["configuration_checksum"] = checksum_arrays(
        {"time_s": arrays["time_s"]},
        {
            "geometry": geometry.__dict__,
            "material": {
                "cortex": material.cortex.__dict__,
                "cytosol": material.cytosol.__dict__,
                "glycocalyx": material.glycocalyx.__dict__,
                "nucleus": material.nucleus.__dict__,
                "poisson_ratio": material.poisson_ratio,
                "parameter_set_id": material.parameter_set_id,
            },
            "config": config.__dict__,
        },
    )
    summary_names = [
        name for name in (
            "membrane_tension_max_principal_n_m",
            "displacement_normal_cell_m",
            "displacement_normal_apical_top_m",
            "displacement_tangential_m",
            "strain_max_principal",
            "strain_energy_density_j_m3",
            "tension_loading_rate_n_m_s",
        ) if name in arrays
    ]
    summaries = {name: field_summary(arrays[name], record.time_s) for name in summary_names}
    return SimulationResult(arrays=arrays, metadata=metadata, summaries=summaries)


def run_registered_comparison(
    record: ArteryRecord,
    geometry: Geometry | None = None,
    material: MaterialState | None = None,
    config: RunConfig | None = None,
) -> dict[str, SimulationResult | np.ndarray | float]:
    geometry = geometry or Geometry()
    material = material or MaterialState()
    base = config or RunConfig()
    wss = run_case(record, geometry, material, RunConfig(**{**base.__dict__, "load_case": "wss_only"}))
    combined = run_case(
        record, geometry, material,
        RunConfig(**{**base.__dict__, "load_case": "wss_plus_lamb_signed"}),
    )
    isotropic = run_case(
        record, geometry, material,
        RunConfig(**{**base.__dict__, "load_case": "isotropic_lamb"}),
    )
    anisotropy_increment = run_case(
        record, geometry, material,
        RunConfig(**{**base.__dict__, "load_case": "anisotropy_increment"}),
    )
    field = "membrane_tension_max_principal_n_m"
    delta = incremental_field(combined.arrays[field], wss.arrays[field])
    iso_delta = incremental_field(isotropic.arrays[field], wss.arrays[field])
    return {
        "wss_only": wss,
        "combined": combined,
        "isotropic_control": isotropic,
        "anisotropy_increment": anisotropy_increment,
        "membrane_tension_increment_lamb_n_m": delta,
        "membrane_tension_increment_isotropic_n_m": iso_delta,
        "incremental_lamb_wss_ratio": norm_ratio(delta, wss.arrays[field]),
        "anisotropy_to_isotropic_increment_ratio": norm_ratio(delta - iso_delta, iso_delta),
    }
