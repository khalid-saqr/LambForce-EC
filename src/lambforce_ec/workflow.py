from __future__ import annotations

from dataclasses import dataclass, replace
from datetime import datetime, timezone
from typing import Any, Mapping
import numpy as np

from .analysis import field_summary, incremental_field, norm_ratio, principal_max
from .constitutive import shear_modulus, sls_complex_modulus
from .harmonics import apply_harmonic_control, forward_real_series, reconstruct_coefficients
from .loads import extract_scalar_loads, select_load_case, spatial_kernel
from .models import ArteryRecord, Geometry, MaterialState, RunConfig
from .protocol import assert_claim_bearing_source
from .solvers import BoundedFDPlateSolver, LumpedFoundationSolver, SpectralPlateSolver
from .structural import foundation_stiffness_map, validate_correlated_glycocalyx
from .validation import assert_conserved, checksum_arrays


@dataclass
class SimulationResult:
    arrays: dict[str, np.ndarray]
    metadata: dict[str, Any]
    summaries: dict[str, dict[str, Any]]


def _solution_values(solution: Any) -> dict[str, np.ndarray]:
    return {
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


def _solve_harmonics(
    q_time_pa: np.ndarray,
    tau_time_pa: np.ndarray,
    time_s: np.ndarray,
    kernel: np.ndarray,
    geometry: Geometry,
    material: MaterialState,
    config: RunConfig,
) -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    active_material = material.for_rheology(config.rheology_mode)
    q_series = apply_harmonic_control(forward_real_series(q_time_pa, time_s), config.harmonic_control)
    tau_series = apply_harmonic_control(
        forward_real_series(tau_time_pa, time_s), config.harmonic_control
    )
    solver: SpectralPlateSolver | BoundedFDPlateSolver
    if config.solver_id == "periodic_spectral_2d":
        solver = SpectralPlateSolver()
    elif config.solver_id == "bounded_fd_2d":
        solver = BoundedFDPlateSolver()
    else:
        raise ValueError("_solve_harmonics supports only spatial solvers.")

    nh = q_series.coefficients.shape[0]
    shape = (nh, config.nx, config.nz)
    names = list(
        _solution_values(
            type(
                "Placeholder",
                (),
                {
                    "load_pa": 0,
                    "displacement_cell_m": 0,
                    "displacement_apical_top_m": 0,
                    "glycocalyx_strain": 0,
                    "glycocalyx_reaction_pa": 0,
                    "curvature_x_m_inv": 0,
                    "curvature_z_m_inv": 0,
                    "curvature_xz_m_inv": 0,
                    "strain_x": 0,
                    "strain_z": 0,
                    "shear_strain_xz": 0,
                    "tension_x_n_m": 0,
                    "tension_z_n_m": 0,
                    "tension_xz_n_m": 0,
                },
            )()
        )
    )
    coeffs = {name: np.zeros(shape, dtype=complex) for name in names}
    residuals: list[float] = []
    conservation_errors: list[float] = []
    work_by_harmonic: list[float] = []
    dissipated_power_by_harmonic: list[float] = []
    matrix_nnz: list[int] = []

    for h in range(nh):
        omega = float(q_series.omega_rad_s[h])
        qh = q_series.coefficients[h] * kernel
        kmap = (
            None
            if config.nuclear_representation == "homogeneous_cell_body"
            else foundation_stiffness_map(geometry, active_material, config, omega)
        )
        if isinstance(solver, SpectralPlateSolver):
            solution = solver.solve(
                qh,
                geometry,
                active_material,
                omega,
                foundation_stiffness_map_n_m3=kmap,
                numerical_tolerance_relative=config.numerical_tolerance_relative,
            )
        else:
            solution = solver.solve(
                qh,
                geometry,
                active_material,
                omega,
                foundation_stiffness_map_n_m3=kmap,
                boundary_condition=config.lateral_support,
                numerical_tolerance_relative=config.numerical_tolerance_relative,
            )
            matrix_nnz.append(solution.matrix_nnz)
        assert_conserved(
            solution.applied_resultant_n,
            solution.reaction_resultant_n,
            config.conservation_tolerance_relative,
            f"normal harmonic {h}",
        )
        for name, value in _solution_values(solution).items():
            coeffs[name][h] = value
        residuals.append(solution.residual_relative_l2)
        denominator = max(abs(solution.applied_resultant_n), 1e-30)
        conservation_errors.append(
            float(abs(solution.reaction_resultant_n - solution.applied_resultant_n) / denominator)
            if abs(solution.applied_resultant_n) >= 1e-30
            else 0.0
        )
        is_nyquist = q_series.n_time % 2 == 0 and h == nh - 1
        single_sided_weight = 1.0 if h == 0 or is_nyquist else 4.0
        work_by_harmonic.append(single_sided_weight * solution.work_measure_j)
        dissipated_power_by_harmonic.append(
            single_sided_weight * solution.average_dissipated_power_w
        )

    arrays = {name: reconstruct_coefficients(value, q_series.n_time) for name, value in coeffs.items()}

    tangential_displacement_coeff = np.zeros(nh, dtype=complex)
    cytosol_shear_strain_coeff = np.zeros(nh, dtype=complex)
    glycocalyx_shear_strain_coeff = np.zeros(nh, dtype=complex)
    wss_tension_resultant_coeff = np.zeros(nh, dtype=complex)
    for h in range(nh):
        omega = float(tau_series.omega_rad_s[h])
        ecyt = sls_complex_modulus(active_material.cytosol, omega)
        eg = sls_complex_modulus(active_material.glycocalyx, omega)
        gcyt = shear_modulus(ecyt, active_material.poisson_ratio)
        gg = shear_modulus(eg, active_material.poisson_ratio)
        tau_h = tau_series.coefficients[h]
        cytosol_shear_strain_coeff[h] = tau_h / gcyt
        glycocalyx_shear_strain_coeff[h] = tau_h / gg
        tangential_displacement_coeff[h] = (
            tau_h * geometry.cell_height_m / gcyt
            + tau_h * geometry.glycocalyx_thickness_m / gg
        )
        if config.membrane_cortex_coupling == "perfectly_bonded":
            wss_tension_resultant_coeff[h] = tau_h * geometry.cortex_thickness_m
        else:
            wss_tension_resultant_coeff[h] = 0.0

    arrays["wall_shear_stress_pa"] = reconstruct_coefficients(
        tau_series.coefficients, tau_series.n_time
    )
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
    arrays["curvature_change_m_inv"] = np.sqrt(
        arrays["curvature_x_m_inv"] ** 2
        + arrays["curvature_z_m_inv"] ** 2
        + 0.5 * arrays["curvature_xz_m_inv"] ** 2
    )

    ec = active_material.cortex.einf_pa
    ecyt = active_material.cytosol.einf_pa
    eg = active_material.glycocalyx.einf_pa
    nu = active_material.poisson_ratio
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
        "matrix_nnz_max": max(matrix_nnz, default=None),
    }
    return arrays, metadata


def _run_lumped(
    q_time_pa: np.ndarray,
    geometry: Geometry,
    material: MaterialState,
    config: RunConfig,
) -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    active = material.for_rheology(config.rheology_mode)
    q0 = complex(np.mean(q_time_pa))
    solution = LumpedFoundationSolver().solve(q0, geometry, active)
    arrays = {
        "displacement_normal_cell_m": np.asarray([solution.displacement_cell_m.real]),
        "displacement_normal_apical_top_m": np.asarray([solution.displacement_apical_top_m.real]),
    }
    metadata = {
        "solver_id": "lumped_0d",
        "solver_version": LumpedFoundationSolver.solver_version,
        "applied_resultant_n": solution.applied_resultant_n.real,
        "reaction_resultant_n": solution.reaction_resultant_n.real,
    }
    return arrays, metadata


def run_case(
    record: ArteryRecord,
    geometry: Geometry | None = None,
    material: MaterialState | None = None,
    config: RunConfig | None = None,
    glycocalyx_thickness_field_m: np.ndarray | None = None,
    source_registry: Mapping[str, Any] | None = None,
) -> SimulationResult:
    geometry = geometry or Geometry()
    material = material or MaterialState()
    config = config or RunConfig()
    record.validate()
    geometry.validate()
    material.validate()
    config.validate()
    validate_correlated_glycocalyx(geometry, material)
    if config.claim_bearing:
        assert_claim_bearing_source(record, source_registry)

    scalar_loads = extract_scalar_loads(record)
    q_time, tau_time = select_load_case(scalar_loads, config.load_case)
    if config.solver_id == "lumped_0d":
        arrays, solver_metadata = _run_lumped(q_time, geometry, material, config)
    else:
        kernel = spatial_kernel(geometry, config, glycocalyx_thickness_field_m)
        arrays, solver_metadata = _solve_harmonics(
            q_time, tau_time, record.time_s, kernel, geometry, material, config
        )
        arrays["spatial_kernel"] = kernel

    arrays["time_s"] = np.asarray(record.time_s, dtype=float)
    arrays["lamb_scalar_selected_pa"] = np.asarray(q_time, dtype=float)
    arrays["wss_scalar_selected_pa"] = np.asarray(tau_time, dtype=float)
    metadata = {
        "artery_id": record.artery_id,
        "artery_name": record.artery_name,
        "protocol_version": config.protocol_version,
        "parameter_freeze_version": config.parameter_freeze_version,
        "parameter_set_id": material.parameter_set_id,
        "structural_model_id": config.resolved_structural_model_id,
        "load_distribution": config.load_distribution,
        "lateral_support": config.lateral_support,
        "membrane_cortex_coupling": config.membrane_cortex_coupling,
        "nuclear_representation": config.nuclear_representation,
        "prestress_state": config.prestress_state,
        "rheology_mode": config.rheology_mode,
        "load_case": config.load_case,
        "harmonic_control": config.harmonic_control,
        "claim_bearing": config.claim_bearing,
        "mesh_or_basis_resolution": {"nx": config.nx, "nz": config.nz},
        "coordinate_system": record.coordinate_convention,
        "source_identifier": record.source_identifier,
        "source_version": record.source_version,
        "source_archive_checksum": record.source_checksum,
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
                "glycocalyx_state_id": material.glycocalyx_state_id,
            },
            "config": config.__dict__,
        },
    )
    summary_names = [
        name
        for name in (
            "membrane_tension_max_principal_n_m",
            "displacement_normal_cell_m",
            "displacement_normal_apical_top_m",
            "displacement_tangential_m",
            "strain_max_principal",
            "curvature_change_m_inv",
            "glycocalyx_strain_normal",
            "glycocalyx_reaction_stress_pa",
            "strain_energy_density_j_m3",
            "tension_loading_rate_n_m_s",
        )
        if name in arrays
    ]
    normal_reference = arrays["lamb_scalar_selected_pa"]
    summaries: dict[str, dict[str, Any]] = {}
    for name in summary_names:
        reference = (
            arrays["wss_scalar_selected_pa"]
            if name in {"displacement_tangential_m"}
            else normal_reference
        )
        summaries[name] = field_summary(arrays[name], record.time_s, reference)
    return SimulationResult(arrays=arrays, metadata=metadata, summaries=summaries)


def run_registered_comparison(
    record: ArteryRecord,
    geometry: Geometry | None = None,
    material: MaterialState | None = None,
    config: RunConfig | None = None,
    source_registry: Mapping[str, Any] | None = None,
) -> dict[str, SimulationResult | np.ndarray | float]:
    geometry = geometry or Geometry()
    material = material or MaterialState()
    base = config or RunConfig()

    def execute(load_case: str) -> SimulationResult:
        return run_case(
            record,
            geometry,
            material,
            replace(base, load_case=load_case),
            source_registry=source_registry,
        )

    wss = execute("wss_only")
    combined = execute("wss_plus_lamb_signed")
    isotropic = execute("isotropic_lamb")
    anisotropy_increment = execute("anisotropy_increment")
    field = "membrane_tension_max_principal_n_m"
    delta_total = incremental_field(combined.arrays[field], wss.arrays[field])
    delta_iso = incremental_field(isotropic.arrays[field], wss.arrays[field])
    delta_aniso = incremental_field(anisotropy_increment.arrays[field], wss.arrays[field])
    return {
        "wss_only": wss,
        "combined": combined,
        "isotropic_control": isotropic,
        "anisotropy_increment": anisotropy_increment,
        "membrane_tension_increment_lamb_n_m": delta_total,
        "membrane_tension_increment_isotropic_n_m": delta_iso,
        "membrane_tension_increment_anisotropy_n_m": delta_aniso,
        "incremental_lamb_wss_ratio": norm_ratio(delta_total, wss.arrays[field]),
        "anisotropy_isotropic_ratio": norm_ratio(delta_aniso, delta_iso),
    }


def required_run_matrix() -> list[dict[str, str]]:
    load_cases = [
        "unloaded",
        "wss_only",
        "lamb_signed_only",
        "wss_plus_lamb_signed",
        "exposure_diagnostic",
        "isotropic_lamb",
        "anisotropy_increment",
        "inward_only",
        "outward_only",
        "zero_normal_load",
    ]
    harmonic_controls = ["fundamental_only", "harmonics_le_2", "full_waveform"]
    rows = [{"load_case": case, "harmonic_control": "full_waveform"} for case in load_cases]
    rows.extend(
        {"load_case": "wss_plus_lamb_signed", "harmonic_control": control}
        for control in harmonic_controls[:-1]
    )
    return rows
