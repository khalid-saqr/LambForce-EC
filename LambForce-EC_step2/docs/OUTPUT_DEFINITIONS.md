# Output Definitions

## Claim-bearing fields

| Output ID | Symbol | SI unit | Definition |
|---|---:|---:|---|
| `membrane_tension_max_principal` | $`T_1(\mathbf{x},t)`$ | N m⁻¹ | Maximum principal membrane–cortex resultant |
| `membrane_tension_increment_lamb` | $`\Delta T_1^L`$ | N m⁻¹ | Combined minus WSS-only tension |
| `strain_max_principal` | $`\varepsilon_1`$ | 1 | Maximum principal strain |
| `displacement_normal` | $`w`$ | m | Wall-normal displacement |
| `displacement_tangential` | $`u_t`$ | m | Tangential displacement |
| `curvature_change` | $`\Delta\kappa`$ | m⁻¹ | Change in apical curvature |
| `glycocalyx_strain` | $`\varepsilon_g`$ | 1 | Glycocalyx compression or extension |
| `glycocalyx_reaction_stress` | $`\sigma_g`$ | Pa | Reaction stress transmitted through the glycocalyx |
| `strain_energy_density` | $`U_v`$ | J m⁻³ | Stored elastic energy density |
| `work_per_cycle` | $`W_{\mathrm{cycle}}`$ | J | Net or dissipated work, defined by model class |
| `tension_loading_rate` | $`\partial T_1/\partial t`$ | N m⁻¹ s⁻¹ | Temporal loading rate |
| `harmonic_gain` | $`G_h`$ | 1 | Output/input amplitude ratio at harmonic $`h`$ |
| `harmonic_phase` | $`\varphi_h`$ | rad | Output–input phase difference |
| `spatial_concentration` | $`C_s`$ | 1 | Spatial peak-to-mean ratio |
| `incremental_lamb_wss_ratio` | $`\mathcal R_L`$ | 1 | Norm of Lamb increment divided by WSS-only norm |

The primary relative metric is

```math
\mathcal R_L
=
\frac{
\left\|
\mathbf T_m^{\mathrm{WSS+Lamb}}
-
\mathbf T_m^{\mathrm{WSS}}
\right\|
}{
\left\|
\mathbf T_m^{\mathrm{WSS}}
\right\|
}.
```

## Required summaries

For every field and artery:

- peak;
- RMS;
- cycle mean;
- 95th percentile;
- peak-to-peak range;
- time of peak;
- spatial location of peak;
- spatial peak-to-mean ratio;
- signed inward/outward asymmetry;
- harmonic magnitude and phase;
- cycle-integrated work or energy.

## Metadata requirements

Every output artifact must include:

```text
artery_id
protocol_version
parameter_freeze_version
parameter_set_id
structural_model_id
solver_id
solver_version
mesh_or_basis_resolution
time_or_harmonic_resolution
coordinate_system
units
source_checksums
configuration_checksum
software_commit
created_utc
```

## Decision relationship

No single output threshold establishes biological meaning. The decision system requires convergent evidence across numerical significance, WSS-present increment, anisotropy attribution, experimental-scale consistency, parameter robustness, structural robustness, cross-artery generality and spectral relevance.
