# Load Definitions and Sign Conventions

## Coordinate convention

- $`\mathbf{n}`$: outward unit normal from the endothelial surface into the fluid.
- $`\mathbf{e}_s`$: local streamwise tangential direction.
- Positive $`q_L`$: load in the $`+\mathbf{n}`$ direction.
- Negative $`q_L`$: load in the $`-\mathbf{n}`$ direction.
- WSS is a signed tangential traction along $`\mathbf{e}_s`$.

The coordinate convention must be recorded in every dataset and must be invariant under mesh orientation changes.

## Hydrodynamic load definitions

```math
q_L(t)=\int_{R-\delta_f}^{R}f_r(r,t)\,dr.
```

```math
q_L^{\mathrm{iso}}(t)=\int_{R-\delta_f}^{R}f_r^{\mathrm{iso}}(r,t)\,dr.
```

```math
\Delta q_L^{\mathrm{aniso}}(t)
=
q_L^{\mathrm{total}}(t)-q_L^{\mathrm{iso}}(t).
```

The exposure diagnostic is

```math
q_L^{\mathrm{exposure}}(t)
=
\int_{R-\delta_f}^{R}|f_r(r,t)|\,dr.
```

Exposure is non-directional and cannot be used as the primary mechanical boundary condition.

## Spatial distribution classes

All classes conserve $`A_{\mathrm{cell}}q_L(t)`$.

### `uniform_apical`

```math
q_L(\mathbf{x},t)=q_L(t).
```

### `glycocalyx_resolved`

The local load is computed through the solved glycocalyx mechanics. Numerical integration must recover the same total resultant.

### `localized_bound`

A predeclared normalized spatial kernel $`K(\mathbf{x})`$ is used only as a conservative concentration bound:

```math
q_L(\mathbf{x},t)
=
A_{\mathrm{cell}}q_L(t)K(\mathbf{x}),
\qquad
\int_{A_{\mathrm{cell}}}K(\mathbf{x})\,dA=1.
```

No kernel width may be fitted after results are inspected.

## Simultaneous loading cases

| Case | Tangential WSS | Normal Lamb load |
|---|---|---|
| unloaded | 0 | 0 |
| WSS only | $`\tau_w(t)`$ | 0 |
| Lamb only | 0 | $`q_L(t)`$ |
| combined | $`\tau_w(t)`$ | $`q_L(t)`$ |
| isotropic control | $`\tau_w(t)`$ | $`q_L^{\mathrm{iso}}(t)`$ |
| anisotropy increment | 0 or registered WSS | $`\Delta q_L^{\mathrm{aniso}}(t)`$ |
| exposure diagnostic | 0 | $`q_L^{\mathrm{exposure}}(t)`$ |
| zero-normal control | $`\tau_w(t)`$ | 0 |

## Generic artery input

The solver must accept additional arteries through a validated record containing either source harmonics or precomputed fields. The input adapter must not contain a fixed list of artery names.

Required metadata:

```text
artery_id
artery_name
radius_m
omega0_rad_s
source_identifier
source_version
source_checksum
coordinate_convention
```

Required precomputed arrays:

```text
radial_coordinate_m
time_s
lamb_density_signed_N_m3
wall_shear_stress_Pa
lamb_density_isotropic_N_m3
```

Optional source-solver inputs:

```text
pressure_gradient_complex_harmonics_Pa_m
harmonic_indices
blood_density_kg_m3
axial_kinematic_viscosity_m2_s
anisotropy_beta
anisotropy_gamma
anisotropy_delta
```
