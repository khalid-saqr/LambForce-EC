# Scientific Contract

## Purpose

`LambForce-EC` will determine whether the anisotropic near-wall Lamb-force field established in the published hydrodynamic study produces a resolved, physically plausible and robust endothelial membrane-loading state while wall shear stress (WSS) is present.

The claim-bearing quantity is the incremental membrane response

```math
\Delta \mathbf{T}^{L}_{m}(\mathbf{x},t)
=
\mathbf{T}^{\mathrm{WSS+Lamb}}_{m}(\mathbf{x},t)
-
\mathbf{T}^{\mathrm{WSS}}_{m}(\mathbf{x},t).
```

## Fixed hydrodynamic ground truth

The following are immutable source definitions from Saqr (2026):

1. anisotropic Womersley governing equations;
2. harmonic velocity reconstruction;
3. vorticity and Lamb-vector definitions;
4. signed radial force-density field;
5. near-wall pillbox integration;
6. published absolute-value exposure integral;
7. WSS;
8. isotropic and anisotropy-specific controls;
9. six arterial waveform definitions;
10. harmonic amplitudes, phases and spectral analysis.

The new model may reproduce these quantities for regression verification, but it must not redefine or tune them.

## New model boundary

The new computational problem starts at the physical interface between the published hydrodynamic fields and endothelial mechanics:

```text
signed near-wall Lamb density + WSS
                 ↓
force-conserving normal and tangential loading
                 ↓
glycocalyx → membrane–cortex → cell body
                 ↓
spatial and temporal endothelial mechanical fields
```

## Input contract

### Required hydrodynamic input

An artery case must provide either:

- source pressure-gradient harmonics and the parameters required to reproduce the anisotropic Womersley solution; or
- precomputed, verified signed Lamb-force density and WSS fields.

Minimum fields:

| Name | Unit | Meaning |
|---|---:|---|
| `artery_id` | 1 | Stable case identifier |
| `radius_m` | m | Vessel radius |
| `omega0_rad_s` | rad s⁻¹ | Fundamental angular frequency |
| `radial_coordinate_m` | m | Near-wall radial coordinate |
| `time_s` | s | Time over one cardiac cycle |
| `lamb_density_signed_N_m3` | N m⁻³ | Signed radial Lamb-force density |
| `wall_shear_stress_Pa` | Pa | Tangential WSS |
| `lamb_density_isotropic_N_m3` | N m⁻³ | Isotropic control |
| `source_identifier` | 1 | DOI, release or dataset identifier |
| `source_checksum` | 1 | SHA-256 checksum |

The six published arteries are mandatory reference cases, not a hard-coded solver limit.

### Required endothelial input

Endothelial geometry, glycocalyx, cortex, cytosol and nucleus values must come from the versioned parameter registry. No claim-bearing parameter may be embedded directly in solver code.

## Output contract

The model must produce complete spatial and temporal fields for:

- membrane tension;
- membrane strain;
- normal and tangential displacement;
- curvature change;
- glycocalyx strain and reaction stress;
- strain-energy density;
- work per cycle;
- loading rate;
- harmonic gain and phase.

Every output must include SI units, coordinate metadata, artery identifier, parameter-set identifier, structural-model identifier, source checksums and protocol version.

## Null and alternative hypotheses

**Null:** the incremental Lamb-induced membrane field is not resolved above combined numerical and registered physiological uncertainty.

**Alternative:** the incremental field is resolved, attributable to the anisotropy-specific input, consistent with documented endothelial mechanical scales and robust across the registered parameter and structural model domain.

## Prohibited mechanisms and parameters

The claim-bearing model must not contain:

- a fitted Lamb-force transfer efficiency;
- an arbitrary equivalent-pressure gain;
- an arbitrary loading area;
- post-result localization tuning;
- outcome-selected prestress;
- artery-specific endothelial material properties in the primary analysis;
- a universal biological threshold;
- Piezo1, calcium or gene-expression parameters.

## Extensibility requirement

The hydrodynamic input layer must accept an arbitrary number of artery cases through a validated schema. The six reference arteries are selected by configuration. Adding an artery must require data and metadata only, not source-code modification.

## Step 1 exit criteria

Step 1 passes only when:

- dimensional closure is documented;
- the loading operator conserves resultant force;
- the minimum mechanics prototype produces stable nonzero normal and tangential responses;
- elastic energy remains non-negative;
- runtime is compatible with subsequent solver-prototype comparison;
- no unresolved claim-bearing gain remains.
