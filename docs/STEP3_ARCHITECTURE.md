# Step 3 production package

## Status

Step 3 implements the first production computational architecture selected by Step 2. It is a mechanics package, not a claim-bearing result set.

## Immutable boundary

The published anisotropic Womersley solution, signed Lamb-force density, pillbox integration, WSS, isotropic control and six arterial waveforms remain hydrodynamic ground truth. The package starts from verified signed near-wall fields:

```text
verified signed Lamb density + WSS
                ↓
radial integration in SI units
                ↓
force-conserving normal distribution + tangential traction
                ↓
registered glycocalyx / cortex / cell-body mechanics
                ↓
spatial, temporal and harmonic mechanical outputs
```

No transfer efficiency, equivalent-pressure gain, arbitrary loading area or post-result localization parameter is present.

## Package layers

- `models.py`: validated artery, geometry, material and run contracts.
- `loads.py`: signed radial integration, isotropic and anisotropy controls, and normalized spatial kernels.
- `constitutive.py`: standard-linear-solid complex moduli and reduced structural properties.
- `solvers/spectral.py`: primary periodic Fourier-spectral plate/foundation solver.
- `solvers/bounded_fd.py`: independent bounded finite-difference verification solver.
- `solvers/lumped.py`: uniform analytical foundation limit.
- `harmonics.py`: endpoint-exclusive FFT decomposition, registered harmonic controls and exact reconstruction.
- `workflow.py`: simultaneous WSS and Lamb loading, output fields, summaries and registered comparisons.
- `io.py`: NPZ input/output, configuration loading, checksums and manifests.
- `cli.py`: validation, synthetic software test generation and execution.

## Reduced structural model

For each spatial and temporal harmonic,

```math
D^*(\omega)\nabla^4w+k_f^*(\omega)w=q_L.
```

The glycocalyx is a series compliance with no traction gain. The transmitted reaction stress equals the applied normal traction. WSS is represented independently by the homogeneous shear limit and contributes a direct apical shear resultant `tau_w * h_c`.

This is the frozen first implementation architecture, not a declaration that a single plate/foundation law is the final biological model. Constitutive and load-distribution interfaces are isolated so a layered or three-dimensional model can replace the reduced law without changing hydrodynamic input records, study definitions, provenance or analysis.

## Load distributions

- `uniform_apical`: exact uniform transmitted traction.
- `localized_bound`: predeclared periodic Gaussian concentration bound, normalized to the same resultant.
- `glycocalyx_resolved`: local stiffness distribution derived from a supplied positive glycocalyx-thickness field and normalized to the same resultant.

All classes satisfy

```math
\int_A q_L(\mathbf{x},t)\,dA=Aq_L(t).
```

## Non-claim-bearing synthetic data

The `synthetic` command creates an input solely to verify software, units, signs, harmonics, controls and provenance. It must never be reported as an arterial or endothelial result.

## Uniform-load consequence

In the linear periodic plate limit, a perfectly uniform normal load produces cell-body and glycocalyx deformation but no curvature-derived in-plane bending tension. This is an analytical consequence of the registered model, not a numerical failure. `localized_bound` and `glycocalyx_resolved` are predeclared structural sensitivity classes; they are not fitted after observing outcomes.
