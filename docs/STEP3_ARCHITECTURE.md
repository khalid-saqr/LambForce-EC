# LambForce-EC production package

## Status

Package version `0.5.0` is the pre-Phase 1 readiness implementation. It is a mechanics and
qualification package, not a claim-bearing result set.

## Immutable boundary

The published anisotropic Womersley solution, signed Lamb-force density, pillbox integration, WSS,
isotropic control, and six arterial waveforms remain hydrodynamic ground truth.

```text
immutable archive + exact artery member
                ↓
checksummed deterministic conversion
                ↓
signed Lamb density + isotropic field + WSS
                ↓
archive harmonic and radial-interval qualification
                ↓
force-conserving mechanics interface
                ↓
registered endothelial outputs and decision gates
```

No transfer efficiency, equivalent-pressure gain, arbitrary loading area, or post-result
localization parameter is present.

## Package layers

- `models.py`: validated artery, geometry, material, and run contracts.
- `provenance.py`: converted-record payload hashing and tamper detection.
- `archive.py`: deterministic archive-member conversion and hydrodynamic qualification.
- `protocol.py`: semantic registry, source, traceability, and claim-bearing gates.
- `loads.py`: signed radial integration, isotropic and anisotropy controls, and normalized kernels.
- `constitutive.py`: elastic and standard-linear-solid reduced properties.
- `solvers/spectral.py`: primary periodic Fourier-spectral plate/foundation solver.
- `solvers/bounded_fd.py`: bounded finite-difference support-class verification solver.
- `solvers/lumped.py`: uniform analytical foundation limit.
- `harmonics.py`: endpoint-exclusive FFT decomposition and registered harmonic controls.
- `workflow.py`: simultaneous WSS and Lamb loading, outputs, summaries, and comparisons.
- `io.py`: integrity-checked NPZ records and checksummed result bundles.
- `cli.py`: repository audit, conversion, qualification, validation, and mechanics execution.

## Provenance layers

The following values are distinct and must not be substituted for one another:

- immutable archive SHA-256;
- exact source-member SHA-256;
- conversion-manifest SHA-256;
- converter Git commit;
- converted-record payload SHA-256;
- source-registry SHA-256 used by a claim-bearing CLI run.

The loader recomputes the converted-record checksum. Claim-bearing execution additionally requires a
unique verified source-registry match with the immutable hydrodynamic-ground-truth role.

## Reduced structural model

For each spatial and temporal harmonic,

```math
D^*(\omega)\nabla^4w+k_f^*(\omega)w=q_L.
```

The glycocalyx is a series compliance with no traction gain. WSS is independent tangential
traction. The structural ensemble contains periodic, simply-supported/compliant, and clamped
bounds, two membrane–cortex coupling limits, and two nuclear representations.

This reduced architecture is replaceable by a later layered or three-dimensional model without
changing the immutable hydrodynamic records or registered study definitions.

## Load distributions

- `uniform_apical`;
- `localized_bound`;
- `glycocalyx_resolved`.

Every class is normalized so that

```math
\int_A q_L(\mathbf{x},t)\,dA=Aq_L(t).
```

## Non-claim-bearing synthetic data

The `synthetic` command verifies software only. Synthetic records have a separate source role and
cannot pass the six-artery claim-bearing contract.

## Phase 1 boundary

Phase 1 may begin only after this package's full CI matrix passes. Archive ingestion then uses
`ingest-archive`; a converted artery becomes eligible only after `verify-hydrodynamics` passes
against a frozen verified source registry.
