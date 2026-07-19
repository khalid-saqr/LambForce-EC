# Phase 0 protocol-surface closure

## Status

`IMPLEMENTATION_SURFACES_CLOSED_ARCHIVE_BLOCKED`

Phase 0 closes the executable protocol surface but does not qualify any six-artery input and does
not generate claim-bearing endothelial results.

## Frozen protocol dimensions

The package validates:

- ten load and sign-control cases;
- three harmonic controls;
- three force-conserving load distributions;
- periodic, compliant-edge, and clamped-edge support classes;
- perfectly bonded and tangential-slip membrane–cortex limits;
- homogeneous and stiff-nuclear-region representations;
- elastic primary rheology and a separate standard-linear-solid extension;
- zero prestress only.

The structural campaign remains the complete 36-state Cartesian product before parameter and
spectral controls.

## Parameter-free structural bounds

- `periodic_monolayer` uses the Fourier-spectral plate/foundation operator.
- `compliant_edge` uses a simply-supported discrete plate bound.
- `clamped_edge` uses a zero-displacement/zero-slope ghost-node closure verified by manufactured
  solutions.
- `tangential_slip_limit` removes direct WSS-to-membrane tension transfer while retaining the shear
  deformation diagnostic.
- `stiff_nuclear_region` uses the frozen nuclear ellipse and registered nuclear modulus.
- nonzero prestress is rejected until a new sourced freeze and impact analysis exist.

## Fail-closed archive boundary

Claim-bearing execution requires all of the following:

- verified immutable archive SHA-256;
- exact artery-member SHA-256;
- conversion-manifest SHA-256;
- converter Git commit;
- recomputed converted-record payload SHA-256;
- complete published hydrodynamic metadata;
- signed, isotropic, and WSS arrays;
- archived harmonic signatures;
- exact radial integration interval.

Synthetic data cannot satisfy the source role or six-artery contract.

## CI gate

The repository audit, complete tests, Ruff, compilation, wheel construction, and installed-wheel
audit must pass on Python 3.10, 3.11, and 3.12 before any readiness change is merged.
