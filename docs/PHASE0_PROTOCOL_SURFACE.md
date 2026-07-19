# Phase 0 protocol-surface closure

## Verdict

`PASS_IMPLEMENTATION_READY`, with the six immutable hydrodynamic archives deliberately unresolved until Phase 1.

Phase 0 closes the gap between the frozen README and the executable package. It does not generate or inspect claim-bearing endothelial results.

## Frozen protocol surface

The package now validates the complete set of required dimensions:

- ten load/control cases, including inward-only, outward-only, and explicit zero-normal-load controls;
- three harmonic controls;
- three force-conserving load distributions;
- periodic, compliant-edge, and clamped-edge support classes;
- perfectly bonded and tangential-slip membrane–cortex limits;
- homogeneous and stiff-nuclear-region representations;
- elastic primary rheology and standard-linear-solid extension;
- zero prestress only.

The structural robustness campaign is frozen as the full Cartesian product: 3 load distributions × 3 lateral supports × 2 coupling limits × 2 nuclear representations = 36 structural states before parameter and spectral controls.

## Parameter-free mathematical bounds

No new fitted biological transfer parameter is introduced.

- `periodic_monolayer` uses the Fourier-spectral plate/foundation operator.
- `compliant_edge` uses a simply supported discrete plate bound.
- `clamped_edge` uses a zero-displacement/zero-slope ghost-node discrete bound.
- `tangential_slip_limit` sets direct WSS-to-membrane tension transfer to zero while retaining the shear deformation diagnostic.
- `stiff_nuclear_region` replaces the cytosol foundation stiffness inside the frozen nuclear ellipse by the registered nuclear modulus divided by frozen nuclear height.
- nonzero prestress is rejected until a new sourced freeze and impact analysis exist.

## Primary versus robustness rheology

The default run is now `rheology_mode: elastic`. In this mode all standard-linear-solid modulus ratios are forced to one without changing the registered relaxed moduli. The SLS extension is executed only after the elastic bounds pass.

## Provenance barrier

A run with `claim_bearing: true` requires a versioned source registry entry whose source identifier, source version, artery ID, and archive SHA-256 all match and whose status is `verified`. The six reference archive records are presently `awaiting_archive`, so scientific execution fails closed rather than substituting synthetic or reconstructed data.

## Remaining blocker

Phase 1 must attach and regression-verify the immutable signed Lamb-density, isotropic-control, and WSS arrays for the six arteries. No other protocol surface remains implicit.
