# Phase 1: published-source hydrodynamic reproduction

## Purpose

Phase 1 reproduces the six published anisotropic Womersley cases from the exact public source used
for the Scientific Reports study. It is not an archive-import stage and it does not run endothelial
mechanics.

## Frozen source chain

The source chain is:

1. paper DOI `10.1038/s41598-026-47474-x`;
2. repository `khalid-saqr/picoNewton` at commit
   `4c3c36db0578373cc4e48d9d8c7e8a85944ed1cb`;
3. notebook `picoNewton_v2.ipynb` with Git blob
   `9d61c237cda75df338ce0383038f7765c886f503`;
4. `registry/published_v2_hydrodynamics.yaml`, containing the frozen scalar and harmonic inputs;
5. the LambForce-EC reproduction commit and numerical profile;
6. the payload checksum of every generated standardized record.

There is no assumed external six-member NPZ archive.

## Two explicit modes

- `historical_v2` preserves the historical differentiation layout and the harmonic-product ordering
  used for regression against the published executable artifact.
- `verified` uses the polynomial-tested differentiation orientation and reconstructs real velocity
  and vorticity fields before nonlinear Lamb-vector multiplication.

The modes are never mixed. A difference between them is reported, not hidden.

## All-six execution

Run:

```bash
lambforce-ec reproduce-hydrodynamics \
  --published-inputs registry/published_v2_hydrodynamics.yaml \
  --reproduction-commit <40-character-LambForce-EC-commit> \
  --profile publication \
  --output <phase1-directory>
```

The command generates twelve records: six arteries in each of the two modes. One frozen algorithm,
one numerical profile, and one tolerance set are used for all arteries. The `verification` profile is
only a software smoke test; the `publication` profile is required for scientific qualification.

## Required historical oracle

Generated records are not claim-bearing merely because reproduction completes. The historical mode
must be regressed against a cold, output-stripped execution of the exact v2 notebook blob. Until that
oracle is generated and reviewed, the reproduction report remains
`REPRODUCED_AWAITING_HISTORICAL_V2_ORACLE`.

The cold execution must preserve the notebook's inputs and executable ordering and export, for every
artery, the signed Lamb load, exposure measure, and WSS waveform. Any discrepancy between the cold
notebook, `historical_v2`, and `verified` modes must be reported with absolute and relative errors.

## Verification

Run:

```bash
lambforce-ec verify-hydrodynamics \
  --reproduction-directory <phase1-directory> \
  --published-inputs registry/published_v2_hydrodynamics.yaml \
  --source-registry registry/source_registry.yaml \
  --output <verification-report.json>
```

Before the historical oracle is attached, the expected result is
`PASS_WORKFLOW_READY_HISTORICAL_ORACLE_PENDING`, with `claim_bearing: false`.

## Promotion gate

The source registry may move to `verified_reproduction` only after:

- both modes execute for all six arteries at the publication profile;
- the historical mode passes the cold-notebook regression;
- the verified mode passes analytical, residual, quadrature, temporal, and isotropic-limit checks;
- all twelve record payload checksums are frozen;
- the reproduction commit and environment are recorded;
- the discrepancy report is reviewed.

Only then may mechanics consume the standardized records as claim-bearing hydrodynamic inputs.
