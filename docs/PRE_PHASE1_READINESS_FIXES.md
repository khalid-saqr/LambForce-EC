# Pre-Phase 1 readiness fixes

## Scope

This change set resolves the repository review findings raised after Phase 0 and before immutable
six-artery archive ingestion. It does not ingest scientific data and does not generate biological
results.

## Corrected failure

PR #4 was merged even though its pull-request CI stopped at `lambforce-ec phase0-check` because
`protocol/readme_traceability.yaml` contained the invalid status `implemed`. The status is corrected,
the validator is strengthened, and CI now uses `fail-fast: false` so every supported Python version
reports independently.

## Integrity model

The former provenance check compared a registry archive checksum with a metadata value. It did not
bind the converted arrays to that archive. The revised interface separates:

- immutable archive checksum;
- exact source-member checksum;
- conversion-manifest checksum;
- converter Git commit;
- converted-record payload checksum.

The payload checksum is recomputed during every load.

## Claim-bearing contract

A claim-bearing six-artery record now requires the complete published hydrodynamic configuration,
the isotropic field, archived harmonic signatures, and the exact radial integration interval. A
missing isotropic field is an error rather than an implicit zero control.

## Governance

The Phase 0 command now validates:

- parameter-registry semantics, not only headers;
- source-registry role, status, checksum-map completeness, and uniqueness;
- exactly 34 README traceability requirements;
- exactly 32 implemented and 2 archive-blocked requirements;
- referenced implementation and test paths when run from a checkout;
- the Git blob SHA of the README used to build the matrix.

Canonical registries, protocols, schemas, and the README are installed with the wheel so the checks
also execute outside the repository.

## Phase 1 commands

Two deterministic commands are added:

- `ingest-archive` converts one exact archive member and reports all candidate checksums;
- `verify-hydrodynamics` evaluates source binding, payload integrity, radial interval, signed and
  isotropic integration, WSS, and archived harmonic reconstruction.

Neither command can turn synthetic data into claim-bearing data.

## Mechanics verification

Manufactured-solution tests are added for the bounded simply-supported (`compliant_edge`) and
clamped-edge fourth-derivative closures. The bounded solver now enforces its configured numerical
residual tolerance instead of discarding it.

## Merge gate

This change must not be merged until the Python 3.10, 3.11, and 3.12 jobs all pass:

- strict repository Phase 0 audit;
- complete test suite;
- Ruff;
- bytecode compilation;
- wheel construction;
- installed-wheel Phase 0 audit from outside the checkout.
