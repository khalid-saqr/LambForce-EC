# Phase 0 validation report

## Status

`SOFTWARE_AND_PROTOCOL_SURFACES_VALIDATED`

PR #5 restored CI credibility and passed the Python 3.10-3.12 matrix, including semantic registries,
traceability, payload tamper detection, mechanics tests, Ruff, compilation, wheel construction, and
installed-wheel validation.

PR #6 corrects the scientific source model used by the next phase. The earlier archive-ingestion
contract is superseded because the publication identifies a public Git notebook artifact rather than
a separate six-member hydrodynamic archive.

## Active pre-scientific gate

The repository now validates:

- the exact paper DOI, picoNewton commit, v2 notebook path and Git blob;
- the frozen six-artery scalar and harmonic input registry;
- explicit `historical_v2` and `verified` reproduction modes;
- repeated deterministic generation of twelve standardized records in the software profile;
- payload tamper detection and fail-closed claim-bearing mechanics;
- the complete existing mechanics and governance test suite.

## Scientific status

The hydrodynamic reproduction has not yet been promoted to claim-bearing status. A cold execution of
the exact v2 notebook, all-six publication-resolution regression, numerical convergence, discrepancy
review, and payload freeze remain required. No claim-bearing endothelial result has been generated.
