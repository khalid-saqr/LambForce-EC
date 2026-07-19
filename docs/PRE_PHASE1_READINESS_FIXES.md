# Pre-Phase 1 readiness corrections

PR #5 correctly restored CI credibility, payload-integrity checks, installed-wheel data access,
semantic governance validation, and bounded-support verification. Its assumption of a separate
immutable hydrodynamic archive was not supported by the publication's actual data-availability path.

Version 0.6.0 retains the useful protections while replacing that source model with:

```text
paper + frozen picoNewton commit + exact v2 notebook blob
+ frozen scalar/harmonic inputs + reproduction commit
+ standardized record payload checksum
```

The archive-import CLI is removed from the active workflow. Phase 1 now reproduces all six arteries
in historical and verified modes and explicitly waits for a cold execution of the exact v2 notebook
before any record can be promoted to claim-bearing status.
