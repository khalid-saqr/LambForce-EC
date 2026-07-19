# Governance source of truth

The runtime and scientific protocol use only these canonical paths:

- `registry/parameter_registry.csv`;
- `registry/parameter_sets.yaml`;
- `registry/source_registry.yaml`;
- `registry/structural_classes.yaml`;
- `configs/reference_arteries.yaml`;
- `protocol/decision_gates.yaml`;
- `protocol/readme_traceability.yaml`;
- `schemas/`.

`LambForce-EC_step2/` is a frozen historical solver-selection package retained for provenance. Its
registries, protocols, and benchmark files are not loaded by package code and must not be edited to
change the active scientific freeze.

The CLI resolves canonical data from an explicit path, the active repository checkout, or the
installed wheel's `share/lambforce_ec` data directory. Every claim-bearing result records the
checksum of the source registry actually used.
 
The original Step 3 root checksum list is retained at `historical/step3/SHA256SUMS` and applies only
to commit `65d4f7252b09b17d466ec13096335129f80465b6`.
