# Data contract

## Scientific boundary

The published anisotropic Womersley solution is immutable hydrodynamic ground truth. This
package does not refit, reconstruct from WSS, or replace the signed Lamb-force field.

A converted artery record is **not claim-bearing** merely because it can be loaded. It becomes
eligible only after all archive, member, conversion, payload, and hydrodynamic-regression checks
pass against the version-controlled source registry.

## Canonical NPZ arrays

Every converted record contains:

- `radial_coordinate_m`: `(n_radial,)`, strictly increasing;
- `time_s`: `(n_time,)`, uniformly sampled and endpoint-exclusive;
- `lamb_density_signed_n_m3`: `(n_radial, n_time)`;
- `lamb_density_isotropic_n_m3`: `(n_radial, n_time)` for every claim-bearing record;
- `wall_shear_stress_pa`: `(n_time,)`;
- `metadata_json`: scalar canonical JSON string.

The outward wall normal is positive. The signed field is integrated without an absolute value for
mechanics. The absolute-value integral is retained only as an exposure diagnostic.

## Four independent provenance checksums

The metadata separates four non-interchangeable checksums:

1. `source_archive_sha256` — checksum of the immutable original archive;
2. `source_member_sha256` — checksum of the exact original artery member used for conversion;
3. `conversion_manifest_sha256` — checksum of the artery conversion manifest;
4. `record_payload_sha256` — checksum of the converted arrays and immutable record metadata.

It also records `converter_commit_sha`, the full Git commit that performed conversion.

`load_artery_npz` recomputes `record_payload_sha256` and rejects altered arrays or immutable
metadata. A copied archive checksum in metadata is therefore insufficient to establish integrity.

## Claim-bearing hydrodynamic contract

For the six registered arteries, `metadata.hydrodynamic_contract` is mandatory and contains:

- `rho_kg_m3`;
- `nu_zz_m2_s`;
- `reference_area_m2`;
- `fluid_integration_depth_m`;
- `womersley_alpha`;
- `radial_collocation_order`;
- `harmonic_rms_tolerance`;
- `source_waveform_identifier`;
- `source_archive_member`;
- archived magnitude and phase arrays for:
  - signed integrated Lamb load;
  - isotropic integrated Lamb load;
  - wall shear stress.

The radial grid must begin at `R - fluid_integration_depth_m`, end at `R`, and contain exactly the
registered radial collocation order. Missing isotropic data are rejected for claim-bearing records;
they are never replaced by a zero field.

## Source registry gate

A claim-bearing run must resolve uniquely to a source-registry record with:

- `claim_bearing: true`;
- `role: immutable_hydrodynamic_ground_truth`;
- `archive_status: verified`;
- matching archive checksum;
- matching artery-member checksum;
- matching conversion-manifest checksum;
- matching converter commit.

Synthetic data are registered separately as `software_validation_only` and can never satisfy this
gate.

## Deterministic conversion

Use:

```text
lambforce-ec ingest-archive \
  --archive <immutable-archive> \
  --member-npz <exact-artery-member.npz> \
  --manifest <artery-conversion-manifest.yaml> \
  --converter-commit <40-character-git-sha> \
  --output <canonical-artery-record.npz>
```

The command does not mark a source verified. It emits candidate checksums for review and registry
freeze.

## Hydrodynamic qualification

After the source registry is frozen as verified, run:

```text
lambforce-ec verify-hydrodynamics \
  --input <canonical-artery-record.npz> \
  --source-registry <frozen-source-registry.yaml> \
  --output <qualification-report.json>
```

The qualification report verifies:

- complete source binding;
- converted-record payload integrity;
- exact radial integration interval;
- signed and isotropic fields;
- WSS waveform;
- harmonic reconstruction error against the archived tolerance;
- explicit anisotropy increment;
- exposure diagnostic;
- `total = isotropic + anisotropy increment`.

## Result bundle

Each mechanics run directory contains:

- `fields.npz`;
- `manifest.json`.

For claim-bearing CLI runs, the manifest records the converted-record checksum and the exact
source-registry checksum in addition to solver, protocol, configuration, and output-file checksums.
