# Hydrodynamic data contract

## Scientific source

The published anisotropic Womersley source is a versioned Git artifact, not a presumed binary data
archive. Its identity is fixed by the paper DOI, picoNewton repository commit, v2 notebook path and
Git blob, frozen input-registry checksum, LambForce-EC reproduction commit, numerical configuration,
and generated-record payload checksum.

## Standardized arrays

Each generated artery record contains:

- `radial_coordinate_m`: strictly increasing near-wall physical radius, shape `(n_radial,)`;
- `time_s`: uniformly sampled, endpoint-exclusive cardiac cycle, shape `(n_time,)`;
- `lamb_density_signed_n_m3`: signed anisotropic radial Lamb-force density,
  shape `(n_radial, n_time)`;
- `lamb_density_isotropic_n_m3`: isotropic-control density with the same shape;
- `wall_shear_stress_pa`: simultaneous tangential traction, shape `(n_time,)`;
- `metadata_json`: canonical provenance and hydrodynamic-contract metadata.

The outward normal is positive. Signed density is integrated without an absolute value for mechanics.
The absolute-value integral remains an exposure diagnostic only.

## Published-source provenance

`metadata.published_source` records:

- source repository and exact commit;
- v2 notebook path and Git blob SHA;
- frozen published-input registry SHA-256;
- LambForce-EC reproduction commit;
- reproduction-configuration SHA-256;
- reproduction mode: `historical_v2` or `verified`.

The current NPZ serializer retains the former checksum field names as compatibility aliases for files
created by package versions 0.5.x and earlier. Their canonical meanings in version 0.6.0 are:

| Compatibility field | Canonical meaning |
|---|---|
| `source_archive_sha256` | canonical source-snapshot identity |
| `source_member_sha256` | exact published-source identity digest |
| `conversion_manifest_sha256` | reproduction-configuration digest |
| `converter_commit_sha` | LambForce-EC reproduction commit |

These aliases are not evidence that an external archive or artery-member NPZ exists. New scientific
logic uses `metadata.published_source` and the source registry.

## Payload integrity

`record_payload_sha256` covers every hydrodynamic array and immutable metadata. It is recomputed on
load. Copying source metadata onto changed arrays cannot authorize a record.

## Hydrodynamic contract

Every six-artery record must include:

- density, axial kinematic viscosity, fundamental frequency, reference area and integration depth;
- Womersley number and radial order;
- source notebook identity and reproduction mode;
- generated harmonic magnitude and phase arrays for signed Lamb load, isotropic Lamb load, and WSS;
- the exact radial interval from `R - delta_f` to `R`.

A missing isotropic field is an error for claim-bearing records and is never replaced by zero.

## Qualification state

A generated record is initially non-claim-bearing. The source registry must remain
`frozen_published_source` until the all-six historical and verified reproduction, cold-notebook
regression, numerical qualification, and payload freeze are complete.
