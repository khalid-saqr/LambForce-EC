# Data contract

## NPZ arrays

A precomputed artery input contains:

- `radial_coordinate_m`: `(n_radial,)`, strictly increasing;
- `time_s`: `(n_time,)`, uniformly sampled and endpoint-exclusive;
- `lamb_density_signed_n_m3`: `(n_radial, n_time)`;
- `wall_shear_stress_pa`: `(n_time,)`;
- `lamb_density_isotropic_n_m3`: optional `(n_radial, n_time)`;
- `metadata_json`: scalar JSON string.

The metadata must contain `artery_id`, `artery_name`, `radius_m`, `omega0_rad_s`, `source_identifier`, `source_version`, `source_checksum` and `coordinate_convention`.

The six published arteries are reference cases selected by configuration. The validator contains no artery-name allowlist, so an additional artery requires only a conforming record and provenance.

## Sign convention

The outward wall normal is positive. The signed radial force-density field is integrated without an absolute value for mechanics. The absolute-value exposure is retained only as a diagnostic.

## Output bundle

Each run directory contains:

- `fields.npz`: complete arrays;
- `manifest.json`: units, solver and protocol versions, source and configuration checksums, mesh/basis resolution, summaries and file checksum.

A result is not claim-bearing unless its source checksum resolves to the immutable published hydrodynamic archive and the registered decision gates are evaluated.
