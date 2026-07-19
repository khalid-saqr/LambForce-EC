# Validation matrix

| Requirement | Verification |
|---|---|
| Signed radial integration | analytical constant and sign-changing fields |
| Exposure remains diagnostic | absolute-field control with no WSS |
| Spatial resultant conservation | uniform, localized, and glycocalyx-resolved kernels |
| Periodic spectral limit | uniform foundation response and single Fourier mode |
| Bounded support classes | manufactured simply-supported and clamped biharmonic fields |
| Harmonic reconstruction | FFT round trip, truncation controls, and archived waveform regression |
| Elastic and viscoelastic behavior | passive SLS limits and zero elastic dissipation |
| Solver residual | configured residual tolerance enforced by both spatial solvers |
| Converted-record integrity | payload checksum recomputed during every load |
| Immutable source binding | archive, member, manifest, converter, and source-registry checks |
| Claim-bearing completeness | isotropic field, hydrodynamic parameters, radial interval, and harmonics required |
| WSS-present comparison | combined and WSS-only runs share all non-load inputs |
| Anisotropy attribution | total, isotropic, and direct anisotropy-increment controls |
| Governance | semantic registries, exact traceability IDs, README blob, code paths, and test symbols |
| Installed package | wheel audit resolves canonical registries and protocols outside the checkout |
| Colab execution | package install, Drive mount, and timestamped runtime directory |

## Stage boundary

The software tests validate implementation and qualification behavior. They do not establish
biological meaning. Phase 1 must ingest the immutable six-artery archive and reproduce each
published signed integral, isotropic limit, WSS waveform, and harmonic signature before mechanics
runs become claim-bearing.
