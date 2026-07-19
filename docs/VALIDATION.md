# Step 3 validation matrix

| Requirement | Test |
|---|---|
| Signed radial integration | constant-field analytical integral |
| Exposure remains non-directional | absolute-field analytical integral |
| Spatial resultant conservation | uniform, localized and glycocalyx-resolved kernels |
| Spectral analytical limit | uniform foundation response and single cosine mode |
| Harmonic reconstruction | FFT round trip and harmonic truncation |
| Complex viscoelastic support | passive SLS harmonic response |
| Solver residual | spectral and bounded finite-difference residual tests |
| Arbitrary artery support | non-reference artery ID accepted by the same schema |
| WSS-present comparison | combined and WSS-only runs share all non-load inputs |
| Anisotropy attribution | total, isotropic and total-minus-isotropic controls generated |
| Provenance | source and configuration checksums written to manifest |
| Colab execution | notebook installs package, mounts Drive and creates timestamped runtime directory |

The tests validate implementation behavior. They do not establish biological meaning or reproduce the six published fields without the immutable source arrays.
