# Step 2 Benchmark Report

## Status

```text
step: 2
status: PASS
primary solver: periodic_spectral_2d
verification solver: bounded_fd_2d
analytical baseline: lumped_0d
full 3D verification: deferred to representative later cases
```

Step 2 compares three numerical candidates under a common load, geometry, material, and harmonic interface. These results select the computational method only; they do not constitute endothelial biological findings.

## Candidate summary

| Candidate | Maximum runtime per benchmark (s) | Maximum residual | Maximum resultant error | Maximum normal displacement (m) |
|---|---:|---:|---:|---:|
| `bounded_fd_2d` | 0.00906049 | 4.082e-16 | 1.732e-16 | 7.004e-10 |
| `lumped_0d` | 1.9153e-05 | 0.000e+00 | 0.000e+00 | 7.000e-10 |
| `periodic_spectral_2d` | 0.000637499 | 3.483e-16 | 4.352e-25 | 7.000e-10 |


## Analytical and convergence results

- Periodic spectral cosine-mode maximum relative error: `8.877e-16`.
- Periodic spectral maximum benchmark runtime: `0.000637499 s`.
- Bounded finite-difference residuals remained below `4.244e-16`.
- Bounded finite-difference center-displacement self-convergence error decreased from `3.768e-02` at `N=8` to `5.857e-04` at `N=40`, relative to the `N=48` reference.

## Selection

The periodic Fourier-spectral solver is selected as the primary method because it:

1. matches the frozen periodic-monolayer structural class;
2. returns complete spatial fields;
3. supports complex harmonic amplitudes directly;
4. preserves the applied resultant;
5. has exact modal analytical solutions;
6. is substantially faster than the bounded finite-difference verification method.

The bounded finite-difference solver is retained as an independent discretization and boundary-condition verification path. The lumped model remains an analytical limit and is rejected as a primary method because it has no spatial field.

## Important scope limit

The numerical method is selected; the final biological constitutive law is not. Step 3 must preserve replaceable interfaces for glycocalyx, membrane–cortex, cytosol, nucleus, and structural boundary classes. No prototype result enables a biological claim.
