# Step 2 Solver Selection

## Decision

```text
primary solver: periodic Fourier-spectral 2D plate/foundation
verification solver: bounded sparse finite-difference 2D plate/foundation
analytical baseline: zero-dimensional uniform foundation
full 3D verification: deferred to representative cases after the primary solver is validated
```

## Why the periodic spectral solver is primary

1. The primary structural class in the frozen README is a periodic endothelial monolayer.
2. The Fourier operator gives exact modal solutions for the linear reduced model.
3. Static and complex harmonic inputs use the same interface.
4. Spatial normal-load distributions remain force conserving.
5. The method exposes displacement, foundation reaction, bending reaction, energy, phase, and harmonic gain fields.
6. Runtime scales well enough for six arteries and parameter robustness.

## Why the other candidates are not primary

The lumped model has no spatial fields. The bounded finite-difference model imposes a finite-domain boundary condition inconsistent with the primary periodic-monolayer assumption and is materially slower, but it provides an independent verification path.

## Scope boundary

This selection freezes the numerical architecture for the first implementation. It does not claim that a single-layer plate/foundation law is the final biological model. Step 3 must expose the constitutive and geometry interfaces so that a more detailed glycocalyx–membrane–cortex formulation can replace the prototype law without changing hydrodynamic inputs, studies, analysis, or provenance.

## Exit status

The machine-readable decision is stored in `benchmarks/solver_selection.json`. Step 2 passes only when every hard gate in that file is true.
