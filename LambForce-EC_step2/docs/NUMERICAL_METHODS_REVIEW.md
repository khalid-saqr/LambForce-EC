# Step 2 Numerical Methods Review

## Purpose

Step 2 selects the numerical method for the LambForce-EC mechanics interface. It does not generate biological claims. Every candidate receives the same normal Lamb loading, tangential WSS loading, material state, and canonical verification cases.

## Common physical interface

The normal problem is represented at this stage by

```math
D^*(\omega)
abla^4 w + k_f^*(\omega)w=q_L,
```

with

```math
D^*(\omega)=rac{E_c^*(\omega)h_c^3}{12(1-
u^2)},
\qquad
k_f^*(\omega)=rac{E_{cyt}^*(\omega)}{h_{cell}}.
```

The tangential verification response uses

```math
u_t=rac{	au_w h_{cell}}{G_{cyt}^*(\omega)}.
```

These equations are method-selection surrogates. Step 2 selects a discretization and software interface; it does not freeze this reduced constitutive model as the final biological interpretation.

## Candidate A: zero-dimensional analytical foundation

The uniform response is $`w=q_L/k_f`$. It is exact for a uniform foundation-only limit and provides a strong unit, sign, force, and complex-modulus benchmark. It cannot return spatial membrane fields and is therefore ineligible as the primary solver.

## Candidate B: periodic Fourier-spectral plate

The periodic solver diagonalizes the plate operator in Fourier space:

```math
\widehat w(\mathbf k,\omega)=
rac{\widehat q_L(\mathbf k,\omega)}{D^*(\omega)|\mathbf k|^4+k_f^*(\omega)}.
```

It directly represents a repeating endothelial monolayer unit cell, supports complex harmonic amplitudes, conserves the spatial resultant, and has an exact analytical benchmark for every Fourier mode.

## Candidate C: bounded finite-difference plate

The verification candidate applies a sparse finite-difference biharmonic operator on a bounded domain. Squaring a Dirichlet Laplacian produces a Navier-type bounded plate. This supplies an independent discretization and a boundary-condition sensitivity check. It is not treated as the primary monolayer model.

## Full three-dimensional finite elements

A full three-dimensional cell is not selected as the primary solver at Step 2. It would make six-artery parameter robustness prohibitively expensive and would introduce meshing and constitutive choices before the reduced interface is validated. Representative 3D verification remains an allowed later milestone.

## Selection logic

The primary solver must provide spatial fields, static and harmonic solutions, force conservation, analytical verification, convergence, and CPU-feasible execution. The periodic spectral solver satisfies these requirements. The bounded finite-difference method is retained as the independent verification solver; the lumped solution remains an analytical limit.
