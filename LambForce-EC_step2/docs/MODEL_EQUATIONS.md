# Model Equations

## 1. Published hydrodynamic definitions

The Lamb vector and vorticity are

```math
\boldsymbol{\ell}=\mathbf{u}\times\boldsymbol{\omega},
\qquad
\boldsymbol{\omega}=\nabla\times\mathbf{u}.
```

For the axisymmetric anisotropic Womersley velocity field,

```math
\ell_r=u_\theta\omega_z-u_z\omega_\theta,
\qquad
f_r(r,t)=\rho\ell_r(r,t).
```

The published non-directional exposure is

```math
F_L^{\mathrm{exposure}}(t)
=
A_{\mathrm{ref}}
\int_{R-\delta_f}^{R}|f_r(r,t)|\,dr.
```

The directional mechanics input is reconstructed from the same archived signed field:

```math
F_L^{\mathrm{signed}}(t)
=
A_{\mathrm{ref}}
\int_{R-\delta_f}^{R}f_r(r,t)\,dr.
```

Area normalization removes dependence on the paper's reference footprint:

```math
q_L(t)
=
\frac{F_L^{\mathrm{signed}}(t)}{A_{\mathrm{ref}}}
=
\int_{R-\delta_f}^{R}f_r(r,t)\,dr
\qquad [\mathrm{Pa}].
```

No empirical gain is introduced.

## 2. Force-conserving endothelial loading

For a cell footprint $`A_{\mathrm{cell}}`$, every admissible spatial distribution must satisfy

```math
\int_{A_{\mathrm{cell}}}q_L(\mathbf{x},t)\,dA
=
A_{\mathrm{cell}}q_L(t).
```

WSS is applied independently as tangential traction:

```math
\mathbf{t}_{\mathrm{WSS}}(\mathbf{x},t)
=
\boldsymbol{\tau}_w(t).
```

The Lamb load is applied in the local wall-normal direction:

```math
\mathbf{t}_{L}(\mathbf{x},t)
=
q_L(\mathbf{x},t)\mathbf{n}.
```

## 3. Minimum Step 1 feasibility model

The Step 1 prototype is not the final biological solver. It is a dimensionally closed lower-complexity test used to determine whether the interface is solvable and numerically measurable.

The apical layer is represented as an equivalent clamped plate supported by a linear cell-body foundation:

```math
D\nabla^4w+k_fw=q_L,
```

with bending rigidity

```math
D=\frac{E_ch_c^3}{12(1-\nu_c^2)}
\qquad [\mathrm{N\,m}]
```

and foundation stiffness

```math
k_f=\frac{E_{\mathrm{cyt}}}{h_{\mathrm{cell}}}
\qquad [\mathrm{N\,m^{-3}}].
```

For an equivalent circular footprint with radius

```math
a=\sqrt{\frac{A_{\mathrm{cell}}}{\pi}},
```

the clamped-plate bending contribution is written as an equivalent normal stiffness

```math
k_b=\frac{64D}{a^4}
\qquad [\mathrm{N\,m^{-3}}].
```

The feasibility estimate is

```math
w_0=\frac{q_L}{k_f+k_b}.
```

The applied and reaction resultants are constrained to agree:

```math
F_{\mathrm{applied}}=q_LA_{\mathrm{cell}},
\qquad
F_{\mathrm{reaction}}=(k_f+k_b)w_0A_{\mathrm{cell}}.
```

The elastic energy estimate is

```math
U=\frac{1}{2}F_{\mathrm{applied}}w_0.
```

Tangential response is checked independently using a homogeneous simple-shear bound:

```math
G_{\mathrm{cyt}}=\frac{E_{\mathrm{cyt}}}{2(1+\nu_{\mathrm{cyt}})},
\qquad
\gamma=\frac{\tau_w}{G_{\mathrm{cyt}}},
\qquad
u_t=\gamma h_{\mathrm{cell}}.
```

## 4. Dimensional closure

| Transformation | Input unit | Operation | Output unit |
|---|---:|---|---:|
| Lamb density to area-normalized normal load | N m⁻³ | integrate over radial depth | N m⁻² = Pa |
| Normal load to resultant | Pa | multiply by cell area | N |
| Plate/foundation stiffness | N m⁻³ | multiply by displacement | Pa |
| Elastic work | N × m | one-half force–displacement product | J |
| WSS to shear strain | Pa / Pa | modulus division | 1 |
| Shear strain to displacement | 1 × m | multiply by cell height | m |

## 5. Step 2 decision left open

Step 1 does not freeze the final PDE or discretization. Step 2 must compare at least:

- a reduced membrane–foundation formulation;
- a layered continuum formulation;
- a higher-fidelity verification formulation when computationally feasible.

The final solver is chosen by conservation, analytical verification, convergence, runtime and qualitative agreement—not by implementation convenience.
