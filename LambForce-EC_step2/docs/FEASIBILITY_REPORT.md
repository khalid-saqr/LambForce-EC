# Step 1 Minimum Feasibility Report

## Status

```text
status: PASS
model role: feasibility only
claim-bearing interpretation: disabled
```

The calculation confirms that the hydrodynamic-to-mechanics interface is dimensionally closed and computationally solvable without a fitted force-transfer gain. It does **not** validate a final endothelial solver or establish biological meaning.

## Prototype definition

The prototype uses:

- the published Lamb-force magnitude envelope of 0.1–7 pN;
- the published 100 µm² reference area to derive an area-normalized normal load of 0.001–0.07 Pa;
- a 36.0 × 32.1 µm endothelial footprint;
- a 5 µm cell height;
- a 0.10 µm cortical thickness;
- cortex modulus states of 1.0 and 5.6 kPa;
- cytosol modulus states of 0.5 and 1.5 kPa;
- Poisson-ratio bounds of 0.45 and 0.49;
- WSS magnitude checks of 1 and 3 Pa.

The normal-load surrogate is a clamped equivalent plate on a linear cell-body foundation. The tangential-load surrogate is a homogeneous simple-shear bound. These are independent feasibility calculations, not the final coupled formulation.

## Numerical results

| Quantity | Minimum | Maximum |
|---|---:|---:|
| Area-normalized Lamb load | 0.001 Pa | 0.07 Pa |
| Normal resultant over reference cell footprint | 1.156 pN | 80.89 pN |
| Normal center displacement | 0.003333 nm | 0.7 nm |
| Positive elastic work | 1.9260e-24 J | 2.8312e-20 J |
| WSS shear strain | 0.001933 | 0.01788 |
| Tangential displacement | 9.667 nm | 89.4 nm |
| Foundation/bending stiffness ratio | 3.443e+05 | 6.07e+06 |
| Maximum force-conservation error | 0.000e+00 | — |
| Prototype execution time | 1.7354e-04 s | — |

## Exit checks

| Check | Result |
|---|---|
| Dimensional closure | PASS |
| Applied/reaction resultant conservation | PASS |
| Non-negative elastic work | PASS |
| Stable nonzero normal response | PASS |
| Stable nonzero tangential response | PASS |
| Runtime below 1 second | PASS |

## Stress-test interpretation

### 1. The interface is numerically resolvable

The predicted normal displacement spans approximately 0.003–0.7 nm in this surrogate. This is nonzero in double precision and far above machine underflow, but a final solver must show that the response remains above its mesh, quadrature and iterative-solver error floor.

### 2. WSS and Lamb loading remain mechanically distinct

The normal and tangential surrogates produce different displacement directions and scales. The design can therefore preserve WSS as tangential traction and the Lamb input as signed normal loading without converting one into the other.

### 3. Resultant conservation does not require a transfer coefficient

The area-normalized load is derived directly from the published signed force-density integral. Applied and reaction resultants agree to the reported numerical tolerance.

### 4. The simple foundation model is too dominant to freeze as the final solver

The cell-body foundation stiffness exceeds the plate-bending contribution by roughly $`3.4\times10^5`$ to $`6.1\times10^6`$. Consequently, this surrogate is controlled almost entirely by the assumed cytosolic foundation. It is suitable as a feasibility and regression bound, but it is not sufficient for the final biological claim.

Step 2 must test whether a layered glycocalyx–cortex–cell-body model changes:

- normal displacement magnitude;
- membrane tension and curvature;
- spatial concentration around the nucleus and boundaries;
- harmonic gain and phase;
- the relative influence of cortical, glycocalyx and cytosolic properties.

### 5. No biological conclusion is enabled

The calculation uses magnitude envelopes rather than artery-specific signed waveforms and does not include spatially resolved glycocalyx mechanics. It cannot determine cross-artery generality, anisotropy attribution or spectral relevance.

## Step 1 feasibility verdict

```text
dimensional closure: passed
force conservation: passed
minimum mechanics response: passed
expected runtime: acceptable
unresolved empirical transfer gains: none
final solver formulation: intentionally not selected
step 2 solver comparison required: yes
```
