<div align="center">

# LambForce-EC

## Endothelial membrane loading from the anisotropic near-wall Lamb-force field

> **Thesis:** A multiscale physical interface reveals whether the integrated anisotropic near-wall inertial force field produces biologically meaningful endothelial membrane loading independently of wall shear stress across six arterial waveforms.

</div>

---

## 1. Scientific scope

`LambForce-EC` is a standalone computational continuation of:

**K. M. Saqr. _A transverse picoNewton force revealed in anisotropic Womersley flow._ Scientific Reports 16, 12584 (2026).**  
[DOI 10.1038/s41598-026-47474-x](https://doi.org/10.1038/s41598-026-47474-x)

The published paper is the hydrodynamic ground truth. It establishes that WSS is a tangential boundary traction, while the Lamb vector is a volumetric velocity–vorticity inertial field. In anisotropic Womersley flow, the near-wall Lamb field is signed, wall-normal, artery-dependent, and spectrally rich.

This repository addresses only the missing mechanobiological bridge:

```text
published anisotropic Womersley fields
                 ↓
signed Lamb-force density + WSS
                 ↓
force-conserving apical loading
                 ↓
glycocalyx → membrane–cortex → cell body
                 ↓
endothelial membrane-loading fields
```

The primary endpoint is endothelial mechanics. Downstream molecular signalling and disease prediction are outside the claim-bearing scope.

---

## 2. Hypothesis and primary comparison

```math
\Delta \mathbf{T}^{L}_{m}(\mathbf{x},t)
=
\mathbf{T}^{\mathrm{WSS+Lamb}}_{m}(\mathbf{x},t)
-
\mathbf{T}^{\mathrm{WSS}}_{m}(\mathbf{x},t).
```

- **Null:** the incremental Lamb-induced membrane field is not resolved above numerical and physiological uncertainty.
- **Alternative:** it is resolved, anisotropy-attributable, experimentally plausible, and robust across the registered model domain.

The project does not retest whether WSS can reconstruct the Lamb field; that distinction is already established by the ground-truth paper.

---

## 3. Hydrodynamic inputs

```math
\boldsymbol{\ell}=\mathbf{u}\times\boldsymbol{\omega},
\qquad
\boldsymbol{\omega}=\nabla\times\mathbf{u},
```

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

The claim-bearing signed input is reconstructed from the archived signed field:

```math
F_L^{\mathrm{signed}}(t)
=
A_{\mathrm{ref}}
\int_{R-\delta_f}^{R}f_r(r,t)\,dr,
```

```math
q_L(t)
=
\frac{F_L^{\mathrm{signed}}(t)}{A_{\mathrm{ref}}}
=
\int_{R-\delta_f}^{R}f_r(r,t)\,dr
\qquad [\mathrm{Pa}].
```

WSS is applied simultaneously as tangential traction, $`\boldsymbol{\tau}_w(t)`$. The area-normalized Lamb input is applied as wall-normal loading. No fitted transfer efficiency is permitted.

Every spatial load model must conserve the resultant:

```math
\int_{A_{\mathrm{cell}}}q_L(\mathbf{x},t)\,dA
=
A_{\mathrm{cell}}q_L(t).
```

---

## 4. Parameter freeze

```text
parameter_freeze_version: 2.0.0
freeze_date: 2026-07-19
```

| Grade | Evidence | Use |
|---|---|---|
| A | Direct endothelial measurement of the same quantity | Primary value/range |
| B | Direct endothelial measurement from another bed, species, or preparation | Sensitivity range |
| C | Published endothelial mechanics model | Reference/secondary value |
| D | Mathematical derivation or structural assumption | Bound or verification only |

No parameter may be changed after inspecting the scientific outcome. Dependent quantities must be derived rather than sampled independently.

---

## 5. Frozen hydrodynamic registry

### Six arterial waveforms

| Artery | Radius `R` (m) | Womersley `alpha` | Six published signed harmonic amplitudes |
|---|---:|---:|---|
| Aortic root | 0.0150 | 22.03 | `[1.00, 0.82, 0.54, 0.33, 0.24, 0.17]` |
| Thoracic aorta | 0.0120 | 17.62 | `[1.00, 0.76, 0.45, 0.28, 0.20, 0.12]` |
| Femoral | 0.0040 | 5.87 | `[1.00, 0.58, 0.10, -0.17, 0.05, 0.04]` |
| Carotid | 0.0035 | 5.14 | `[1.00, 0.63, 0.31, 0.15, 0.10, 0.06]` |
| Iliac | 0.0045 | 6.61 | `[1.00, 0.51, 0.12, -0.11, 0.05, 0.03]` |
| Brachial | 0.0020 | 2.94 | `[1.00, 0.49, 0.16, -0.05, 0.02, 0.01]` |

### Constitutive and numerical rules

| Parameter | Symbol | Frozen value/rule | Unit | Source |
|---|---:|---:|---:|---|
| Anisotropy ratio | $`\beta=\nu_{z\theta}/\nu_{zz}`$ | `[-0.1, 0.1]` | 1 | Ground-truth paper |
| Anisotropy ratio | $`\gamma=\nu_{\theta z}/\nu_{zz}`$ | `[-0.1, 0.1]` | 1 | Ground-truth paper |
| Azimuthal viscosity ratio | $`\delta=\nu_{\theta\theta}/\nu_{zz}`$ | `[0.9, 1.1]` | 1 | Ground-truth paper |
| Isotropic limit | $`(\beta,\gamma,\delta)`$ | `(0,0,1)` | 1 | Exact limit |
| Radial collocation order | $`N`$ | 150 | nodes | Ground-truth paper |
| Harmonic truncation | $`H`$ | RMS reconstruction error `<1e-3` | 1 | Ground-truth paper |
| Fundamental frequency | $`\omega_0`$ | artery-specific source value | rad s⁻¹ | Waveform source |
| Harmonic frequency | $`\omega_h=h\omega_0`$ | derived | rad s⁻¹ | Definition |
| Blood density | $`\rho`$ | exact archived paper value | kg m⁻³ | Paper configuration |
| Axial viscosity | $`\nu_{zz}`$ | exact archived paper value | m² s⁻¹ | Paper configuration |
| Reference area | $`A_{\mathrm{ref}}`$ | exact archived paper value | m² | Paper configuration |
| Fluid integration depth | $`\delta_f`$ | exact archived paper value | m | Paper configuration |

The coupled model's time grid, quadrature order, and mesh density are selected only by convergence. They are not inherited constants.

### Required fields

| Field | Symbol | Unit |
|---|---:|---:|
| Signed Lamb-force density | $`f_r(r,t)`$ | N m⁻³ |
| Signed area-normalized load | $`q_L(t)`$ | Pa |
| Exposure force | $`F_L^{\mathrm{exposure}}(t)`$ | N |
| WSS | $`\tau_w(t)`$ | Pa |
| Isotropic field | $`f_r^{\mathrm{iso}}(r,t)`$ | N m⁻³ |
| Anisotropy increment | $`\Delta f_r^{\mathrm{aniso}}`$ | N m⁻³ |
| Harmonic amplitudes/phases | $`\widehat q_{L,h},\widehat\tau_{w,h},\varphi_h`$ | Pa, rad |

---

## 6. Frozen endothelial geometry

The fluid integration depth $`\delta_f`$ is not endothelial cell height.

| Parameter | Symbol | Frozen value/range | Unit | Grade | Source |
|---|---:|---:|---:|---|---|
| Streamwise footprint | $`L_x`$ | 36.0 | µm | C | Dabagh et al. 2014 |
| Transverse footprint | $`L_z`$ | 32.1 | µm | C | Dabagh et al. 2014 |
| Footprint area | $`A_{\mathrm{cell}}=L_xL_z`$ | 1155.6 | µm² | D | Derived |
| Surface-height reference | $`h_{\mathrm{cell}}`$ | 5.0 | µm | C | Dabagh et al. 2014 |
| Surface-height range | $`h_{\mathrm{cell}}`$ | 3.71–5.11 | µm | C | `4.41 ± 0.70 µm` source range |
| Reference aspect ratio | $`q`$ | 1.12 | 1 | C | Dabagh et al. 2014 |
| Aspect-ratio range | $`q`$ | 0.81–1.43 | 1 | C | `1.12 ± 0.31` source range |
| Native arterial aspect ratio | $`q`$ | 2.81 | 1 | B | Han et al. |
| Stretched arterial aspect ratio | $`q`$ | 3.65 | 1 | B | Han et al. |
| Apical cortex thickness | $`h_c`$ | 0.10 | µm | C | Dabagh et al. |
| Nuclear in-plane axes | $`a_n,b_n`$ | 8.0, 6.0 | µm | C | Dabagh et al. |
| Nuclear height | $`h_n`$ | 2.50 | µm | C | Dabagh et al. |

For an area-preserving elliptical state with $`q=a/b`$:

```math
a=\sqrt{\frac{A_{\mathrm{cell}}q}{\pi}},
\qquad
b=\sqrt{\frac{A_{\mathrm{cell}}}{\pi q}}.
```

---

## 7. Frozen material registry

### Glycocalyx

| Parameter | Symbol | Frozen set | Unit | Grade | Source |
|---|---:|---:|---:|---|---|
| Thickness | $`h_g`$ | `{0.11, 0.50, 1.00}` | µm | A–C | AFM/confocal studies and EC models |
| Young modulus | $`E_g`$ | `{25, 390, 1000}` | Pa | A–C | AFM studies and EC models |
| Poisson ratio | $`\nu_g`$ | `{0.45, 0.49}` | 1 | D | Near-incompressible bounds |
| Primary rheology | — | elastic | — | D | No unsupported relaxation constant |

Thickness and modulus from the same experimental interpretation are correlated and must not be independently recombined.

### Cell body

```math
G=\frac{E}{2(1+\nu)},
\qquad
K=\frac{E}{3(1-2\nu)}.
```

| Component | Symbol | Frozen set/range | Unit | Grade | Source |
|---|---:|---:|---:|---|---|
| Apical cortex | $`E_c`$ | `{1000, 5600}` | Pa | B/C | EC model and HAEC AFM |
| Cytosol | $`E_{\mathrm{cyt}}`$ | `500–1500` | Pa | A–C | Compression, AFM, EC models |
| Nucleus in cell | $`E_n`$ | `5000–6000` | Pa | A/C | Compression and EC models |
| Isolated-nucleus bound | $`E_n`$ | 8000 | Pa | A | Compression experiment |
| Soft-component Poisson ratio | $`\nu`$ | `{0.45, 0.49}` | 1 | D | Near-incompressible bounds |

`E`, `G`, and `K` must not be sampled independently.

### Viscoelastic robustness

```math
E^*(\omega)=E_\infty+(E_0-E_\infty)\frac{i\omega\tau}{1+i\omega\tau}.
```

| Component | Relaxation time | Unit | Grade |
|---|---:|---:|---|
| Apical cortex | 0.01–0.10 | s | C |
| Cytosol | 1.0–5.0 | s | C |
| Nucleus | 0.10–0.50 | s | C |

The modulus-ratio set is `E0/Einf ∈ {1, 2, 5}`. The primary claim must first pass the elastic bounds.

---

## 8. Structural model classes

- **Load distribution:** `uniform_apical`, `glycocalyx_resolved`, `localized_bound`.
- **Lateral support:** `periodic_monolayer`, `compliant_edge`, `clamped_edge`.
- **Membrane–cortex coupling:** `perfectly_bonded`, `tangential_slip_limit`.
- **Nuclear representation:** `homogeneous_cell_body`, `stiff_nuclear_region`.
- **Prestress:** zero in the primary model; nonzero prestress requires a separately sourced freeze.

The main conclusion may not depend only on a concentrated load, clamped boundary, or single nuclear representation.

---

## 9. Prohibited free parameters

The claim-bearing model cannot contain:

- Lamb-force transfer efficiency;
- arbitrary equivalent-pressure conversion;
- arbitrary loading area;
- continuously fitted localization;
- outcome-selected prestress;
- independent sampling of `E`, `G`, and `K`;
- artery-specific endothelial properties in the primary analysis;
- a universal biological threshold;
- any parameter adjusted after results are inspected.

---

## 10. Required simulations

For every artery, solve the same endothelial model for:

1. unloaded reference;
2. WSS only;
3. signed Lamb only;
4. WSS plus signed Lamb;
5. exposure-only diagnostic;
6. isotropic Lamb field;
7. anisotropy-specific increment;
8. fundamental only;
9. harmonics $`h\le2`$;
10. full waveform;
11. inward-only and outward-only controls;
12. zero-normal-load control;
13. elastic bounds and viscoelastic extension;
14. every structural model class.

---

## 11. Claim-bearing outputs

| Output | Symbol | SI unit |
|---|---:|---:|
| Maximum principal membrane tension | $`T_1(\mathbf{x},t)`$ | N m⁻¹ |
| Lamb-induced incremental tension | $`\Delta T_1^L`$ | N m⁻¹ |
| Maximum principal strain | $`\varepsilon_1`$ | 1 |
| Normal displacement | $`w`$ | m |
| Curvature change | $`\Delta\kappa`$ | m⁻¹ |
| Glycocalyx strain/reaction stress | $`\varepsilon_g,\sigma_g`$ | 1, Pa |
| Strain-energy density | $`U`$ | J m⁻³ |
| Work per cycle | $`W_{\mathrm{cycle}}`$ | J |
| Tension loading rate | $`\partial T_1/\partial t`$ | N m⁻¹ s⁻¹ |
| Harmonic gain/phase | $`G_h,\varphi_h`$ | 1, rad |
| Spatial concentration | $`C_s`$ | 1 |
| Incremental Lamb/WSS ratio | $`\mathcal R_L`$ | 1 |

```math
\mathcal R_L=
\frac{\left\|\mathbf T_m^{\mathrm{WSS+Lamb}}-\mathbf T_m^{\mathrm{WSS}}\right\|}
{\left\|\mathbf T_m^{\mathrm{WSS}}\right\|}.
```

Each output requires peak, RMS, cycle mean, 95th percentile, peak-to-peak range, peak time/location, spatial concentration, signed asymmetry, harmonic magnitude/phase, and cycle-integrated work.

---

## 12. Decision gates

The Lamb field is classified as computationally biologically meaningful only when all gates pass:

1. incremental loading exceeds numerical error;
2. the effect persists with WSS present;
3. the effect is attributable to the anisotropy increment;
4. predicted mechanics are consistent with documented endothelial scales;
5. the result survives all registered parameter sets;
6. the result survives all structural classes;
7. the effect occurs in at least four of six arteries;
8. the full waveform differs from the fundamental and $`h\le2`$ controls;
9. the frozen protocol predates the scientific run.

Allowed classifications:

```text
robustly_biologically_meaningful
mechanically_present_but_parameter_conditional
mechanically_present_but_below_experimental_scale
not_resolved_above_numerical_error
negative_under_registered_model_domain
```

---

## 13. Verification requirements

Before scientific runs:

- verify source checksums;
- reproduce the published six-artery hydrodynamics;
- recover the isotropic Womersley limit;
- verify signed integration and units;
- verify `total = isotropic + anisotropy increment`;
- conserve the applied resultant;
- pass zero-load, rigid, and zero-stiffness limits;
- demonstrate positive elastic energy and non-negative dissipation;
- demonstrate mesh, time, harmonic, and quadrature convergence;
- verify coordinate consistency and deterministic archives.

---

## 14. Machine-readable registry

The implementation must mirror this README in version-controlled CSV or YAML with:

```text
parameter_id, symbol, description, si_unit,
frozen_value, frozen_lower, frozen_upper, frozen_set,
source_doi, source_url, source_cell_type, source_vascular_bed,
measurement_method, source_strength, independent_or_derived,
correlation_group, primary_or_secondary, claim_bearing,
transformation_rule, notes
```

---

## 15. Verifiable sources

### Hydrodynamic basis

1. Saqr KM. _A transverse picoNewton force revealed in anisotropic Womersley flow._ [DOI 10.1038/s41598-026-47474-x](https://doi.org/10.1038/s41598-026-47474-x).
2. Willemet M, Chowienczyk P, Alastruey J. Arterial waveform database. [DOI 10.1152/ajpheart.00175.2015](https://doi.org/10.1152/ajpheart.00175.2015).

### Endothelial geometry and mechanics

3. Dabagh M et al. _Shear-induced force transmission in a multicomponent, multicell model of the endothelium._ [DOI 10.1098/rsif.2014.0431](https://doi.org/10.1098/rsif.2014.0431); [PMCID PMC4233691](https://pmc.ncbi.nlm.nih.gov/articles/PMC4233691/).
4. Dabagh M et al. _Mechanotransmission in endothelial cells subjected to oscillatory and multi-directional shear flow._ [DOI 10.1098/rsif.2017.0185](https://doi.org/10.1098/rsif.2017.0185); [PMCID PMC5454307](https://pmc.ncbi.nlm.nih.gov/articles/PMC5454307/).
5. Han et al. Arterial endothelial morphology under longitudinal stretch. [PMID 18922530](https://pubmed.ncbi.nlm.nih.gov/18922530/); [PMCID PMC2823635](https://pmc.ncbi.nlm.nih.gov/articles/PMC2823635/).
6. Ferko et al. Endothelial stress amplification around cellular structures. [PMID 17160699](https://pubmed.ncbi.nlm.nih.gov/17160699/).

### Glycocalyx and apical transmission

7. Bai K, Wang W. Glycocalyx development and mechanics. [DOI 10.1098/rsif.2011.0901](https://doi.org/10.1098/rsif.2011.0901).
8. Marsh G, Waugh RE. Endothelial glycocalyx AFM protocol. [DOI 10.3791/50163](https://doi.org/10.3791/50163).
9. Direct glycocalyx thickness/modulus measurement. [PMID 32135082](https://pubmed.ncbi.nlm.nih.gov/32135082/).
10. Weinbaum S et al. Glycocalyx-mediated mechanotransduction. [DOI 10.1073/pnas.1332808100](https://doi.org/10.1073/pnas.1332808100).
11. Glycocalyx–spectrin transmission into membrane tension. [DOI 10.1038/s41556-022-00953-5](https://doi.org/10.1038/s41556-022-00953-5).

### Cortex, cytoplasm, and nucleus

12. Caille N et al. _Contribution of the nucleus to the mechanical properties of endothelial cells._ [DOI 10.1016/S0021-9290(01)00201-9](https://doi.org/10.1016/S0021-9290(01)00201-9); [PMID 11784536](https://pubmed.ncbi.nlm.nih.gov/11784536/).
13. Costa KD et al. Human aortic endothelial AFM stiffness. [PMID 16524328](https://pubmed.ncbi.nlm.nih.gov/16524328/).
14. Sato M et al. _Micromechanical Architecture of the Endothelial Cell Cortex._ [PMCID PMC1305044](https://pmc.ncbi.nlm.nih.gov/articles/PMC1305044/).
15. Direction-dependent endothelial rheology. [PMCID PMC2563098](https://pmc.ncbi.nlm.nih.gov/articles/PMC2563098/).

---

## 16. Change control and status

Any change to a frozen value, source, structural class, or gate requires a new freeze version and an impact analysis. Silent changes and post-result tuning are prohibited.

```text
scope: endothelial membrane mechanics
protocol: frozen at version 2.0.0
implementation: not yet released
claim-bearing results: not yet generated
```

The first implementation milestone is the machine-readable parameter registry, followed by analytical and numerical verification before any six-artery scientific run.

---

## Citation

```bibtex
@article{Saqr2026LambForce,
  author  = {Saqr, Khalid M.},
  title   = {A transverse picoNewton force revealed in anisotropic Womersley flow},
  journal = {Scientific Reports},
  year    = {2026},
  volume  = {16},
  pages   = {12584},
  doi     = {10.1038/s41598-026-47474-x}
}
```

**Author:** Khalid M. Saqr  
**Contact:** `k.saqr@aast.edu`
