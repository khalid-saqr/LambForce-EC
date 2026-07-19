<div align="center">

# LambForce-EC

## Multiscale endothelial membrane loading driven by the anisotropic near-wall Lamb-force field

**Frozen model thesis**

> A multiscale physical interface reveals whether the integrated anisotropic near-wall inertial force field produces biologically meaningful endothelial membrane loading independently of wall shear stress across six arterial waveforms.

</div>

---

## 1. Scope and ground truth

This repository will implement a new computational continuation of:

**K. M. Saqr, “A transverse picoNewton force revealed in anisotropic Womersley flow,” Scientific Reports 16, 12584 (2026).**  
DOI: [10.1038/s41598-026-47474-x](https://doi.org/10.1038/s41598-026-47474-x)

The following are treated as **established hydrodynamic ground truth** and will not be retested as the new paper's central claim:

- the near-wall Lamb vector is a volumetric inertial field, not wall shear stress;
- the anisotropic Lamb-force waveform is fundamentally distinct from WSS;
- the integrated endothelial-scale Lamb force is artery-dependent, signed, directional, and spectrally rich;
- the six published arterial waveforms and their anisotropic Womersley solutions are immutable source inputs.

The new model addresses only the missing mechanobiological bridge:

```text
published anisotropic near-wall Lamb field
                 ↓
signed endothelial-scale normal resultant
                 ↓
glycocalyx → apical membrane–cortex → cell body
                 ↓
spatial and temporal endothelial membrane-loading fields
```

Piezo1, calcium, gene expression, disease prediction, and post-hoc biological thresholds are outside the primary scope.

---

## 2. Meaning of “frozen”

The parameter freeze is versioned as:

```text
parameter_freeze_version: 1.0.0
freeze_date: 2026-07-19
```

“Frozen” means:

1. no parameter may be selected after seeing whether the Lamb-force hypothesis passes;
2. every claim-bearing quantity has a unit, source, source-strength grade, and allowed value or range;
3. dependent quantities are calculated from independent quantities rather than sampled separately;
4. unsupported continuous tuning parameters are replaced by predeclared structural limits;
5. any future change requires a new freeze version and documented scientific reason.

This protocol does **not** claim zero biological uncertainty. It eliminates undocumented, unconstrained, and outcome-tuned parameter choices.

### Source-strength grades

| Grade | Definition | Permitted use |
|---|---|---|
| A | Direct endothelial measurement of the same parameter | Primary value or range |
| B | Direct endothelial measurement from another vascular bed, phenotype, or preparation | Primary uncertainty range |
| C | Value adopted in a published endothelial computational model | Reference value or secondary range |
| D | Derived mathematical quantity or explicit structural assumption | Cannot independently determine the conclusion |

---

## 3. Frozen physical interpretation

### 3.1 Signed mechanical input

The mechanical solver must use the signed wall-normal resultant:

\[
F_{L}^{\mathrm{signed}}(t)
=
A_{\mathrm{EC}}
\int_{R-\delta_{\mathrm{EC}}}^{R}
\rho\,\ell_r(r,t)\,dr.
\]

The absolute-value quantity remains a non-directional exposure measure:

\[
F_{L}^{\mathrm{exposure}}(t)
=
A_{\mathrm{EC}}
\int_{R-\delta_{\mathrm{EC}}}^{R}
\left|\rho\,\ell_r(r,t)\right|\,dr.
\]

`F_L_signed` is the claim-bearing load. `F_L_exposure` is descriptive only.

### 3.2 No fitted transfer efficiency

The model must not contain a free “Lamb-force transfer efficiency.” The applied normal load must conserve the published resultant:

\[
\int_{A_{\mathrm{EC}}} q_L(\mathbf{x},t)\,dA
=
F_{L}^{\mathrm{signed}}(t).
\]

Any admissible spatial load-distribution model must satisfy this identity numerically.

### 3.3 WSS and Lamb force are simultaneous, orthogonal inputs

- WSS acts as tangential traction, `tau_w(t)` in Pa.
- Lamb forcing acts as a signed wall-normal resultant, `F_L_signed(t)` in N.
- The principal comparison is the incremental membrane response:

\[
\Delta \mathbf{T}^{L}_{m}(\mathbf{x},t)
=
\mathbf{T}^{\mathrm{WSS+Lamb}}_{m}(\mathbf{x},t)
-
\mathbf{T}^{\mathrm{WSS}}_{m}(\mathbf{x},t).
\]

The model does not attempt to reconstruct the Lamb field from WSS; that physical nonredundancy is already ground truth.

---

## 4. Immutable hydrodynamic registry

All source arrays must be imported from an immutable release of the published `picoNewton` repository, with commit SHA and SHA-256 checksums recorded in every run manifest.

### 4.1 Endothelial-scale control volume

| ID | Symbol | Frozen value | SI unit | Status | Source |
|---|---:|---:|---:|---|---|
| `cell_footprint_area` | \(A_{\mathrm{EC}}\) | \(100\times10^{-12}\) | m² | Immutable | Published picoNewton baseline geometry |
| `cell_control_volume` | \(V_{\mathrm{EC}}\) | \(1.0\times10^{-15}\) | m³ | Immutable | Published picoNewton baseline geometry |
| `control_volume_depth` | \(\delta_{\mathrm{EC}}=V/A\) | \(1.0\times10^{-5}\) | m | Derived | Exact derivation from the two values above |

Equivalent human-readable values:

- `A_EC = 100 µm²`
- `V_EC = 1000 µm³`
- `delta_EC = 10 µm`

### 4.2 Six arterial inputs

| Artery | Radius \(R\) (m) | Womersley \(\alpha\) | Six published signed harmonic amplitudes |
|---|---:|---:|---|
| Aortic root | 0.0150 | 22.03 | `[1.00, 0.82, 0.54, 0.33, 0.24, 0.17]` |
| Thoracic aorta | 0.0120 | 17.62 | `[1.00, 0.76, 0.45, 0.28, 0.20, 0.12]` |
| Femoral | 0.0040 | 5.87 | `[1.00, 0.58, 0.10, -0.17, 0.05, 0.04]` |
| Carotid | 0.0035 | 5.14 | `[1.00, 0.63, 0.31, 0.15, 0.10, 0.06]` |
| Iliac | 0.0045 | 6.61 | `[1.00, 0.51, 0.12, -0.11, 0.05, 0.03]` |
| Brachial | 0.0020 | 2.94 | `[1.00, 0.49, 0.16, -0.05, 0.02, 0.01]` |

### 4.3 Hydrodynamic quantities required from the source archive

| ID | Symbol | Unit | Frozen rule |
|---|---:|---:|---|
| `lamb_density_signed` | \(f_r(r,t)=\rho\ell_r\) | N m⁻³ | Full signed radial field |
| `lamb_force_signed` | \(F_L^{signed}(t)\) | N | Reintegrated without absolute value |
| `lamb_force_exposure` | \(F_L^{exposure}(t)\) | N | Published absolute-value exposure |
| `wall_shear_stress` | \(\tau_w(t)\) | Pa | Published tangential wall traction |
| `lamb_force_isotropic` | \(F_L^{iso}(t)\) | N | Isotropic control |
| `lamb_force_anisotropy_increment` | \(\Delta F_L^{aniso}\) | N | `total - isotropic` |
| `force_ratio` | \(\chi(t)=|F_r|/|F_z|\) | 1 | Published directional descriptor |
| `force_angle` | \(\phi(t)\) | rad | Published resultant angle |
| `harmonic_frequency` | \(\omega_h=h\omega_0\) | rad s⁻¹ | Derived |
| `harmonic_complex_amplitudes` | \(\widehat G_h\) | Pa m⁻¹ | Immutable source data |
| `blood_density` | \(\rho\) | kg m⁻³ | Exact published configuration value |
| `axial_kinematic_viscosity` | \(\nu_{zz}\) | m² s⁻¹ | Exact published configuration value |

### 4.4 Anisotropy and numerical values

| ID | Symbol | Frozen value/range | Unit | Source |
|---|---:|---:|---:|---|
| `anisotropy_beta` | \(\beta=\nu_{z\theta}/\nu_{zz}\) | `[-0.1, 0.1]` | 1 | Published study |
| `anisotropy_gamma` | \(\gamma=\nu_{\theta z}/\nu_{zz}\) | `[-0.1, 0.1]` | 1 | Published study |
| `anisotropy_delta` | \(\delta=\nu_{\theta\theta}/\nu_{zz}\) | `[0.9, 1.1]` | 1 | Published study |
| `isotropic_beta` | \(\beta\) | 0 | 1 | Exact isotropic limit |
| `isotropic_gamma` | \(\gamma\) | 0 | 1 | Exact isotropic limit |
| `isotropic_delta` | \(\delta\) | 1 | 1 | Exact isotropic limit |
| `radial_spectral_order` | \(N\) | 150 | nodes | Published production resolution |
| `time_points_per_cycle` | \(N_t\) | 2048 | samples | Frozen publication workflow |
| `near_wall_quadrature_points` | \(N_q\) | 256 | points | Frozen publication workflow |
| `harmonic_reconstruction_error` | — | `< 1e-3` | normalized RMS | Published truncation rule |

---

## 5. Frozen cell geometry

The mechanical model must preserve `A_EC = 100 µm²` and `V_EC = 1000 µm³` for consistency with the hydrodynamic integration.

### 5.1 Independent geometry parameters

Only these geometry quantities are independent:

| ID | Symbol | Frozen values | Unit | Grade | Source |
|---|---:|---:|---:|---|---|
| `cell_area` | \(A_{EC}\) | 100 | µm² | Ground truth | Published picoNewton geometry |
| `cell_volume` | \(V_{EC}\) | 1000 | µm³ | Ground truth | Published picoNewton geometry |
| `cell_aspect_ratio_native` | \(q\) | 2.81 | 1 | B | Porcine carotid endothelial measurement |
| `cell_aspect_ratio_stretched` | \(q\) | 3.65 | 1 | B | Same study; elevated axial stretch |
| `cortex_thickness` | \(h_c\) | 0.10 | µm | C | Published endothelial multicomponent models |
| `glycocalyx_thickness_low` | \(h_g\) | 0.11 | µm | A | Direct AFM thin-layer interpretation |
| `glycocalyx_thickness_reference` | \(h_g\) | 0.50 | µm | C | Published mechanotransmission model |
| `glycocalyx_thickness_high` | \(h_g\) | 1.00 | µm | A/B | Direct AFM/confocal experimental envelope |
| `nucleus_height` | \(h_n\) | 2.50 | µm | C | Published endothelial geometry model |
| `nucleus_center_height` | \(z_n\) | 1.25 | µm | C | Published endothelial geometry model |

### 5.2 Derived geometry

The cell footprint is elliptical. For aspect ratio \(q=a/b\):

\[
a=\sqrt{\frac{A_{EC}q}{\pi}},
\qquad
b=\sqrt{\frac{A_{EC}}{\pi q}}.
\]

The mechanical cell height is derived, not fitted:

\[
h_{cell}=\frac{V_{EC}}{A_{EC}}=10\ \mu\mathrm m.
\]

Nuclear in-plane dimensions are scaled from the published endothelial model ratios rather than introduced as independent parameters:

\[
\frac{a_n}{a}=\frac{8}{18},
\qquad
\frac{b_n}{b}=\frac{6}{16.05}.
\]

This preserves the source-model nucleus-to-cell proportions while maintaining the ground-truth cell area.

### 5.3 Geometry sources

- Han et al., arterial endothelial aspect ratios `2.81 ± 0.25` and `3.65 ± 0.38`: [PMCID PMC2823635](https://pmc.ncbi.nlm.nih.gov/articles/PMC2823635/), [PMID 18922530](https://pubmed.ncbi.nlm.nih.gov/18922530/).
- Dabagh et al., seven-cell endothelial geometry, 500 nm glycocalyx, 100 nm cortex, and nucleus dimensions: [PMCID PMC5454307](https://pmc.ncbi.nlm.nih.gov/articles/PMC5454307/).
- Marsh and Waugh, direct endothelial glycocalyx AFM protocol: [DOI 10.3791/50163](https://doi.org/10.3791/50163).
- Gao and Lipowsky, endothelial glycocalyx composition and thickness: [DOI 10.1016/j.mvr.2010.06.005](https://doi.org/10.1016/j.mvr.2010.06.005).

---

## 6. Frozen elastic material parameters

### 6.1 Constitutive rule

The independent material inputs are Young's modulus `E` and Poisson ratio `nu`. Shear and bulk moduli must be derived:

\[
G=\frac{E}{2(1+\nu)},
\qquad
K=\frac{E}{3(1-2\nu)}.
\]

`E`, `G`, and `K` must never be sampled independently.

### 6.2 Claim-bearing elastic registry

| Component | ID | Symbol | Frozen values | SI unit | Grade | Verifiable source |
|---|---|---:|---:|---:|---|---|
| Glycocalyx | `glycocalyx_modulus_low` | \(E_g\) | 25 | Pa | A | Direct AFM, 110 nm layer |
| Glycocalyx | `glycocalyx_modulus_reference` | \(E_g\) | 390 | Pa | A/C | AFM envelope and endothelial FE reference |
| Cortex | `cortex_modulus_reference` | \(E_c\) | 1000 | Pa | C | Endothelial multicomponent FE model |
| Cortex | `cortex_modulus_high` | \(E_c\) | 5600 | Pa | B | Human aortic EC AFM over cytoskeletal structures |
| Cytosol | `cytosol_modulus_low` | \(E_{cyt}\) | 500 | Pa | A/B | Endothelial cell compression and FE fit |
| Cytosol | `cytosol_modulus_high` | \(E_{cyt}\) | 1500 | Pa | B | Human aortic EC AFM adjacent cytoplasm |
| Nucleus | `nucleus_modulus_reference` | \(E_n\) | 5000 | Pa | A/B | Endothelial cell compression and FE fit |
| Nucleus | `nucleus_modulus_high` | \(E_n\) | 8000 | Pa | A/B | Isolated endothelial nuclei |
| All soft components | `poisson_ratio_primary` | \(\nu\) | 0.45 | 1 | D | Near-incompressible lower structural limit |
| All soft components | `poisson_ratio_upper` | \(\nu\) | 0.49 | 1 | D | Near-incompressible upper structural limit |

### 6.3 Elastic sources

- Endothelial glycocalyx thickness `110 nm` and modulus `0.025 kPa`, plus whole-cell modulus `3.0–6.5 kPa`: [PMID 32135082](https://pubmed.ncbi.nlm.nih.gov/32135082/).
- HUVEC glycocalyx development and mechanical properties: [DOI 10.1098/rsif.2011.0901](https://doi.org/10.1098/rsif.2011.0901).
- Cytoplasm approximately `500 Pa`, in-cell nucleus approximately `5000 Pa`, isolated nucleus approximately `8000 Pa`: [PMID 11784536](https://pubmed.ncbi.nlm.nih.gov/11784536/).
- Endothelial FE reference values `500 Pa` cytosol, `6000 Pa` nucleus, `1000 Pa` cortex, and `390 Pa` glycocalyx: [PMCID PMC4233691](https://pmc.ncbi.nlm.nih.gov/articles/PMC4233691/).
- Human aortic endothelial cytoskeletal-region and adjacent-cytoplasm stiffness measurements: [PMID 16524328](https://pubmed.ncbi.nlm.nih.gov/16524328/).

---

## 7. Frozen viscoelastic extension

The **primary claim must pass both elastic limits** before the viscoelastic extension is interpreted.

Each viscoelastic component uses a standard-linear-solid complex modulus:

\[
E^*(\omega)
=
E_{\infty}
+
(E_0-E_{\infty})
\frac{i\omega\tau}{1+i\omega\tau}.
\]

Only `E0`, `Einf`, and `tau` are independent.

| Component | ID | Frozen range | Unit | Grade | Use |
|---|---|---:|---:|---|---|
| Apical cortex | `cortex_relaxation_time` | 0.01–0.10 | s | C | Secondary frequency-response analysis |
| Cytosol | `cytosol_relaxation_time` | 1.0–5.0 | s | C | Secondary low-frequency support |
| Nucleus | `nucleus_relaxation_time` | 0.10–0.50 | s | C | Secondary nuclear-support analysis |
| Glycocalyx | — | none | — | — | Elastic in primary model; no unsupported relaxation constant |

Frozen modulus-ratio brackets:

```text
E0 / Einf ∈ {1, 2, 5}
```

These are three discrete structural brackets, not a continuously fitted parameter.

Sources:

- Endothelial multicomponent dynamic model and material architecture: [PMCID PMC5454307](https://pmc.ncbi.nlm.nih.gov/articles/PMC5454307/).
- Dynamic endothelial nuclear/cytoplasmic stiffness modeling under physiological shear: [PMCID PMC10627551](https://pmc.ncbi.nlm.nih.gov/articles/PMC10627551/).
- Direction-dependent endothelial cytoplasmic rheology: [PMCID PMC2563098](https://pmc.ncbi.nlm.nih.gov/articles/PMC2563098/).

---

## 8. Parameters explicitly excluded from tuning

The following are prohibited as free parameters in the claim-bearing model:

| Excluded parameter | Frozen treatment |
|---|---|
| Lamb-force transfer efficiency | Removed; exact resultant conservation enforced |
| Arbitrary loading area | Removed; `A_EC = 100 µm²` |
| Independent cell length and width | Derived from area and aspect ratio |
| Independent cell height | Derived from `V_EC/A_EC` |
| Independent `E`, `G`, and `K` | `G` and `K` derived from `E` and `nu` |
| Arbitrary localization factor | Replaced by three force-conserving structural models |
| Baseline membrane prestress | Primary model fixed at zero; nonzero values prohibited until a direct endothelial source is registered |
| Membrane pressure ceiling | Removed |
| Piezo1 channel count | Outside scope |
| Piezo1 spatial fraction | Outside scope |
| Calcium conversion gain | Outside scope |
| Active contractile force | Outside primary scope |
| Focal-adhesion molecular count | Outside primary scope |
| Adherens-junction molecular count | Outside primary scope |
| Universal biological threshold | Prohibited; biological meaning requires convergent evidence |

---

## 9. Frozen structural model classes

These are categorical models, not continuous fit parameters.

### 9.1 Normal Lamb-load distribution

All cases conserve the same signed resultant.

1. `uniform_apical` — uniform normal load over the complete `A_EC`.
2. `glycocalyx_continuum` — distribution obtained from the solved glycocalyx layer.
3. `peripheral_upper_bound` — peripheral concentration normalized to conserve total force.

The main conclusion must not depend only on `peripheral_upper_bound`.

### 9.2 Boundary support

1. `periodic_monolayer` — primary representation.
2. `compliant_edge` — secondary finite-neighbour support.
3. `clamped_edge` — stiff upper bound only.

The main conclusion must not depend only on `clamped_edge`.

### 9.3 Membrane–cortex attachment

1. `perfectly_bonded`.
2. `tangential_slip_limit`.

No unmeasured adhesion coefficient may be introduced before a direct source is registered.

### 9.4 Nuclear representation

1. `homogeneous_foundation` — no explicit nucleus.
2. `stiff_nuclear_patch` — source-scaled explicit nucleus.

The main conclusion must remain qualitatively consistent under both.

---

## 10. Frozen simulation matrix

Every artery must be solved under the same endothelial model for:

1. unloaded cell;
2. WSS only;
3. signed Lamb force only;
4. WSS plus signed Lamb force;
5. exposure-only diagnostic;
6. isotropic Lamb component;
7. anisotropy-specific Lamb increment;
8. fundamental harmonic only;
9. harmonics `h <= 2`;
10. complete harmonic waveform;
11. inward-only control;
12. outward-only control;
13. elastic instantaneous limit;
14. elastic equilibrium limit;
15. full viscoelastic model;
16. zero-normal-load control;
17. each load-distribution class;
18. each boundary-support class;
19. each membrane–cortex attachment class.

No artery-specific endothelial parameters are permitted in the primary analysis.

---

## 11. Claim-bearing outputs

The package must archive complete spatial and temporal fields, not only aggregate RMS values.

| Output | Symbol | SI unit | Meaning |
|---|---:|---:|---|
| Maximum principal membrane tension | \(T_1(\mathbf{x},t)\) | N m⁻¹ | Primary membrane mechanical state |
| Lamb-induced incremental tension | \(\Delta T_1^L\) | N m⁻¹ | Independent Lamb contribution with WSS present |
| Maximum principal strain | \(\varepsilon_1\) | 1 | Local deformation |
| Normal displacement | \(w\) | m | Wall-normal deformation |
| Curvature change | \(\Delta\kappa\) | m⁻¹ | Local membrane geometry change |
| Glycocalyx strain | \(\varepsilon_g\) | 1 | Glycocalyx compression/extension |
| Glycocalyx reaction stress | \(\sigma_g\) | Pa | Transmitted apical load |
| Strain-energy density | \(U\) | J m⁻³ | Stored mechanical energy |
| Work per cycle | \(W_{cycle}\) | J | Repeated mechanical exposure |
| Tension loading rate | \(\partial T_1/\partial t\) | N m⁻¹ s⁻¹ | High-frequency mechanical stimulus |
| Harmonic gain | \(G_h\) | 1 | Mechanical transmission by harmonic |
| Harmonic phase | \(\varphi_h\) | rad | Phase lag |
| Spatial peak-to-mean ratio | \(C_s\) | 1 | Mechanical concentration |
| Incremental Lamb/WSS norm ratio | \(\mathcal{R}_L\) | 1 | Relative mechanical contribution |

The primary incremental ratio is:

\[
\mathcal{R}_L
=
\frac{
\left\|
\mathbf{T}_m^{WSS+Lamb}
-
\mathbf{T}_m^{WSS}
\right\|
}{
\left\|
\mathbf{T}_m^{WSS}
\right\|
}.
\]

Required summaries for every field:

- peak;
- RMS;
- cycle mean;
- 95th percentile;
- peak-to-peak range;
- time and location of peak;
- spatial maximum-to-mean ratio;
- harmonic magnitude and phase for `h = 1...6`;
- signed inward/outward asymmetry;
- cycle-integrated work.

---

## 12. Biological-meaningfulness decision rule

No single arbitrary force or tension threshold is permitted.

The Lamb field is classified as **computationally biologically meaningful** only when all gates pass:

1. `numerical_significance` — incremental membrane loading exceeds the combined discretization and solver error floor;
2. `wss_present` — the incremental effect persists in `WSS + Lamb`, not only in Lamb-only simulations;
3. `anisotropy_attribution` — the effect is present in the anisotropy-specific increment and is not explained by the isotropic component;
4. `experimental_scale_consistency` — predicted deformation, tension, stress, or local reaction force is consistent with experimentally observed endothelial mechanical scales;
5. `parameter_robustness` — the result survives the complete registered glycocalyx, cortex, cytosol, and nuclear ranges;
6. `structural_robustness` — the result survives the predeclared load-distribution, support, attachment, and nuclear model classes;
7. `cross_artery_generality` — the predeclared effect passes in at least four of six arteries;
8. `spectral_relevance` — the full waveform produces a resolved change in loading rate, phase, work, or peak tension relative to the `h=1` and `h<=2` controls;
9. `no_post_hoc_tuning` — the exact freeze version predates the publication run.

Allowed final classifications:

```text
robustly_biologically_meaningful
mechanically_present_but_parameter_conditional
mechanically_present_but_below_experimental_scale
not_resolved_above_numerical_error
negative_under_registered_model_domain
```

---

## 13. Verification requirements

Before scientific runs are enabled, the package must pass:

- source-data checksum verification;
- reproduction of the published six-artery hydrodynamic outputs;
- signed-force integration check;
- isotropic-limit check;
- total = isotropic + anisotropy-increment identity;
- load-resultant conservation for every spatial distribution;
- zero-load response;
- rigid and zero-stiffness limits;
- elastic instantaneous and equilibrium limits;
- positive strain energy and viscoelastic dissipation;
- mesh convergence;
- temporal/harmonic convergence;
- coordinate-rotation consistency;
- deterministic rerun and archive checksums.

---

## 14. Machine-readable registry schema

The implementation must reproduce this freeze in `parameters/parameter_registry.csv` or YAML with these fields:

```text
parameter_id
symbol
description
si_unit
frozen_value
frozen_lower
frozen_upper
frozen_set
source_doi
source_url
source_cell_type
source_vascular_bed
measurement_method
source_strength
independent_or_derived
correlation_group
primary_or_secondary
claim_bearing
transformation_rule
notes
```

Mandatory correlation groups:

- `A_EC`, `V_EC`, and `h_cell`;
- `A_EC`, aspect ratio, semi-major axis, and semi-minor axis;
- `E`, `nu`, `G`, and `K`;
- `E0`, `Einf`, and `tau`;
- glycocalyx thickness and modulus when obtained from the same AFM interpretation.

---

## 15. Verifiable source register

### Hydrodynamic ground truth

1. Saqr KM. *A transverse picoNewton force revealed in anisotropic Womersley flow.* Scientific Reports 16, 12584 (2026). [DOI 10.1038/s41598-026-47474-x](https://doi.org/10.1038/s41598-026-47474-x).
2. Willemet M, Chowienczyk P, Alastruey J. Virtual healthy-subject arterial waveform database. [DOI 10.1152/ajpheart.00175.2015](https://doi.org/10.1152/ajpheart.00175.2015).

### Endothelial glycocalyx

3. Marsh G, Waugh RE. Direct AFM measurement protocol. [DOI 10.3791/50163](https://doi.org/10.3791/50163).
4. Bai K, Wang W. Spatio-temporal glycocalyx development and mechanical properties. [DOI 10.1098/rsif.2011.0901](https://doi.org/10.1098/rsif.2011.0901).
5. Endothelial glycocalyx thickness `110 nm` and modulus `0.025 kPa`. [PMID 32135082](https://pubmed.ncbi.nlm.nih.gov/32135082/).
6. Weinbaum S, Tarbell JM, Damiano ER. Glycocalyx structure and function. [DOI 10.1146/annurev.bioeng.9.060906.151959](https://doi.org/10.1146/annurev.bioeng.9.060906.151959).
7. Weinbaum et al. Glycocalyx mechanotransduction and cortical linkage. [DOI 10.1073/pnas.1332808100](https://doi.org/10.1073/pnas.1332808100).

### Endothelial geometry and mechanics

8. Han et al. Arterial endothelial morphology under longitudinal stretch. [PMCID PMC2823635](https://pmc.ncbi.nlm.nih.gov/articles/PMC2823635/).
9. Caille et al. Cytoplasmic and nuclear elastic moduli. [PMID 11784536](https://pubmed.ncbi.nlm.nih.gov/11784536/).
10. Dabagh et al. Multicomponent endothelial force transmission. [PMCID PMC4233691](https://pmc.ncbi.nlm.nih.gov/articles/PMC4233691/).
11. Dabagh et al. Oscillatory and multidirectional endothelial mechanotransmission. [PMCID PMC5454307](https://pmc.ncbi.nlm.nih.gov/articles/PMC5454307/).
12. Ferko et al. Endothelial FE stress amplification around nucleus and focal adhesions. [PMID 17160699](https://pubmed.ncbi.nlm.nih.gov/17160699/).
13. Human aortic endothelial AFM stiffness measurements. [PMID 16524328](https://pubmed.ncbi.nlm.nih.gov/16524328/).
14. Dynamic nuclear/cytoplasmic stiffness under shear. [PMCID PMC10627551](https://pmc.ncbi.nlm.nih.gov/articles/PMC10627551/).
15. Endothelial anisotropic rheology. [PMCID PMC2563098](https://pmc.ncbi.nlm.nih.gov/articles/PMC2563098/).

### Direct membrane-mechanics comparison scales

16. Endothelial membrane tether extraction mechanics. [PMCID PMC2990408](https://pmc.ncbi.nlm.nih.gov/articles/PMC2990408/).
17. Spectrin–glycocalyx transmission into endothelial membrane tension. [DOI 10.1038/s41556-022-00953-5](https://doi.org/10.1038/s41556-022-00953-5).

---

## 16. Change-control rule

The parameter registry may change only through a pull request containing:

- the proposed value or range;
- the old and new source;
- the source's cell type, vascular bed, and measurement method;
- the reason the existing source is inadequate;
- an impact analysis showing whether prior conclusions change;
- a new `parameter_freeze_version`.

Silent parameter changes are prohibited.

---

## 17. Current project status

```text
repository state: protocol and parameter freeze
solver state: not yet implemented
scientific outcome: not yet evaluated
claims enabled: false
```

The next implementation artifact should be the machine-readable parameter registry corresponding exactly to this README, followed by verification tests before any claim-bearing simulation is run.
