# Auto-generated results summary

All numbers below come from this run (seed=0); interpret with domain knowledge.

## 0. Cohort
n = 1842 spectra after QC; grouping unit for CV / clustering: `formula`.

| step | n |
|---|---|
| All samples | 2047 |
| Spectrum quality in {excellent,fair} | 1905 |
| SNR >= 20 | 1880 |
| Detected peaks >= 3 | 1845 |
| Orientation = unoriented | 1842 |

## Key findings (auto-generated from this run's numbers)

Hypothesis-driven pairs fixed in code (HEADLINE_PAIRS); report as exploratory if chosen after viewing results.

| structure | spectrum | partial rho [95% CI] | q |
|---|---|---|---|
| ln Bond reduced mass | Low-band fraction | +0.33 [+0.29, +0.37] | <0.001 |
| ln Bond reduced mass | ln N_peak | -0.33 [-0.37, -0.29] | <0.001 |
| ln Bond reduced mass | High-band fraction | -0.32 [-0.37, -0.28] | <0.001 |
| Bond electroneg. diff. | ln W1 distance | -0.20 [-0.24, -0.15] | <0.001 |
| Site-mixing entropy | ln Mean linewidth gamma | +0.28 [+0.23, +0.32] | <0.001 |
| Electroneg. dispersion | High-band fraction | +0.26 [+0.21, +0.30] | <0.001 |
| N anion-group types | ln Median peak spacing | -0.04 [-0.08, +0.01] | 0.130 |
| N anion-group types | ln N_peak | +0.28 [+0.23, +0.32] | <0.001 |

- Group-theoretical N_Raman and observed N_peak show a weak monotonic association (Spearman rho = 0.33).
- Structure/chemistry adds out-of-fold information beyond family + crystal system (paired dR2 95% CI > 0) for: ln N_peak (+0.07 [+0.04, +0.10]); Peak-intensity evenness (+0.03 [+0.01, +0.06]); ln Mean linewidth gamma (+0.24 [+0.20, +0.29]); ln Median peak spacing (+0.06 [+0.02, +0.10]); Peak overlap ratio (+0.08 [+0.05, +0.12]); High-band fraction (+0.04 [+0.01, +0.07]); ln W1 distance (+0.05 [+0.02, +0.08]).
- No detectable gain over family + crystal system (CI includes or is below 0) for: Low-band fraction, Mid-band fraction.
- Measurement conditions alone out-predict all structure/chemistry descriptors for: Peak-intensity evenness (treat these endpoints as instrument-dominated).
- Endpoints with out-of-fold R2 whose 95% CI excludes 0: ln N_peak (R2=0.38), Peak-intensity evenness (R2=0.09), ln Mean linewidth gamma (R2=0.26), ln Median peak spacing (R2=0.18), Peak overlap ratio (R2=0.17), Low-band fraction (R2=0.39), Mid-band fraction (R2=0.18), High-band fraction (R2=0.34), ln W1 distance (R2=0.29).
- Largest joint-permutation group per endpoint (dR2 >= 0.02): ln N_peak: Bond geometry (dR2=0.18); Peak-intensity evenness: Bond geometry (dR2=0.10); ln Mean linewidth gamma: Chemistry (dR2=0.24); ln Median peak spacing: Symmetry / group theory (dR2=0.17); Peak overlap ratio: Bond geometry (dR2=0.13); Low-band fraction: Bond geometry (dR2=0.29); Mid-band fraction: Bond geometry (dR2=0.23); High-band fraction: Chemistry (dR2=0.32); ln W1 distance: Bond geometry (dR2=0.21).
- QC sensitivity across 4 alternative settings: sign agreement of main strong partial correlations 1.00-1.00; sign agreement with q<0.05 retained 1.00-1.00; rank agreement of the whole partial-correlation matrix 0.99-1.00.

## 1. Strongest partial correlations (adjusted for ln N_atom, ln SNR)

| structure | spectrum | partial rho [95% CI] | q | n |
|---|---|---|---|---|
| Raman-active fraction | Mapping efficiency eta | -0.48 [-0.52, -0.44] | <0.001 | 1833 |
| Centrosymmetric (0/1) | Mapping efficiency eta | +0.37 [+0.33, +0.41] | <0.001 | 1833 |
| ln Bond reduced mass | Low-band fraction | +0.33 [+0.29, +0.37] | <0.001 | 1807 |
| ln Bond reduced mass | ln N_peak | -0.33 [-0.37, -0.29] | <0.001 | 1807 |
| ln Bond reduced mass | High-band fraction | -0.32 [-0.37, -0.28] | <0.001 | 1807 |
| Degeneracy fraction | Mapping efficiency eta | +0.32 [+0.28, +0.37] | <0.001 | 1833 |
| N anion-group types | ln W1 distance | -0.32 [-0.37, -0.28] | <0.001 | 1842 |
| N anion-group types | High-band fraction | +0.32 [+0.28, +0.36] | <0.001 | 1842 |
| Electroneg. dispersion | Low-band fraction | -0.32 [-0.36, -0.28] | <0.001 | 1842 |
| Max Raman irrep mult. | Mapping efficiency eta | -0.32 [-0.36, -0.27] | <0.001 | 1833 |
| ln Wyckoff orbits | Mapping efficiency eta | -0.30 [-0.34, -0.26] | <0.001 | 1833 |
| Bond electroneg. diff. | Low-band fraction | -0.28 [-0.33, -0.24] | <0.001 | 1807 |
| Site-mixing entropy | ln Mean linewidth gamma | +0.28 [+0.23, +0.32] | <0.001 | 1842 |
| N anion-group types | ln N_peak | +0.28 [+0.23, +0.32] | <0.001 | 1842 |
| ln SOAP N environments | Mapping efficiency eta | -0.27 [-0.31, -0.23] | <0.001 | 1833 |

## 2. Categorical factors (top 10 eta2_H)

| factor | spectrum | eta2_H | q |
|---|---|---|---|
| family | ln N_peak | 0.27 | <0.001 |
| family | Low-band fraction | 0.25 | <0.001 |
| family | High-band fraction | 0.24 | <0.001 |
| family | ln W1 distance | 0.23 | <0.001 |
| family | Mid-band fraction | 0.15 | <0.001 |
| n_anion_group_types_cat | High-band fraction | 0.14 | <0.001 |
| n_anion_group_types_cat | ln W1 distance | 0.13 | <0.001 |
| crystal_system | Mapping efficiency eta | 0.12 | <0.001 |
| n_anion_group_types_cat | ln N_peak | 0.11 | <0.001 |
| family | Peak overlap ratio | 0.09 | <0.001 |

## 3. Regression: significant standardised coefficients (BH q<0.05; M2 vs M3 tells whether an effect survives within-family)

| tier | response | term | beta [95% CI] | q | R2 |
|---|---|---|---|---|---|
| M2 extended | Low-band fraction | ln Bond reduced mass | +0.40 [+0.33, +0.48] | <0.001 | 0.36 |
| M2 extended | Mapping efficiency eta | ln N_atom (primitive) | -0.38 [-0.47, -0.29] | <0.001 | 0.30 |
| M3 extended + family FE | Mapping efficiency eta | ln N_atom (primitive) | -0.37 [-0.46, -0.28] | <0.001 | 0.31 |
| M2 extended | High-band fraction | ln Bond reduced mass | -0.31 [-0.39, -0.24] | <0.001 | 0.24 |
| M2 extended | ln W1 distance | Bond electroneg. diff. | -0.30 [-0.43, -0.17] | <0.001 | 0.22 |
| M3 extended + family FE | ln Mean linewidth gamma | Electroneg. dispersion | -0.30 [-0.44, -0.15] | <0.001 | 0.12 |
| M2 extended | Mapping efficiency eta | Degeneracy fraction | +0.27 [+0.18, +0.36] | <0.001 | 0.30 |
| M3 extended + family FE | Mapping efficiency eta | Degeneracy fraction | +0.27 [+0.18, +0.35] | <0.001 | 0.31 |
| M2 extended | ln N_peak | ln Bond reduced mass | -0.26 [-0.33, -0.19] | <0.001 | 0.29 |
| M2 extended | ln W1 distance | ln Bond reduced mass | +0.25 [+0.16, +0.35] | <0.001 | 0.22 |
| M2 extended | ln Median peak spacing | Electroneg. dispersion | -0.25 [-0.41, -0.10] | 0.004 | 0.06 |
| M3 extended + family FE | ln W1 distance | Bond electroneg. diff. | -0.25 [-0.37, -0.13] | <0.001 | 0.29 |
| M2 extended | Peak overlap ratio | ln Bond reduced mass | +0.25 [+0.16, +0.34] | <0.001 | 0.07 |
| M3 extended + family FE | ln Mean linewidth gamma | ln Bond reduced mass | +0.24 [+0.13, +0.34] | <0.001 | 0.12 |
| M2 extended | ln W1 distance | Electroneg. dispersion | +0.23 [+0.08, +0.38] | 0.009 | 0.22 |
| M3 extended + family FE | ln N_peak | ln Bond reduced mass | -0.22 [-0.30, -0.14] | <0.001 | 0.37 |
| M3 extended + family FE | Low-band fraction | ln Bond reduced mass | +0.22 [+0.14, +0.30] | <0.001 | 0.45 |
| M2 extended | ln Mean linewidth gamma | Bond electroneg. diff. | +0.22 [+0.09, +0.34] | 0.003 | 0.07 |
| M2 extended | ln Median peak spacing | Bond electroneg. diff. | +0.21 [+0.08, +0.34] | 0.004 | 0.06 |
| M2 extended | ln N_peak | Electroneg. dispersion | +0.21 [+0.07, +0.35] | 0.008 | 0.29 |

## 4. Out-of-fold R2 [95% cluster-bootstrap CI], grouped CV

| spectrum | Measurement (SNR, wavelength) | Categorical (family, crystal system) | Size / complexity | Symmetry / group theory | Bond geometry | Chemistry | Local environment (SOAP) | All structure/chemistry |
|---|---|---|---|---|---|---|---|---|
| ln N_peak | 0.16 [0.12, 0.19] | 0.31 [0.27, 0.35] | 0.17 [0.13, 0.21] | 0.21 [0.17, 0.25] | 0.31 [0.27, 0.35] | 0.30 [0.25, 0.34] | 0.18 [0.14, 0.22] | 0.38 [0.35, 0.42] |
| Peak-intensity evenness | 0.17 [0.14, 0.20] | 0.06 [0.03, 0.08] | 0.03 [0.00, 0.05] | 0.03 [-0.00, 0.06] | 0.07 [0.04, 0.10] | 0.09 [0.07, 0.12] | 0.03 [0.01, 0.06] | 0.09 [0.05, 0.12] |
| ln Mean linewidth gamma | 0.08 [0.06, 0.11] | 0.01 [-0.00, 0.04] | 0.01 [-0.01, 0.03] | 0.06 [0.03, 0.08] | 0.11 [0.08, 0.14] | 0.19 [0.15, 0.23] | 0.02 [-0.01, 0.04] | 0.26 [0.22, 0.30] |
| ln Median peak spacing | 0.04 [0.01, 0.08] | 0.12 [0.07, 0.15] | 0.03 [-0.01, 0.07] | 0.11 [0.06, 0.16] | 0.10 [0.06, 0.15] | 0.07 [0.05, 0.10] | 0.02 [-0.01, 0.05] | 0.18 [0.13, 0.23] |
| Peak overlap ratio | 0.07 [0.04, 0.11] | 0.08 [0.06, 0.12] | 0.03 [0.00, 0.06] | 0.05 [0.02, 0.09] | 0.08 [0.05, 0.12] | 0.12 [0.08, 0.15] | 0.01 [-0.02, 0.04] | 0.17 [0.13, 0.20] |
| Low-band fraction | 0.05 [0.01, 0.09] | 0.38 [0.33, 0.44] | 0.07 [0.04, 0.10] | 0.10 [0.06, 0.14] | 0.38 [0.32, 0.44] | 0.38 [0.32, 0.43] | 0.08 [0.04, 0.12] | 0.39 [0.33, 0.45] |
| Mid-band fraction | -0.00 [-0.02, 0.02] | 0.18 [0.14, 0.21] | -0.02 [-0.04, 0.00] | 0.05 [0.02, 0.08] | 0.14 [0.10, 0.18] | 0.11 [0.07, 0.14] | 0.01 [-0.02, 0.04] | 0.18 [0.14, 0.21] |
| High-band fraction | 0.02 [-0.01, 0.04] | 0.30 [0.25, 0.35] | 0.08 [0.04, 0.12] | 0.08 [0.04, 0.13] | 0.30 [0.25, 0.34] | 0.27 [0.22, 0.32] | 0.08 [0.04, 0.11] | 0.34 [0.29, 0.39] |
| ln W1 distance | 0.11 [0.08, 0.14] | 0.24 [0.20, 0.28] | 0.09 [0.05, 0.12] | 0.09 [0.05, 0.13] | 0.27 [0.22, 0.31] | 0.27 [0.23, 0.31] | 0.09 [0.05, 0.12] | 0.29 [0.24, 0.33] |

### Added value of structure/chemistry over baselines (paired dR2)

| spectrum | All - Categorical | All - Measurement |
|---|---|---|
| ln N_peak | +0.07 [+0.04, +0.10] | +0.23 [+0.18, +0.28] |
| Peak-intensity evenness | +0.03 [+0.01, +0.06] | -0.08 [-0.13, -0.03] |
| ln Mean linewidth gamma | +0.24 [+0.20, +0.29] | +0.18 [+0.13, +0.22] |
| ln Median peak spacing | +0.06 [+0.02, +0.10] | +0.14 [+0.08, +0.19] |
| Peak overlap ratio | +0.08 [+0.05, +0.12] | +0.10 [+0.05, +0.14] |
| Low-band fraction | +0.01 [-0.03, +0.04] | +0.34 [+0.28, +0.40] |
| Mid-band fraction | +0.00 [-0.04, +0.04] | +0.18 [+0.13, +0.22] |
| High-band fraction | +0.04 [+0.01, +0.07] | +0.32 [+0.27, +0.38] |
| ln W1 distance | +0.05 [+0.02, +0.08] | +0.19 [+0.13, +0.23] |

### Leave-one-family-out extrapolation (all features)

| target | lofo_skill_pooled | lofo_r2_pooled | grouped_cv_r2_same_rows | generalisation gap |
|---|---|---|---|---|
| ln N_peak | 0.30 | 0.21 | 0.37 | 0.16 |
| Peak-intensity evenness | 0.04 | 0.01 | 0.10 | 0.09 |
| ln Mean linewidth gamma | 0.17 | 0.17 | 0.25 | 0.09 |
| ln Median peak spacing | -0.00 | -0.04 | 0.18 | 0.22 |
| Peak overlap ratio | 0.04 | 0.02 | 0.17 | 0.15 |
| Low-band fraction | 0.35 | 0.27 | 0.40 | 0.12 |
| Mid-band fraction | 0.03 | -0.01 | 0.18 | 0.19 |
| High-band fraction | 0.29 | 0.22 | 0.34 | 0.13 |
| ln W1 distance | 0.22 | 0.16 | 0.30 | 0.15 |

## 5. Sparse relations (nested-CV Lasso; percentages = bootstrap selection frequency)

- **ln N_peak** (R2_cv=0.27 [0.24, 0.32], 17 terms): `z(ln N_peak) ~ +0.14*[Electroneg. dispersion](97%) -0.11*[Mean bond length](100%) -0.11*[Raman irrep entropy](100%) +0.09*[ln N_Raman (group theory)](89%) -0.09*[ln Bond reduced mass](95%) +0.07*[ln N_atom (primitive)](78%)`
- **Peak-intensity evenness** (R2_cv=0.02 [0.01, 0.03], 5 terms): `z(Peak-intensity evenness) ~ -0.08*[ln Bond reduced mass](98%) +0.08*[ln Wyckoff orbits](77%) -0.02*[Degeneracy fraction](66%) +0.01*[Polyhedral polymerisation](62%) +0.00*[ln SOAP N environments](44%)`
- **ln Mean linewidth gamma** (R2_cv=0.15 [0.13, 0.18], 5 terms): `z(ln Mean linewidth gamma) ~ +0.20*[Site-mixing entropy](100%) +0.11*[ln Wyckoff orbits](100%) -0.11*[Oxyanion fraction](100%) +0.05*[ln Bond reduced mass](96%) +0.04*[Mean coordination no.](92%)`
- **ln Median peak spacing** (R2_cv=0.03 [0.02, 0.04], 16 terms): `z(ln Median peak spacing) ~ -0.13*[ln N_Raman (group theory)](99%) +0.11*[Raman irrep entropy](100%) +0.09*[Missingness: ln SOAP min env. separation](99%) +0.08*[Max Raman irrep mult.](98%) +0.07*[Degeneracy fraction](100%) -0.07*[Missingness: Oxyanion fraction](98%)`
- **Peak overlap ratio** (R2_cv=0.03 [0.02, 0.04], 4 terms): `z(Peak overlap ratio) ~ +0.05*[ln Wyckoff orbits](85%) +0.02*[ln SOAP N environments](61%) +0.01*[Site-mixing entropy](69%) +0.00*[ln Bond reduced mass](50%)`
- **Low-band fraction** (R2_cv=0.35 [0.30, 0.39], 7 terms): `z(Low-band fraction) ~ +0.24*[Missingness: Oxyanion fraction](100%) +0.18*[ln Bond reduced mass](100%) -0.14*[Electroneg. dispersion](100%) -0.06*[N anion-group types](93%) -0.02*[ln Wyckoff orbits](61%) -0.00*[Bond electroneg. diff.](51%)`
- **Mid-band fraction** (R2_cv=0.11 [0.08, 0.14], 12 terms): `z(Mid-band fraction) ~ -0.13*[Oxyanion fraction](100%) -0.13*[Mean bond length](100%) +0.12*[Electroneg. dispersion](100%) +0.12*[Mean coordination no.](100%) -0.08*[Missingness: Oxyanion fraction](97%) +0.07*[Polyhedral polymerisation](100%)`
- **High-band fraction** (R2_cv=0.27 [0.22, 0.31], 10 terms): `z(High-band fraction) ~ +0.14*[Bond electroneg. diff.](100%) -0.12*[ln Bond reduced mass](94%) -0.12*[Missingness: Polyhedral polymerisation](100%) +0.11*[N anion-group types](99%) -0.11*[Missingness: Oxyanion fraction](98%) +0.07*[ln Wyckoff orbits](83%)`
- **ln W1 distance** (R2_cv=0.22 [0.17, 0.25], 10 terms): `z(ln W1 distance) ~ +0.19*[Missingness: Polyhedral polymerisation](100%) -0.15*[Bond electroneg. diff.](100%) +0.09*[Mean bond length](99%) +0.07*[Missingness: Oxyanion fraction](98%) -0.06*[ln SOAP N environments](78%) -0.05*[ln Wyckoff orbits](84%)`

## 6. Joint group-permutation importance (dR2)

| group | ln N_peak | Peak-intensity evenness | ln Mean linewidth gamma | ln Median peak spacing | Peak overlap ratio | Low-band fraction | Mid-band fraction | High-band fraction | ln W1 distance |
|---|---|---|---|---|---|---|---|---|---|
| Size / complexity | 0.06 | 0.00 | 0.01 | 0.01 | 0.00 | 0.00 | 0.01 | 0.01 | 0.01 |
| Symmetry / group theory | 0.08 | 0.03 | 0.04 | 0.17 | 0.08 | 0.01 | 0.04 | 0.01 | 0.01 |
| Bond geometry | 0.18 | 0.10 | 0.16 | 0.08 | 0.13 | 0.29 | 0.23 | 0.17 | 0.21 |
| Chemistry | 0.12 | 0.06 | 0.24 | 0.11 | 0.12 | 0.15 | 0.09 | 0.32 | 0.15 |
| Local environment (SOAP) | 0.01 | 0.01 | 0.02 | 0.01 | 0.01 | 0.01 | 0.03 | 0.00 | 0.02 |

## 7. PLS cumulative Q2 (grouped CV)
1 LV = 0.125, 2 LV = 0.149, 3 LV = 0.167, 4 LV = 0.173, 5 LV = 0.175, 6 LV = 0.182

## 8. Group-theory ceiling
Spearman rho(N_Raman, N_peak) = 0.33.

## 9. Sensitivity of partial correlations to QC thresholds

| spec | n | n_strong_pairs_main | frac_same_sign | frac_same_sign_and_q_lt_0p05 | spearman_of_all_rho_vs_main |
|---|---|---|---|---|---|
| SNR>=10 | 1864 | 32 | 1.00 | 1.00 | 1.00 |
| SNR>=40 | 1781 | 32 | 1.00 | 1.00 | 1.00 |
| peaks>=2 | 1872 | 32 | 1.00 | 1.00 | 1.00 |
| peaks>=5 | 1720 | 32 | 1.00 | 1.00 | 0.99 |
