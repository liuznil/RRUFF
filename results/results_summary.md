# Auto-generated results summary

All numbers below come from this run (seed=0); interpret with domain knowledge.

## 0. Cohort
n = 1838 spectra after QC; grouping unit for CV / clustering: `formula`.

| step | n |
|---|---|
| All samples | 2043 |
| Spectrum quality in {excellent,fair} | 1901 |
| SNR >= 20 | 1876 |
| Detected peaks >= 3 | 1841 |
| Orientation = unoriented | 1838 |

## Key findings (auto-generated from this run's numbers)

Hypothesis-driven pairs fixed in code (HEADLINE_PAIRS); report as exploratory if chosen after viewing results.

| structure | spectrum | partial rho [95% CI] | q |
|---|---|---|---|
| ln Bond reduced mass | Low-band fraction | +0.33 [+0.29, +0.37] | <0.001 |
| ln Bond reduced mass | ln N_peak | -0.33 [-0.37, -0.29] | <0.001 |
| ln Bond reduced mass | High-band fraction | -0.32 [-0.37, -0.28] | <0.001 |
| Bond electroneg. diff. | ln W1 distance | -0.20 [-0.24, -0.15] | <0.001 |
| Site-mixing entropy | ln Median linewidth gamma | +0.28 [+0.23, +0.32] | <0.001 |
| Electroneg. dispersion | High-band fraction | +0.26 [+0.21, +0.30] | <0.001 |
| N anion-group types | ln N_peak | +0.28 [+0.23, +0.32] | <0.001 |
| Polyhedral polymerisation | Spectral centroid | -0.37 [-0.42, -0.32] | <0.001 |

- Group-theoretical N_Raman and observed N_peak show a weak monotonic association (Spearman rho = 0.33).
- Structural and composition-derived descriptors add out-of-fold information beyond family + crystal system (paired dR2 95% CI > 0) for: ln N_peak (+0.07 [+0.04, +0.10]); ln Median linewidth gamma (+0.23 [+0.20, +0.27]); ln Median peak spacing (+0.05 [+0.01, +0.08]); Peak overlap ratio (+0.07 [+0.03, +0.11]); High-band fraction (+0.04 [+0.01, +0.07]); ln W1 distance (+0.05 [+0.02, +0.08]); Spectral centroid (+0.09 [+0.05, +0.12]); Spectral skewness (+0.09 [+0.06, +0.12]).
- No detectable gain over family + crystal system (CI includes or is below 0) for: Peak-intensity evenness, Low-band fraction, Mid-band fraction.
- Measurement conditions alone out-predict all structural and composition-derived descriptors for: Peak-intensity evenness (treat these endpoints as instrument-dominated).
- Endpoints with out-of-fold R2 whose 95% CI excludes 0: ln N_peak (R2=0.38), Peak-intensity evenness (R2=0.08), ln Median linewidth gamma (R2=0.25), ln Median peak spacing (R2=0.16), Peak overlap ratio (R2=0.16), Low-band fraction (R2=0.39), Mid-band fraction (R2=0.17), High-band fraction (R2=0.34), ln W1 distance (R2=0.29), Spectral centroid (R2=0.43), Spectral skewness (R2=0.44).
- Largest joint-permutation group per endpoint (dR2 >= 0.02): ln N_peak: Bond geometry (dR2=0.19); Peak-intensity evenness: Bond geometry (dR2=0.10); ln Median linewidth gamma: Composition-derived chemistry (dR2=0.25); ln Median peak spacing: Symmetry / group theory (dR2=0.14); Peak overlap ratio: Bond geometry (dR2=0.13); Low-band fraction: Bond geometry (dR2=0.35); Mid-band fraction: Bond geometry (dR2=0.19); High-band fraction: Composition-derived chemistry (dR2=0.32); ln W1 distance: Bond geometry (dR2=0.21); Spectral centroid: Composition-derived chemistry (dR2=0.37); Spectral skewness: Composition-derived chemistry (dR2=0.31).
- QC sensitivity across 4 alternative settings: sign agreement of main strong partial correlations 1.00-1.00; sign agreement with q<0.05 retained 1.00-1.00; rank agreement of the whole partial-correlation matrix 0.99-1.00.

## 1. Strongest partial correlations (adjusted for ln N_atom, ln SNR)

| structure | spectrum | partial rho [95% CI] | q | n |
|---|---|---|---|---|
| Raman-active fraction | Mapping efficiency eta | -0.48 [-0.51, -0.44] | <0.001 | 1829 |
| Polyhedral polymerisation | Spectral skewness | +0.41 [+0.36, +0.46] | <0.001 | 1244 |
| Centrosymmetric (0/1) | Mapping efficiency eta | +0.38 [+0.33, +0.42] | <0.001 | 1829 |
| Polyhedral polymerisation | Spectral centroid | -0.37 [-0.42, -0.32] | <0.001 | 1244 |
| ln Bond reduced mass | ln N_peak | -0.33 [-0.37, -0.29] | <0.001 | 1804 |
| ln Bond reduced mass | Low-band fraction | +0.33 [+0.29, +0.37] | <0.001 | 1804 |
| N anion-group types | High-band fraction | +0.33 [+0.28, +0.37] | <0.001 | 1838 |
| N anion-group types | ln W1 distance | -0.32 [-0.37, -0.28] | <0.001 | 1838 |
| ln Bond reduced mass | High-band fraction | -0.32 [-0.37, -0.28] | <0.001 | 1804 |
| Degeneracy fraction | Mapping efficiency eta | +0.32 [+0.28, +0.36] | <0.001 | 1829 |
| Electroneg. dispersion | Low-band fraction | -0.32 [-0.36, -0.28] | <0.001 | 1838 |
| Max Raman irrep mult. | Mapping efficiency eta | -0.31 [-0.36, -0.27] | <0.001 | 1829 |
| N anion-group types | Spectral centroid | +0.30 [+0.26, +0.34] | <0.001 | 1838 |
| ln Wyckoff orbits | Mapping efficiency eta | -0.30 [-0.34, -0.25] | <0.001 | 1829 |
| Bond electroneg. diff. | Low-band fraction | -0.28 [-0.33, -0.24] | <0.001 | 1804 |

## 2. Categorical factors (top 10 eta2_H)

| factor | spectrum | eta2_H | q |
|---|---|---|---|
| family | Spectral skewness | 0.33 | <0.001 |
| family | Spectral centroid | 0.33 | <0.001 |
| family | ln N_peak | 0.27 | <0.001 |
| family | High-band fraction | 0.25 | <0.001 |
| family | Low-band fraction | 0.25 | <0.001 |
| family | ln W1 distance | 0.23 | <0.001 |
| family | Mid-band fraction | 0.15 | <0.001 |
| n_anion_group_types_cat | Spectral centroid | 0.15 | <0.001 |
| n_anion_group_types_cat | High-band fraction | 0.14 | <0.001 |
| n_anion_group_types_cat | ln W1 distance | 0.13 | <0.001 |

## 3. Regression: significant standardised coefficients (BH q<0.05; M2 vs M3 tells whether an effect survives within-family)

| tier | response | term | beta [95% CI] | q | R2 |
|---|---|---|---|---|---|
| M2 extended | Low-band fraction | ln Bond reduced mass | +0.40 [+0.33, +0.48] | <0.001 | 0.36 |
| M2 extended | Mapping efficiency eta | ln N_atom (primitive) | -0.38 [-0.47, -0.29] | <0.001 | 0.30 |
| M3 extended + family FE | Mapping efficiency eta | ln N_atom (primitive) | -0.37 [-0.46, -0.27] | <0.001 | 0.31 |
| M2 extended | High-band fraction | ln Bond reduced mass | -0.31 [-0.38, -0.23] | <0.001 | 0.24 |
| M3 extended + family FE | Spectral skewness | Electroneg. dispersion | -0.30 [-0.45, -0.15] | <0.001 | 0.42 |
| M3 extended + family FE | ln Median linewidth gamma | Electroneg. dispersion | -0.30 [-0.44, -0.15] | <0.001 | 0.12 |
| M2 extended | ln W1 distance | Bond electroneg. diff. | -0.30 [-0.43, -0.17] | <0.001 | 0.22 |
| M2 extended | Spectral centroid | ln Bond reduced mass | -0.28 [-0.35, -0.20] | <0.001 | 0.25 |
| M2 extended | Mapping efficiency eta | Degeneracy fraction | +0.27 [+0.18, +0.36] | <0.001 | 0.30 |
| M3 extended + family FE | Mapping efficiency eta | Degeneracy fraction | +0.27 [+0.18, +0.35] | <0.001 | 0.31 |
| M2 extended | ln N_peak | ln Bond reduced mass | -0.26 [-0.33, -0.19] | <0.001 | 0.29 |
| M2 extended | ln W1 distance | ln Bond reduced mass | +0.25 [+0.15, +0.35] | <0.001 | 0.22 |
| M2 extended | ln Median peak spacing | Electroneg. dispersion | -0.25 [-0.41, -0.10] | 0.004 | 0.06 |
| M3 extended + family FE | ln W1 distance | Bond electroneg. diff. | -0.25 [-0.37, -0.13] | <0.001 | 0.29 |
| M2 extended | Peak overlap ratio | ln Bond reduced mass | +0.25 [+0.15, +0.34] | <0.001 | 0.07 |
| M2 extended | Spectral skewness | ln Bond reduced mass | +0.24 [+0.17, +0.32] | <0.001 | 0.26 |
| M3 extended + family FE | ln N_peak | ln Bond reduced mass | -0.23 [-0.31, -0.15] | <0.001 | 0.37 |
| M3 extended + family FE | ln Median linewidth gamma | ln Bond reduced mass | +0.23 [+0.12, +0.33] | <0.001 | 0.12 |
| M2 extended | ln W1 distance | Electroneg. dispersion | +0.23 [+0.08, +0.38] | 0.008 | 0.22 |
| M3 extended + family FE | Low-band fraction | ln Bond reduced mass | +0.22 [+0.14, +0.30] | <0.001 | 0.45 |

## 4. Out-of-fold R2 [95% cluster-bootstrap CI], grouped CV

| spectrum | Measurement (SNR, wavelength) | Categorical (family, crystal system) | Size / complexity | Symmetry / group theory | Bond geometry | Composition-derived chemistry | Local environment (SOAP) | All structure/chemistry |
|---|---|---|---|---|---|---|---|---|
| ln N_peak | 0.16 [0.12, 0.20] | 0.31 [0.27, 0.35] | 0.17 [0.13, 0.20] | 0.20 [0.16, 0.24] | 0.31 [0.27, 0.35] | 0.31 [0.26, 0.35] | 0.18 [0.14, 0.22] | 0.38 [0.35, 0.42] |
| Peak-intensity evenness | 0.16 [0.13, 0.19] | 0.06 [0.03, 0.09] | 0.03 [0.00, 0.06] | 0.02 [-0.02, 0.05] | 0.07 [0.04, 0.09] | 0.09 [0.06, 0.12] | 0.04 [0.01, 0.07] | 0.08 [0.05, 0.12] |
| ln Median linewidth gamma | 0.09 [0.06, 0.12] | 0.02 [0.00, 0.04] | 0.00 [-0.02, 0.02] | 0.05 [0.02, 0.08] | 0.10 [0.07, 0.13] | 0.19 [0.15, 0.23] | 0.02 [-0.01, 0.05] | 0.25 [0.22, 0.29] |
| ln Median peak spacing | 0.04 [0.01, 0.08] | 0.12 [0.08, 0.16] | 0.03 [-0.00, 0.07] | 0.10 [0.05, 0.14] | 0.10 [0.06, 0.14] | 0.08 [0.06, 0.11] | 0.01 [-0.01, 0.03] | 0.16 [0.12, 0.21] |
| Peak overlap ratio | 0.07 [0.04, 0.10] | 0.09 [0.05, 0.12] | 0.02 [-0.01, 0.05] | 0.04 [0.01, 0.07] | 0.09 [0.05, 0.12] | 0.12 [0.09, 0.15] | 0.02 [-0.01, 0.04] | 0.16 [0.12, 0.20] |
| Low-band fraction | 0.04 [0.01, 0.07] | 0.38 [0.32, 0.43] | 0.08 [0.05, 0.11] | 0.10 [0.06, 0.14] | 0.38 [0.32, 0.44] | 0.38 [0.31, 0.43] | 0.08 [0.05, 0.12] | 0.39 [0.33, 0.45] |
| Mid-band fraction | -0.01 [-0.03, 0.01] | 0.17 [0.13, 0.21] | -0.02 [-0.04, 0.00] | 0.05 [0.02, 0.08] | 0.14 [0.10, 0.18] | 0.11 [0.08, 0.15] | 0.01 [-0.02, 0.04] | 0.17 [0.14, 0.22] |
| High-band fraction | 0.02 [-0.01, 0.04] | 0.30 [0.25, 0.34] | 0.08 [0.05, 0.12] | 0.08 [0.04, 0.12] | 0.30 [0.24, 0.34] | 0.28 [0.23, 0.33] | 0.09 [0.05, 0.12] | 0.34 [0.28, 0.38] |
| ln W1 distance | 0.11 [0.08, 0.14] | 0.24 [0.20, 0.29] | 0.09 [0.06, 0.12] | 0.09 [0.06, 0.13] | 0.26 [0.22, 0.30] | 0.27 [0.23, 0.32] | 0.09 [0.06, 0.13] | 0.29 [0.25, 0.34] |
| Spectral centroid | 0.03 [0.00, 0.06] | 0.34 [0.30, 0.37] | 0.08 [0.05, 0.12] | 0.11 [0.08, 0.15] | 0.32 [0.28, 0.36] | 0.36 [0.31, 0.40] | 0.10 [0.06, 0.13] | 0.43 [0.39, 0.47] |
| Spectral skewness | 0.05 [0.01, 0.08] | 0.35 [0.31, 0.38] | 0.09 [0.05, 0.12] | 0.13 [0.09, 0.16] | 0.33 [0.29, 0.37] | 0.39 [0.35, 0.43] | 0.10 [0.06, 0.12] | 0.44 [0.40, 0.47] |

### Added value of structure/chemistry over baselines (paired dR2)

| spectrum | All - Categorical | All - Measurement |
|---|---|---|
| ln N_peak | +0.07 [+0.04, +0.10] | +0.23 [+0.17, +0.29] |
| Peak-intensity evenness | +0.02 [-0.00, +0.05] | -0.08 [-0.13, -0.04] |
| ln Median linewidth gamma | +0.23 [+0.20, +0.27] | +0.17 [+0.12, +0.21] |
| ln Median peak spacing | +0.05 [+0.01, +0.08] | +0.12 [+0.07, +0.17] |
| Peak overlap ratio | +0.07 [+0.03, +0.11] | +0.09 [+0.05, +0.13] |
| Low-band fraction | +0.01 [-0.03, +0.05] | +0.35 [+0.29, +0.41] |
| Mid-band fraction | -0.00 [-0.04, +0.04] | +0.18 [+0.14, +0.23] |
| High-band fraction | +0.04 [+0.01, +0.07] | +0.32 [+0.26, +0.37] |
| ln W1 distance | +0.05 [+0.02, +0.08] | +0.19 [+0.13, +0.24] |
| Spectral centroid | +0.09 [+0.05, +0.12] | +0.40 [+0.35, +0.44] |
| Spectral skewness | +0.09 [+0.06, +0.12] | +0.39 [+0.33, +0.44] |

### Leave-one-family-out extrapolation (all features)

| target | lofo_skill_pooled | lofo_r2_pooled | grouped_cv_r2_same_rows | generalisation gap |
|---|---|---|---|---|
| ln N_peak | 0.30 | 0.21 | 0.37 | 0.16 |
| Peak-intensity evenness | 0.05 | 0.02 | 0.10 | 0.08 |
| ln Median linewidth gamma | 0.17 | 0.16 | 0.25 | 0.09 |
| ln Median peak spacing | -0.00 | -0.04 | 0.17 | 0.21 |
| Peak overlap ratio | 0.02 | -0.01 | 0.16 | 0.16 |
| Low-band fraction | 0.35 | 0.28 | 0.40 | 0.12 |
| Mid-band fraction | 0.02 | -0.02 | 0.18 | 0.20 |
| High-band fraction | 0.28 | 0.21 | 0.34 | 0.13 |
| ln W1 distance | 0.22 | 0.15 | 0.31 | 0.15 |
| Spectral centroid | 0.33 | 0.26 | 0.43 | 0.17 |
| Spectral skewness | 0.32 | 0.25 | 0.44 | 0.19 |

## 5. Sparse relations (nested-CV Lasso; percentages = bootstrap selection frequency)

- **ln N_peak** (R2_cv=0.27 [0.23, 0.30], 13 terms): `z(ln N_peak) ~ -0.14*[ln Bond reduced mass](99%) +0.13*[Electroneg. dispersion](99%) +0.08*[ln SOAP N environments](81%) +0.07*[N anion-group types](100%) -0.07*[Mean bond length](98%) +0.07*[ln N_Raman (group theory)](78%)`
- **Peak-intensity evenness** (R2_cv=0.03 [0.01, 0.03], 2 terms): `z(Peak-intensity evenness) ~ +0.02*[ln Wyckoff orbits](54%) -0.02*[ln Bond reduced mass](70%)`
- **ln Median linewidth gamma** (R2_cv=0.16 [0.12, 0.18], 8 terms): `z(ln Median linewidth gamma) ~ +0.22*[Site-mixing entropy](100%) +0.17*[ln Wyckoff orbits](98%) -0.13*[Oxyanion fraction](100%) +0.13*[ln Bond reduced mass](100%) +0.07*[Mean coordination no.](100%) -0.05*[Missingness: Oxyanion fraction](91%)`
- **ln Median peak spacing** (R2_cv=0.03 [0.02, 0.04], 2 terms): `z(ln Median peak spacing) ~ +0.05*[Degeneracy fraction](95%) +0.05*[Missingness: ln SOAP min env. separation](86%)`
- **Peak overlap ratio** (R2_cv=0.02 [0.01, 0.03], 4 terms): `z(Peak overlap ratio) ~ +0.14*[ln SOAP N environments](91%) +0.07*[ln Bond reduced mass](97%) +0.05*[Site-mixing entropy](98%) +0.00*[Missingness: Polyhedral polymerisation](50%)`
- **Low-band fraction** (R2_cv=0.36 [0.31, 0.40], 10 terms): `z(Low-band fraction) ~ +0.24*[Missingness: Oxyanion fraction](100%) +0.18*[ln Bond reduced mass](100%) -0.13*[Electroneg. dispersion](99%) -0.06*[N anion-group types](98%) -0.02*[ln Wyckoff orbits](49%) -0.01*[Bond electroneg. diff.](62%)`
- **Mid-band fraction** (R2_cv=0.09 [0.06, 0.11], 5 terms): `z(Mid-band fraction) ~ +0.10*[Electroneg. dispersion](100%) -0.08*[Missingness: Oxyanion fraction](97%) -0.06*[Mean bond length](95%) +0.05*[Mean coordination no.](97%) -0.04*[Oxyanion fraction](94%)`
- **High-band fraction** (R2_cv=0.27 [0.23, 0.31], 7 terms): `z(High-band fraction) ~ +0.12*[Bond electroneg. diff.](100%) -0.11*[Missingness: Oxyanion fraction](100%) -0.11*[Missingness: Polyhedral polymerisation](100%) -0.11*[ln Bond reduced mass](95%) +0.08*[N anion-group types](98%) +0.06*[ln SOAP N environments](78%)`
- **ln W1 distance** (R2_cv=0.24 [0.20, 0.27], 9 terms): `z(ln W1 distance) ~ +0.20*[Missingness: Polyhedral polymerisation](100%) -0.17*[Bond electroneg. diff.](100%) -0.13*[ln SOAP N environments](90%) +0.11*[Mean bond length](100%) +0.06*[Missingness: Oxyanion fraction](98%) +0.05*[Polyhedral polymerisation](98%)`
- **Spectral centroid** (R2_cv=0.35 [0.31, 0.38], 13 terms): `z(Spectral centroid) ~ -0.27*[Polyhedral polymerisation](100%) -0.25*[Missingness: Polyhedral polymerisation](100%) +0.16*[Bond electroneg. diff.](100%) -0.14*[ln Bond reduced mass](94%) -0.11*[Missingness: Oxyanion fraction](100%) +0.08*[ln SOAP N environments](97%)`
- **Spectral skewness** (R2_cv=0.35 [0.31, 0.38], 7 terms): `z(Spectral skewness) ~ +0.30*[Missingness: Oxyanion fraction](100%) +0.18*[Missingness: Polyhedral polymerisation](100%) +0.16*[Polyhedral polymerisation](100%) -0.06*[Bond electroneg. diff.](87%) +0.02*[ln Bond reduced mass](56%) +0.00*[Degeneracy fraction](49%)`

## 6. Joint group-permutation importance (dR2)

| group | ln N_peak | Peak-intensity evenness | ln Median linewidth gamma | ln Median peak spacing | Peak overlap ratio | Low-band fraction | Mid-band fraction | High-band fraction | ln W1 distance | Spectral centroid | Spectral skewness |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Size / complexity | 0.03 | 0.01 | 0.01 | 0.01 | -0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | -0.00 |
| Symmetry / group theory | 0.07 | 0.03 | 0.04 | 0.14 | 0.07 | 0.02 | 0.02 | 0.01 | 0.01 | 0.02 | 0.02 |
| Bond geometry | 0.18 | 0.10 | 0.14 | 0.06 | 0.13 | 0.35 | 0.19 | 0.15 | 0.21 | 0.23 | 0.22 |
| Composition-derived chemistry | 0.10 | 0.05 | 0.26 | 0.09 | 0.11 | 0.12 | 0.08 | 0.32 | 0.15 | 0.37 | 0.31 |
| Local environment (SOAP) | 0.02 | 0.01 | 0.02 | 0.01 | 0.02 | 0.01 | 0.02 | 0.00 | 0.03 | 0.02 | 0.01 |

## 7. PLS cumulative Q2 (grouped CV)
1 LV = 0.147, 2 LV = 0.176, 3 LV = 0.189, 4 LV = 0.205, 5 LV = 0.210, 6 LV = 0.216

## 8. Group-theory ceiling
Spearman rho(N_Raman, N_peak) = 0.33.

## 9. Sensitivity of partial correlations to QC thresholds

| spec | n | n_strong_pairs_main | frac_same_sign | frac_same_sign_and_q_lt_0p05 | spearman_of_all_rho_vs_main |
|---|---|---|---|---|---|
| SNR>=10 | 1860 | 41 | 1.00 | 1.00 | 1.00 |
| SNR>=40 | 1777 | 41 | 1.00 | 1.00 | 1.00 |
| peaks>=2 | 1868 | 41 | 1.00 | 1.00 | 1.00 |
| peaks>=5 | 1716 | 41 | 1.00 | 1.00 | 0.99 |
