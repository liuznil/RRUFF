# Figure captions (descriptive; result statements are in paper_results.md)

## Figure 1
Study design. Raman spectra from the RRUFF Raman directory are linked to crystallographic and composition information from the RRUFF AMCSD directory; structural and composition-derived descriptors are then constructed for prediction. (a) Analysis framework; (b) quality-control funnel; (c) median standardised spectral profile per chemical family.

## Figure 2
Group-theoretical mode availability versus observed peaks. (a) N_Raman against N_peak with the 1:1 line (all active modes resolved). (b, c) Partial Spearman association of symmetry descriptors with the peak-count residual given ln N_Raman, i.e. resolved peaks beyond what mode availability implies; the ratio eta = N_peak/N_Raman is deliberately not used because it is coupled to N_Raman by construction.

## Figure 3
Strongest bond-geometry associations with peak count, band fractions, and spectral distance. Points show observations; lines show decile-binned medians. Annotations report partial Spearman rho, 95% CI, and BH-FDR q.

## Figure 4
Strongest composition-derived associations with linewidth, spectral shape, and peak count. Weak or null exploratory pairs are retained in the full correlation table rather than the main figure.

## Figure 5
Out-of-fold R2 (grouped 5-fold CV by chemical formula) across the informative spectral endpoints. Measurement = SNR + excitation wavelength; categorical = chemical family + crystal system; eta is excluded because it is definitionally coupled to peak count and Raman-active mode count.

## Figure 6
Joint held-out permutation importance (change in R2 when an entire descriptor group is permuted jointly). Values quantify predictive information, not causal effect.

## Figure 7
Integrated findings: (a) incremental out-of-fold value beyond family + crystal system; (b) descriptor group most often ranked first by joint permutation; (c) measurement-only versus structure/chemistry predictive performance.

## Figure S1
Quality-control funnel and median standardised spectral profile per chemical family.

## Figure S2
Complete partial Spearman correlation matrix (adjusted for ln N_atom, ln SNR). * q<0.05, ** q<0.01, *** q<0.001 (BH-FDR over all cells).

## Figure S3
Kruskal-Wallis effect size eta^2_H of categorical factors on each spectral endpoint.

## Figure S4
Standardised OLS coefficients (95% CI, formula-cluster-robust SE) for the extended model (M2) and the extended model with chemical-family fixed effects (M3). Filled: BH q<0.05.

## Figure S5
Leave-one-family-out skill score relative to the training-set mean for each held-out family.

## Figure S6
(a) Lasso coefficients (top 18 by magnitude; 'Missingness' rows are indicator variables). (b) Cumulative grouped-CV Q2 of PLS. (c) Rank agreement of the partial-correlation matrix under alternative QC thresholds.
