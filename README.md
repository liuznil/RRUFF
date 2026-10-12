# Prediction of global Raman spectral features in minerals from structure and composition

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22899440.svg)](https://doi.org/10.5281/zenodo.23318733)

This repository contains the complete, reproducible data science workflow for the study: **"Prediction of global Raman spectral features in minerals from structure and composition"**, submitted to *Physical Chemistry Chemical Physics*.

The workflow covers structural descriptor engineering, Raman endpoint extraction,
grouped machine-learning validation, leave-one-family-out (LOFO) extrapolation
tests, nested Lasso selection, permutation/drop-group attribution, and
statistical analysis of structure--spectrum relationships across 1,838
quality-controlled Raman spectra from the RRUFF database.

## Reproducibility

Install the dependencies listed in `requirements.txt`, place the RRUFF Raman
files and AMCSD CIF files under the documented `data/` layout, and run:

```bash
python scripts/build_dataset.py --root . --jobs 1
python scripts/data_analysis.py --csv results/dataset.csv --out results
```

For a quick validation run:

```bash
python scripts/data_analysis.py --csv results/dataset.csv --out results_fast --fast --skip-sensitivity
```

The dataset builder records the selected AMCSD and Raman identifiers and rejects
ambiguous RRUFF identifiers that occur under more than one normalized mineral
name. Rejected candidates are recorded in `results/dataset_failures.csv` rather
than being silently paired. This protects the merged dataset from unverified
polymorph or composition matches, including duplicate-ID cases such as R230009.

To reproduce the additional analyses addressing the reviewer comments, run:

```bash
python scripts/structure_features.py reviewer-supplement --csv results/dataset.csv --out results
```

This uses repeated five-fold formula-grouped validation and 300 formula-cluster
bootstrap samples by default. It writes endpoint equations, representative
mineral mode-assignment examples, a complete drop-group ablation table, model
comparisons (including combined categorical/numeric and nearest-neighbour
baselines), descriptor-applicability and missingness sensitivities, richer SOAP
and multi-output comparisons, alternative preprocessing results, and a
training-fold empirical-reference sensitivity for the W1 endpoint under
`results/tables/`. The submission-ready Word supporting information is saved as
`results/supplementary.docx`.

Raman features are calculated after sorting and restricting spectra to
100--1300 cm^-1. The workflow applies rolling-minimum baseline correction,
Savitzky--Golay smoothing, robust noise thresholds, peak detection, and
half-prominence FWHM estimation. In addition to peak and band descriptors, it
calculates Wasserstein distance to a fixed uniform reference over the full
fingerprint window, spectral centroid, and spectral skewness.

All formula-grouped predictive splits keep repeated chemical formulas within a
single fold. Nested Lasso preprocessing is fitted independently within each
inner training fold: imputation, missingness-indicator selection, and scaling
are never learned from an inner validation fold. The analysis also writes a
complete-case grouped-CV sensitivity table to quantify the dependence of
performance on descriptor availability. Structure extraction records block-level
availability and extraction failures, including the explicitly non-applicable
SOAP case for non-idealised cells.

## 📂 Directory Structure

```text
├── data/                       # Placeholder for raw data
│   ├── AMCSD/                  # AMCSD CIF files
│   │   └── cif/
│   ├── raman/                  # Processed Raman text files
│   │   ├── excellent_unoriented/
│   │   ├── excellent_oriented/
│   │   ├── fair_unoriented/
│   │   └── ...
├── results/                    # Output directory for models, metrics, and figures
│   ├── figures/
│   │   ├── main/
│   │   └── supplementary/
│   ├── regression_reports/
│   ├── tables/
│   ├── dataset_failures.csv    # Dropped or rejected matches
│   ├── dataset.csv             # Merged structure--Raman dataset
│   └── supplementary.docx     # PCCP supporting information
├── scripts/                    # Executable scripts to reproduce main results
│   ├── build_dataset.py
│   ├── data_analysis.py
│   ├── reviewer_supplement.py  # Compatibility entry point
│   ├── raman_features.py
│   └── structure_features.py   # Structure extraction and reviewer analyses
├── README.md                    # This file
└── requirements.txt             # Environment specification
```

## Main outputs

`results/tables/` contains the reproducible tabular outputs, including:

- `table_cv_performance.csv`: grouped-CV performance for measurement,
  categorical, descriptor-group, and full models.
- `table_cv_paired_differences.csv`: paired incremental performance and
  drop-group ablation estimates.
- `table_perm_importance_group.csv` and
  `table_perm_importance_feature.csv`: held-out permutation importance.
- `table_lasso_performance.csv`, `table_lasso_coefficients.csv`, and
  `table_lasso_selection_frequency.csv`: nested Lasso results and stability
  selection.
- `table_sensitivity_qc.csv`: quality-control threshold sensitivity.
- `table_sensitivity_missingness.csv`: complete-case versus native-missingness
  grouped-CV sensitivity.
- `table_group_theory_ceiling.csv`: group-theoretical Raman mode counts versus
  experimentally detected peaks in the 100--1300 cm^-1 window.
- `table_endpoint_equations.csv`, `table_drop_group_complete.csv`, and
  `table_mode_assignment_examples.csv`: endpoint definitions, ablation
  estimates with confidence intervals, and illustrative mode-window cases.
- `table_reviewer_predictor_baselines.csv`,
  `table_reviewer_baseline_differences.csv`,
  `table_sensitivity_applicability.csv`,
  `table_feature_missingness_status.csv`,
  `table_soap_richer_comparison.csv`,
  `table_multioutput_comparison.csv`,
  `table_multioutput_differences.csv`,
  `table_sensitivity_preprocessing.csv`, and
  `table_w1_reference_sensitivity.csv`: supplementary reviewer-response
  comparisons and sensitivities.

The publication manuscript is in `results/main.tex`. It documents the
identifier-matching audit, spectral-window limitation, endpoint definitions,
fold-local preprocessing, attribution analyses, and the distinction between
measurement-only and structure/composition models.
