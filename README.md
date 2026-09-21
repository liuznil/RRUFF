# Structural Determinants of Raman Spectral Resolvability in Crystalline Solids

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22871576.svg)](https://doi.org/10.5281/zenodo.22871576) <!-- 发布后替换为真实 Zenodo DOI -->

This repository contains the complete, reproducible data science workflow for the study: **"How much of a mineral's Raman spectrum can crystal structure and chemistry predict? A grouped-validation study"**, submitted to *Patterns*. 

We provide the code for structural descriptor engineering, rigorous machine learning validation (Nested Grouped CV, LOFO-CV), and statistical analysis of structure-spectrum relationships across 1,842 high-quality Raman spectra from the RRUFF database.

## 📂 Directory Structure

```text
├── data/                       # Placeholder for raw data
│   ├── AMCSD                   # Placeholder for AMCSD
│   ├── raman                   # Placeholder for raman data 
├── results/                    # Output directory for models, metrics, and figures
│   ├── figures/                
│   ├── figures_audit/          
│   ├── regression_reports/     
│   ├── tables/                 
│   ├── dataset_failures.csv    # Data parse failed
│   └── dataset.csv             # Data parse success
├── scripts/                    # Executable scripts to reproduce main results
│   ├── build_dataset.py
│   ├── data_analysis.py
│   ├── raman_features.py
│   └── structure_features.py
├── README.md                    # This file
└── requirements.txt             # Environment specification
