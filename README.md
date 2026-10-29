# PDAC FDG-PET Molecular Analysis

## Overview

A reproducible Python-based bioinformatics analysis of molecular differences between two pancreatic ductal adenocarcinoma (PDAC) cell lines, PSN1 and HPAC, in the context of [18F]FDG PET imaging research.

The project integrates public cancer cell-line transcriptomics, proteomics, and proliferation-associated transcriptional signatures.

## Research objectives

- Compare glucose transport and metabolism-related RNA expression between PSN1 and HPAC.
- Compare selected protein abundances using Sanger proteomics.
- Calculate Hallmark E2F Targets and G2/M Checkpoint RNA signature scores.
- Evaluate SLC2A1 and HK2 RNA expression relative to proliferation-associated transcriptional signatures.
- Explore corresponding protein-level associations.
- Generate reproducible scientific figures.

## Data sources

- DepMap Public 26Q3: cell-line metadata and RNA expression.
- Sanger mass spectrometry proteomics (2022).
- MSigDB Hallmark gene sets, version 2026.2.Hs.

Raw datasets are not redistributed in this repository. Users should obtain the source data from their respective providers.

## Methods

Analyses were implemented in Python using pandas, NumPy, SciPy, and Matplotlib.

Hallmark E2F Targets and G2/M Checkpoint scores were calculated from RNA expression using gene-wise z-score normalization relative to reference PDAC cell lines, followed by averaging across matched signature genes.

RNA residuals were calculated by fitting separate linear regressions between each Hallmark score and SLC2A1 or HK2 RNA expression in reference PDAC cell lines, excluding PSN1 and HPAC.

Additional exploratory analyses evaluate GLUT1 and HK2 protein abundance relative to RNA-derived Hallmark scores.

## Repository structure

- `scripts/`: Python scripts for data processing, analysis, and visualization.
- `data/`: Local source datasets (excluded from GitHub).
- `results/`: Generated tables and figures (excluded from GitHub).

## Reproducibility

The analysis scripts are provided for transparency and reproducibility. Instructions for obtaining and preparing the source datasets will be expanded in a future update.

## Publication

Associated manuscript: in preparation.



Research project developed as part of a pancreatic cancer molecular imaging and bioinformatics investigation.

