# Euro Area Inflation Forecasting: Machine Learning vs. Econometric Models

Replication code for my master's thesis:  
**"Can Machine Learning Outperform Traditional Econometric Models? Evidence from Euro Area Inflation Forecasting"**

## Overview

This repository contains R code to forecast euro area **core HICP inflation** using a structured hierarchy of seven models, designed to decompose forecasting gains into contributions from *data richness* vs. *nonlinearity*. Forecast horizons are **h = 1, 3, and 12 months** over a recursive out-of-sample evaluation starting January 1999.

## Models

| # | Model | Type |
|---|-------|------|
| 1 | Random Walk | Benchmark |
| 2 | UCSV (Chan 2018 reparameterisation) | Benchmark |
| 3 | Dynamic Factor Model (DFM) | Traditional multivariate |
| 4 | Ridge Regression | Machine learning |
| 5 | LASSO | Machine learning |
| 6 | Random Forest | Machine learning |
| 7 | XGBoost | Machine learning |

## Data

- **EA-MD-QD dataset** (February 2026 version, [Zenodo record 18804061](https://zenodo.org/record/18804061)) — 47 monthly series, preprocessed via the EA-MD-QD MATLAB routine (`routine_data.m`, settings: `country=EA`, `frequency=M`, `transformation=TR2`)
- Supplementary series: Brent crude (FRED), TTF gas, CISS (ECB SDW), FAO food prices, NY Fed GSCPI, ECB SPF inflation expectations, ECB Wage Tracker

## Repository Structure
```
Master_Thesis/
├── data/               # Preprocessed input data (eadataM_NA_TR2.xlsx)
├── models/
│   ├── ucsv/           # UCSV Gibbs sampler (translated from Chan's MATLAB)
│   ├── dfm/            # Dynamic Factor Model
│   └── ml/             # Ridge, LASSO, Random Forest, XGBoost
├── evaluation/         # RMSFE computation and model comparison
├── figures/            # Output plots
└── README.md
```

## Requirements

- **R** (≥ 4.2) with packages: `here`, `stochvol`, `glmnet`, `randomForest`, `xgboost`
- **MATLAB** (optional, for data preprocessing via EA-MD-QD scripts)

## Usage

1. Download and preprocess the EA-MD-QD dataset using `routine_data.m`
2. Place the output `eadataM_NA_TR2.xlsx` in `data/`
3. Run each model script from the project root (`.Rproj`)

## References

- Bańbura & Bobeica (2023), *Does the Phillips curve help forecast euro area inflation?*
- Chan (2018), *Bayesian estimation of time-series models*
- Stock & Watson (2007), *Why has US inflation become harder to forecast?*
- Goulet Coulombe et al. (2022), *How is machine learning useful for macroeconomic forecasting?*
- Naghi et al. (2024), *Forecasting Dutch inflation using ML methods*

## Status

Work in progress — master's thesis, submission 2025/26.
