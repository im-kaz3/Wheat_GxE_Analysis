# Wheat GxE, Micro-Phenology, and Critical Window Alignment Pipeline

## Overview
This repository contains the complete "A to Z" analytical pipeline utilized in the manuscript to evaluate 200 wheat genotypes under chronic vegetative thermal stress. It encompasses everything from the raw data ingestion and variance partitioning to the novel **Critical Window Alignment (CWA)** framework in a single, unified codebase.

## Code Architecture (`01_Master_Pipeline_Public.R`)
This script executes the comprehensive biophysical and statistical analyses, moving progressively from raw spatial data to advanced predictive modeling.

### Core Modules:
1. **Spatial Normalization:** Extracts Best Linear Unbiased Estimators (BLUEs) via Restricted Maximum Likelihood (REML) mixed models (Alpha-Lattice & RCBD).
2. **GxE Variance Partitioning:** Evaluates temporal resilience and exact variance structures (Genotypic Repeatability) across optimal and heat-stressed environments.
3. **Categorical Trait Tagging & Phenological Clustering:** Calculates Peak Phase frequencies, AUTPC (Area Under the Tiller Progress Curve), and utilizes K-Means and PCA to isolate morphological strategies ("Fast" vs. "Slow" Establishers).
4. **Tiller Stability Index (TSI):** Mathematically de-risks spatial deployment to identify universal agronomic elites.
5. **Structural Equation Modeling (SEM):** Path analysis of yield components.
6. **Critical Window Alignment (CWA):** 
   - Transforms cumulative trait counts into discrete, non-overlapping interval gains (Δ).
   - Decomposes the $R^2$ of yield explained by the intervals using **Shapley values (LMG method)** to avoid order bias.
   - Executes **1,000 genotype-cluster bootstrap resamples** to map the "Decision Window" for terminal yield.

### Usage
- Ensure the raw datasets (`FinalDataYear1.xlsx` and `Weather_Master.xlsx`) are placed in the same directory as the script.
- The script uses `getwd()` to define the project root. Ensure your R session is working in this directory.
- Simply source or run `01_Master_Pipeline_Public.R`.
- All outputs (figures, matrices, and statistical tables) will be automatically generated inside a newly created `Results_Output` subdirectory.

## Dependencies
The pipeline is written in R (v4.0+) and requires the following libraries:
- **Data Manipulation**: `dplyr`, `tidyr`, `lubridate`, `readxl`
- **Statistical Modeling**: `lme4`, `emmeans`, `lavaan`, `factoextra`, `cluster`
- **Visualization**: `ggplot2`, `cowplot`, `patchwork`, `ggpubr`, `ggalluvial`

## References
* Bates, D., Mächler, M., Bolker, B., & Walker, S. (2015). Fitting Linear Mixed-Effects Models Using lme4. *Journal of Statistical Software*, 67(1), 1-48.
* Shapley, L. S. (1953). A value for n-person games. *Contributions to the Theory of Games*, 2(28), 307-317.
* Lindeman, R. H., Merenda, P. F., & Gold, R. Z. (1980). *Introduction to Bivariate and Multivariate Analysis*. Scott, Foresman.
