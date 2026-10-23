# Wheat GxE, Micro-Phenology, and Critical Window Alignment Pipeline

## Overview
This repository contains the complete "A to Z" analytical pipeline utilized in the manuscript to evaluate 200 wheat genotypes under chronic vegetative thermal stress. It encompasses everything from the raw data ingestion and variance partitioning to the novel **Critical Window Alignment (CWA)** framework in a single, unified codebase.

## Code Architecture (`01_Master_Pipeline_Public.R`)
This script executes the comprehensive biophysical and statistical analyses, moving progressively from raw spatial data to advanced predictive modeling.

# Dissecting Vegetative Micro-Phenology and Thermal Plasticity in Spring Wheat Across Contrasting Thermal Regimes

**Authors:** Apoorva Ashu (PhD) 
**Institution:** Division of Plant Physiology, ICAR–Indian Agricultural Research Institute (IARI), New Delhi, India  

---

## 1. Overview & Repository Architecture

This repository contains the complete, self-contained, end-to-end R analytical pipeline (`Wheat_GxE_Master_Analytical_Pipeline.R`), raw field trial and meteorological datasheets, GWAS-ready phenotype matrices (`Plot_ID` = Genotype 1 to 200), and the locked publication assets (Figures 1–8, Supplementary Figures S1–S6, Tables 1–8, and Supplementary Tables S1–S16).

### Experimental Regimes ($n = 200$ Spring Wheat Genotypes $\times 4$ Thermal Environments)
| Field Sheet Code | Manuscript Environment Label | Agro-Climatic Zone & Location | Sowing Date | Vegetative Window (`E` to `T4` DAS) | Stem Elongation (`SE` DAS) |
| :--- | :--- | :--- | :--- | :--- | :--- |
| `Delhi_S2` | **NWPZ Baseline** | North Western Plains Zone (ICAR-IARI, New Delhi; $28.64^\circ\text{N}, 77.16^\circ\text{E}$) | 21 Nov 2025 | 13–85 DAS (72 d) | 91 DAS |
| `Delhi_S1` | **NWPZ Warm** | North Western Plains Zone (ICAR-IARI, New Delhi; Early Vegetative Heat) | 03 Nov 2025 | 9–80 DAS (71 d) | 94 DAS |
| `Dharwad_S1` | **PZ Warm I** | Peninsular Zone (UAS Dharwad, Karnataka; $15.49^\circ\text{N}, 74.98^\circ\text{E}$) | 10 Nov 2025 | 5–45 DAS (40 d) | 65 DAS |
| `Dharwad_S2` | **PZ Warm II** | Peninsular Zone (UAS Dharwad, Karnataka; Chronic Vegetative Heat) | 20 Nov 2025 | 7–49 DAS (42 d) | 70 DAS |

---

## 2. Input Datasheets Required by the Pipeline

All four experimental input workbooks are included in the repository root (and archived in `01_Experimental_Datasheets_Weather_and_Parameters/`):
1. **`Master_Data.xlsx`**: Plot-level phenotypic measurements (`T1_Count`, `T2_Count`, `T3_Count`, `Tmax_Count` [`T4_Count`], `PT_Count`, `GrainYld`, `TKW`) across `Delhi_S2`, `Delhi_S1`, `Dharwad_S1`, and `Dharwad_S2` ($400\text{ plots/environment} = 1,600\text{ plots}$).
2. **`Weather_Master.xlsx`**: Daily meteorological records (`Tmax`, `Tmin`, `RH_Max`, `RH_Min`, `Rainfall`, `Sunshine`) and diurnal crown-soil temperatures (`9am`, `1pm`, `5pm` at 5 cm depth) for New Delhi and Dharwad.
3. **`GS population information.xlsx`**: Canonical mapping of `Genotype No.` (`Plot_ID` 1 to 200) to `Genotype Name` and physical field serial numbers (`Sr. No.` 1–200 in Rep 1 and 201–400 in Rep 2).
4. **`Data_Dharwad.xlsx`**: Original Dharwad field book reference sheets (`EarlyR1`, `EarlyR2`, `NormalR1`, `NormalR2`) used for spatial ordering verification.

---

## 3. End-to-End Reproduction Instructions

### System & Package Requirements
- **R Version:** $\ge 4.4.0$ (tested on R 4.5.1, Windows/Linux/macOS)
- **Required CRAN Packages (auto-installed if missing):**
  `lme4`, `car`, `emmeans`, `nlme`, `readxl`, `openxlsx`, `dplyr`, `tidyr`, `tibble`, `purrr`, `stringr`, `lubridate`, `metan`, `broom`, `scales`, `cluster`, `ggplot2`, `ggrepel`, `ggpubr`, `patchwork`, `cowplot`, `ggdist`, `ggsci`, `mgcv`, `digest`

### Single-Command Execution
From the repository root containing `Wheat_GxE_Master_Analytical_Pipeline.R` and the four `.xlsx` input files:
```bash
Rscript Wheat_GxE_Master_Analytical_Pipeline.R
```

---

## 4. Complete Locked Numbering Manifest (Sections, Figures, and Tables)

### 4.1 Analytical Sections
| Pipeline Section | Asset Folder | Manuscript Scope |
| :--- | :--- | :--- |
| **Section 1** | `Section_1_Weather_and_Treatments` | Meteorological characterization, diurnal thermal anomalies ($\Delta T_{\max}, \Delta T_{\min}$), crown-soil buffering, $\text{GDD}_{\text{std}}$, and Nighttime Heat Units ($\text{NHU}$) |
| **Section 2.0** | `Section_2_BLUE_Extraction` | Non-destructive Delhi Rep-2 block mapping (`Block_true`), spatial layout diagnostics, and Stage-1 Linear Mixed Model BLUE extraction |
| **Section 2.1** | `Section_2.1_Phenotypic_Distributions` | Multi-environment phenotypic distributions, paired genotype reaction norms, and repeated-measures LMM ANOVA |
| **Section 2.2** | `Section_2.2_Phenology_Classifications` | Unsupervised K-Means tillering vigor & trajectory classifications, silhouette validation, and McNemar regime transitions |
| **Section 2.3** | `Section_2.3_Phenology_Indices` | Micro-phenology architecture ($\text{AUTPC}_{\text{norm}}$ and Tiller Survival Index [$\text{TSI}$]) and out-of-sample cross-regime yield validation |
| **Section 2.4** | `Section_2.4_Thermal_Modeling` | Non-linear thermal and chronological tillering trajectories via REML-penalized Generalized Additive Models (GAMs) |
| **Section 2.5** | `Section_2.5_Categorical_Classifications` | Three-factor categorical tillering strategies ($2 \times 2 \times 2 = 8$ profiles), multi-regime yield dynamics, and composite strategy rankings |
| **Section 2.6** | `Section_2.6_Phenotypic_Plasticity` | Finlay–Wilkinson OLS and empirical Bayes LMM (BLUP) reaction norms, LOEO stability diagnostics, and adaptation quadrants |
| **Section 2.7** | `Section_2.7_Variance_and_Heritability` | Three-tier variance partitioning (within-trial, Smith-weighted two-stage, and plot-level RCBD) and Cullis/Piepho broad-sense heritability ($H^2$) |
| **Section 2.8** | `Section_2.8_GGE_Stability` | Multi-environment GGE biplot suite, zone-specific (NWPZ vs. PZ) adaptation trade-offs, and Unified Ideotype Z-Score selection |
| **Section 10** | `GWAS_Datasheets` | Complete 4-environment Stage-1 BLUE and cross-environment plasticity/stability GWAS phenotype datasheets ordered strictly by `Plot_ID` (Genotype 1 to 200) |

### 4.2 Main Text Figures (`Figure 1` – `Figure 8`)
| Figure ID | Section Folder | Primary Filename | Description |
| :--- | :--- | :--- | :--- |
| **Figure 1** | `Section_1_Weather_and_Treatments` | `Figure_Thermal_Disparities_and_Stage_Dynamics_v2.png` | Diurnal air thermal anomalies, crown-soil thermal buffering, and stage-resolved $\text{GDD}_{\text{std}}$ vs. $\text{NHU}$ (split parts: `Figure_1_Part1_...png`, `Figure_1_Part2_...png`) |
| **Figure 2** | `Section_2.1_Phenotypic_Distributions` | `Global_Parameter_Paired_Genotype_Profiles.png` | 9-panel multi-environment phenotypic distributions and paired genotype reaction norms across all 4 regimes |
| **Figure 3** | `Section_2.2_Phenology_Classifications` | `Compound_Phenological_Transitions_Vigor_Trajectory_TCE.png` | Unsupervised K-Means tillering vigor & trajectory classifications, McNemar regime shifts, and TCE yield architecture |
| **Figure 4** | `Section_2.3_Phenology_Indices` | `TSI_Compound_Yield_Validation_v2.png` | Micro-phenology index architecture ($\text{AUTPC}_{\text{norm}}$ and $\text{TSI}$), multi-stress transitions, and grain yield validation |
| **Figure 5** | `Section_2.4_Thermal_Modeling` | `Thermal_GAM_Trajectories.png` | REML-GAM non-linear tillering response trajectories across $\text{GDD}_{\text{std}}$ and chronological DAE |
| **Figure 6** | `Section_2.5_Categorical_Classifications` | `Figure_Strategy_Yield_Dynamics_Complete_v2.png` | Three-factor tillering strategy classification, regime plasticity transitions, and 8-profile multi-regime yield performance |
| **Figure 7** | `Section_2.6_Phenotypic_Plasticity` | `Finlay_Wilkinson_Biplot.png` | Finlay–Wilkinson joint regression stability biplot (plasticity slope $b_i$ vs. 4-environment mean grain yield) |
| **Figure 8** | `Section_2.8_GGE_Stability` | `Figure_Top_Genotypes_Ideotype_Profiles_v2.png` | Trans-zonal adaptation trade-off map (A), 4-environment yield ranking (B), and reaction norms / standardized profiles of Top 10 Ideotypes vs. Popular Released Varieties (C1–C2) |

### 4.3 Supplementary Figures (`Figure S1` – `Figure S6`)
| Figure ID | Section Folder | Primary Filename | Description |
| :--- | :--- | :--- | :--- |
| **Figure S1** | `Section_1_Weather_and_Treatments` | `Figure_S1_Integrated_ThermoPhenology_v2.png` | Continuous daily meteorological trajectories, canopy-soil microclimate, and stage windows across all 4 environments |
| **Figure S2** | `Section_2.4_Thermal_Modeling` | `Thermal_GAM_Genotype_Faceted_Supplement.png` | Genotype-faceted REML-GAM tillering trajectories across thermal accumulation metrics |
| **Figure S3** | `Section_2.5_Categorical_Classifications` | `Figure_Strategy_Regime_Rankings_Comparison.png` | Multi-regime yield rank progression (bump chart) and composite strategy score across the 8 categorical tillering profiles |
| **Figure S4** | `Section_2.6_Phenotypic_Plasticity` | `Finlay_Wilkinson_Reaction_Norms.png` | Finlay–Wilkinson linear reaction norms of the 200 genotypes across the environmental index ($I_j$), stratified by quadrant |
| **Figure S5** | `Section_2.6_Phenotypic_Plasticity` | `Finlay_Wilkinson_BLUP_Shrinkage.png` | Empirical Bayes shrinkage comparison of fixed-effects OLS vs. random-slopes LMM (BLUP) plasticity slopes |
| **Figure S6** | `Section_2.8_GGE_Stability` | `Figure_S6_GGE_Composite.png` | Two-panel GGE biplot suite: Mega-Environment Which-Won-Where polygon (A) and Mean vs. Stability AEC view (B) |

### 4.4 Main Text Tables (`Table 1` – `Table 8`) & Supplementary Tables (`Table S1` – `Table S16`)
| Table ID | Section Folder | Primary CSV / Workbook | Summary |
| :--- | :--- | :--- | :--- |
| **Table 1** | `Section_1_Weather_and_Treatments` | `Table_1_Vegetative_Thermal_and_Crown_Soil_Regimes_v2.csv` | Trial coordinates, sowing dates, stage-resolved $\text{GDD}_{\text{std}}$, $\text{NHU}$, and diurnal thermal anomalies |
| **Table 2** | `Section_2.1_Phenotypic_Distributions` | `Table_2_Multi_Environment_Phenotypic_and_ANOVA_Summary.csv` | Environment-wise means $\pm$ SD, ranges, and LMM regime contrasts across 9 phenotypic traits |
| **Table 3** | `Section_2.2_Phenology_Classifications` | `Table_3_Phenological_Typology_and_Regime_Transitions.csv` | Phenological typology classifications (Vigor, Trajectory, TCE) and McNemar transition tests |
| **Table 4** | `Section_2.3_Phenology_Indices` | `Table_4_MicroPhenology_AUTPC_and_TSI_Yield_Correlations.csv` | Within-environment $\text{AUTPC}_{\text{norm}}$ and $\text{TSI}$ correlations with grain yield and productive tillers |
| **Table 5** | `Section_2.4_Thermal_Modeling` | `Table_5_GAM_Tillering_Trajectories_and_Stage_Means.csv` | REML-GAM effective degrees of freedom ($\text{edf}$), $F$-statistics, deviance explained, and stage means |
| **Table 6** | `Section_2.5_Categorical_Classifications` | `Table_6_Strategy_Profile_Regime_Rankings.csv` | Multi-regime grain yield rankings, composite strategy scores, and omnibus tests across the 8 strategy profiles |
| **Table 7** | `Section_2.7_Variance_and_Heritability` | `Table_7_MAIN_TEXT_Stratified_Variance_and_Heritability.xlsx` | Within-trial RCBD/spatial variance components and two-stage Cullis/Piepho broad-sense heritability ($H^2$) |
| **Table 8** | `Section_2.8_GGE_Stability` | `Table_Top10_vs_Checks_Ideotype_Summary_v2.csv` | Composite Ideotype Z-Score, zone-specific indices, GGE stability penalty, and strategy profiles of Top 10 lines vs. Checks |
| **Table S0** | `Section_2_BLUE_Extraction` | `Table_Genotype_BLUEs_and_SEs.csv` (+ 9 diagnostic CSVs) | Complete Stage-1 Genotype $\times$ Environment BLUEs & SEs ($7,200$ estimates) and block/spatial diagnostics |
| **Tables S1–S2** | `Section_1_Weather_and_Treatments` | `Section_1_Manuscript_and_Supplementary_Tables_v2.xlsx` | Stage-specific thermal disparity dynamics (`Table S1`) and crown-soil diurnal temperature observations (`Table S2`) |
| **Table S3** | `Section_2.1_Phenotypic_Distributions` | `Table_S3_Pairwise_Regime_Contrasts_and_Heteroskedasticity.csv` | Tukey-adjusted pairwise regime contrasts and Brown–Forsythe / Fligner–Killeen variance heterogeneity tests |
| **Tables S4–S5** | `Section_2.2_Phenology_Classifications` | `Section_2.2_Manuscript_and_Supplementary_Tables.xlsx` | K-Means cluster centroids & silhouette validity (`Table S4`) and 200-genotype class assignments (`Table S5`) |
| **Tables S6–S7** | `Section_2.3_Phenology_Indices` | `Section_2.3_Manuscript_and_Supplementary_Tables.xlsx` | Multi-stress TSI transitions, 200-genotype TSI scores, out-of-sample validation (`Table S6`), and AUTPC anchor sensitivity (`Table S7`) |
| **Tables S8–S10** | `Section_2.4_Thermal_Modeling` | `Section_2.4_Manuscript_and_Supplementary_Tables.xlsx` | GAM model selection AIC (`Table S8`), smooth/parametric coefficients (`Table S9`), and empirical thermal growth rates (`Table S10`) |
| **Table S11** | `Section_2.5_Categorical_Classifications` | `Section_2.5_Manuscript_and_Supplementary_Tables.xlsx` | Side-by-side omnibus tests (`S11A`), cluster-free shape regression (`S11B`), and secondary tiller/TCE class tables (`S11C–D`) |
| **Tables S12–S14**| `Section_2.6_Phenotypic_Plasticity` | `Table_S12_to_S14_SUPPLEMENTARY_Finlay_Wilkinson_Plasticity.xlsx` | FW quadrant summary & LOEO diagnostics (`Table S12`), developmental cluster cross-tabulation (`Table S13`), and all 200 FW slopes (`Table S14`) |
| **Table S15** | `Section_2.7_Variance_and_Heritability` | `Table_S15_SUPPLEMENTARY_Stage2_and_Delhi_Plot_ANOVA.xlsx` | Two-Stage Smith-weighted multi-environment ANOVA (`S15A`) and Delhi plot-level RCBD ANOVA (`S15B`) |
| **Table S16** | `Section_2.8_GGE_Stability` | `Table_8_MAIN_TEXT_and_Table_S16_SUPPLEMENTARY_Ideotype_Selection.xlsx` | All 200 genotype ideotype rankings (`S16A`), named cultivars (`S16B`), NWPZ zone ranking (`S16C`), and PZ zone ranking (`S16D`) |
da, P. F., & Gold, R. Z. (1980). *Introduction to Bivariate and Multivariate Analysis*. Scott, Foresman.
