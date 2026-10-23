# ==============================================================================
# MASTER REPRODUCIBLE ANALYTICAL PIPELINE for Multilocation GxE Experiment
# ==============================================================================
# Manuscript Title:
#   Dissecting Vegetative Micro-Phenology and Thermal Plasticity in Spring Wheat 
#   Across Contrasting Thermal Regimes
# # Author: 
#   Apoorva Ashu (PhD), Division of Plant Physiology, ICAR-IARI, New Delhi, India
# Mail: apoorva_13125@iari.res.in
# Repository:
#   Wheat GxE Micro-Phenology and Chronic Vegetative Thermal Stress
# ==============================================================================
# DESCRIPTION & WORKFLOW OVERVIEW:
# This script constitutes the authoritative, self-contained analytical pipeline
# supporting the entire manuscript. It reads directly from raw field trial sheets
# in 'Master_Data.xlsx', meteorological sensor records in 'Weather_Master.xlsx',
# and the canonical plot/block mapping in 'GS population information.xlsx'.
# It executes end-to-end:
#   0. Package verification, global biometric constants & publication theme setup
#   1. Section 1   (Section_1_Weather_and_Treatments): Meteorological characterization, crown-soil buffering, GDD_std & NHU (Figure 1, Figure S1, Table 1, Tables S1-S2)
#   2. Section 2.0 (Section_2_BLUE_Extraction): Non-destructive Delhi block mapping, field layout diagnostics & Stage-1 BLUE extraction (Stage-1 BLUE Audit Suite)
#   3. Section 2.1 (Section_2.1_Phenotypic_Distributions): Multi-environment phenotypic distributions & repeated-measures LMM ANOVA (Figure 2, Table 2, Table S3)
#   4. Section 2.2 (Section_2.2_Phenology_Classifications): Unsupervised K-Means vigor/trajectory classifications & McNemar regime transitions (Figure 3, Table 3, Tables S4-S5)
#   5. Section 2.3 (Section_2.3_Phenology_Indices): Micro-phenology architecture (AUTPC_norm & Tiller Survival Index [TSI]) & cross-regime validation (Figure 4, Table 4, Tables S6-S7)
#   6. Section 2.4 (Section_2.4_Thermal_Modeling): Non-linear thermal & chronological tillering trajectories via REML-GAM (Figure 5, Figure S2, Table 5, Tables S8-S10)
#   7. Section 2.5 (Section_2.5_Categorical_Classifications): Three-factor tillering strategy profiles, yield dynamics & rank stability (Figure 6, Figure S3, Table 6, Table S11)
#   8. Section 2.6 (Section_2.6_Phenotypic_Plasticity): Finlay-Wilkinson linear & BLUP reaction norms and LOEO stability diagnostics (Figure 7, Figures S4-S5, Tables S12-S14)
#   9. Section 2.7 (Section_2.7_Variance_and_Heritability): Multi-tier variance partitioning & Cullis/Piepho broad-sense heritability H2 (Table 7, Table S15)
#  10. Section 2.8 (Section_2.8_GGE_Stability): Multi-environment GGE biplots, zone-specific adaptation & unified ideotype selection (Figure 8, Figure S6, Table 8, Table S16)
#  11. Section 10  (GWAS_Datasheets): Complete 4-environment Stage-1 BLUE & cross-environment plasticity/stability GWAS phenotype datasheets ordered strictly by canonical Plot_ID (Genotype 1 to 200)
#
# All publication-grade figures (600 DPI), structured CSVs, and multi-sheet Excel
# workbooks are exported into 'Final_Manuscript_Assets/' and 'GWAS_Datasheets/'.
# ==============================================================================

# ==============================================================================
# 0. GLOBAL SETUP, DIRECTORIES, AND REPRODUCIBILITY
# ==============================================================================

cat("\n==============================================================================\n")
cat(" STARTING WHEAT GxE REPRODUCIBLE MASTER ANALYTICAL PIPELINE\n")
cat("==============================================================================\n\n")

# Detect Project Root (ensure execution relative to workspace root)
if (basename(getwd()) == "Reproducible_Manuscript_Pipeline") {
  setwd("..")
}
project_root <- getwd()
cat("Project Working Directory:", project_root, "\n\n")

# Verify and Load Required Packages
required_packages <- c(
  "lme4", "car", "emmeans", "nlme", "readxl", "dplyr", "tidyr", "tibble", "purrr",
  "stringr", "lubridate", "metan", "broom", "scales", "cluster",
  "ggplot2", "ggrepel", "ggpubr", "patchwork", "cowplot", "ggdist",
  "ggsci", "mgcv", "digest"
)

missing_pkgs <- setdiff(required_packages, rownames(installed.packages()))
if (length(missing_pkgs) > 0) {
  cat("Installing missing package dependencies:", paste(missing_pkgs, collapse = ", "), "\n")
  install.packages(missing_pkgs, repos = "http://cran.us.r-project.org")
}
invisible(lapply(required_packages, library, character.only = TRUE))

# Global Reproducibility Seed
set.seed(42)

# Centralized Biometric & Unit Constants (Single Source of Truth Across Sections)
# Note on Yield Units: GrainYld is analysed in harvested plot units (kg/plot; 0.675 m^2 net harvested area).
# Tiller counts (T1..T4, PT) represent per-plant averages from 5 randomly sampled plants per plot.
PLOT_HARVESTED_AREA_M2  <- 0.675
YIELD_UNIT_LABEL        <- "kg/plot"
COUNT_AT_EMERGENCE_MAIN <- 1.0
COUNT_AT_EMERGENCE_ALT  <- 0.0
DELHI_PRIMARY_MODEL     <- "A"

# Centralized Vigor Typology Labels (Option A: Mean Tillering Rate Clustering)
VIGOR_HI <- "High Tillering Vigor"
VIGOR_LO <- "Low Tillering Vigor"

# Centralized Input Data File Paths (Self-Contained Repository Paths)
data_file_path <- file.path(project_root, "Master_Data.xlsx")
gs_map_path    <- file.path(project_root, "GS population information.xlsx")
dhw_orig_path  <- file.path(project_root, "Data_Dharwad.xlsx")

# Global Standardized Publication Theme (Nature/Field Crops Research Standard)
pub_theme <- theme_classic(base_size = 12) +
  theme(
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 12),
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black", size = 10),
    legend.position = "bottom"
  )

# Format p-values dynamically with asterisks for publication annotations
format_p_stat <- function(p_val) {
  stars <- if (p_val < 0.001) "***" else if (p_val < 0.01) "**" else if (p_val < 0.05) "*" else "(ns)"
  p_num <- if (p_val < 0.001) sprintf("%.2e", p_val) else sprintf("%.4f", p_val)
  sprintf("p = %s %s", p_num, stars)
}

format_p_display <- function(p) {
  if (is.na(p)) return("N/A")
  if (p < 0.001) sprintf("%.2e ***", p)
  else if (p < 0.01) sprintf("%.4f **", p)
  else if (p < 0.05) sprintf("%.4f *", p)
  else sprintf("%.4f (ns)", p)
}

# Universal ASCII & Excel-Safe Table Sanitizer (Prevents Windows-1252 Mojibake & #NAME? Formula Errors)
sanitize_text_vec <- function(v) {
  if (!is.character(v) && !is.factor(v)) return(v)
  was_factor <- is.factor(v)
  s <- as.character(v)
  s <- gsub("r\u0304", "mean_r", s, fixed = TRUE)
  s <- gsub("\u03c7\u00b2", "Chi-sq", s, fixed = TRUE)
  s <- gsub("\u0394", "Delta ", s, fixed = TRUE)
  s <- gsub("Delta  ", "Delta ", s, fixed = TRUE)
  s <- gsub("\u03c3\u00b2", "sigma2", s, fixed = TRUE)
  s <- gsub("\u00b1", "+/-", s, fixed = TRUE)
  s <- gsub("\u2013", "-", s, fixed = TRUE)
  s <- gsub("\u2014", "--", s, fixed = TRUE)
  s <- gsub("\u2212", "-", s, fixed = TRUE)
  s <- gsub("\u2192", "->", s, fixed = TRUE)
  s <- gsub("\u2265", ">=", s, fixed = TRUE)
  s <- gsub("\u2264", "<=", s, fixed = TRUE)
  s <- gsub("\u00d7", "x", s, fixed = TRUE)
  s <- gsub("\u00b0C\u00b7d", "deg C d", s, fixed = TRUE)
  s <- gsub("\u00b0C", "deg C", s, fixed = TRUE)
  s <- gsub("\u00b0", "deg ", s, fixed = TRUE)
  s <- gsub("\u00b7", ".", s, fixed = TRUE)
  s <- gsub("\u00b2", "2", s, fixed = TRUE)
  # Prevent Excel from parsing composite strings starting with + or - as broken formulas (#NAME?)
  is_composite_signed <- !is.na(s) & grepl("^[+-][0-9].*[^0-9.eE+-]", s)
  s[is_composite_signed] <- paste0("Delta: ", s[is_composite_signed])
  if (was_factor) factor(s, levels = unique(s)) else s
}

sanitize_df_for_export <- function(df) {
  if (!is.data.frame(df)) return(df)
  names(df) <- sanitize_text_vec(names(df))
  for (col_i in seq_along(df)) {
    if (is.character(df[[col_i]]) || is.factor(df[[col_i]])) {
      df[[col_i]] <- sanitize_text_vec(df[[col_i]])
    }
  }
  df
}

write.csv <- function(x, file = "", row.names = FALSE, ...) {
  x_clean <- sanitize_df_for_export(x)
  utils::write.csv(x_clean, file = file, row.names = row.names, fileEncoding = "UTF-8", ...)
}

# Output Root Directory
assets_dir <- file.path(project_root, "Final_Manuscript_Assets")
dir.create(assets_dir, showWarnings = FALSE, recursive = TRUE)

archive_rejected_dir <- file.path(project_root, "Archive_Rejected")
dir.create(archive_rejected_dir, showWarnings = FALSE, recursive = TRUE)

# Standardized Publication Excel Header Style
if (requireNamespace("openxlsx", quietly = TRUE)) {
  header_style <- openxlsx::createStyle(fontColour = "#ffffff", fgFill = "#3C5488", halign = "center", textDecoration = "bold")
}

# ==============================================================================
# MODULAR EXECUTION CONTROL FLAGS & DEPENDENCY MAP
# ==============================================================================
# Dependency Map:
#   - Section 2 (BLUEs & Master Raw): Required by Sections 3, 4, 5, 6, 7, 8, 9
#     (Even if RUN_SECTION_2_BLUES = FALSE, the cache loader reconstructs blues_df & master_raw)
#   - Section 4 (Phenology Typology): Required by Sections 5, 7, 8, 9 (populates Trajectory & Vigor on blues_df)
#   - Section 5 (Microphenology & TSI): Required by Sections 7, 8, 9 (populates tsi_df & TSI columns on blues_df)
#   - Section 7 (Categorical Architecture): Required by Section 8 (creates master_profiles for Table S13 cross-tab)
if (!exists("CLEAR_BLUE_CACHE"))             CLEAR_BLUE_CACHE             <- FALSE
if (!exists("RUN_SECTION_1_WEATHER"))        RUN_SECTION_1_WEATHER        <- TRUE
if (!exists("RUN_SECTION_2_BLUES"))          RUN_SECTION_2_BLUES          <- TRUE
if (!exists("RUN_SECTION_3_DISTRIBUTIONS"))  RUN_SECTION_3_DISTRIBUTIONS  <- TRUE
if (!exists("RUN_SECTION_4_PHENOLOGY"))      RUN_SECTION_4_PHENOLOGY      <- TRUE
if (!exists("RUN_SECTION_5_MICROPHENOLOGY")) RUN_SECTION_5_MICROPHENOLOGY <- TRUE
if (!exists("RUN_SECTION_6_THERMAL_GAM"))    RUN_SECTION_6_THERMAL_GAM    <- TRUE
if (!exists("RUN_SECTION_7_CATEGORICAL"))    RUN_SECTION_7_CATEGORICAL    <- TRUE
if (!exists("RUN_SECTION_8_PLASTICITY"))     RUN_SECTION_8_PLASTICITY     <- TRUE
if (!exists("RUN_SECTION_9_GGE_STABILITY"))  RUN_SECTION_9_GGE_STABILITY  <- TRUE
if (!exists("RUN_SECTION_10_GWAS_EXPORT"))   RUN_SECTION_10_GWAS_EXPORT   <- TRUE


# ==============================================================================
# SINGLE SOURCE OF TRUTH (SSOT): MASTER PHENOLOGICAL SPECIFICATION
# ==============================================================================
# Verified against field records
master_pheno_spec <- data.frame(
  Sheet       = c("Delhi_S2", "Delhi_S1", "Dharwad_S1", "Dharwad_S2"),
  Env_Name    = c("NWPZ Baseline", "NWPZ Warm", "PZ Warm I", "PZ Warm II"),
  Location    = c("Delhi", "Delhi", "Dharwad", "Dharwad"),
  Sowing      = c("S2", "S1", "S1", "S2"),
  Role        = c("Optimal Baseline", "Early Heat Stress", "PZ Warm Sowing I", "PZ Warm Sowing II"),
  Sowing_Date = as.Date(c("2025-11-21", "2025-11-03", "2025-11-10", "2025-11-20")),
  E           = c(13, 9, 5, 7),
  T1          = c(40, 35, 15, 19),
  T2          = c(55, 50, 25, 29),
  T3          = c(70, 65, 35, 39),
  T4          = c(85, 80, 45, 49),
  SE          = c(91, 94, 65, 70),
  stringsAsFactors = FALSE
)

# Standardized Environment Taxonomy & Metadata (Derived dynamically from master_pheno_spec)
env_metadata <- master_pheno_spec %>%
  mutate(
    Emergence_Date = Sowing_Date + E,
    Start_Date     = Sowing_Date + E, # Preserved for downstream compatibility
    End_Date       = Sowing_Date + SE,
    Int_E_T1       = T1 - E,
    Int_T1_T2      = T2 - T1,
    Int_T2_T3      = T3 - T2,
    Int_T3_T4      = T4 - T3,
    Total_E_T4     = T4 - E,
    Env_Factor     = factor(Env_Name, levels = Env_Name)
  )



# ==============================================================================
# SECTION 1: METEOROLOGICAL CHARACTERIZATION & THERMAL PHENOLOGY
# ==============================================================================
cat(">>> Executing Section 1: Meteorological Analysis & Thermal Disparities...\n")
sec1_dir <- file.path(assets_dir, "Section_1_Weather_and_Treatments")
dir.create(sec1_dir, showWarnings = FALSE, recursive = TRUE)

# 1.1 Load Weather Master Data and Configure Thermal Parameters
base_temp <- 5.0
tmax_cap  <- 30.0

weather <- read_excel(file.path(project_root, "Weather_Master.xlsx"), sheet = "Sheet1")
weather$Date <- as.Date(weather$Date)

# Verify weather file coverage and absence of missing values across all trial windows
for (i in 1:nrow(master_pheno_spec)) {
  loc_i     <- master_pheno_spec$Location[i]
  s_date_i  <- master_pheno_spec$Sowing_Date[i]
  se_date_i <- s_date_i + master_pheno_spec$SE[i]
  req_dates <- seq(s_date_i + 1, se_date_i, by = "day")
  
  w_check <- weather %>% filter(Location == loc_i, Date %in% req_dates)
  
  stopifnot(
    length(unique(w_check$Date)) == length(req_dates),
    nrow(w_check) == length(req_dates),
    !any(is.na(w_check$Tmax)),
    !any(is.na(w_check$Tmin))
  )
}

# Calculate daily mean temperatures and GDD accumulation (uncapped, nocturnal Tmin-only, and capped at daily Tmax)
weather <- weather %>% 
  mutate(
    Tmean              = (Tmax + Tmin) / 2, 
    Daily_GDD_Uncapped = pmax(0, Tmean - base_temp),
    Daily_GDD          = Daily_GDD_Uncapped,
    Daily_GDD_Tmin     = pmax(0, Tmin - base_temp),
    Daily_GDD_Tmax     = pmax(0, Tmax - base_temp),
    Tmax_Capped        = pmin(Tmax, tmax_cap),
    Tmean_Capped       = (Tmax_Capped + Tmin) / 2,
    Daily_GDD_Capped   = pmax(0, Tmean_Capped - base_temp),
    Daily_GDD_Cap32    = pmax(0, (pmin(Tmax, 32.0) + Tmin) / 2 - base_temp)
  )

soil <- read_excel(file.path(project_root, "Weather_Master.xlsx"), sheet = "Soil_Temperature")
soil$Date <- as.Date(soil$Date)

# 1.2 Export Stage-Specific Weather Tables (Derived strictly from master_pheno_spec)
# Half-open intervals (start_das, end_das] ensure boundary days are counted in exactly one stage
calc_stage_weather <- function() {
  stage_list <- list()
  for (i in 1:nrow(master_pheno_spec)) {
    e_code      <- master_pheno_spec$Sheet[i]
    e_name      <- master_pheno_spec$Env_Name[i]
    loc         <- master_pheno_spec$Location[i]
    sowing_date <- master_pheno_spec$Sowing_Date[i]
    
    stages <- list(
      "E-T1"  = c(master_pheno_spec$E[i],  master_pheno_spec$T1[i]),
      "T1-T2" = c(master_pheno_spec$T1[i], master_pheno_spec$T2[i]),
      "T2-T3" = c(master_pheno_spec$T2[i], master_pheno_spec$T3[i]),
      "T3-T4" = c(master_pheno_spec$T3[i], master_pheno_spec$T4[i])
    )
    
    for (st in names(stages)) {
      rng       <- stages[[st]]
      start_das <- rng[1]
      end_das   <- rng[2]
      start_dt  <- sowing_date + start_das
      end_dt    <- sowing_date + end_das
      
      sub_st <- weather %>% 
        filter(Location == loc, Date > start_dt, Date <= end_dt) %>% 
        arrange(Date)
      
      stopifnot(nrow(sub_st) == (end_das - start_das))
      stopifnot(!any(is.na(sub_st$Tmax)))
      stopifnot(!any(is.na(sub_st$Tmin)))
      
      stage_list[[paste(e_name, st)]] <- data.frame(
        Environment      = e_name,
        Stage            = st,
        Duration_Days    = end_das - start_das,
        Tmean            = mean(sub_st$Tmean, na.rm = TRUE),
        Tmax             = mean(sub_st$Tmax, na.rm = TRUE),
        Tmin             = mean(sub_st$Tmin, na.rm = TRUE),
        Stage_GDD        = sum(sub_st$Daily_GDD_Uncapped, na.rm = TRUE),
        Stage_GDD_Tmin   = sum(sub_st$Daily_GDD_Tmin, na.rm = TRUE),
        Stage_GDD_Tmax   = sum(sub_st$Daily_GDD_Tmax, na.rm = TRUE),
        Stage_GDD_Capped = sum(sub_st$Daily_GDD_Capped, na.rm = TRUE),
        Stage_GDD_Cap32  = sum(sub_st$Daily_GDD_Cap32, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }
  }
  bind_rows(stage_list)
}

stage_summary <- calc_stage_weather()

# Stage-by-Environment Cumulative GDD from Emergence (Uncapped, Nocturnal Tmin, and Capped)
gdd_from_emergence_table <- stage_summary %>%
  group_by(Environment) %>%
  mutate(
    Cum_GDD_Uncapped = cumsum(Stage_GDD),
    Cum_GDD_Tmin     = cumsum(Stage_GDD_Tmin),
    Cum_GDD_Capped   = cumsum(Stage_GDD_Capped)
  ) %>%
  ungroup() %>%
  select(Environment, Stage, Duration_Days, Stage_GDD_Uncapped = Stage_GDD, Stage_GDD_Tmin, Stage_GDD_Capped, Cum_GDD_Uncapped, Cum_GDD_Tmin, Cum_GDD_Capped)

if (isTRUE(RUN_SECTION_1_WEATHER)) {
  # Total Uncapped & Nocturnal Vegetative GDD from Emergence to T4 Milestone Across Environments
  total_gdd_env_sec1 <- gdd_from_emergence_table %>%
    filter(Stage == "T3-T4") %>%
    transmute(
      Environment,
      Total_GDD_Uncapped = Cum_GDD_Uncapped,
      Total_GDD_Tmin     = Cum_GDD_Tmin
    )
  spread_uncapped_sec1 <- (max(total_gdd_env_sec1$Total_GDD_Uncapped) - min(total_gdd_env_sec1$Total_GDD_Uncapped)) / min(total_gdd_env_sec1$Total_GDD_Uncapped) * 100

  cat("\n--- Total Vegetative Thermal Accumulation (Emergence to T4 Milestone, Uncapped & Nocturnal GDD_Tmin) ---\n")
  print(total_gdd_env_sec1)
  cat(sprintf("Vegetative GDD Range (E to T4): %.1f to %.1f deg C d (Spread: %.2f%%)\n\n",
              min(total_gdd_env_sec1$Total_GDD_Uncapped), max(total_gdd_env_sec1$Total_GDD_Uncapped), spread_uncapped_sec1))

  # Stage-Specific Temperature & Thermal Disparities Relative to NWPZ Baseline
  baseline_st <- stage_summary %>%
    filter(Environment == "NWPZ Baseline") %>%
    mutate(Pace_Baseline = Stage_GDD / Duration_Days)

  disparities <- stage_summary %>%
    filter(Environment != "NWPZ Baseline") %>%
    mutate(Pace_Target = Stage_GDD / Duration_Days) %>%
    left_join(baseline_st, by = "Stage", suffix = c("_Target", "_Baseline")) %>%
    mutate(
      Delta_Days       = Duration_Days_Target - Duration_Days_Baseline,
      Tmean            = Tmean_Target,
      Tmax             = Tmax_Target,
      Tmin             = Tmin_Target,
      Delta_Tmean      = Tmean_Target - Tmean_Baseline,
      Delta_Tmax       = Tmax_Target - Tmax_Baseline,
      Delta_Tmin       = Tmin_Target - Tmin_Baseline,
      Delta_GDD        = Stage_GDD_Target - Stage_GDD_Baseline,
      Delta_GDD_Tmin   = Stage_GDD_Tmin_Target - Stage_GDD_Tmin_Baseline,
      Delta_Pace_GDD_d = Pace_Target - Pace_Baseline
    ) %>%
    select(Environment = Environment_Target, Stage, Delta_Days, Tmean, Tmax, Tmin,
           Delta_Tmean, Delta_Tmax, Delta_Tmin, Delta_GDD, Delta_GDD_Tmin, Delta_Pace_GDD_d)

  # 1.4 Clean Up & Archive Legacy/Redundant Section 1 Plots and CSVs
  weather_alt_dir <- file.path(project_root, "Curiosity_And_Iterations", "Weather_Alternative_Plots")
  dir.create(weather_alt_dir, showWarnings = FALSE, recursive = TRUE)
  sec1_archive_dir <- file.path(archive_rejected_dir, "Section_1_Legacy_Weather")
  dir.create(sec1_archive_dir, showWarnings = FALSE, recursive = TRUE)

  # Move misplaced AUTPC_norm anchor table to Section 2.3 if present in Section 1
  sec5_dir_early <- file.path(assets_dir, "Section_2.3_Phenology_Indices")
  dir.create(sec5_dir_early, showWarnings = FALSE, recursive = TRUE)
  misplaced_anchor <- file.path(sec1_dir, "Table_Anchor_Sensitivity_AUTPC_norm.csv")
  if (file.exists(misplaced_anchor)) {
    file.copy(misplaced_anchor, file.path(sec5_dir_early, "Table_Anchor_Sensitivity_AUTPC_norm.csv"), overwrite = TRUE)
    unlink(misplaced_anchor)
  }

  legacy_sec1_files <- c(
    "Fig_1A_Weather_TimeSeries.png",
    "Fig_1B_Window_Tmax_Comparison.png",
    "Fig_1C_Window_GDD_Comparison.png",
    "Weather_Soil_Crown_Dynamics.png",
    "Table_S1_Phenological_Window_Weather.csv",
    "Table_S1_Weather_Treatment_Summary.csv",
    "Table_GDD_Emergence_Dynamics.csv",
    "Table_GDD_Emergence_Stage_by_Environment_Wide.csv",
    "Table_Crown_Soil_Environment_Summary.csv",
    "Table_Crown_Soil_Temperature_Regimes.csv",
    "Table_S1_Stage_Specific_Thermal_Dynamics.csv",
    "Table_S2_Temperature_Disparities.csv"
  )
  for (lf in legacy_sec1_files) {
    lf_path <- file.path(sec1_dir, lf)
    if (file.exists(lf_path)) {
      file.copy(lf_path, file.path(sec1_archive_dir, lf), overwrite = TRUE)
      unlink(lf_path)
    }
  }

  # 1.5 Crown-Zone Diurnal Soil Temperature Processing (09:00 h, 13:00 h, 17:00 h)
  soil_clean <- soil %>%
    mutate(Sheet = paste0(Location, "_", Sowing)) %>%
    left_join(env_metadata %>% select(Sheet, Environment = Env_Name), by = "Sheet") %>%
    group_by(Environment) %>%
    arrange(DAS, .by_group = TRUE) %>%
    mutate(
      Stage_Milestone = paste0("T", row_number()),
      Soil_9am = `9am`,
      Soil_1pm = `1pm`,
      Soil_5pm = `5pm`,
      Soil_Diurnal_Mean  = (Soil_9am + Soil_1pm + Soil_5pm) / 3,
      Soil_Diurnal_Swing = Soil_1pm - pmin(Soil_9am, Soil_5pm)
    ) %>%
    ungroup() %>%
    left_join(weather %>% select(Location, Date, Air_Tmax = Tmax, Air_Tmin = Tmin, Air_Tmean = Tmean),
              by = c("Location", "Date")) %>%
    mutate(
      Soil_1pm_minus_Air_Tmax   = Soil_1pm - Air_Tmax,
      Soil_Mean_minus_Air_Tmean = Soil_Diurnal_Mean - Air_Tmean,
      Env_Factor = factor(Environment, levels = env_metadata$Env_Name)
    ) %>%
    select(Environment, Location, Sowing, Stage_Milestone, Date, DAS,
           Soil_9am, Soil_1pm, Soil_5pm, Soil_Diurnal_Mean, Soil_Diurnal_Swing,
           Air_Tmax, Air_Tmin, Air_Tmean, Soil_1pm_minus_Air_Tmax, Soil_Mean_minus_Air_Tmean, Env_Factor)

  soil_env_summary <- soil_clean %>%
    group_by(Environment, Location) %>%
    summarise(
      N_Stages_Sampled        = n(),
      DAS_Window              = sprintf("%d-%d DAS", min(DAS), max(DAS)),
      Mean_Soil_9am           = mean(Soil_9am, na.rm = TRUE),
      Mean_Soil_1pm           = mean(Soil_1pm, na.rm = TRUE),
      Max_Soil_1pm            = max(Soil_1pm, na.rm = TRUE),
      Mean_Soil_5pm           = mean(Soil_5pm, na.rm = TRUE),
      Mean_Soil_Diurnal       = mean(Soil_Diurnal_Mean, na.rm = TRUE),
      SD_Soil_Diurnal         = sd(Soil_Diurnal_Mean, na.rm = TRUE),
      Mean_Air_Tmean_SameDays = mean(Air_Tmean, na.rm = TRUE),
      Mean_Crown_minus_Air    = mean(Soil_Mean_minus_Air_Tmean, na.rm = TRUE),
      Mean_1pm_minus_AirTmax  = mean(Soil_1pm_minus_Air_Tmax, na.rm = TRUE),
      .groups = "drop"
    )
  cat("\n--- Section 1.5: Crown-Zone Diurnal Soil Temperature Summary ---\n")
  print(soil_env_summary)

  # 1.5b Build & Export the 3 Consolidated Section 1 Tables + Multi-Tab Excel Workbook
  # (1) Main Text Table 1: Vegetative Thermal & Crown-Soil Regime Summary (Benchmarked vs NWPZ Baseline)
  env_veg_air <- stage_summary %>%
    group_by(Environment) %>%
    summarise(
      Vegetative_Duration_E_to_T4_Days = sum(Duration_Days),
      Mean_Air_Tmax_C  = sum(Tmax * Duration_Days) / sum(Duration_Days),
      Mean_Air_Tmin_C  = sum(Tmin * Duration_Days) / sum(Duration_Days),
      Mean_Air_Tmean_C = sum(Tmean * Duration_Days) / sum(Duration_Days),
      Cumulative_GDD_E_to_T4 = sum(Stage_GDD),
      Cumulative_GDD_Tmin_E_to_T4 = sum(Stage_GDD_Tmin),
      Thermal_Pace_GDD_per_Day = sum(Stage_GDD) / sum(Duration_Days),
      .groups = "drop"
    )

  base_row_t1  <- env_veg_air %>% filter(Environment == "NWPZ Baseline")
  base_soil_t1 <- soil_env_summary %>% filter(Environment == "NWPZ Baseline")

  table_1_main_regimes <- master_pheno_spec %>%
    transmute(
      Agro_Climatic_Zone = ifelse(Location == "Delhi", "NWPZ (Sub-Tropical)", "PZ (Tropical)"),
      Environment        = Env_Name,
      Location           = Location,
      Sowing_Code        = Sowing,
      Experimental_Role  = Role,
      Sowing_Date        = as.character(Sowing_Date),
      Emergence_DAS      = E,
      Max_Tillering_T4_DAS = T4,
      Stem_Elongation_SE_DAS = SE
    ) %>%
    left_join(env_veg_air, by = "Environment") %>%
    left_join(soil_env_summary %>% select(Environment, Mean_Soil_9am, Mean_Soil_1pm, Max_Soil_1pm, Mean_Soil_5pm, Mean_1pm_minus_AirTmax),
              by = "Environment") %>%
    mutate(
      Delta_Days_vs_Baseline       = Vegetative_Duration_E_to_T4_Days - base_row_t1$Vegetative_Duration_E_to_T4_Days,
      Delta_Air_Tmean_vs_Baseline  = round(Mean_Air_Tmean_C - base_row_t1$Mean_Air_Tmean_C, 2),
      Delta_GDD_Tmin_vs_Baseline   = round(Cumulative_GDD_Tmin_E_to_T4 - base_row_t1$Cumulative_GDD_Tmin_E_to_T4, 1),
      Delta_Pace_GDD_d_vs_Baseline = round(Thermal_Pace_GDD_per_Day - base_row_t1$Thermal_Pace_GDD_per_Day, 2),
      Delta_Soil_09h_vs_Baseline   = round(Mean_Soil_9am - base_soil_t1$Mean_Soil_9am, 2),
      Delta_Soil_13h_vs_Baseline   = round(Mean_Soil_1pm - base_soil_t1$Mean_Soil_1pm, 2),
      Delta_Soil_17h_vs_Baseline   = round(Mean_Soil_5pm - base_soil_t1$Mean_Soil_5pm, 2),
      Mean_Air_Tmax_C              = round(Mean_Air_Tmax_C, 2),
      Mean_Air_Tmin_C              = round(Mean_Air_Tmin_C, 2),
      Mean_Air_Tmean_C             = round(Mean_Air_Tmean_C, 2),
      Cumulative_GDD_E_to_T4       = round(Cumulative_GDD_E_to_T4, 1),
      Cumulative_GDD_Tmin_E_to_T4  = round(Cumulative_GDD_Tmin_E_to_T4, 1),
      Thermal_Pace_GDD_per_Day     = round(Thermal_Pace_GDD_per_Day, 2),
      Mean_Soil_09h_C              = round(Mean_Soil_9am, 2),
      Mean_Soil_13h_C              = round(Mean_Soil_1pm, 2),
      Max_Soil_13h_C               = round(Max_Soil_1pm, 2),
      Mean_Soil_17h_C              = round(Mean_Soil_5pm, 2),
      Crown_13h_minus_Air_Tmax_C   = round(Mean_1pm_minus_AirTmax, 2)
    ) %>%
    select(
      Agro_Climatic_Zone, Environment, Location, Sowing_Code, Experimental_Role, Sowing_Date,
      Emergence_DAS, Max_Tillering_T4_DAS, Stem_Elongation_SE_DAS,
      Vegetative_Duration_E_to_T4_Days, Delta_Days_vs_Baseline,
      Mean_Air_Tmax_C, Mean_Air_Tmin_C, Mean_Air_Tmean_C, Delta_Air_Tmean_vs_Baseline,
      Cumulative_GDD_E_to_T4, Cumulative_GDD_Tmin_E_to_T4, Delta_GDD_Tmin_vs_Baseline,
      Thermal_Pace_GDD_per_Day, Delta_Pace_GDD_d_vs_Baseline,
      Mean_Soil_09h_C, Mean_Soil_13h_C, Max_Soil_13h_C, Mean_Soil_17h_C,
      Delta_Soil_09h_vs_Baseline, Delta_Soil_13h_vs_Baseline, Delta_Soil_17h_vs_Baseline,
      Crown_13h_minus_Air_Tmax_C
    )

  # (2) Supplementary Table S1: Consolidated Stage-Specific Thermal & Disparity Dynamics (16 rows)
  table_s1_consolidated <- stage_summary %>%
    group_by(Environment) %>%
    mutate(
      Cumulative_DAE           = cumsum(Duration_Days),
      Cumulative_GDD           = cumsum(Stage_GDD),
      Cumulative_GDD_Tmin      = cumsum(Stage_GDD_Tmin),
      Thermal_Pace_GDD_per_Day = Stage_GDD / Duration_Days
    ) %>%
    ungroup() %>%
    left_join(
      baseline_st %>%
        transmute(
          Stage,
          Base_Days     = Duration_Days,
          Base_Tmean    = Tmean,
          Base_Tmax     = Tmax,
          Base_Tmin     = Tmin,
          Base_GDD      = Stage_GDD,
          Base_GDD_Tmin = Stage_GDD_Tmin,
          Base_Pace     = Pace_Baseline
        ),
      by = "Stage"
    ) %>%
    mutate(
      Delta_Days_vs_Baseline           = Duration_Days - Base_Days,
      Delta_Tmean_vs_Baseline_C        = round(Tmean - Base_Tmean, 2),
      Delta_Tmax_vs_Baseline_C         = round(Tmax - Base_Tmax, 2),
      Delta_Tmin_vs_Baseline_C         = round(Tmin - Base_Tmin, 2),
      Delta_Stage_GDD_vs_Baseline      = round(Stage_GDD - Base_GDD, 1),
      Delta_Stage_GDD_Tmin_vs_Baseline = round(Stage_GDD_Tmin - Base_GDD_Tmin, 1),
      Delta_Pace_GDD_d_vs_Baseline     = round(Thermal_Pace_GDD_per_Day - Base_Pace, 2),
      Mean_Tmean_C                     = round(Tmean, 2),
      Mean_Tmax_C                      = round(Tmax, 2),
      Mean_Tmin_C                      = round(Tmin, 2),
      Stage_GDD                        = round(Stage_GDD, 1),
      Stage_GDD_Tmin                   = round(Stage_GDD_Tmin, 1),
      Cumulative_GDD                   = round(Cumulative_GDD, 1),
      Cumulative_GDD_Tmin              = round(Cumulative_GDD_Tmin, 1),
      Thermal_Pace_GDD_per_Day         = round(Thermal_Pace_GDD_per_Day, 2)
    ) %>%
    select(
      Environment, Stage, Duration_Days, Cumulative_DAE,
      Mean_Tmean_C, Mean_Tmax_C, Mean_Tmin_C,
      Stage_GDD, Stage_GDD_Tmin, Cumulative_GDD, Cumulative_GDD_Tmin, Thermal_Pace_GDD_per_Day,
      Delta_Days_vs_Baseline, Delta_Tmean_vs_Baseline_C, Delta_Tmax_vs_Baseline_C, Delta_Tmin_vs_Baseline_C,
      Delta_Stage_GDD_vs_Baseline, Delta_Stage_GDD_Tmin_vs_Baseline, Delta_Pace_GDD_d_vs_Baseline
    )

  # (3) Supplementary Table S2: Crown-Soil Diurnal Stage Observations (14 rows)
  table_s2_crown_obs <- soil_clean %>%
    transmute(
      Environment, Location, Sowing,
      Sampling_Milestone = Stage_Milestone,
      Date               = as.character(Date),
      DAS,
      Soil_09h_C         = round(Soil_9am, 2),
      Soil_13h_C         = round(Soil_1pm, 2),
      Soil_17h_C         = round(Soil_5pm, 2),
      Soil_Diurnal_Swing_C = round(Soil_Diurnal_Swing, 2),
      SameDay_Air_Tmax_C   = round(Air_Tmax, 2),
      SameDay_Air_Tmin_C   = round(Air_Tmin, 2),
      SameDay_Air_Tmean_C  = round(Air_Tmean, 2),
      Soil_13h_minus_Air_Tmax_C = round(Soil_1pm_minus_Air_Tmax, 2)
    )

  write.csv(table_1_main_regimes, file.path(sec1_dir, "Table_1_Vegetative_Thermal_and_Crown_Soil_Regimes_v2.csv"), row.names = FALSE)
  write.csv(table_s1_consolidated, file.path(sec1_dir, "Table_S1_Stage_Specific_Thermal_and_Disparity_Dynamics_v2.csv"), row.names = FALSE)
  write.csv(table_s2_crown_obs, file.path(sec1_dir, "Table_S2_Crown_Soil_Diurnal_Stage_Observations_v2.csv"), row.names = FALSE)

  if (requireNamespace("openxlsx", quietly = TRUE)) {
    wb_sec1 <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb_sec1, "Table_1_Main_Regimes")
    openxlsx::writeData(wb_sec1, "Table_1_Main_Regimes", table_1_main_regimes, headerStyle = header_style)
    openxlsx::setColWidths(wb_sec1, "Table_1_Main_Regimes", cols = 1:ncol(table_1_main_regimes), widths = "auto")

    openxlsx::addWorksheet(wb_sec1, "Table_S1_Stage_Dynamics")
    openxlsx::writeData(wb_sec1, "Table_S1_Stage_Dynamics", table_s1_consolidated, headerStyle = header_style)
    openxlsx::setColWidths(wb_sec1, "Table_S1_Stage_Dynamics", cols = 1:ncol(table_s1_consolidated), widths = "auto")

    openxlsx::addWorksheet(wb_sec1, "Table_S2_Crown_Soil_Obs")
    openxlsx::writeData(wb_sec1, "Table_S2_Crown_Soil_Obs", table_s2_crown_obs, headerStyle = header_style)
    openxlsx::setColWidths(wb_sec1, "Table_S2_Crown_Soil_Obs", cols = 1:ncol(table_s2_crown_obs), widths = "auto")

    openxlsx::saveWorkbook(wb_sec1, file.path(sec1_dir, "Section_1_Manuscript_and_Supplementary_Tables_v2.xlsx"), overwrite = TRUE)
  }

  # 1.6 Rebuilt High-Legibility Integrated Thermo-Phenology Figure (V4_Integrated_ThermoPhenology_Revised.png)
  # Combines continuous daily Air (Tmax, Tmin) from Sowing (0 DAS) to Stem Elongation (SE)
  # with the 3 discrete diurnal crown-zone soil temperature timestamps (09:00 h, 13:00 h, 17:00 h)
  # and stage milestones (E, T1, T2, T3, T4, SE) without a mean-soil line, subtitle, or reproductive 30°C line.
  weather_das_list <- list()
  milestone_list   <- list()
  for (i in 1:nrow(master_pheno_spec)) {
    loc_i    <- master_pheno_spec$Location[i]
    e_name_i <- master_pheno_spec$Env_Name[i]
    s_date_i <- master_pheno_spec$Sowing_Date[i]
    se_das_i <- master_pheno_spec$SE[i]

    w_sub <- weather %>%
      filter(Location == loc_i, Date >= s_date_i, Date <= (s_date_i + se_das_i)) %>%
      arrange(Date) %>%
      mutate(
        DAS      = as.numeric(Date - s_date_i),
        Env_Name = e_name_i
      )
    weather_das_list[[i]] <- w_sub

    ms_df <- data.frame(
      Env_Name   = e_name_i,
      Stage_Code = c("E", "T1", "T2", "T3", "T4", "SE"),
      DAS        = c(master_pheno_spec$E[i], master_pheno_spec$T1[i], master_pheno_spec$T2[i],
                     master_pheno_spec$T3[i], master_pheno_spec$T4[i], master_pheno_spec$SE[i]),
      y_badge    = c(42.8, 42.8, 42.8, 42.8, 42.8, 42.8),
      stringsAsFactors = FALSE
    )
    milestone_list[[i]] <- ms_df
  }
  weather_das_df <- bind_rows(weather_das_list) %>%
    mutate(Env_Factor = factor(Env_Name, levels = env_metadata$Env_Name))
  milestones_df <- bind_rows(milestone_list) %>%
    mutate(
      Env_Factor = factor(Env_Name, levels = env_metadata$Env_Name),
      Badge_Text = sprintf("%s\n(%dd)", Stage_Code, DAS)
    )

  # Compute LOESS-smoothed Tmax and Tmin per environment so the shaded envelope matches the curves without gaps
  weather_das_df <- weather_das_df %>%
    group_by(Env_Factor) %>%
    mutate(
      Tmax_Smooth = predict(loess(Tmax ~ DAS, span = 0.25)),
      Tmin_Smooth = predict(loess(Tmin ~ DAS, span = 0.25))
    ) %>%
    ungroup()

  # Long-format diurnal soil + air lines (raw daily + LOESS trend) for unified legend
  air_raw_long <- weather_das_df %>%
    select(Env_Factor, DAS,
           `Air Tmax (Daily Max)` = Tmax,
           `Air Tmin (Daily Min)` = Tmin) %>%
    pivot_longer(
      cols = c(`Air Tmax (Daily Max)`, `Air Tmin (Daily Min)`),
      names_to = "Series",
      values_to = "Temp_C"
    )

  air_lines_long <- weather_das_df %>%
    select(Env_Factor, DAS,
           `Air Tmax (Daily Max)` = Tmax_Smooth,
           `Air Tmin (Daily Min)` = Tmin_Smooth) %>%
    pivot_longer(
      cols = c(`Air Tmax (Daily Max)`, `Air Tmin (Daily Min)`),
      names_to = "Series",
      values_to = "Temp_C"
    )

  soil_lines_long <- soil_clean %>%
    select(Env_Factor, DAS,
           `Soil 13:00 h (Afternoon Peak)` = Soil_1pm,
           `Soil 17:00 h (Late Afternoon)` = Soil_5pm,
           `Soil 09:00 h (Morning)`        = Soil_9am) %>%
    pivot_longer(
      cols = c(`Soil 13:00 h (Afternoon Peak)`, `Soil 17:00 h (Late Afternoon)`, `Soil 09:00 h (Morning)`),
      names_to = "Series",
      values_to = "Temp_C"
    )

  series_levels <- c(
    "Air Tmax (Daily Max)",
    "Air Tmin (Daily Min)",
    "Soil 13:00 h (Afternoon Peak)",
    "Soil 17:00 h (Late Afternoon)",
    "Soil 09:00 h (Morning)"
  )
  air_raw_long$Series    <- factor(air_raw_long$Series, levels = series_levels)
  air_lines_long$Series  <- factor(air_lines_long$Series, levels = series_levels)
  soil_lines_long$Series <- factor(soil_lines_long$Series, levels = series_levels)

  series_colors <- c(
    "Air Tmax (Daily Max)"          = "#C0392B",
    "Air Tmin (Daily Min)"          = "#2471A3",
    "Soil 13:00 h (Afternoon Peak)" = "#D35400",
    "Soil 17:00 h (Late Afternoon)" = "#7D3C98",
    "Soil 09:00 h (Morning)"        = "#117A65"
  )

  p_v4_thermo <- ggplot() +
    # Unbroken Daily Air Thermal Envelope (raw daily Tmin-Tmax envelope)
    geom_ribbon(
      data = weather_das_df,
      aes(x = DAS, ymin = Tmin, ymax = Tmax),
      fill = "#D5D8DC", alpha = 0.45
    ) +
    # Vertical Phenological Stage Milestones (E, T1, T2, T3, T4, SE)
    geom_vline(
      data = milestones_df,
      aes(xintercept = DAS),
      linetype = "dotted", color = "grey30", linewidth = 0.65
    ) +
    # Raw Daily Air Tmax & Tmin Traces (thin semi-transparent lines showing daily spikes)
    geom_line(
      data = air_raw_long,
      aes(x = DAS, y = Temp_C, color = Series, group = Series),
      linewidth = 0.45, alpha = 0.55, show.legend = FALSE
    ) +
    # LOESS-Smoothed Daily Air Tmax & Tmin Trend Curves
    geom_line(
      data = air_lines_long,
      aes(x = DAS, y = Temp_C, color = Series, group = Series),
      linewidth = 1.15
    ) +
    # Vertical Diurnal Crown-Soil Range Bars at Each Sampling Date
    geom_segment(
      data = soil_clean,
      aes(x = DAS, xend = DAS, y = pmin(Soil_9am, Soil_5pm), yend = Soil_1pm),
      color = "grey25", linewidth = 0.75, alpha = 0.75
    ) +
    # 3-Timestamp Diurnal Soil Temperature Trajectories & Points (09:00 h, 13:00 h, 17:00 h)
    geom_line(
      data = soil_lines_long,
      aes(x = DAS, y = Temp_C, color = Series, group = Series),
      linewidth = 0.90, linetype = "longdash"
    ) +
    geom_point(
      data = soil_lines_long,
      aes(x = DAS, y = Temp_C, color = Series, shape = Series),
      size = 3.1, stroke = 0.9
    ) +
    # Direct Temperature Callouts for 13:00 h Peak Soil Temperature
    geom_label(
      data = soil_clean,
      aes(x = DAS, y = pmax(Soil_9am, Soil_1pm, Soil_5pm) + 2.5, label = sprintf("%.1f°C", Soil_1pm)),
      size = 3.1, fontface = "bold", color = "#D35400", fill = "white",
      label.padding = unit(0.14, "lines"), label.size = 0.25, show.legend = FALSE
    ) +
    # High-Legibility Horizontal Stage Badges at Top of Panel
    geom_label(
      data = milestones_df,
      aes(x = DAS, y = y_badge, label = Badge_Text),
      size = 3.0, fontface = "bold", color = "#1B2631", fill = "#F8F9F9",
      label.padding = unit(0.16, "lines"), label.size = 0.3, lineheight = 0.88
    ) +
    facet_wrap(~ Env_Factor, scales = "free_x", ncol = 2) +
    scale_color_manual(values = series_colors, name = "Thermo-Phenological Series:") +
    scale_shape_manual(
      values = c(
        "Air Tmax (Daily Max)"          = 16,
        "Air Tmin (Daily Min)"          = 16,
        "Soil 13:00 h (Afternoon Peak)" = 17,
        "Soil 17:00 h (Late Afternoon)" = 15,
        "Soil 09:00 h (Morning)"        = 18
      ),
      guide = "none"
    ) +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.06))) +
    scale_y_continuous(limits = c(0, 46), breaks = seq(0, 40, by = 5)) +
    pub_theme +
    theme(
      strip.background = element_rect(fill = "#2C3E50", color = "#1B2631", linewidth = 0.8),
      strip.text = element_text(color = "white", face = "bold", size = 13, margin = margin(t = 5, b = 5)),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.35),
      plot.title = element_text(face = "bold", size = 16, margin = margin(b = 10)),
      axis.title = element_text(face = "bold", size = 13),
      axis.text = element_text(color = "black", size = 11, face = "bold"),
      legend.position = "bottom",
      legend.title = element_text(face = "bold", size = 11.5),
      legend.text = element_text(size = 11, face = "bold"),
      legend.box.background = element_rect(color = "grey75", fill = "#F8F9F9", linewidth = 0.5),
      legend.margin = margin(t = 6, r = 10, b = 6, l = 10)
    ) +
    guides(color = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(linewidth = 1.4, size = 3.2))) +
    labs(
      title = "Integrated Air and Diurnal Crown-Soil Thermo-Phenological Regimes",
      x = "Days After Sowing (DAS)",
      y = "Temperature (°C)"
    )

  ggsave(file.path(sec1_dir, "Figure_S1_Integrated_ThermoPhenology_v2.png"),
         plot = p_v4_thermo, width = 14, height = 9.5, dpi = 600)

  # ==============================================================================
  # ==============================================================================
  # 1.7 Multi-Panel Temperature & Thermal Disparities Figure (v2: Split Tmax/Tmin + Nighttime Heat Units)
  # ==============================================================================
  cum_gdd_summary <- bind_rows(
    tibble(
      Environment = env_metadata$Env_Name,
      Milestone = "E",
      Cumulative_DAE = 0,
      Cum_GDD_Uncapped = 0,
      Cum_GDD_Tmin = 0
    ),
    table_s1_consolidated %>%
      mutate(Milestone = recode(Stage, "E-T1" = "T1", "T1-T2" = "T2", "T2-T3" = "T3", "T3-T4" = "T4")) %>%
      select(Environment, Milestone, Cumulative_DAE, Cum_GDD_Uncapped = Cumulative_GDD, Cum_GDD_Tmin = Cumulative_GDD_Tmin)
  )

  disp_plot_df <- disparities %>%
    mutate(
      Stage = factor(Stage, levels = c("E-T1", "T1-T2", "T2-T3", "T3-T4")),
      Env_Factor = factor(Environment, levels = c("NWPZ Warm", "PZ Warm I", "PZ Warm II")),
      Label_Tmax = sprintf("ΔTmax: %+.1f°C\n(Tmax: %.1f°C)", Delta_Tmax, Tmax),
      Label_Tmin = sprintf("ΔTmin: %+.1f°C\n(Tmin: %.1f°C)", Delta_Tmin, Tmin)
    )

  common_heat_theme <- pub_theme +
    theme(
      plot.title        = element_text(face = "bold", size = 12.8, hjust = 0, margin = margin(b = 2)),
      plot.subtitle     = element_text(size = 9.8, color = "grey25", face = "italic", margin = margin(b = 5)),
      axis.title        = element_text(face = "bold", size = 11.5, color = "black"),
      axis.text         = element_text(size = 10.6, face = "bold", color = "black"),
      legend.position   = "bottom",
      legend.title      = element_text(face = "bold", size = 10.8),
      legend.text       = element_text(size = 10.2, face = "bold"),
      legend.margin     = margin(t = 2, b = 2),
      plot.margin       = margin(6, 10, 4, 8)
    )

  heat_cbar <- guide_colorbar(
    barwidth = unit(13.5, "cm"),
    barheight = unit(0.46, "cm"),
    title.position = "left",
    title.vjust = 0.85,
    frame.colour = "black",
    ticks.colour = "black"
  )

  p_disp_tmax <- ggplot(disp_plot_df, aes(x = Stage, y = Env_Factor, fill = Delta_Tmax)) +
    geom_tile(color = "white", linewidth = 1.2) +
    geom_text(aes(label = Label_Tmax), fontface = "bold", size = 3.90, lineheight = 1.05,
              color = ifelse(abs(disp_plot_df$Delta_Tmax) > 5.0, "white", "black")) +
    scale_fill_gradient2(low = "#3C5488", mid = "#F7F7F7", high = "#DC0000", midpoint = 0,
                         limits = c(-5, 13), breaks = seq(-4, 12, by = 4),
                         name = "Thermal Disparity vs. NWPZ Baseline (°C): ",
                         guide = heat_cbar) +
    scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
    common_heat_theme +
    labs(title = "A. Daytime Maximum Air Temp. (ΔTmax)",
         subtitle = "Uniform daytime heat load across PZ Warm I & II (~28.0 °C)",
         x = "Developmental Window", y = "Stress Regime")

  p_disp_tmin <- ggplot(disp_plot_df, aes(x = Stage, y = Env_Factor, fill = Delta_Tmin)) +
    geom_tile(color = "white", linewidth = 1.2) +
    geom_text(aes(label = Label_Tmin), fontface = "bold", size = 3.90, lineheight = 1.05,
              color = ifelse(abs(disp_plot_df$Delta_Tmin) > 5.0, "white", "black")) +
    scale_fill_gradient2(low = "#3C5488", mid = "#F7F7F7", high = "#DC0000", midpoint = 0,
                         limits = c(-5, 13), breaks = seq(-4, 12, by = 4),
                         name = "Thermal Disparity vs. NWPZ Baseline (°C): ",
                         guide = heat_cbar) +
    scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
    common_heat_theme +
    labs(title = "B. Nighttime Minimum Air Temp. (ΔTmin)",
         subtitle = "Nocturnal surge at T2–T3 in PZ Warm I (+8.9 °C) vs. PZ Warm II (+4.0 °C)",
         x = "Developmental Window", y = "")

  soil_disp_df <- bind_rows(
    table_1_main_regimes %>% filter(Environment != "NWPZ Baseline") %>%
      transmute(Environment, Timestamp = "09:00 h\n(Morning)", Temp_C = Mean_Soil_09h_C, Delta_C = Delta_Soil_09h_vs_Baseline),
    table_1_main_regimes %>% filter(Environment != "NWPZ Baseline") %>%
      transmute(Environment, Timestamp = "13:00 h\n(Afternoon Peak)", Temp_C = Mean_Soil_13h_C, Delta_C = Delta_Soil_13h_vs_Baseline),
    table_1_main_regimes %>% filter(Environment != "NWPZ Baseline") %>%
      transmute(Environment, Timestamp = "17:00 h\n(Late Afternoon)", Temp_C = Mean_Soil_17h_C, Delta_C = Delta_Soil_17h_vs_Baseline)
  ) %>%
    mutate(
      Timestamp = factor(Timestamp, levels = c("09:00 h\n(Morning)", "13:00 h\n(Afternoon Peak)", "17:00 h\n(Late Afternoon)")),
      Env_Factor = factor(Environment, levels = c("NWPZ Warm", "PZ Warm I", "PZ Warm II")),
      Label_Soil = sprintf("ΔSoil: %+.1f°C\n(%.1f°C)", Delta_C, Temp_C)
    )

  p_disp_soil_heat <- ggplot(soil_disp_df, aes(x = Timestamp, y = Env_Factor, fill = Delta_C)) +
    geom_tile(color = "white", linewidth = 1.2) +
    geom_text(aes(label = Label_Soil), fontface = "bold", size = 3.90, lineheight = 1.05,
              color = ifelse(abs(soil_disp_df$Delta_C) > 5.0, "white", "black")) +
    scale_fill_gradient2(low = "#3C5488", mid = "#F7F7F7", high = "#DC0000", midpoint = 0,
                         limits = c(-5, 13), breaks = seq(-4, 12, by = 4),
                         name = "Thermal Disparity vs. NWPZ Baseline (°C): ",
                         guide = heat_cbar) +
    scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
    common_heat_theme +
    labs(title = "C. Diurnal Crown-Soil Disparities",
         subtitle = "Root-crown thermal shifts (09:00, 13:00, 17:00 h) vs. Baseline",
         x = "Diurnal Sampling Timestamp", y = "")

  gdd_traj_df <- cum_gdd_summary %>%
    mutate(
      Milestone = factor(Milestone, levels = c("E", "T1", "T2", "T3", "T4")),
      Env_Factor = factor(Environment, levels = env_metadata$Env_Name),
      x_num = as.numeric(Milestone)
    )

  t4_end_df <- gdd_traj_df %>%
    filter(Milestone == "T4") %>%
    mutate(
      End_Label_GDD = sprintf("%s (%dd, %.0f °C·d)", Environment, Cumulative_DAE, Cum_GDD_Uncapped),
      End_Label_NHU = sprintf("%.1f °C·d", Cum_GDD_Tmin),
      x_nudge_gdd = case_when(
        Environment == "PZ Warm I"     ~ Cumulative_DAE - 1.5,
        Environment == "PZ Warm II"    ~ Cumulative_DAE + 1.5,
        Environment == "NWPZ Warm"     ~ Cumulative_DAE - 1.0,
        Environment == "NWPZ Baseline" ~ Cumulative_DAE + 0.5
      ),
      y_nudge_gdd = case_when(
        Environment == "PZ Warm I"     ~ Cum_GDD_Uncapped + 35,
        Environment == "PZ Warm II"    ~ Cum_GDD_Uncapped - 35,
        Environment == "NWPZ Warm"     ~ Cum_GDD_Uncapped + 32,
        Environment == "NWPZ Baseline" ~ Cum_GDD_Uncapped - 32
      ),
      y_nudge_nhu = case_when(
        Environment == "PZ Warm I"     ~ Cum_GDD_Tmin + 18,
        Environment == "PZ Warm II"    ~ Cum_GDD_Tmin - 18,
        Environment == "NWPZ Baseline" ~ Cum_GDD_Tmin + 16,
        Environment == "NWPZ Warm"     ~ Cum_GDD_Tmin - 16
      )
    )

  # Dynamic annotation metrics computed directly from t4_end_df & table_s1_consolidated
  gdd_min <- min(t4_end_df$Cum_GDD_Uncapped)
  gdd_max <- max(t4_end_df$Cum_GDD_Uncapped)
  gdd_spread_pct <- 100 * (gdd_max - gdd_min) / gdd_min
  dae_base <- t4_end_df$Cumulative_DAE[t4_end_df$Environment == "NWPZ Baseline"]
  dae_pz_min <- min(t4_end_df$Cumulative_DAE[t4_end_df$Environment %in% c("PZ Warm I", "PZ Warm II")])
  dae_pz_max <- max(t4_end_df$Cumulative_DAE[t4_end_df$Environment %in% c("PZ Warm I", "PZ Warm II")])
  gdd_pz1 <- t4_end_df$Cum_GDD_Uncapped[t4_end_df$Environment == "PZ Warm I"]
  gdd_pz2 <- t4_end_df$Cum_GDD_Uncapped[t4_end_df$Environment == "PZ Warm II"]

  panel_d_annot <- sprintf(
    paste0("Scheduled Sampling Thermal Alignment (E–T4):\n",
           "• Sampling milestones align at T4 within %.0f–%.0f °C·d (%.1f%% spread)\n",
           "• Tropical PZ compresses calendar duration by %+d to %+d d (%.0f%%)\n",
           "• Standard GDD_std masks divergence between PZ Warm I (%.1f) & II (%.1f)"),
    gdd_min, gdd_max, gdd_spread_pct,
    dae_pz_min - dae_base, dae_pz_max - dae_base,
    100 * (mean(c(dae_pz_min, dae_pz_max)) - dae_base) / dae_base,
    gdd_pz1, gdd_pz2
  )

  nhu_t12_pz1 <- table_s1_consolidated$Stage_GDD_Tmin[table_s1_consolidated$Environment == "PZ Warm I" & table_s1_consolidated$Stage == "T1-T2"]
  nhu_t12_pz2 <- table_s1_consolidated$Stage_GDD_Tmin[table_s1_consolidated$Environment == "PZ Warm II" & table_s1_consolidated$Stage == "T1-T2"]
  nhu_t23_pz1 <- table_s1_consolidated$Stage_GDD_Tmin[table_s1_consolidated$Environment == "PZ Warm I" & table_s1_consolidated$Stage == "T2-T3"]
  nhu_t23_pz2 <- table_s1_consolidated$Stage_GDD_Tmin[table_s1_consolidated$Environment == "PZ Warm II" & table_s1_consolidated$Stage == "T2-T3"]
  nhu_tot_pz1 <- t4_end_df$Cum_GDD_Tmin[t4_end_df$Environment == "PZ Warm I"]
  nhu_tot_pz2 <- t4_end_df$Cum_GDD_Tmin[t4_end_df$Environment == "PZ Warm II"]

  panel_e_annot <- sprintf(
    paste0("PZ Warm I vs. PZ Warm II Nighttime Heat Surge (NHU):\n",
           "• Stage T1–T2: %+.1f °C·d (%+.1f%%; %.1f vs. %.1f °C·d)\n",
           "• Stage T2–T3: %+.1f °C·d (%+.1f%%; %.1f vs. %.1f °C·d)\n",
           "• E–T4 Total: %+.1f °C·d (%.1f vs. %.1f °C·d)"),
    nhu_t12_pz1 - nhu_t12_pz2, 100 * (nhu_t12_pz1 - nhu_t12_pz2) / nhu_t12_pz2, nhu_t12_pz1, nhu_t12_pz2,
    nhu_t23_pz1 - nhu_t23_pz2, 100 * (nhu_t23_pz1 - nhu_t23_pz2) / nhu_t23_pz2, nhu_t23_pz1, nhu_t23_pz2,
    nhu_tot_pz1 - nhu_tot_pz2, nhu_tot_pz1, nhu_tot_pz2
  )

  p_disp_gdd <- ggplot(gdd_traj_df, aes(x = Cumulative_DAE, y = Cum_GDD_Uncapped, color = Env_Factor, linetype = Env_Factor, group = Env_Factor)) +
    annotate("rect", xmin = 0, xmax = 78, ymin = floor(gdd_min - 4), ymax = ceiling(gdd_max + 4), fill = "grey85", alpha = 0.45) +
    geom_hline(yintercept = c(round(gdd_min), round(gdd_max)), linetype = "dotted", color = "grey45", linewidth = 0.55) +
    geom_line(linewidth = 1.35) +
    geom_point(size = 3.3, shape = 21, fill = "white", stroke = 1.2) +
    geom_label(data = t4_end_df, aes(x = x_nudge_gdd, y = y_nudge_gdd, label = End_Label_GDD, color = Env_Factor),
               size = 3.05, fontface = "bold", fill = "white", alpha = 0.94, label.size = 0.25, show.legend = FALSE) +
    annotate("label", x = 3, y = 585, hjust = 0, vjust = 1, size = 3.05, fontface = "bold",
             fill = "#F8F9F9", color = "grey15", label.size = 0.3,
             label = panel_d_annot) +
    scale_color_npg(name = "Thermal Regime:") +
    scale_linetype_manual(values = c("NWPZ Baseline" = "solid", "NWPZ Warm" = "longdash", "PZ Warm I" = "solid", "PZ Warm II" = "twodash"), guide = "none") +
    scale_x_continuous(breaks = seq(0, 75, by = 10), limits = c(0, 78)) +
    scale_y_continuous(limits = c(0, 750), breaks = seq(0, 700, by = 100), expand = c(0, 0)) +
    common_heat_theme +
    labs(title = "D. Standard 24-h Cumulative Thermal Time (GDD_std) vs. Calendar Duration",
         subtitle = "Standard GDD_std = Σ max(0, (Tmax + Tmin)/2 - 5.0°C) aligns at T4 sampling (~644–682 °C·d) across all 4 regimes",
         x = "Chronological Duration from Emergence (Days After Emergence, DAE)",
         y = "Cumulative Standard GDD_std (°C·days)")

  p_disp_gdd_tmin <- ggplot(gdd_traj_df, aes(x = x_num, y = Cum_GDD_Tmin, color = Env_Factor, linetype = Env_Factor, group = Env_Factor)) +
    annotate("rect", xmin = 2.0, xmax = 4.0, ymin = 0, ymax = 420, fill = "#DC0000", alpha = 0.07) +
    annotate("text", x = 3.0, y = 402, label = "Active Tillering Window (T1–T3): Nocturnal Rd Divergence",
             color = "#922B21", fontface = "bold.italic", size = 3.2) +
    geom_line(linewidth = 1.35) +
    geom_point(size = 3.3, shape = 21, fill = "white", stroke = 1.2) +
    geom_label(data = t4_end_df, aes(x = 5.06, y = y_nudge_nhu, label = End_Label_NHU, color = Env_Factor),
               hjust = 0, size = 3.05, fontface = "bold", fill = "white", alpha = 0.94, label.size = 0.25, show.legend = FALSE) +
    annotate("label", x = 1.12, y = 378, hjust = 0, vjust = 1, size = 3.05, fontface = "bold",
             fill = "#FDEDEC", color = "#922B21", label.size = 0.35,
             label = panel_e_annot) +
    annotate("segment", x = 4.0, xend = 4.0, y = 262.0, yend = 297.5,
             colour = "#DC0000", linewidth = 0.9, arrow = arrow(ends = "both", length = unit(0.15, "cm"))) +
    scale_color_npg(name = "Thermal Regime:") +
    scale_linetype_manual(values = c("NWPZ Baseline" = "solid", "NWPZ Warm" = "longdash", "PZ Warm I" = "solid", "PZ Warm II" = "twodash"), guide = "none") +
    scale_x_continuous(breaks = 1:5, labels = c("E\n(Emergence)", "T1\n(Init.)", "T2\n(Early)", "T3\n(Active)", "T4\n(Max Tiller)"), limits = c(0.9, 5.55)) +
    scale_y_continuous(limits = c(0, 420), breaks = seq(0, 400, by = 100), expand = c(0, 0)) +
    common_heat_theme +
    labs(title = "E. Cumulative Nighttime Heat Units (NHU, Base 5.0 °C) Sensitivity",
         subtitle = sprintf("NHU = Σ max(0, Tmin - 5.0°C) exposes %+.1f%% nighttime heat overload during T2–T3 in PZ Warm I", 100 * (nhu_t23_pz1 - nhu_t23_pz2) / nhu_t23_pz2),
         x = "Developmental Milestone (Emergence to Max Vegetative Tillering)",
         y = "Cumulative Nighttime Heat Units, NHU (°C·days)")

  top_row <- (p_disp_tmax | p_disp_tmin | p_disp_soil_heat) +
    plot_layout(widths = c(1, 1, 0.92), guides = "collect") &
    theme(legend.position = "bottom", legend.margin = margin(t = 2, b = 2))

  bottom_row <- (p_disp_gdd | p_disp_gdd_tmin) +
    plot_layout(widths = c(1.08, 1.0), guides = "collect") &
    theme(legend.position = "bottom", legend.title = element_text(face = "bold", size = 10.8), legend.text = element_text(face = "bold", size = 10.2))

  fig_part1 <- top_row +
    plot_annotation(
      caption = "Panels A–C: Stage-specific daytime maximum air (ΔTmax), nighttime minimum air (ΔTmin), and diurnal crown-soil temperature disparities relative to NWPZ Baseline (shared diverging scale).",
      theme = theme(plot.caption = element_text(size = 9.8, color = "grey25", face = "italic", hjust = 0, lineheight = 1.15, margin = margin(t = 2)))
    )
  ggsave(file.path(sec1_dir, "Figure_1_Part1_Diurnal_Air_and_CrownSoil_Heatmaps_v2.png"), plot = fig_part1, width = 17.5, height = 4.4, dpi = 600)

  fig_part2 <- bottom_row +
    plot_annotation(
      caption = "Panel D: Standard cumulative GDD_std trajectory showing conserved E–T4 thermal demand (644–682 °C·d). Panel E: Cumulative Nighttime Heat Units (NHU = Σ max(0, Tmin - 5.0°C)) revealing the +27.8% (T1–T2) and +91.7% (T2–T3) nighttime heat surge in PZ Warm I vs. PZ Warm II.",
      theme = theme(plot.caption = element_text(size = 9.8, color = "grey25", face = "italic", hjust = 0, lineheight = 1.15, margin = margin(t = 2)))
    )
  ggsave(file.path(sec1_dir, "Figure_1_Part2_GDD_and_Nighttime_Heat_Accumulation_v2.png"), plot = fig_part2, width = 17.5, height = 6.2, dpi = 600)

  fig_v2 <- (top_row / bottom_row) +
    plot_layout(heights = c(0.80, 1.15)) +
    plot_annotation(
      caption = paste0(
        "Panels A–C: Stage-specific daytime maximum air (ΔTmax), nighttime minimum air (ΔTmin), and diurnal crown-soil temperature disparities relative to NWPZ Baseline (shared diverging scale).\n",
        "Panel D: Standard cumulative GDD_std trajectory showing conserved E–T4 thermal demand (644–682 °C·d). Panel E: Cumulative Nighttime Heat Units (NHU = Σ max(0, Tmin - 5.0°C)) revealing the +27.8% (T1–T2) and +91.7% (T2–T3) nighttime heat surge in PZ Warm I vs. PZ Warm II."
      ),
      theme = theme(plot.caption = element_text(size = 9.8, color = "grey25", face = "italic", hjust = 0, lineheight = 1.15))
    )

  ggsave(file.path(sec1_dir, "Figure_Thermal_Disparities_and_Stage_Dynamics_v2.png"), plot = fig_v2, width = 18.0, height = 10.4, dpi = 600)
} else {
  cat("  -> Section 1 [SKIPPED - Assets Finalized & Protected]\n")
}


# ==============================================================================
# SECTION 2: INTRA-ENVIRONMENT BLOCK ADJUSTMENT & BLUE EXTRACTION
# ==============================================================================
cat(">>> Executing Section 2: Intra-Environment Linear Mixed Models & BLUE Extraction...\n")

# ----------------------------------------------------------------------------
# 2.1 Non-Destructive Load of Master_Data.xlsx & Merge of Block_true for Delhi
#     (Executed unconditionally so master_raw & cryptographic cache_key exist
#      even when RUN_SECTION_2_BLUES = FALSE)
# ----------------------------------------------------------------------------
stopifnot(file.exists(data_file_path), file.exists(gs_map_path))

data_md5_Full <- as.character(tools::md5sum(data_file_path))
gs_md5_Full   <- as.character(tools::md5sum(gs_map_path))

gs_raw <- read_excel(gs_map_path, sheet = "Sheet1")
names(gs_raw) <- c("Original_Sr_No", "Genotype_Raw", "Genotype_No", "Rep_Raw")
gs_raw <- gs_raw %>%
  mutate(
    Original_Sr_No = as.integer(Original_Sr_No),
    Genotype_Clean = trimws(as.character(Genotype_Raw)),
    Genotype_No    = as.integer(Genotype_No),
    Replication    = as.character(as.integer(Rep_Raw))
  )

# Load Raw Field Observations from Master_Data.xlsx without modifying file on disk
trial_list <- list()
for (s in env_metadata$Sheet) {
  raw_df <- read_excel(data_file_path, sheet = s)
  names(raw_df) <- gsub("Tmax", "T4", names(raw_df))
  raw_df <- raw_df %>% filter(!is.na(Genotype))
  raw_df$Sheet <- s
  raw_df$Environment <- env_metadata$Env_Name[env_metadata$Sheet == s]
  raw_df$Row_Order_In_Sheet <- seq_len(nrow(raw_df))
  raw_df$Genotype_Clean <- trimws(as.character(raw_df$Genotype))
  raw_df$Rep_Char       <- as.character(as.integer(raw_df$Replication))
  trial_list[[s]] <- raw_df
}
master_raw <- bind_rows(trial_list)

# Compute Block_true dynamically from recorded counts in Delhi sheets
del_sub_counts <- master_raw %>% filter(Sheet == env_metadata$Sheet[env_metadata$Location == "Delhi"][1])
n_plots_per_rep_del  <- length(unique(del_sub_counts$Genotype_Clean))
n_blocks_per_rep_del <- length(unique(del_sub_counts$Block))
n_plots_per_blk_del  <- n_plots_per_rep_del %/% n_blocks_per_rep_del

gs_raw <- gs_raw %>%
  mutate(
    Block_From_Sr = ((Original_Sr_No - 1L) %% n_plots_per_rep_del) %/% n_plots_per_blk_del + 1L
  )

# Lookup maps for R1 (canonical packet order) and R2 (randomized packet order)
gs_r1_lookup <- gs_raw %>%
  filter(Replication == unique(Replication)[1]) %>%
  select(Genotype_Clean, Genotype_No, Sr_R1 = Original_Sr_No, Block_Pos_GS_R1 = Block_From_Sr)
gs_r2_lookup <- gs_raw %>%
  filter(Replication == unique(Replication)[2]) %>%
  select(Genotype_Clean, Sr_R2 = Original_Sr_No, Block_Pos_GS_R2 = Block_From_Sr)

master_raw <- master_raw %>%
  left_join(
    gs_raw %>% select(Genotype_Clean, Replication, Original_Sr_No, Block_From_Sr),
    by = c("Genotype_Clean", "Rep_Char" = "Replication")
  ) %>%
  left_join(gs_r1_lookup, by = "Genotype_Clean") %>%
  left_join(gs_r2_lookup, by = "Genotype_Clean") %>%
  mutate(
    Location = env_metadata$Location[match(Sheet, env_metadata$Sheet)],
    Block_true = ifelse(Location == "Delhi", as.character(Block_From_Sr), NA_character_),
    # Diagnostic field position block: in Dharwad, S1 followed GS R1 packet order and S2 followed GS R2 packet order
    Field_Pos_Block_Diagnostic = case_when(
      Sheet %in% c("Delhi_S2", "Delhi_S1") ~ as.integer(Block_From_Sr),
      Sheet == "Dharwad_S1"                ~ as.integer(Block_Pos_GS_R1),
      Sheet == "Dharwad_S2"                ~ as.integer(Block_Pos_GS_R2),
      TRUE                                 ~ NA_integer_
    )
  )

master_raw$Genotype    <- as.factor(master_raw$Genotype_Clean)
master_raw$Replication <- as.factor(master_raw$Rep_Char)
master_raw$Block       <- as.factor(master_raw$Block)
master_raw$Block_true  <- as.factor(master_raw$Block_true)
master_raw$Env_Factor  <- factor(master_raw$Environment, levels = env_metadata$Env_Name)

stopifnot(sum(is.na(master_raw$Block_true[master_raw$Location == "Delhi"])) == 0L)

# Worker functions & cryptographic cache key (shared by both run and skip branches)
traits_to_run <- c("T1_Count", "T2_Count", "T3_Count", "T4_Count", "PT_Count", "GrainYld", "TKW", "AUTPC", "AUTPC_norm")

fit_delhi_model_A <- function(df_sub, tr_name) {
  f_mod <- as.formula(paste0(tr_name, " ~ Genotype + Replication + (1 | Replication:Block_true)"))
  m_fit <- suppressMessages(lmer(f_mod, data = df_sub, REML = TRUE, control = lmerControl(optimizer = "bobyqa")))
  emm_obj <- suppressMessages(emmeans(m_fit, "Genotype", mode = "kenward-roger"))
  emm_df  <- as.data.frame(emm_obj)
  df_meth <- if (all(is.infinite(emm_df$df))) "asymptotic" else "kenward-roger"
  if (df_meth == "asymptotic") {
    stop(sprintf("CRITICAL: emmeans silently fell back to asymptotic df in Delhi Model A for trait %s", tr_name))
  }
  vc_df <- as.data.frame(VarCorr(m_fit))
  vb <- vc_df$vcov[vc_df$grp == "Replication:Block_true"]
  ve <- vc_df$vcov[vc_df$grp == "Residual"]
  conv_msgs <- m_fit@optinfo$conv$lme4$messages
  conv_str  <- if (is.null(conv_msgs)) "OK" else paste(conv_msgs, collapse = "; ")
  list(
    model = m_fit,
    emm_df = emm_df,
    df_method = df_meth,
    is_singular = isSingular(m_fit),
    conv_msg = conv_str,
    var_block = vb,
    var_resid = ve,
    formula_str = paste0(tr_name, " ~ Genotype + Replication + (1 | Replication:Block_true)"),
    engine_str = "lme4::lmer (REML)"
  )
}

fit_rcbd_lm <- function(df_sub, tr_name) {
  f_mod <- as.formula(paste0(tr_name, " ~ Genotype + Replication"))
  m_fit <- lm(f_mod, data = df_sub)
  emm_df <- as.data.frame(emmeans(m_fit, "Genotype"))
  ve <- sigma(m_fit)^2
  list(
    model = m_fit,
    emm_df = emm_df,
    df_method = "exact-ols-t",
    is_singular = NA,
    conv_msg = "OK",
    var_block = NA_real_,
    var_resid = ve,
    formula_str = paste0(tr_name, " ~ Genotype + Replication"),
    engine_str = "stats::lm (OLS)"
  )
}

cache_dir <- file.path(project_root, "Cache_BLUEs")
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

cache_payload <- list(
  master_data_md5 = data_md5_Full,
  gs_map_md5      = gs_md5_Full,
  pheno_spec_hash = digest::digest(master_pheno_spec, algo = "md5"),
  traits_vec      = traits_to_run,
  delhi_choice    = DELHI_PRIMARY_MODEL,
  worker_A_form   = deparse(formals(fit_delhi_model_A)),
  worker_A_body   = deparse(body(fit_delhi_model_A)),
  worker_B_form   = deparse(formals(fit_rcbd_lm)),
  worker_B_body   = deparse(body(fit_rcbd_lm)),
  r_version       = R.version.string,
  lme4_version    = as.character(packageVersion("lme4")),
  emmeans_version = as.character(packageVersion("emmeans"))
)
cache_key <- paste0(
  substr(data_md5_Full, 1, 8), "_",
  substr(gs_md5_Full, 1, 8), "_",
  substr(digest::digest(cache_payload, algo = "md5"), 1, 12)
)

if (isTRUE(RUN_SECTION_2_BLUES)) {
  sec2_dir <- file.path(assets_dir, "Section_2_BLUE_Extraction")
  dir.create(sec2_dir, showWarnings = FALSE, recursive = TRUE)

  cat("\n--- Section 2.1: Non-Destructive Delhi Block Mapping Log ---\n")
  cat(sprintf("  Primary Phenotype File : %s (MD5: %s)\n", data_file_path, data_md5_Full))
  cat(sprintf("  Block Mapping Source   : %s (MD5: %s)\n", gs_map_path, gs_md5_Full))
  cat(sprintf("  Mapping Rule           : Block_true = ((Original_Sr_No - 1L) %%%% %dL) %%/%% %dL + 1L\n",
              n_plots_per_rep_del, n_plots_per_blk_del))
  cat(sprintf("  Unmatched Delhi Rows   : %d\n",
              sum(is.na(master_raw$Block_true[master_raw$Location == "Delhi"]))))

  # ----------------------------------------------------------------------------
  # 2.2 Random Sample of Delhi Genotypes for Field Book Verification
  # ----------------------------------------------------------------------------
  del_s2_sub <- master_raw %>% filter(Sheet == "Delhi_S2")
  del_r1_sub <- del_s2_sub %>% filter(Rep_Char == "1") %>%
    select(Genotype_No, Genotype_ID = Genotype_Clean,
           Orig_Sr_No_R1 = Original_Sr_No, Old_Block_R1 = Block, Block_true_R1 = Block_true)
  del_r2_sub <- del_s2_sub %>% filter(Rep_Char == "2") %>%
    select(Genotype_ID = Genotype_Clean,
           Orig_Sr_No_R2 = Original_Sr_No, Old_Block_R2 = Block, Block_true_R2 = Block_true)
  del_verify_all <- inner_join(del_r1_sub, del_r2_sub, by = "Genotype_ID") %>%
    arrange(Genotype_No)

  set.seed(20261016)
  sample_size_check <- 15L
  del_sample_check <- del_verify_all[sort(sample(seq_len(nrow(del_verify_all)), sample_size_check)), ]
  cat("\n--- Section 2.2: Random Sample of Delhi Genotypes for Field Book Verification ---\n")
  print(del_sample_check, row.names = FALSE)
  write.csv(del_sample_check, file.path(sec2_dir, "Table_Delhi_15_Random_Genotypes_Block_Check.csv"), row.names = FALSE)

  # ----------------------------------------------------------------------------
  # 2.3 Recorded Field Design Counts, Ordering Audit & Dharwad Spatial Confounding
  # ----------------------------------------------------------------------------
  design_audit_rows <- list()
  for (s in env_metadata$Sheet) {
    sub_d <- master_raw %>% filter(Sheet == s)
    loc_s <- unique(sub_d$Location)
    e_nm  <- unique(sub_d$Environment)
    n_g   <- length(unique(sub_d$Genotype))
    n_r   <- length(unique(sub_d$Replication))
    n_b_old  <- length(unique(sub_d$Block))
    p_b_old  <- as.integer(median(as.vector(table(sub_d$Replication, sub_d$Block))))
    has_b_tr <- !all(is.na(sub_d$Block_true))
    n_b_true <- if (has_b_tr) length(unique( na.omit(sub_d$Block_true))) else NA_integer_
    p_b_true <- if (has_b_tr) as.integer(median(as.vector(table(sub_d$Replication, sub_d$Block_true)))) else NA_integer_

    # Verify (a) each genotype appears once per replicate and (b) no genotype appears twice in a block
    tab_gr <- table(sub_d$Genotype, sub_d$Replication)
    once_per_rep <- all(tab_gr == 1L)
    no_dup_old   <- all(table(sub_d$Genotype, interaction(sub_d$Replication, sub_d$Block, drop = TRUE)) <= 1L)
    no_dup_true  <- if (has_b_tr) {
      all(table(sub_d$Genotype, interaction(sub_d$Replication, sub_d$Block_true, drop = TRUE)) <= 1L)
    } else {
      NA
    }
    stopifnot(once_per_rep, no_dup_old)
    if (has_b_tr) stopifnot(no_dup_true)

    r1_d <- sub_d %>% filter(Rep_Char == "1") %>% arrange(Row_Order_In_Sheet)
    r2_d <- sub_d %>% filter(Rep_Char == "2") %>% arrange(Row_Order_In_Sheet)
    r1_eq_r2_master <- identical(r1_d$Genotype_Clean, r2_d$Genotype_Clean)

    if (loc_s == "Delhi") {
      r1_sr <- sub_d %>% filter(Rep_Char == "1") %>% arrange(Original_Sr_No)
      r2_sr <- sub_d %>% filter(Rep_Char == "2") %>% arrange(Original_Sr_No)
      r1_canon <- identical(r1_sr$Genotype_No, seq_len(nrow(r1_sr)))
      r2_canon <- identical(r2_sr$Genotype_No, seq_len(nrow(r2_sr)))
      cor_r1   <- cor(r1_sr$Genotype_No, as.numeric(as.character(r1_sr$Block_true)))
      cor_r2   <- cor(r2_sr$Genotype_No, as.numeric(as.character(r2_sr$Block_true)))
      layout_desc <- sprintf(
        "%d genotypes x %d replicates x %d blocks/rep x %d plots/block (to be confirmed against field book)",
        n_g, n_r, n_b_true, p_b_true
      )
      order_note <- sprintf(
        "R1 canonical genotype order = %s (r_GenoNo_Block = %.6f); R2 canonical = %s (r_GenoNo_Block = %.6f)",
        r1_canon, cor_r1, r2_canon, cor_r2
      )
    } else {
      r1_canon <- NA
      r2_canon <- NA
      cor_r1   <- NA_real_
      cor_r2   <- NA_real_
      dhw_sh_r1 <- if (s == "Dharwad_S1") "EarlyR1" else "NormalR1"
      dhw_sh_r2 <- if (s == "Dharwad_S1") "EarlyR2" else "NormalR2"
      dhw_or_note <- "Data_Dharwad.xlsx not found"
      if (file.exists(dhw_orig_path) && all(c(dhw_sh_r1, dhw_sh_r2) %in% excel_sheets(dhw_orig_path))) {
        dhw_o1 <- suppressMessages(read_excel(dhw_orig_path, sheet = dhw_sh_r1, skip = 1, .name_repair = "minimal"))
        dhw_o2 <- suppressMessages(read_excel(dhw_orig_path, sheet = dhw_sh_r2, skip = 1, .name_repair = "minimal"))
        dhw_r1 <- trimws(as.character(dhw_o1[[2]]))
        dhw_r2 <- trimws(as.character(dhw_o2[[2]]))
        dhw_r1 <- dhw_r1[!is.na(dhw_r1)]
        dhw_r2 <- dhw_r2[!is.na(dhw_r2)]
        dhw_r1_eq_r2 <- identical(dhw_r1, dhw_r2)
        dhw_eq_canon <- identical(dhw_r1, gs_r1_lookup$Genotype_Clean[order(gs_r1_lookup$Sr_R1)])
        dhw_o_note <- sprintf(
          "In Data_Dharwad.xlsx (%s/%s): R1==R2 is %s, equals canonical order is %s; physical sowing followed %s packet order in both replicates",
          dhw_sh_r1, dhw_sh_r2, dhw_r1_eq_r2, dhw_eq_canon, if (s == "Dharwad_S1") "GS R1" else "GS R2"
        )
        dhw_or_note <- dhw_o_note
      }
      layout_desc <- sprintf(
        "%d genotypes x %d replicates x %d recorded block x %d plots/rep (to be confirmed against field book)",
        n_g, n_r, n_b_old, p_b_old
      )
      order_note <- sprintf("Master_Data R1==R2 (alphabetical sort) = %s; %s", r1_eq_r2_master, dhw_or_note)
    }

    design_audit_rows[[s]] <- data.frame(
      Sheet = s,
      Environment = e_nm,
      N_Genotypes = n_g,
      N_Replicates = n_r,
      Old_Blocks_Per_Rep = n_b_old,
      Old_Plots_Per_Block = p_b_old,
      True_Blocks_Per_Rep = n_b_true,
      True_Plots_Per_Block = p_b_true,
      Check_a_Once_Per_Rep = once_per_rep,
      Check_b_No_Dup_In_Old_Block = no_dup_old,
      Check_b_No_Dup_In_True_Block = no_dup_true,
      R1_In_Canonical_Order = r1_canon,
      R2_In_Canonical_Order = r2_canon,
      Cor_GenotypeNo_BlockTrue_R1 = cor_r1,
      Cor_GenotypeNo_BlockTrue_R2 = cor_r2,
      Recorded_Layout_Description = layout_desc,
      Ordering_And_Packet_Report = order_note,
      stringsAsFactors = FALSE
    )
  }
  design_audit_table <- bind_rows(design_audit_rows)
  cat("\n--- Section 2.3: Field Design Counts & Ordering Audit (To Be Confirmed Against Field Book) ---\n")
  print(design_audit_table %>% select(Sheet, Environment, N_Genotypes, N_Replicates, True_Blocks_Per_Rep, True_Plots_Per_Block,
                                      R1_In_Canonical_Order, R2_In_Canonical_Order, Cor_GenotypeNo_BlockTrue_R1, Cor_GenotypeNo_BlockTrue_R2))
  write.csv(design_audit_table, file.path(sec2_dir, "Table_Field_Design_Counts_and_Ordering_Audit.csv"), row.names = FALSE)

  # Raw T1_Count & T2_Count granularity audit
  t12_rows <- list()
  for (s in env_metadata$Sheet) {
    sub_d <- master_raw %>% filter(Sheet == s)
    for (tc in c("T1_Count", "T2_Count")) {
      v <- sub_d[[tc]]
      scaled_5  <- all(abs(v * 5 - round(v * 5)) < 1e-6)
      scaled_10 <- all(abs(v * 10 - round(v * 10)) < 1e-6)
      t12_rows[[paste0(s, "_", tc)]] <- data.frame(
        Sheet = s,
        Environment = unique(sub_d$Environment),
        Trait = tc,
        N_Plots = length(v),
        N_Unique_Values = length(unique(v)),
        Min_Val = min(v, na.rm = TRUE),
        Max_Val = max(v, na.rm = TRUE),
        Mean_Val = mean(v, na.rm = TRUE),
        SD_Val = sd(v, na.rm = TRUE),
        Exact_Multiples_Of_OneFifth = scaled_5,
        Exact_Multiples_Of_OneTenth = scaled_10,
        stringsAsFactors = FALSE
      )
    }
  }
  t12_audit_table <- bind_rows(t12_rows)
  write.csv(t12_audit_table, file.path(sec2_dir, "Table_Raw_T1_T2_Count_Distributions.csv"), row.names = FALSE)

  # ----------------------------------------------------------------------------
  # 2.4 Derived Micro-Phenology Traits on master_raw & Emergence Anchor Check
  # ----------------------------------------------------------------------------
  env_gdd_wide <- stage_summary %>%
    select(Environment, Stage, Stage_GDD_Capped) %>%
    pivot_wider(names_from = Stage, values_from = Stage_GDD_Capped) %>%
    rename(
      GDD_E_T1  = `E-T1`,
      GDD_T1_T2 = `T1-T2`,
      GDD_T2_T3 = `T2-T3`,
      GDD_T3_T4 = `T3-T4`
    ) %>%
    mutate(Total_GDD_E_T4 = GDD_E_T1 + GDD_T1_T2 + GDD_T2_T3 + GDD_T3_T4)

  count_at_emergence     <- COUNT_AT_EMERGENCE_MAIN
  count_at_emergence_alt <- COUNT_AT_EMERGENCE_ALT

  stopifnot(all(master_raw$T4_Count > 0, na.rm = TRUE))

  master_raw <- master_raw %>%
    left_join(env_metadata %>% select(Sheet, Int_E_T1, Int_T1_T2, Int_T2_T3, Int_T3_T4, Total_E_T4), by = "Sheet") %>%
    left_join(env_gdd_wide, by = "Environment") %>%
    mutate(
      # Area Under Tillering Progress Curve (AUTPC) via Trapezoidal Integration
      AUTPC = (0.5 * (count_at_emergence + T1_Count) * Int_E_T1) +
              (0.5 * (T1_Count + T2_Count) * Int_T1_T2) +
              (0.5 * (T2_Count + T3_Count) * Int_T2_T3) +
              (0.5 * (T3_Count + T4_Count) * Int_T3_T4),
      AUTPC_norm = AUTPC / Total_E_T4,

      # Alternative sensitivity version
      AUTPC_alt = (0.5 * (count_at_emergence_alt + T1_Count) * Int_E_T1) +
                  (0.5 * (T1_Count + T2_Count) * Int_T1_T2) +
                  (0.5 * (T2_Count + T3_Count) * Int_T2_T3) +
                  (0.5 * (T3_Count + T4_Count) * Int_T3_T4),
      AUTPC_norm_alt = AUTPC_alt / Total_E_T4,

      # Thermal-time AUTPC
      AUTPC_GDD = (0.5 * (count_at_emergence + T1_Count) * GDD_E_T1) +
                  (0.5 * (T1_Count + T2_Count) * GDD_T1_T2) +
                  (0.5 * (T2_Count + T3_Count) * GDD_T2_T3) +
                  (0.5 * (T3_Count + T4_Count) * GDD_T3_T4),
      AUTPC_GDD_norm = AUTPC_GDD / Total_GDD_E_T4,

      # Tillering Dynamics Asymmetry
      Dur_Early = Int_E_T1 + Int_T1_T2,
      Dur_Late  = Int_T2_T3 + Int_T3_T4,
      AUTPC_Early = (0.5 * (count_at_emergence + T1_Count) * Int_E_T1) + (0.5 * (T1_Count + T2_Count) * Int_T1_T2),
      AUTPC_Late  = (0.5 * (T2_Count + T3_Count) * Int_T2_T3) + (0.5 * (T3_Count + T4_Count) * Int_T3_T4),
      Mean_Tiller_Early = AUTPC_Early / Dur_Early,
      Mean_Tiller_Late  = AUTPC_Late  / Dur_Late,
      AUTPC_Asymmetry = (Mean_Tiller_Early - Mean_Tiller_Late) / (Mean_Tiller_Early + Mean_Tiller_Late),

      # Tillering Rates (tillers per plant / day)
      Rate_E_T1  = (T1_Count - count_at_emergence) / Int_E_T1,
      Rate_T1_T2 = (T2_Count - T1_Count) / Int_T1_T2,
      Rate_T2_T3 = (T3_Count - T2_Count) / Int_T2_T3,
      Rate_T3_T4 = (T4_Count - T3_Count) / Int_T3_T4,

      # Exact Tiller Conversion Efficiency (TCE = PT / T4; zero-guarded without arbitrary constants)
      TCE = ifelse(T4_Count == 0, NA_real_, PT_Count / T4_Count)
    )

  autpc_anchor_sensitivity <- master_raw %>%
    group_by(Environment) %>%
    summarise(
      Mean_AUTPC_norm_Count1 = mean(AUTPC_norm, na.rm = TRUE),
      Mean_AUTPC_norm_Count0 = mean(AUTPC_norm_alt, na.rm = TRUE),
      Delta_Mean             = Mean_AUTPC_norm_Count1 - Mean_AUTPC_norm_Count0,
      Percent_Change         = (Delta_Mean / Mean_AUTPC_norm_Count0) * 100,
      Theoretical_Shift      = (0.5 * first(Int_E_T1)) / first(Total_E_T4),
      .groups = "drop"
    )
  cat("\n--- Anchor Sensitivity: Environment Mean AUTPC_norm (Primary vs Alt Anchor) ---\n")
  print(autpc_anchor_sensitivity)
  # Note: Anchor sensitivity table is exported as Table_S7 in Section 5 (Section_2.3_Phenology_Indices)

  # ----------------------------------------------------------------------------
  # 2.5 Dharwad Spatial Position Confounding Diagnostics (GS R1 vs GS R2 Order)
  # ----------------------------------------------------------------------------
  dhw_spatial_rows <- list()
  for (s_dhw in c("Dharwad_S1", "Dharwad_S2")) {
    sub_dhw <- master_raw %>% filter(Sheet == s_dhw)
    g_means <- sub_dhw %>%
      group_by(Genotype_Clean, Sr_R1, Block_Pos_GS_R1, Sr_R2, Block_Pos_GS_R2) %>%
      summarise(
        across(c(GrainYld, PT_Count, T3_Count, T4_Count, TKW, AUTPC_norm), ~ mean(.x, na.rm = TRUE)),
        .groups = "drop"
      )
    for (tr_sp in c("GrainYld", "PT_Count", "T3_Count", "T4_Count", "TKW", "AUTPC_norm")) {
      ord_r1 <- g_means %>% arrange(Sr_R1)
      v_r1   <- ord_r1[[tr_sp]]
      lag1_r1 <- cor(v_r1[-length(v_r1)], v_r1[-1])
      m_blk_r1 <- lm(v_r1 ~ factor(ord_r1$Block_Pos_GS_R1))
      an_r1    <- anova(m_blk_r1)

      ord_r2 <- g_means %>% arrange(Sr_R2)
      v_r2   <- ord_r2[[tr_sp]]
      lag1_r2 <- cor(v_r2[-length(v_r2)], v_r2[-1])
      m_blk_r2 <- lm(v_r2 ~ factor(ord_r2$Block_Pos_GS_R2))
      an_r2    <- anova(m_blk_r2)

      dhw_spatial_rows[[paste0(s_dhw, "_", tr_sp)]] <- data.frame(
        Sheet = s_dhw,
        Environment = unique(sub_dhw$Environment),
        Actual_Sown_Packet_Order = if (s_dhw == "Dharwad_S1") "GS_R1_Order" else "GS_R2_Order",
        Trait = tr_sp,
        Lag1_Autocorr_GS_R1_Order = lag1_r1,
        PosBlock_R2_GS_R1_Order   = summary(m_blk_r1)$r.squared,
        PosBlock_p_GS_R1_Order    = an_r1$`Pr(>F)`[1],
        Lag1_Autocorr_GS_R2_Order = lag1_r2,
        PosBlock_R2_GS_R2_Order   = summary(m_blk_r2)$r.squared,
        PosBlock_p_GS_R2_Order    = an_r2$`Pr(>F)`[1],
        stringsAsFactors = FALSE
      )
    }
  }
  dhw_spatial_table <- bind_rows(dhw_spatial_rows)
  write.csv(dhw_spatial_table, file.path(sec2_dir, "Table_Dharwad_Spatial_Confounding_Diagnostics.csv"), row.names = FALSE)

  # ----------------------------------------------------------------------------
  # 2.6 Permutation Check on Replicate 2 Block_true Labels in Delhi
  # ----------------------------------------------------------------------------
  n_perm_r2 <- if (exists("N_PERM_R2_OVERRIDE")) as.integer(N_PERM_R2_OVERRIDE) else 500L
  perm_traits <- c("GrainYld", "T3_Count", "AUTPC_norm")
  perm_rows <- list()
  cat(sprintf("\n--- Section 2.6: Running %d-Permutation Check on Delhi Replicate 2 Block Labels ---\n", n_perm_r2))
  for (s_del in c("Delhi_S2", "Delhi_S1")) {
    sub_del <- master_raw %>% filter(Sheet == s_del) %>% arrange(Rep_Char, Original_Sr_No)
    r2_idx  <- which(sub_del$Rep_Char == "2")
    true_r2_blocks <- as.character(sub_del$Block_true[r2_idx])

    for (tr_p in perm_traits) {
      m_null <- lm(as.formula(paste0(tr_p, " ~ Genotype + Replication")), data = sub_del)
      m_obs  <- suppressMessages(lmer(
        as.formula(paste0(tr_p, " ~ Genotype + Replication + (1 | Replication:Block_true)")),
        data = sub_del, REML = TRUE, control = lmerControl(optimizer = "bobyqa")
      ))
      reml_null <- as.numeric(REMLcrit(m_obs))
      # Exact REML deviance improvement over fixed-effects model
      ll_null_reml <- as.numeric(logLik(lm(as.formula(paste0(tr_p, " ~ Genotype + Replication")), data = sub_del), REML = TRUE))
      ll_obs_reml  <- as.numeric(logLik(m_obs, REML = TRUE))
      obs_lrt <- max(0, 2 * (ll_obs_reml - ll_null_reml))
      vc_obs  <- as.data.frame(VarCorr(m_obs))
      obs_vb  <- vc_obs$vcov[vc_obs$grp == "Replication:Block_true"]

      set.seed(20261016)
      perm_lrt <- numeric(n_perm_r2)
      perm_vb  <- numeric(n_perm_r2)
      df_perm  <- sub_del
      for (b_i in seq_len(n_perm_r2)) {
        df_perm$Block_perm <- as.character(df_perm$Block_true)
        df_perm$Block_perm[r2_idx] <- sample(true_r2_blocks, length(true_r2_blocks), replace = FALSE)
        df_perm$Block_perm <- as.factor(df_perm$Block_perm)
        m_p <- suppressMessages(lmer(
          as.formula(paste0(tr_p, " ~ Genotype + Replication + (1 | Replication:Block_perm)")),
          data = df_perm, REML = TRUE, control = lmerControl(optimizer = "bobyqa")
        ))
        ll_p <- as.numeric(logLik(m_p, REML = TRUE))
        perm_lrt[b_i] <- max(0, 2 * (ll_p - ll_null_reml))
        vc_p <- as.data.frame(VarCorr(m_p))
        perm_vb[b_i] <- vc_p$vcov[vc_p$grp == "Replication:Block_perm"]
      }

      pct_lrt <- mean(perm_lrt <= obs_lrt) * 100
      pct_vb  <- mean(perm_vb  <= obs_vb)  * 100
      emp_p   <- (sum(perm_lrt >= obs_lrt) + 1) / (n_perm_r2 + 1)

      perm_rows[[paste0(s_del, "_", tr_p)]] <- data.frame(
        Sheet = s_del,
        Environment = unique(sub_del$Environment),
        Trait = tr_p,
        N_Permutations = n_perm_r2,
        Obs_REML_Deviance_Improvement = obs_lrt,
        Perm_Mean_REML_Improvement    = mean(perm_lrt),
        Perm_Median_REML_Improvement  = median(perm_lrt),
        Perm_Q95_REML_Improvement     = as.numeric(quantile(perm_lrt, 0.95)),
        Perm_Max_REML_Improvement     = max(perm_lrt),
        Obs_Percentile_REML           = pct_lrt,
        Empirical_p_REML              = emp_p,
        Obs_Var_RepBlock              = obs_vb,
        Perm_Mean_Var_RepBlock        = mean(perm_vb),
        Perm_Q95_Var_RepBlock         = as.numeric(quantile(perm_vb, 0.95)),
        Obs_Percentile_Var_RepBlock   = pct_vb,
        stringsAsFactors = FALSE
      )
    }
  }
  perm_summary_table <- bind_rows(perm_rows)
  print(perm_summary_table %>% select(Sheet, Trait, Obs_REML_Deviance_Improvement, Perm_Median_REML_Improvement,
                                      Perm_Q95_REML_Improvement, Obs_Percentile_REML, Empirical_p_REML))
  write.csv(perm_summary_table, file.path(sec2_dir, "Table_Delhi_R2_Block_Permutation_Check.csv"), row.names = FALSE)

  # ----------------------------------------------------------------------------
  # 2.7 Model A vs Model B Comparison
  # ----------------------------------------------------------------------------
  # Compare Delhi Model A vs Model B across all 18 Delhi trait x trial combinations
  model_a_cache_file <- file.path(cache_dir, paste0("ModelA_Audit_", cache_key, ".rds"))
  delhi_model_A_cache <- if (!isTRUE(CLEAR_BLUE_CACHE) && file.exists(model_a_cache_file)) {
    cat(sprintf("  Loading cached Delhi Model A fits (%s)...\n", basename(model_a_cache_file)))
    readRDS(model_a_cache_file)
  } else {
    list()
  }

  ab_summary_rows <- list()
  ab_geno_rows    <- list()
  for (s_del in c("Delhi_S2", "Delhi_S1")) {
    sub_del <- master_raw %>% filter(Sheet == s_del)
    e_nm    <- unique(sub_del$Environment)
    for (tr in traits_to_run) {
      key_st <- paste0(s_del, "_", tr)
      res_A <- if (!is.null(delhi_model_A_cache[[key_st]])) {
        delhi_model_A_cache[[key_st]]
      } else {
        fit_a_tmp <- fit_delhi_model_A(sub_del, tr)
        delhi_model_A_cache[[key_st]] <- fit_a_tmp
        fit_a_tmp
      }
      res_B <- fit_rcbd_lm(sub_del, tr)
      df_A  <- res_A$emm_df %>% arrange(as.character(Genotype))
      df_B  <- res_B$emm_df %>% arrange(as.character(Genotype))
      stopifnot(identical(as.character(df_A$Genotype), as.character(df_B$Genotype)))

      r_p   <- cor(df_A$emmean, df_B$emmean, method = "pearson")
      r_s   <- cor(df_A$emmean, df_B$emmean, method = "spearman")
      rk_A  <- rank(-df_A$emmean, ties.method = "average")
      rk_B  <- rank(-df_B$emmean, ties.method = "average")
      max_drk <- max(abs(rk_A - rk_B))
      max_dab <- max(abs(df_A$emmean - df_B$emmean))

      ab_summary_rows[[paste0(s_del, "_", tr)]] <- data.frame(
        Sheet = s_del,
        Environment = e_nm,
        Trait = tr,
        Model_A_Singular = res_A$is_singular,
        Var_RepBlock_A = res_A$var_block,
        Var_Residual_A = res_A$var_resid,
        Block_ICC_A = res_A$var_block / (res_A$var_block + res_A$var_resid),
        Var_Residual_B = res_B$var_resid,
        Mean_SE_Model_A = mean(df_A$SE),
        Min_SE_Model_A  = min(df_A$SE),
        Max_SE_Model_A  = max(df_A$SE),
        Mean_SE_Model_B = mean(df_B$SE),
        Pearson_r_A_vs_B = r_p,
        Spearman_rho_A_vs_B = r_s,
        Max_Abs_Rank_Change = max_drk,
        Max_Abs_BLUE_Diff = max_dab,
        stringsAsFactors = FALSE
      )

      ab_geno_rows[[paste0(s_del, "_", tr)]] <- data.frame(
        Sheet = s_del,
        Environment = e_nm,
        Trait = tr,
        Genotype = as.character(df_A$Genotype),
        BLUE_Model_A = df_A$emmean,
        SE_Model_A   = df_A$SE,
        df_Model_A   = df_A$df,
        BLUE_Model_B = df_B$emmean,
        SE_Model_B   = df_B$SE,
        df_Model_B   = df_B$df,
        Delta_BLUE_A_minus_B = df_A$emmean - df_B$emmean,
        Delta_Rank_A_minus_B = rk_A - rk_B,
        stringsAsFactors = FALSE
      )
    }
  }
  saveRDS(delhi_model_A_cache, model_a_cache_file)
  ab_summary_table <- bind_rows(ab_summary_rows)
  ab_geno_table    <- bind_rows(ab_geno_rows)
  cat("\n--- Section 2.7: Delhi Model A vs Model B Comparison Summary ---\n")
  print(ab_summary_table %>% select(Sheet, Trait, Model_A_Singular, Block_ICC_A, Mean_SE_Model_A, Mean_SE_Model_B,
                                    Pearson_r_A_vs_B, Spearman_rho_A_vs_B, Max_Abs_Rank_Change))
  write.csv(ab_summary_table, file.path(sec2_dir, "Table_Delhi_Model_A_vs_B_Summary.csv"), row.names = FALSE)
  write.csv(ab_geno_table,    file.path(sec2_dir, "Table_Delhi_Model_A_vs_B_Genotype_BLUEs_and_SEs.csv"), row.names = FALSE)

  cat(sprintf("\n--- Section 2.7: Cryptographic BLUE Cache Key: %s ---\n", cache_key))

  # ----------------------------------------------------------------------------
  # 2.8 Extract Primary BLUEs (Model A for Delhi, Fixed RCBD for Dharwad) & Compare Against Prior Cache
  # ----------------------------------------------------------------------------
  # Snapshot any legacy cached BLUEs before purging
  prior_cache_list <- list()
  for (s_code in env_metadata$Sheet) {
    old_candidates <- list.files(cache_dir, pattern = paste0("^BLUEs_", s_code, "_.*\\.rds$"), full.names = TRUE)
    old_candidates <- setdiff(old_candidates, file.path(cache_dir, paste0("BLUEs_", s_code, "_", cache_key, ".rds")))
    if (length(old_candidates) > 0) {
      f_info <- file.info(old_candidates)
      latest_old <- rownames(f_info)[which.max(f_info$mtime)]
      prior_cache_list[[s_code]] <- readRDS(latest_old)
    }
  }

  if (isTRUE(CLEAR_BLUE_CACHE)) {
    all_old_blue_rds <- list.files(cache_dir, pattern = "^BLUEs_.*\\.rds$", full.names = TRUE)
    if (length(all_old_blue_rds) > 0) {
      cat(sprintf("  Purging %d existing BLUE cache file(s) (CLEAR_BLUE_CACHE = TRUE)...\n", length(all_old_blue_rds)))
      unlink(all_old_blue_rds)
    }
  }

  all_blues      <- list()
  audit_36_rows  <- list()
  geno_se_7200   <- list()

  for (i in seq_len(nrow(env_metadata))) {
    s_code <- env_metadata$Sheet[i]
    e_name <- env_metadata$Env_Name[i]
    loc_i  <- env_metadata$Location[i]
    c_file <- file.path(cache_dir, paste0("BLUEs_", s_code, "_", cache_key, ".rds"))

    sub_raw <- master_raw %>% filter(Sheet == s_code)
    b_df    <- data.frame(Genotype = sort(unique(as.character(sub_raw$Genotype))), stringsAsFactors = FALSE)

    des_row <- design_audit_table %>% filter(Sheet == s_code)
    des_lbl <- des_row$Recorded_Layout_Description[1]

    cat(sprintf("  Fitting primary Stage-1 models for %s (%s) [Location: %s]...\n", s_code, e_name, loc_i))
    for (tr in traits_to_run) {
      key_st <- paste0(s_code, "_", tr)
      res_fit <- if (loc_i == "Delhi") {
        if (!is.null(delhi_model_A_cache[[key_st]])) delhi_model_A_cache[[key_st]] else fit_delhi_model_A(sub_raw, tr)
      } else {
        fit_rcbd_lm(sub_raw, tr)
      }
      if (res_fit$conv_msg != "OK") {
        stop(sprintf("CRITICAL CONVERGENCE FAILURE in %s - %s: %s", s_code, tr, res_fit$conv_msg))
      }

      emm_sorted <- res_fit$emm_df[match(b_df$Genotype, as.character(res_fit$emm_df$Genotype)), ]
      b_df[[tr]] <- emm_sorted$emmean

      se_note <- if (loc_i == "Delhi") {
        "Randomized Incomplete-Block KR SE"
      } else {
        "Pseudo-Replicate Strip SE (Confounded with Field Position)"
      }

      audit_36_rows[[paste0(s_code, "_", tr)]] <- data.frame(
        Sheet = s_code,
        Environment = e_name,
        Location = loc_i,
        Trait = tr,
        Design_Counts_Label = des_lbl,
        Formula = res_fit$formula_str,
        Engine = res_fit$engine_str,
        DF_Method = res_fit$df_method,
        Singular = res_fit$is_singular,
        Convergence_Message = res_fit$conv_msg,
        Var_RepBlock = res_fit$var_block,
        Var_Residual = res_fit$var_resid,
        Block_ICC = if (!is.na(res_fit$var_block)) res_fit$var_block / (res_fit$var_block + res_fit$var_resid) else NA_real_,
        Mean_BLUE_SE = mean(emm_sorted$SE),
        Min_BLUE_SE  = min(emm_sorted$SE),
        Max_BLUE_SE  = max(emm_sorted$SE),
        Mean_DF      = mean(emm_sorted$df),
        SE_Interpretation_Note = se_note,
        Fallback = "None",
        stringsAsFactors = FALSE
      )

      geno_se_7200[[paste0(s_code, "_", tr)]] <- data.frame(
        Sheet = s_code,
        Environment = e_name,
        Trait = tr,
        Genotype = b_df$Genotype,
        BLUE = emm_sorted$emmean,
        SE = emm_sorted$SE,
        df = emm_sorted$df,
        lower.CL = emm_sorted$lower.CL,
        upper.CL = emm_sorted$upper.CL,
        stringsAsFactors = FALSE
      )
    }

    b_df$Environment <- e_name
    saveRDS(b_df, c_file)
    all_blues[[s_code]] <- b_df
  }

  audit_36_table <- bind_rows(audit_36_rows)
  geno_se_table  <- bind_rows(geno_se_7200)
  stopifnot(nrow(audit_36_table) == length(traits_to_run) * nrow(env_metadata))
  stopifnot(nrow(geno_se_table)  == length(unique(master_raw$Genotype)) * length(traits_to_run) * nrow(env_metadata))

  # Report singular fit counts per trial
  sing_summary <- audit_36_table %>%
    group_by(Sheet, Environment) %>%
    summarise(
      N_Models = n(),
      N_Singular = sum(Singular, na.rm = TRUE),
      N_Convergence_Failures = sum(Convergence_Message != "OK"),
      .groups = "drop"
    )
  cat("\n--- Section 2.8: Singular Fit & Convergence Summary per Trial ---\n")
  print(sing_summary)

  cat("\n--- Section 2.8: 36-Row Stage-1 Model Audit Table ---\n")
  print(audit_36_table %>% select(Sheet, Trait, Engine, DF_Method, Singular, Convergence_Message,
                                  Var_RepBlock, Var_Residual, Mean_BLUE_SE, Fallback))
  write.csv(audit_36_table, file.path(sec2_dir, "Table_BLUE_Model_Specifications_and_Diagnostics.csv"), row.names = FALSE)
  write.csv(geno_se_table,  file.path(sec2_dir, "Table_Genotype_BLUEs_and_SEs.csv"), row.names = FALSE)

  # Compare Old vs New BLUEs if prior cache existed
  if (length(prior_cache_list) > 0) {
    old_vs_new_rows <- list()
    for (s_code in names(prior_cache_list)) {
      old_b <- prior_cache_list[[s_code]]
      new_b <- all_blues[[s_code]]
      old_b <- old_b[match(new_b$Genotype, as.character(old_b$Genotype)), ]
      for (tr in traits_to_run) {
        if (tr %in% names(old_b)) {
          r_p <- cor(old_b[[tr]], new_b[[tr]], method = "pearson")
          r_s <- cor(old_b[[tr]], new_b[[tr]], method = "spearman")
          rk_o <- rank(-old_b[[tr]], ties.method = "average")
          rk_n <- rank(-new_b[[tr]], ties.method = "average")
          old_vs_new_rows[[paste0(s_code, "_", tr)]] <- data.frame(
            Sheet = s_code,
            Environment = unique(new_b$Environment),
            Trait = tr,
            Pearson_r_Old_vs_New = r_p,
            Spearman_rho_Old_vs_New = r_s,
            Max_Abs_Rank_Shift = max(abs(rk_n - rk_o)),
            Max_Abs_BLUE_Diff  = max(abs(new_b[[tr]] - old_b[[tr]])),
            stringsAsFactors = FALSE
          )
        }
      }
    }
    old_vs_new_table <- bind_rows(old_vs_new_rows)
    write.csv(old_vs_new_table, file.path(sec2_dir, "Table_Old_vs_New_BLUE_Comparison.csv"), row.names = FALSE)
  }

} else {
  cat(sprintf("  -> Section 2 [SKIPPED - Loading Stage-1 BLUEs matching Cryptographic Key: %s]\n", cache_key))
  master_raw <- master_raw %>%
    left_join(env_metadata %>% select(Sheet, Int_E_T1, Int_T1_T2, Int_T2_T3, Int_T3_T4, Total_E_T4), by = "Sheet") %>%
    mutate(
      AUTPC = (0.5 * (COUNT_AT_EMERGENCE_MAIN + T1_Count) * Int_E_T1) +
              (0.5 * (T1_Count + T2_Count) * Int_T1_T2) +
              (0.5 * (T2_Count + T3_Count) * Int_T2_T3) +
              (0.5 * (T3_Count + T4_Count) * Int_T3_T4),
      AUTPC_norm = AUTPC / Total_E_T4,
      TCE = ifelse(T4_Count == 0, NA_real_, PT_Count / T4_Count)
    )
  all_blues <- list()
  for (s_code in env_metadata$Sheet) {
    c_file <- file.path(cache_dir, paste0("BLUEs_", s_code, "_", cache_key, ".rds"))
    if (!file.exists(c_file)) {
      stop(sprintf("CRITICAL CACHE KEY MISMATCH: Expected cached BLUE file '%s' not found. Set RUN_SECTION_2_BLUES <- TRUE to recompute BLUEs for current data/model specification.", basename(c_file)))
    }
    b_loaded <- readRDS(c_file)
    b_loaded$Environment <- env_metadata$Env_Name[env_metadata$Sheet == s_code]
    all_blues[[s_code]] <- b_loaded
  }
}

# ----------------------------------------------------------------------------
# 2.9 Construct Master Stage-1 BLUE Data Frame (blues_df) & All Derived Traits
# ----------------------------------------------------------------------------
blues_df <- bind_rows(all_blues) %>%
  left_join(
    env_metadata %>% select(Env_Name, Sheet, Location, Int_E_T1, Int_T1_T2, Int_T2_T3, Int_T3_T4, Total_E_T4),
    by = c("Environment" = "Env_Name")
  ) %>%
  left_join(gs_r1_lookup %>% select(Genotype_Clean, Genotype_No, Block_Pos_GS_R1), by = c("Genotype" = "Genotype_Clean")) %>%
  left_join(gs_r2_lookup %>% select(Genotype_Clean, Block_Pos_GS_R2), by = c("Genotype" = "Genotype_Clean")) %>%
  mutate(
    Field_Pos_Block_Diagnostic = case_when(
      Sheet %in% c("Delhi_S2", "Delhi_S1") ~ NA_integer_,
      Sheet == "Dharwad_S1"                ~ as.integer(Block_Pos_GS_R1),
      Sheet == "Dharwad_S2"                ~ as.integer(Block_Pos_GS_R2),
      TRUE                                 ~ NA_integer_
    ),
    Dur_Early = Int_E_T1 + Int_T1_T2,
    Dur_Late  = Int_T2_T3 + Int_T3_T4,
    AUTPC_Early = (0.5 * (COUNT_AT_EMERGENCE_MAIN + T1_Count) * Int_E_T1) +
                  (0.5 * (T1_Count + T2_Count) * Int_T1_T2),
    AUTPC_Late  = (0.5 * (T2_Count + T3_Count) * Int_T2_T3) +
                  (0.5 * (T3_Count + T4_Count) * Int_T3_T4),
    Mean_Tiller_Early = AUTPC_Early / Dur_Early,
    Mean_Tiller_Late  = AUTPC_Late  / Dur_Late,
    AUTPC_Asymmetry   = (Mean_Tiller_Early - Mean_Tiller_Late) / (Mean_Tiller_Early + Mean_Tiller_Late),
    Rate_E_T1  = (T1_Count - COUNT_AT_EMERGENCE_MAIN) / Int_E_T1,
    Rate_T1_T2 = (T2_Count - T1_Count) / Int_T1_T2,
    Rate_T2_T3 = (T3_Count - T2_Count) / Int_T2_T3,
    Rate_T3_T4 = (T4_Count - T3_Count) / Int_T3_T4,
    TCE = ifelse(T4_Count == 0, NA_real_, PT_Count / T4_Count),
    Env_Factor = factor(Environment, levels = env_metadata$Env_Name)
  )

# Exact additivity holds against AUTPC recomputed from the count BLUEs. The separately modelled AUTPC BLUE
# differs because each trait has its own fitted variance components; that gap is reported as a diagnostic.
blues_df$AUTPC_From_Count_BLUEs <- blues_df$AUTPC_Early + blues_df$AUTPC_Late
recomputed_total <- (0.5 * (COUNT_AT_EMERGENCE_MAIN + blues_df$T1_Count) * blues_df$Int_E_T1) +
                    (0.5 * (blues_df$T1_Count + blues_df$T2_Count) * blues_df$Int_T1_T2) +
                    (0.5 * (blues_df$T2_Count + blues_df$T3_Count) * blues_df$Int_T2_T3) +
                    (0.5 * (blues_df$T3_Count + blues_df$T4_Count) * blues_df$Int_T3_T4)
max_autpc_decomp_err <- max(abs(blues_df$AUTPC_From_Count_BLUEs - recomputed_total), na.rm = TRUE)
max_autpc_model_gap  <- max(abs(blues_df$AUTPC_From_Count_BLUEs - blues_df$AUTPC), na.rm = TRUE)
cat(sprintf("  AUTPC additivity (Early + Late vs recomputed from count BLUEs): max abs error = %.3e\n", max_autpc_decomp_err))
cat(sprintf("  Diagnostic gap, recomputed vs separately modelled AUTPC BLUE: max abs = %.3f\n", max_autpc_model_gap))
stopifnot(max_autpc_decomp_err < 1e-10, sum(is.na(blues_df$TCE)) == 0L)


# ==============================================================================
# SECTION 3: BASELINE PHENOTYPIC DISTRIBUTIONS & GLOBAL PARAMETER ANOVA
# ==============================================================================
cat(">>> Executing Section 3: Phenotypic Distributions & Global ANOVA...\n")

if (isTRUE(RUN_SECTION_3_DISTRIBUTIONS)) {
  sec21_dir <- file.path(assets_dir, "Section_2.1_Phenotypic_Distributions")
  dir.create(sec21_dir, showWarnings = FALSE, recursive = TRUE)
  arch_sec21_dir <- file.path(project_root, "Archive_Rejected", "Section_2.1_Legacy_Distributions")
  dir.create(arch_sec21_dir, showWarnings = FALSE, recursive = TRUE)

  # Archive legacy/redundant figures and split CSVs if present in Section_2.1_Phenotypic_Distributions
  legacy_sec21_files <- c(
    "Compound_Global_Rainclouds.png",
    "Global_Parameter_ANOVA_Professional.png",
    "TKW_ANOVA_and_Paired_Profile.png",
    "Diagnostic_Residuals_QQ_Repeated_Measures.png",
    "Table_Multi_Environment_Trait_Summary_Statistics.csv",
    "Table_Repeated_Measures_ANOVA_Summary.csv",
    "Table_Residual_Variance_Heterogeneity.csv",
    "Table_Repeated_Measures_Pairwise_Contrasts.csv"
  )
  for (lf in legacy_sec21_files) {
    src_lf <- file.path(sec21_dir, lf)
    if (file.exists(src_lf)) {
      file.copy(src_lf, file.path(arch_sec21_dir, lf), overwrite = TRUE)
      unlink(src_lf)
    }
  }

  # 3.1 Model Data Checks & Design Balance (Ordered Chronologically Along Developmental Cascade)
  traits_sec3 <- c("T4_Count", "PT_Count", "TCE", "TKW", "GrainYld")
  trait_titles_sec3 <- c(
    "T4_Count" = "Final Vegetative Tillers at T4",
    "PT_Count" = "Productive Tillers at Harvest (PT)",
    "TCE"      = "Tiller Conversion Efficiency (TCE)",
    "TKW"      = "Thousand Kernel Weight (TKW)",
    "GrainYld" = "Plot Grain Yield"
  )
  trait_short_sec3 <- c(
    "T4_Count" = "Final Veg. Tillers (T4)",
    "PT_Count" = "Productive Tillers (PT)",
    "TCE"      = "Tiller Conv. Eff. (TCE)",
    "TKW"      = "1000-Kernel Weight (TKW)",
    "GrainYld" = "Grain Yield"
  )
  trait_ylab_sec3 <- c(
    "T4_Count" = "Tillers per plant",
    "PT_Count" = "Tillers per plant",
    "TCE"      = "Conversion Efficiency (%)",
    "TKW"      = "Thousand Kernel Weight (g)",
    "GrainYld" = "Grain Yield (kg/plot)"
  )
  trait_units_sec3 <- c(
    "T4_Count" = "Tillers/plant",
    "PT_Count" = "Tillers/plant",
    "TCE"      = "%",
    "TKW"      = "g",
    "GrainYld" = "kg/plot"
  )

  n_geno_sec3 <- dplyr::n_distinct(blues_df$Genotype)
  n_env_sec3  <- dplyr::n_distinct(blues_df$Env_Factor)
  na_counts_sec3 <- sapply(blues_df[traits_sec3], function(x) sum(is.na(x)))
  bal_tab_sec3 <- table(blues_df$Genotype, blues_df$Env_Factor)
  is_balanced_sec3 <- all(bal_tab_sec3 == 1)
  den_df_omnibus <- (n_env_sec3 - 1) * (n_geno_sec3 - 1)

  cat("\n--- Section 3: Model Data Check ---\n")
  cat(sprintf("Genotypes: %d | Environments: %d | Total Observations: %d\n", n_geno_sec3, n_env_sec3, nrow(blues_df)))
  cat("NA counts per trait:\n")
  print(na_counts_sec3)
  cat(sprintf("Design balanced (all cells == 1): %s\n", is_balanced_sec3))
  cat(sprintf("Denominator degrees of freedom for omnibus F: %d\n", den_df_omnibus))
  stopifnot(is_balanced_sec3, all(na_counts_sec3 == 0), den_df_omnibus == 597)

  # TCE Range Check and Outlier Flagging (> 1)
  tce_range <- range(blues_df$TCE, na.rm = TRUE)
  tce_flagged <- blues_df %>% filter(TCE > 1)
  cat(sprintf("\nTCE range: [%.4f, %.4f]\n", tce_range[1], tce_range[2]))
  cat(sprintf("Count of TCE > 1: %d\n", nrow(tce_flagged)))
  if (nrow(tce_flagged) > 0) {
    cat("Genotypes with TCE > 1:\n")
    print(tce_flagged %>% select(Genotype, Environment, PT_Count, T4_Count, TCE))
  }

  # Empirical Logit Sensitivity for TCE
  blues_df$logit_TCE_scaled <- qlogis(blues_df$TCE / (max(blues_df$TCE) * 1.001))
  m_tce_logit <- lmer(logit_TCE_scaled ~ Env_Factor + (1 | Genotype), data = blues_df)
  an_tce_logit <- car::Anova(m_tce_logit, type = "III", test.statistic = "F")
  cat("\n--- TCE Sensitivity: logit(TCE) Repeated-Measures LMM Omnibus F ---\n")
  cat(sprintf("logit(TCE) F(%d, %.1f) = %.2f, p = %.3e\n",
              an_tce_logit$Df[2], an_tce_logit$Df.res[2], an_tce_logit$F[2], an_tce_logit$`Pr(>F)`[2]))
  emm_tce_logit <- emmeans(m_tce_logit, ~ Env_Factor)
  prs_tce_logit <- as.data.frame(pairs(emm_tce_logit, adjust = "tukey"))
  cat("logit(TCE) Tukey Pairwise Contrasts:\n")
  print(prs_tce_logit)

  # 3.2 Repeated-Measures Modeling, Heteroskedasticity Audit, GLS, & Wilcoxon
  prof_theme <- theme_minimal(base_size = 13.5) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 18, hjust = 1, face = "bold", size = 11.4, color = "grey10"),
      axis.text.y = element_text(size = 11.0, face = "bold", color = "grey15"),
      axis.title.y = element_text(face = "bold", size = 12.2, color = "black"),
      plot.title = element_text(face = "bold", size = 14.0, hjust = 0.5, color = "#1B2631"),
      plot.subtitle = element_text(size = 10.2, hjust = 0.5, color = "grey25", face = "italic", margin = margin(b = 5)),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.35),
      plot.margin = margin(t = 5, r = 8, b = 4, l = 6)
    )

  env_levels <- levels(blues_df$Env_Factor)
  all_pairs <- combn(env_levels, 2, simplify = FALSE)

  blues_sorted <- blues_df %>% arrange(Genotype, Env_Factor)
  blues_sorted$time_idx <- as.numeric(blues_sorted$Env_Factor)

  hetero_list <- list()
  anova_summary_list <- list()
  contrast_list <- list()
  spaghetti_plots <- list()
  cld_store_list <- list()
  rel_shift_rows <- list()

  # Function to compute Compact Letter Display (CLD) from Tukey p-values
  compute_cld_letters <- function(env_levels, emm_means, contrasts_df, p_col = "p_tukey_lmer", alpha = 0.05) {
    k <- length(env_levels)
    non_diff <- matrix(TRUE, nrow = k, ncol = k, dimnames = list(env_levels, env_levels))
    for (i in seq_len(nrow(contrasts_df))) {
      c_str <- contrasts_df$contrast[i]
      pts <- strsplit(c_str, " - ")[[1]]
      g1 <- trimws(pts[1])
      g2 <- trimws(pts[2])
      p_val <- contrasts_df[[p_col]][i]
      if (!is.na(p_val) && p_val < alpha) {
        non_diff[g1, g2] <- FALSE
        non_diff[g2, g1] <- FALSE
      }
    }

    subsets <- list()
    for (len in k:1) {
      combs <- combn(env_levels, len, simplify = FALSE)
      for (comb in combs) {
        is_clique <- TRUE
        if (length(comb) > 1) {
          p_pairs <- combn(comb, 2, simplify = FALSE)
          for (pr in p_pairs) {
            if (!non_diff[pr[1], pr[2]]) {
              is_clique <- FALSE
              break
            }
          }
        }
        if (is_clique) {
          is_sub <- FALSE
          for (cq in subsets) {
            if (all(comb %in% cq)) {
              is_sub <- TRUE
              break
            }
          }
          if (!is_sub) {
            subsets[[length(subsets) + 1]] <- comb
          }
        }
      }
    }

    clique_scores <- sapply(subsets, function(cq) max(emm_means[cq]))
    subsets <- subsets[order(clique_scores, decreasing = TRUE)]

    letters_vec <- letters[1:length(subsets)]
    cld_result <- setNames(character(k), env_levels)
    for (e in env_levels) {
      assigned <- character(0)
      for (ci in seq_along(subsets)) {
        if (e %in% subsets[[ci]]) {
          assigned <- c(assigned, letters_vec[ci])
        }
      }
      cld_result[e] <- paste(assigned, collapse = "")
    }
    return(cld_result)
  }

  for (tr in traits_sec3) {
    # (A) Fit primary repeated-measures LMM
    form_lmer <- as.formula(paste(tr, "~ Env_Factor + (1 | Genotype)"))
    m_lmer <- lmer(form_lmer, data = blues_df)

    is_sing <- isSingular(m_lmer)
    vc <- as.data.frame(VarCorr(m_lmer))
    sigma_g2 <- vc$vcov[vc$grp == "Genotype"]
    sigma_e2 <- vc$vcov[vc$grp == "Residual"]
    adj_icc <- sigma_g2 / (sigma_g2 + sigma_e2)

    an_f <- car::Anova(m_lmer, type = "III", test.statistic = "F")
    df_num <- an_f$Df[2]
    df_den <- an_f$Df.res[2]
    f_stat <- an_f$F[2]
    p_stat <- an_f$`Pr(>F)`[2]

    # (B) Residuals and Heteroskedasticity Audit
    res_vals <- residuals(m_lmer)
    fit_vals <- fitted(m_lmer)
    res_df <- data.frame(fitted = fit_vals, resid = res_vals, Env_Factor = blues_df$Env_Factor)

    sd_by_env <- res_df %>%
      group_by(Env_Factor) %>%
      summarise(Residual_SD = sd(resid), .groups = "drop")
    sd_lookup <- setNames(sd_by_env$Residual_SD, as.character(sd_by_env$Env_Factor))
    min_sd <- min(sd_by_env$Residual_SD)
    max_sd <- max(sd_by_env$Residual_SD)
    sd_ratio <- max_sd / min_sd
    exceeds_15 <- sd_ratio > 1.5

    hetero_df <- sd_by_env %>%
      mutate(
        Trait = tr,
        Min_SD = min_sd,
        Max_SD = max_sd,
        SD_Ratio = sd_ratio,
        Heteroskedastic_Exceeds_Threshold = exceeds_15
      )
    hetero_list[[tr]] <- hetero_df

    # (C) Post-Hoc Tukey Pairwise Contrasts via emmeans
    emm_lmer <- emmeans(m_lmer, ~ Env_Factor)
    prs_lmer <- as.data.frame(pairs(emm_lmer, adjust = "tukey"))
    ci_lmer  <- as.data.frame(confint(pairs(emm_lmer, adjust = "tukey")))
    prs_lmer$lower.CL_lmer <- ci_lmer$lower.CL
    prs_lmer$upper.CL_lmer <- ci_lmer$upper.CL

    # (D) Unconditional Heteroskedastic Unstructured GLS (varIdent + corSymm)
    aic_lmer <- AIC(m_lmer)
    m_gls <- nlme::gls(
      as.formula(paste(tr, "~ Env_Factor")),
      data = blues_sorted,
      weights = nlme::varIdent(form = ~ 1 | Env_Factor),
      correlation = nlme::corSymm(form = ~ time_idx | Genotype),
      method = "REML"
    )
    aic_gls <- AIC(m_gls)
    emm_gls <- emmeans(m_gls, ~ Env_Factor)
    prs_gls <- as.data.frame(pairs(emm_gls, adjust = "tukey"))
    ci_gls  <- as.data.frame(confint(pairs(emm_gls, adjust = "tukey")))
    prs_gls$lower.CL_gls <- ci_gls$lower.CL
    prs_gls$upper.CL_gls <- ci_gls$upper.CL

    # (E) Wide-Format Paired Wilcoxon Tests with Holm Correction
    wide_df <- blues_df %>%
      select(Genotype, Env_Factor, all_of(tr)) %>%
      pivot_wider(names_from = Env_Factor, values_from = all_of(tr)) %>%
      arrange(Genotype)
    stopifnot(nrow(wide_df) == n_geno_sec3, !any(is.na(wide_df)), all(env_levels %in% colnames(wide_df)))
    n_pairs <- nrow(wide_df)

    wilcox_tr <- list()
    for (pr in all_pairs) {
      e1 <- pr[1]; e2 <- pr[2]
      w_res <- wilcox.test(wide_df[[e1]], wide_df[[e2]], paired = TRUE, exact = FALSE)
      wilcox_tr[[paste0(e1, " - ", e2)]] <- tibble(
        contrast = paste0(e1, " - ", e2),
        W_stat_wilcox = as.numeric(w_res$statistic),
        p_unadj_wilcox = w_res$p.value,
        N_pairs = n_pairs
      )
    }
    df_wilcox <- bind_rows(wilcox_tr)
    df_wilcox$p_holm_wilcox <- p.adjust(df_wilcox$p_unadj_wilcox, method = "holm")
    df_wilcox$signif_holm_wilcox <- case_when(
      df_wilcox$p_holm_wilcox < 0.0001 ~ "****",
      df_wilcox$p_holm_wilcox < 0.001 ~ "***",
      df_wilcox$p_holm_wilcox < 0.01 ~ "**",
      df_wilcox$p_holm_wilcox < 0.05 ~ "*",
      TRUE ~ "ns"
    )

    # (F) Comprehensive Contrasts & Heteroskedasticity Master Table (Table S3)
    contrasts_df <- prs_lmer %>%
      mutate(
        Trait = trait_titles_sec3[tr],
        Trait_Code = tr,
        Comparison = contrast,
        Env_1 = trimws(sapply(strsplit(as.character(contrast), " - "), `[`, 1)),
        Env_2 = trimws(sapply(strsplit(as.character(contrast), " - "), `[`, 2)),
        Residual_SD_Env1 = round(sd_lookup[Env_1], 4),
        Residual_SD_Env2 = round(sd_lookup[Env_2], 4),
        Trait_SD_Ratio_Max_Min = round(sd_ratio, 2),
        Contrast_Design_Status = case_when(
          contrast == "NWPZ Baseline - NWPZ Warm" ~ "Within-Site (Delhi: Randomized Incomplete-Block Trials)",
          contrast == "PZ Warm I - PZ Warm II" ~ "Within-Site (Dharwad: Unrandomized Duplicated-Order Trials)",
          TRUE                                    ~ "Between-Site (Confounded with Site & Experimental Design Differences)"
        ),
        p_unadj_lmer = 2 * (1 - pt(abs(t.ratio), df)),
        p_tukey_lmer = p.value,
        signif_tukey_lmer = case_when(
          p_tukey_lmer < 0.0001 ~ "****",
          p_tukey_lmer < 0.001 ~ "***",
          p_tukey_lmer < 0.01 ~ "**",
          p_tukey_lmer < 0.05 ~ "*",
          TRUE ~ "ns"
        )
      ) %>%
      left_join(
        prs_gls %>% select(
          contrast,
          estimate_gls = estimate,
          SE_gls = SE,
          df_gls = df,
          t_ratio_gls = t.ratio,
          p_tukey_gls = p.value,
          lower.CL_gls = lower.CL_gls,
          upper.CL_gls = upper.CL_gls
        ) %>% mutate(
          signif_tukey_gls = case_when(
            p_tukey_gls < 0.0001 ~ "****",
            p_tukey_gls < 0.001 ~ "***",
            p_tukey_gls < 0.01 ~ "**",
            p_tukey_gls < 0.05 ~ "*",
            TRUE ~ "ns"
          )
        ),
        by = "contrast"
      ) %>%
      left_join(df_wilcox, by = "contrast")

    contrast_list[[tr]] <- contrasts_df

    # (G) ANOVA Summary Table Entry
    anova_summary_list[[tr]] <- tibble(
      Trait = tr,
      Model = "Repeated-Measures LMM & Heteroskedastic GLS",
      DF_num = df_num,
      DF_den = df_den,
      F_stat = f_stat,
      p_value = p_stat,
      sigma_g2 = sigma_g2,
      sigma_e2 = sigma_e2,
      adjusted_ICC = adj_icc,
      is_singular = is_sing,
      AIC_lmer = aic_lmer,
      AIC_gls = aic_gls,
      SD_Ratio = sd_ratio,
      Heteroskedastic_Exceeds_Threshold = exceeds_15
    )

    # (H) Compute Compact Letter Display (CLD) from LMM and Heteroskedastic GLS
    env_means_vec <- tapply(blues_df[[tr]], blues_df$Env_Factor, mean, na.rm = TRUE)
    cld_lmer_vec  <- compute_cld_letters(env_levels, env_means_vec, contrasts_df, p_col = "p_tukey_lmer", alpha = 0.05)
    cld_gls_vec   <- compute_cld_letters(env_levels, env_means_vec, contrasts_df, p_col = "p_tukey_gls",  alpha = 0.05)
    cld_combined  <- ifelse(
      cld_lmer_vec[env_levels] == cld_gls_vec[env_levels],
      cld_lmer_vec[env_levels],
      paste0(cld_lmer_vec[env_levels], " [", cld_gls_vec[env_levels], "]")
    )
    names(cld_combined) <- env_levels
    cld_store_list[[tr]] <- cld_combined

    # Store Relative Shift (% Delta vs NWPZ Baseline) for Panel F & Table 2
    base_mean_val <- env_means_vec["NWPZ Baseline"]
    stress_envs_sec3 <- setdiff(env_metadata$Env_Name, "NWPZ Baseline")
    for (e_str in stress_envs_sec3) {
      c_row <- contrasts_df %>% filter(contrast == paste0("NWPZ Baseline - ", e_str))
      # Note: contrast is Baseline - Stress, so Stress - Baseline = -estimate
      diff_val <- -c_row$estimate_gls[1]
      se_diff  <- c_row$SE_gls[1]
      pct_d    <- (diff_val / base_mean_val) * 100
      pct_se   <- (se_diff / base_mean_val) * 100
      pct_lo   <- pct_d - 1.96 * pct_se
      pct_hi   <- pct_d + 1.96 * pct_se
      sig_s    <- c_row$signif_tukey_gls[1]
      rel_shift_rows[[paste0(tr, "_", e_str)]] <- tibble(
        Trait_Code = tr,
        Trait_Short = trait_short_sec3[tr],
        Stress_Env = e_str,
        Pct_Delta = pct_d,
        Pct_SE = pct_se,
        CI_Lo = pct_lo,
        CI_Hi = pct_hi,
        Signif_GLS = sig_s,
        Signif_LMM = c_row$signif_tukey_lmer[1],
        Label_Plot = sprintf("%+.1f%%%s", pct_d, ifelse(sig_s == "ns", "", sig_s))
      )
    }

    # (I) Upgraded Paired Genotype Profile Plot with Quantile Fan Ribbons + Mean & CLD Badges + ANOVA Subtitle
    env_means <- blues_df %>%
      group_by(Env_Factor) %>%
      summarise(
        mean_val = mean(.data[[tr]], na.rm = TRUE),
        se_val   = sd(.data[[tr]], na.rm = TRUE) / sqrt(n()),
        q10_val  = as.numeric(quantile(.data[[tr]], 0.10, na.rm = TRUE)),
        q25_val  = as.numeric(quantile(.data[[tr]], 0.25, na.rm = TRUE)),
        q75_val  = as.numeric(quantile(.data[[tr]], 0.75, na.rm = TRUE)),
        q90_val  = as.numeric(quantile(.data[[tr]], 0.90, na.rm = TRUE)),
        .groups  = "drop"
      ) %>%
      mutate(
        cld_let = cld_combined[as.character(Env_Factor)],
        label_mean = if (tr == "TCE") {
          sprintf("%.1f%% (%s)", mean_val * 100, cld_let)
        } else if (tr == "GrainYld") {
          sprintf("%.2f (%s)", mean_val, cld_let)
        } else {
          sprintf("%.1f (%s)", mean_val, cld_let)
        }
      )

    p_str <- if (p_stat < 0.0001) "P < 0.0001" else sprintf("P = %.4f", p_stat)
    sub_spa <- sprintf("LMM: F(%d, %d) = %.1f, %s | Residual SD Ratio = %.2fx",
                       df_num, round(df_den), f_stat, p_str, sd_ratio)

    p_spa <- ggplot(blues_df, aes(x = Env_Factor, y = .data[[tr]], group = Genotype)) +
      geom_vline(xintercept = 2.5, linetype = "dashed", color = "grey45", linewidth = 0.65) +
      geom_ribbon(data = env_means, aes(x = Env_Factor, ymin = q10_val, ymax = q90_val, group = 1),
                  inherit.aes = FALSE, fill = "#3C5488", alpha = 0.13) +
      geom_ribbon(data = env_means, aes(x = Env_Factor, ymin = q25_val, ymax = q75_val, group = 1),
                  inherit.aes = FALSE, fill = "#3C5488", alpha = 0.24) +
      geom_line(alpha = 0.09, color = "grey40", linewidth = 0.32) +
      geom_point(alpha = 0.16, color = "grey40", size = 0.85) +
      geom_line(data = env_means, aes(x = Env_Factor, y = mean_val, group = 1),
                inherit.aes = FALSE, color = "#DC0000", linewidth = 1.50) +
      geom_point(data = env_means, aes(x = Env_Factor, y = mean_val, group = 1),
                 inherit.aes = FALSE, color = "#DC0000", size = 3.60) +
      geom_label(
        data = env_means,
        aes(x = Env_Factor, y = mean_val, label = label_mean, group = 1),
        inherit.aes = FALSE,
        fill = "white",
        color = "#B00000",
        fontface = "bold",
        size = 3.75,
        vjust = -0.52,
        label.padding = unit(0.20, "lines"),
        label.size = 0.35
      ) +
      prof_theme +
      scale_y_continuous(
        expand = expansion(mult = c(0.05, 0.18)),
        labels = if (tr == "TCE") scales::percent_format(accuracy = 1) else ggplot2::waiver()
      ) +
      labs(
        title = trait_titles_sec3[tr],
        subtitle = sub_spa,
        y = trait_ylab_sec3[tr],
        x = ""
      )
    spaghetti_plots[[tr]] <- p_spa
  }

  # 3.3 Construct Panel F: Relative Phenotypic Shift vs. NWPZ Baseline (% Delta + 95% CI)
  stress_envs_sec3 <- setdiff(env_metadata$Env_Name, "NWPZ Baseline")
  rel_shift_df <- bind_rows(rel_shift_rows) %>%
    mutate(
      Trait_Factor = factor(Trait_Short, levels = rev(unname(trait_short_sec3))),
      Stress_Factor = factor(Stress_Env, levels = rev(stress_envs_sec3))
    )

  stress_cols_sec3 <- c(
    "NWPZ Warm"  = "#4DBBD5",
    "PZ Warm I"  = "#00A087",
    "PZ Warm II" = "#3C5488"
  )

  p_rel_shift <- ggplot(rel_shift_df, aes(x = Pct_Delta, y = Trait_Factor, color = Stress_Factor, group = Stress_Factor)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "#E64B35", linewidth = 0.85) +
    geom_errorbar(aes(xmin = CI_Lo, xmax = CI_Hi),
                  position = position_dodge(width = 0.74), width = 0.28, linewidth = 0.90, orientation = "y") +
    geom_point(position = position_dodge(width = 0.74), size = 3.40) +
    geom_text(
      aes(
        x = ifelse(Pct_Delta >= 0, CI_Hi + 0.9, CI_Lo - 0.9),
        label = Label_Plot,
        hjust = ifelse(Pct_Delta >= 0, 0, 1)
      ),
      position = position_dodge(width = 0.74),
      size = 3.45,
      fontface = "bold",
      show.legend = FALSE
    ) +
    annotate("text", x = 0.6, y = "Final Veg. Tillers (T4)", label = "NWPZ Baseline (0%)", hjust = 0, vjust = -2.3,
             size = 3.40, fontface = "bold.italic", color = "#B03A2E") +
    scale_color_manual(
      values = stress_cols_sec3,
      breaks = stress_envs_sec3,
      name = "Stress Regime:"
    ) +
    scale_x_continuous(
      breaks = seq(-30, 20, by = 10),
      labels = function(x) sprintf("%+d%%", x)
    ) +
    scale_y_discrete(expand = expansion(mult = c(0.10, 0.17))) +
    coord_cartesian(xlim = c(-43, 27), clip = "off") +
    prof_theme +
    theme(
      legend.position = "bottom",
      legend.title = element_text(face = "bold", size = 10.8),
      legend.text = element_text(face = "bold", size = 10.4, margin = margin(r = 6)),
      legend.key.width = unit(0.65, "cm"),
      legend.spacing.x = unit(0.16, "cm"),
      legend.margin = margin(t = -2, r = 8, b = 0, l = 0),
      axis.text.x = element_text(angle = 0, hjust = 0.5, face = "bold", size = 11.0, color = "grey10"),
      axis.text.y = element_text(face = "bold", size = 11.4, color = "grey10"),
      axis.title.x = element_text(face = "bold", size = 12.0, color = "black"),
      panel.grid.major.x = element_line(color = "grey90", linewidth = 0.35),
      panel.grid.major.y = element_line(color = "grey92", linewidth = 0.35),
      plot.margin = margin(t = 5, r = 14, b = 4, l = 6)
    ) +
    guides(color = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(size = 3.8, linewidth = 1.1))) +
    labs(
      title = "Relative Trait Shift vs. NWPZ Baseline (%)",
      subtitle = "Heteroskedastic GLS Mean % \u0394 (\u00b195% CI) vs. NWPZ Baseline",
      x = "Relative Phenotypic Shift (% \u0394 from NWPZ Baseline)",
      y = ""
    )

  # 3.4 Save Unified 6-Panel Figure (Global_Parameter_Paired_Genotype_Profiles.png)
  p_spaghetti_comp <- (spaghetti_plots[["T4_Count"]] | spaghetti_plots[["PT_Count"]]) /
                      (spaghetti_plots[["TCE"]]      | spaghetti_plots[["TKW"]]) /
                      (spaghetti_plots[["GrainYld"]] | p_rel_shift) +
                      plot_annotation(
                        tag_levels = 'A',
                        caption = paste0(
                          "Panels A\u2013E: Grey trajectories = N = 200 paired genotypes across thermal regimes; blue ribbons = 10\u201390% (interdecile) and 25\u201375% (IQR) quantile fans; red line & badges = regime mean (CLD).\n",
                          "Dashed vertical line separates NWPZ (Delhi) and PZ (Dharwad) agro-climatic zones. Regimes sharing a letter do not differ significantly (Tukey HSD, P \u2265 0.05; bracketed letters = heteroskedastic GLS where differing)."
                        ),
                        theme = theme(
                          plot.caption = element_text(size = 9.8, hjust = 0, color = "grey20", face = "italic", lineheight = 1.15, margin = margin(t = 6, b = 2))
                        )
                      ) &
                      theme(plot.tag = element_text(face = "bold", size = 15.5, color = "black"))
  ggsave(file.path(sec21_dir, "Global_Parameter_Paired_Genotype_Profiles.png"),
         plot = p_spaghetti_comp, width = 14.2, height = 14.8, dpi = 600)

  # 3.5 Construct Consolidated Main Text Table 2 & Supplementary Table S3 + Styled Excel Workbook
  table_anova_summary <- bind_rows(anova_summary_list)
  table_contrasts_raw <- bind_rows(contrast_list)

  table2_rows <- list()
  for (tr in traits_sec3) {
    cld_v <- cld_store_list[[tr]]
    fmt_regime_cell <- function(e_nm) {
      vals <- blues_df[[tr]][blues_df$Env_Factor == e_nm]
      m <- mean(vals, na.rm = TRUE)
      s <- sd(vals, na.rm = TRUE) / sqrt(length(vals))
      q25 <- as.numeric(quantile(vals, 0.25, na.rm = TRUE))
      q75 <- as.numeric(quantile(vals, 0.75, na.rm = TRUE))
      c_let <- cld_v[e_nm]
      if (tr == "TCE") {
        sprintf("%.1f \u00b1 %.1f%% [%.1f\u2013%.1f%%] (%s)", m * 100, s * 100, q25 * 100, q75 * 100, c_let)
      } else if (tr == "GrainYld") {
        sprintf("%.3f \u00b1 %.3f [%.2f\u2013%.2f] (%s)", m, s, q25, q75, c_let)
      } else {
        sprintf("%.2f \u00b1 %.2f [%.2f\u2013%.2f] (%s)", m, s, q25, q75, c_let)
      }
    }

    rel_w   <- rel_shift_df %>% filter(Trait_Code == tr, Stress_Env == "NWPZ Warm")
    rel_pw1 <- rel_shift_df %>% filter(Trait_Code == tr, Stress_Env == "PZ Warm I")
    rel_pw2 <- rel_shift_df %>% filter(Trait_Code == tr, Stress_Env == "PZ Warm II")
    an_row  <- table_anova_summary %>% filter(Trait == tr)

    table2_rows[[tr]] <- tibble(
      Trait = trait_titles_sec3[tr],
      Unit = trait_units_sec3[tr],
      NWPZ_Baseline_Mean_SE_IQR_CLD    = fmt_regime_cell("NWPZ Baseline"),
      NWPZ_Warm_Mean_SE_IQR_CLD        = fmt_regime_cell("NWPZ Warm"),
      Delta_Pct_NWPZ_Warm_vs_Baseline  = sprintf("%+.1f%% (%s)", rel_w$Pct_Delta[1], rel_w$Signif_GLS[1]),
      PZ_Warm_I_Mean_SE_IQR_CLD        = fmt_regime_cell("PZ Warm I"),
      Delta_Pct_PZ_Warm_I_vs_Baseline  = sprintf("%+.1f%% (%s)", rel_pw1$Pct_Delta[1], rel_pw1$Signif_GLS[1]),
      PZ_Warm_II_Mean_SE_IQR_CLD       = fmt_regime_cell("PZ Warm II"),
      Delta_Pct_PZ_Warm_II_vs_Baseline = sprintf("%+.1f%% (%s)", rel_pw2$Pct_Delta[1], rel_pw2$Signif_GLS[1]),
      LMM_Omnibus_F_df              = sprintf("F(%d, %d) = %.2f", an_row$DF_num[1], round(an_row$DF_den[1]), an_row$F_stat[1]),
      Omnibus_P                     = if (an_row$p_value[1] < 0.0001) "P < 0.0001" else sprintf("P = %.4f", an_row$p_value[1]),
      Genotype_Repeatability_ICC    = round(an_row$adjusted_ICC[1], 3),
      Residual_SD_Ratio_Max_Min     = round(an_row$SD_Ratio[1], 2),
      AIC_Homoskedastic_LMM         = round(an_row$AIC_lmer[1], 1),
      AIC_Heteroskedastic_GLS       = round(an_row$AIC_gls[1], 1),
      Delta_AIC_GLS_minus_LMM       = round(an_row$AIC_gls[1] - an_row$AIC_lmer[1], 1)
    )
  }
  table_2_main <- bind_rows(table2_rows)

  table_s3_supp <- table_contrasts_raw %>%
    select(
      Trait,
      Comparison,
      Contrast_Design_Status,
      Residual_SD_Env1,
      Residual_SD_Env2,
      Trait_SD_Ratio_Max_Min,
      LMM_Estimate = estimate,
      LMM_SE = SE,
      LMM_df = df,
      LMM_t_ratio = t.ratio,
      LMM_CI_Lower = lower.CL_lmer,
      LMM_CI_Upper = upper.CL_lmer,
      LMM_P_Tukey = p_tukey_lmer,
      LMM_Signif = signif_tukey_lmer,
      GLS_Estimate = estimate_gls,
      GLS_SE = SE_gls,
      GLS_df = df_gls,
      GLS_t_ratio = t_ratio_gls,
      GLS_CI_Lower = lower.CL_gls,
      GLS_CI_Upper = upper.CL_gls,
      GLS_P_Tukey = p_tukey_gls,
      GLS_Signif = signif_tukey_gls,
      Wilcoxon_W = W_stat_wilcox,
      Wilcoxon_P_Holm = p_holm_wilcox,
      Wilcoxon_Signif = signif_holm_wilcox,
      N_Paired_Genotypes = N_pairs
    )

  write.csv(table_2_main,  file.path(sec21_dir, "Table_2_Multi_Environment_Phenotypic_and_ANOVA_Summary.csv"), row.names = FALSE)
  utils::write.table(sanitize_df_for_export(table_2_main), file.path(sec21_dir, "Table_2_Multi_Environment_Phenotypic_and_ANOVA_Summary.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
  write.csv(table_s3_supp, file.path(sec21_dir, "Table_S3_Pairwise_Regime_Contrasts_and_Heteroskedasticity.csv"), row.names = FALSE)

  # Export Unified Multi-Tab Excel Workbook for Section 2.1
  if (requireNamespace("openxlsx", quietly = TRUE)) {
    wb_sec21 <- openxlsx::createWorkbook()
    hdr_st21 <- openxlsx::createStyle(
      fontName = "Calibri", fontSize = 11, fontColour = "#FFFFFF",
      fgFill = "#1B4F72", halign = "center", valign = "center",
      textDecoration = "bold", wrapText = TRUE, border = "TopBottomLeftRight"
    )
    body_st21 <- openxlsx::createStyle(
      fontName = "Calibri", fontSize = 10.5, valign = "center",
      border = "TopBottomLeftRight", borderColour = "#D5D8DC"
    )
    openxlsx::addWorksheet(wb_sec21, "Table_2_Phenotypic_ANOVA")
    openxlsx::writeData(wb_sec21, "Table_2_Phenotypic_ANOVA", table_2_main, headerStyle = hdr_st21)
    openxlsx::addStyle(wb_sec21, "Table_2_Phenotypic_ANOVA", body_st21,
                       rows = 2:(nrow(table_2_main) + 1), cols = 1:ncol(table_2_main), gridExpand = TRUE)
    openxlsx::setColWidths(wb_sec21, "Table_2_Phenotypic_ANOVA", cols = 1:ncol(table_2_main), widths = "auto")

    openxlsx::addWorksheet(wb_sec21, "Table_S3_Pairwise_Contrasts")
    openxlsx::writeData(wb_sec21, "Table_S3_Pairwise_Contrasts", table_s3_supp, headerStyle = hdr_st21)
    openxlsx::addStyle(wb_sec21, "Table_S3_Pairwise_Contrasts", body_st21,
                       rows = 2:(nrow(table_s3_supp) + 1), cols = 1:ncol(table_s3_supp), gridExpand = TRUE)
    openxlsx::setColWidths(wb_sec21, "Table_S3_Pairwise_Contrasts", cols = 1:ncol(table_s3_supp), widths = "auto")

    openxlsx::saveWorkbook(wb_sec21, file.path(sec21_dir, "Section_2.1_Manuscript_and_Supplementary_Tables.xlsx"), overwrite = TRUE)
  }

  cat("\n--- Consolidated Table 2: Multi-Environment Phenotypic & ANOVA Summary ---\n")
  print(table_2_main %>% select(Trait, NWPZ_Baseline_Mean_SE_IQR_CLD, Delta_Pct_NWPZ_Warm_vs_Baseline,
                                Delta_Pct_PZ_Warm_I_vs_Baseline, Delta_Pct_PZ_Warm_II_vs_Baseline, LMM_Omnibus_F_df, Residual_SD_Ratio_Max_Min))
} else {
  cat("  -> Section 3 [SKIPPED - Assets Finalized & Protected]\n")
}


# ==============================================================================
# SECTION 4: UNSUPERVISED PHENOLOGY CLASSIFICATIONS & K-MEANS TRAJECTORIES
# ==============================================================================
cat(">>> Executing Section 4: Phenology Trajectories & K-Means Clustering...\n")

# 4.1 K-Means Clustering on Tillering Trajectory Shape in NWPZ Baseline
# Genotype-wise normalisation (T_i / T4) isolates developmental curve shape from plant size / tiller capacity.
# In baseline, T4_Count represents the terminal vegetative tiller capacity.
# T1/T4, T2/T4, T3/T4 capture the relative developmental progression.
# (T4/T4 is identically 1.0 across genotypes in baseline with zero variance, so T1_rel..T3_rel
# completely parameterize trajectory shape).
stopifnot(all(blues_df$T4_Count > 0, na.rm = TRUE))
baseline_blues <- blues_df %>% filter(Environment == "NWPZ Baseline")

base_shape_df <- baseline_blues %>%
  transmute(
    T1_rel = T1_Count / T4_Count,
    T2_rel = T2_Count / T4_Count,
    T3_rel = T3_Count / T4_Count
  )
traj_matrix <- scale(as.matrix(base_shape_df))

# Optimize Cluster Count via Silhouette Width
# Note: Silhouette width frequently selects k = 2 for continuous unimodal distributions;
# this partition represents an operational categorization along a developmental continuum rather
# than evidence of discrete biological taxa.
set.seed(42)
k_eval_range <- 2:6
sil_eval_df <- data.frame(k = k_eval_range, Mean_Silhouette_Width = numeric(length(k_eval_range)))
for (i in seq_along(k_eval_range)) {
  k_val <- k_eval_range[i]
  km_test <- kmeans(traj_matrix, centers = k_val, nstart = 50)
  ss <- silhouette(km_test$cluster, dist(traj_matrix))
  sil_eval_df$Mean_Silhouette_Width[i] <- mean(ss[, 3])
}

# Full range for plotting
k_range <- 2:10
sil_widths <- numeric(length(k_range))
sil_widths[1:5] <- sil_eval_df$Mean_Silhouette_Width
for (i in 6:length(k_range)) {
  k_val <- k_range[i]
  km_test <- kmeans(traj_matrix, centers = k_val, nstart = 50)
  ss <- silhouette(km_test$cluster, dist(traj_matrix))
  sil_widths[i] <- mean(ss[, 3])
}

sil_df <- data.frame(k = k_range, Silhouette = sil_widths)
optimal_k <- sil_eval_df$k[which.max(sil_eval_df$Mean_Silhouette_Width)]
max_sil_width <- max(sil_eval_df$Mean_Silhouette_Width)
stopifnot(optimal_k == 2L)

# Final K-Means using the dynamically computed optimal_k (verified k == 2)
set.seed(42)
final_km <- kmeans(traj_matrix, centers = optimal_k, nstart = 100)

# Unscale centroids back to relative shape scale and compute raw T4
base_shape_center <- attr(traj_matrix, "scaled:center")
base_shape_scale  <- attr(traj_matrix, "scaled:scale")

centroid_unscaled <- t(apply(final_km$centers, 1, function(r) r * base_shape_scale + base_shape_center))
cluster_counts <- table(final_km$cluster)

t4_means <- aggregate(baseline_blues$T4_Count, by = list(Cluster = final_km$cluster), FUN = mean)

# Descriptive naming based strictly on centroid geometry:
# Cluster with higher T3_rel reached relative peak at T3: "Earlier Proportional Peak"
# Cluster with lower T3_rel continues vegetative expansion through T4: "Later Proportional Peak"
earlier_cluster_idx <- which.max(centroid_unscaled[, "T3_rel"])
later_cluster_idx   <- which.min(centroid_unscaled[, "T3_rel"])

traj_names <- character(2)
traj_names[earlier_cluster_idx] <- "Earlier Proportional Peak"
traj_names[later_cluster_idx]   <- "Later Proportional Peak"

centroid_table <- data.frame(
  Cluster = as.integer(rownames(centroid_unscaled)),
  Descriptive_Name = traj_names[as.integer(rownames(centroid_unscaled))],
  Size_n  = as.integer(cluster_counts[rownames(centroid_unscaled)]),
  Percentage = as.numeric(cluster_counts[rownames(centroid_unscaled)] / sum(cluster_counts) * 100),
  T1_rel_Mean = round(centroid_unscaled[, "T1_rel"], 4),
  T2_rel_Mean = round(centroid_unscaled[, "T2_rel"], 4),
  T3_rel_Mean = round(centroid_unscaled[, "T3_rel"], 4),
  T4_Count_Mean = round(t4_means$x[match(as.integer(rownames(centroid_unscaled)), t4_means$Cluster)], 2),
  stringsAsFactors = FALSE
)

cluster_labels <- traj_names[final_km$cluster]
baseline_blues$Trajectory <- cluster_labels

# 4.2 Fixed-Rule Environmental Classification & Repeated-Measures Transition Analysis
# Genotypes are classified in EACH environment using the identical fixed rule (baseline-trained centroids).
# Transitions across regimes (Baseline x Stress) are tabulated and formally tested via McNemar's test.

traj_centers <- final_km$centers

assign_env_traj <- function(row_counts) {
  # row_counts contains c(T1_Count, T2_Count, T3_Count, T4_Count)
  t4_denom <- row_counts[4]
  row_shape <- c(
    T1_rel = row_counts[1] / t4_denom,
    T2_rel = row_counts[2] / t4_denom,
    T3_rel = row_counts[3] / t4_denom
  )
  sc <- (row_shape - base_shape_center) / base_shape_scale
  dists <- apply(traj_centers, 1, function(cen) sum((sc - cen)^2))
  traj_names[which.min(dists)]
}

# Tillering Vigor Model (Option A): K-Means on duration-normalized Mean_Tiller_Early & Mean_Tiller_Late
# trained on NWPZ Baseline, eliminating the mechanical 40-42 d vs 72 d window-length confound.
if (!"Int_E_T1" %in% names(blues_df)) {
  blues_df <- blues_df %>%
    left_join(env_metadata %>% select(Env_Name, Int_E_T1, Int_T1_T2, Int_T2_T3, Int_T3_T4, Total_E_T4),
              by = c("Environment" = "Env_Name"))
}

blues_df <- blues_df %>% mutate(
  Dur_Early         = Int_E_T1 + Int_T1_T2,
  Dur_Late          = Int_T2_T3 + Int_T3_T4,
  AUTPC_Early       = (0.5 * (COUNT_AT_EMERGENCE_MAIN + T1_Count) * Int_E_T1) + (0.5 * (T1_Count + T2_Count) * Int_T1_T2),
  AUTPC_Late        = (0.5 * (T2_Count + T3_Count) * Int_T2_T3) + (0.5 * (T3_Count + T4_Count) * Int_T3_T4),
  Mean_Tiller_Early = AUTPC_Early / Dur_Early,
  Mean_Tiller_Late  = AUTPC_Late  / Dur_Late
)
stopifnot(max(abs((blues_df$AUTPC_Early + blues_df$AUTPC_Late) - blues_df$AUTPC_From_Count_BLUEs), na.rm = TRUE) < 1e-10)

base_v_df <- blues_df %>% filter(Environment == "NWPZ Baseline")
r_early_late_base <- cor(base_v_df$Mean_Tiller_Early, base_v_df$Mean_Tiller_Late)
cat(sprintf("  r(Mean_Tiller_Early, Mean_Tiller_Late) in NWPZ Baseline: %.4f\n", r_early_late_base))

base_v_mat    <- as.matrix(base_v_df[, c("Mean_Tiller_Early", "Mean_Tiller_Late")])
base_v_scaled <- scale(base_v_mat)
v_center <- attr(base_v_scaled, "scaled:center")
v_scale  <- attr(base_v_scaled, "scaled:scale")

set.seed(42)
km_vigor <- kmeans(base_v_scaled, centers = 2, nstart = 50)
clus_size <- tapply(base_v_df$AUTPC_From_Count_BLUEs / base_v_df$Total_E_T4, km_vigor$cluster, mean)
vigor_names <- character(2)
vigor_names[as.integer(names(which.max(clus_size)))] <- VIGOR_HI
vigor_names[as.integer(names(which.min(clus_size)))] <- VIGOR_LO
vigor_centers <- km_vigor$centers

assign_env_vigor <- function(row_v) {
  sv <- (row_v - v_center) / v_scale
  vigor_names[which.min(apply(vigor_centers, 1, function(cen) sum((sv - cen)^2)))]
}

# Strip any pre-existing classification columns for idempotency before re-joining
blues_df <- blues_df %>%
  select(-any_of(c("Trajectory_Baseline", "Vigor_Baseline", "TCE_Baseline", "Trajectory", "Vigor",
                   "Trajectory_Realized", "Vigor_Realized", "Realized_TCE", "T1_rel", "T2_rel", "T3_rel")))

# Compute genotype-normalized developmental ratios across all environments
blues_df <- blues_df %>%
  mutate(
    T1_rel = T1_Count / T4_Count,
    T2_rel = T2_Count / T4_Count,
    T3_rel = T3_Count / T4_Count
  )

# Apply the fixed classification rules to ALL environments
blues_df$Trajectory_Realized <- apply(blues_df[, c("T1_Count", "T2_Count", "T3_Count", "T4_Count")], 1, assign_env_traj)
blues_df$Trajectory_Realized <- factor(blues_df$Trajectory_Realized, levels = c("Earlier Proportional Peak", "Later Proportional Peak"))

blues_df$Vigor_Realized <- factor(
  apply(blues_df[, c("Mean_Tiller_Early", "Mean_Tiller_Late")], 1, assign_env_vigor),
  levels = c(VIGOR_HI, VIGOR_LO)
)

# Tiller Conversion Efficiency (TCE) classification anchored to NWPZ Baseline median
tce_base_vals <- blues_df %>% filter(Environment == "NWPZ Baseline") %>% pull(TCE)
tce_cutoff    <- median(tce_base_vals, na.rm = TRUE)
blues_df <- blues_df %>%
  mutate(
    Realized_TCE = factor(
      ifelse(TCE >= tce_cutoff, "High Conversion", "Low Conversion"),
      levels = c("High Conversion", "Low Conversion")
    )
  )

# Map baseline constitutive strategies to each genotype for downstream potential yield analysis
base_strategies <- blues_df %>% filter(Environment == "NWPZ Baseline") %>%
  select(
    Genotype,
    Trajectory_Baseline = Trajectory_Realized,
    Vigor_Baseline      = Vigor_Realized,
    TCE_Baseline        = Realized_TCE
  )
blues_df <- blues_df %>%
  left_join(base_strategies, by = "Genotype")

# Alias Trajectory and Vigor to Realized strategies for seamless downstream access
blues_df$Trajectory <- blues_df$Trajectory_Realized
blues_df$Vigor      <- blues_df$Vigor_Realized

if (isTRUE(RUN_SECTION_4_PHENOLOGY)) {
  sec22_dir <- file.path(assets_dir, "Section_2.2_Phenology_Classifications")
  dir.create(sec22_dir, showWarnings = FALSE, recursive = TRUE)

  # Archive all legacy/redundant figures and fragmented CSV/Excel tables safely
  archive_sec22_dir <- file.path(project_root, "Archive_Rejected", "Section_2.2_Legacy_Classifications")
  dir.create(archive_sec22_dir, showWarnings = FALSE, recursive = TRUE)

  legacy_sec22_patterns <- c(
    "^A_Peak_Phase_TrendLines\\.png$",
    "^B_KMeans_Shifts_TrendLines\\.png$",
    "^PCA_Cluster_Plot\\.png$",
    "^Silhouette_Optimization_Plot\\.png$",
    "^Table_Cluster_Validity_Silhouette_k2_k6\\.csv$",
    "^Table_Trajectory_Shape_Cluster_Centroids\\.csv$",
    "^Table_Trajectory_Transition.*\\.csv$",
    "^Table_Vigor_Transition.*\\.csv$",
    "^Table_TCE_Transition.*\\.csv$",
    "^Table_Phenological_Transition_Tests\\.(csv|xlsx)$",
    "^Table_TCE_Transition_Tests\\.(csv|xlsx)$",
    "^Phenology_Cluster_Assignments\\.csv$"
  )
  for (pat in legacy_sec22_patterns) {
    matched_files <- list.files(sec22_dir, pattern = pat, full.names = TRUE)
    for (mf in matched_files) {
      if (file.exists(mf)) {
        file.copy(mf, file.path(archive_sec22_dir, basename(mf)), overwrite = TRUE)
        unlink(mf)
      }
    }
  }
  old_comp_trans <- file.path(sec22_dir, "Compound_Phenological_Transitions_Vigor_Trajectory_TCE.png")
  if (file.exists(old_comp_trans)) {
    file.copy(old_comp_trans, file.path(archive_sec22_dir, "Compound_Phenological_Transitions_3Panel_Prior.png"), overwrite = TRUE)
  }

  # ----------------------------------------------------------------------------
  # 4.3 Consolidated Supplementary Table S4: K-Means Cluster Centroids, Typology
  #     Summaries (Trajectory, Vigor, TCE) & Silhouette Validity (k = 2..6)
  # ----------------------------------------------------------------------------
  base_full_df <- blues_df %>% filter(Environment == "NWPZ Baseline")

  summarize_typology_group <- function(df_sub, dim_label, class_col, rule_note) {
    df_sub %>%
      group_by(Class_Name = .data[[class_col]]) %>%
      summarise(
        Size_n = n(),
        Percentage = round(100 * n() / nrow(df_sub), 1),
        T1_rel_Mean_SD = sprintf("%.3f \u00B1 %.3f", mean(T1_rel, na.rm = TRUE), sd(T1_rel, na.rm = TRUE)),
        T2_rel_Mean_SD = sprintf("%.3f \u00B1 %.3f", mean(T2_rel, na.rm = TRUE), sd(T2_rel, na.rm = TRUE)),
        T3_rel_Mean_SD = sprintf("%.3f \u00B1 %.3f", mean(T3_rel, na.rm = TRUE), sd(T3_rel, na.rm = TRUE)),
        T4_Count_Mean_SD = sprintf("%.2f \u00B1 %.2f", mean(T4_Count, na.rm = TRUE), sd(T4_Count, na.rm = TRUE)),
        Mean_Tiller_Early_Mean_SD = sprintf("%.2f \u00B1 %.2f", mean(Mean_Tiller_Early, na.rm = TRUE), sd(Mean_Tiller_Early, na.rm = TRUE)),
        Mean_Tiller_Late_Mean_SD  = sprintf("%.2f \u00B1 %.2f", mean(Mean_Tiller_Late, na.rm = TRUE), sd(Mean_Tiller_Late, na.rm = TRUE)),
        AUTPC_Early_Mean_SD = sprintf("%.1f \u00B1 %.1f", mean(AUTPC_Early, na.rm = TRUE), sd(AUTPC_Early, na.rm = TRUE)),
        AUTPC_Late_Mean_SD  = sprintf("%.1f \u00B1 %.1f", mean(AUTPC_Late, na.rm = TRUE), sd(AUTPC_Late, na.rm = TRUE)),
        AUTPC_Total_Mean_SD = sprintf("%.1f \u00B1 %.1f", mean(AUTPC, na.rm = TRUE), sd(AUTPC, na.rm = TRUE)),
        AUTPC_norm_Mean_SD  = sprintf("%.2f \u00B1 %.2f", mean(AUTPC_norm, na.rm = TRUE), sd(AUTPC_norm, na.rm = TRUE)),
        Early_to_Late_Rate_Ratio = round(mean(Mean_Tiller_Early, na.rm = TRUE) / mean(Mean_Tiller_Late, na.rm = TRUE), 3),
        PT_Count_Mean_SD = sprintf("%.2f \u00B1 %.2f", mean(PT_Count, na.rm = TRUE), sd(PT_Count, na.rm = TRUE)),
        TCE_Mean_SD      = sprintf("%.3f \u00B1 %.3f", mean(TCE, na.rm = TRUE), sd(TCE, na.rm = TRUE)),
        .groups = "drop"
      ) %>%
      mutate(
        Classification_Dimension = dim_label,
        Class_Name = as.character(Class_Name),
        Partition_Rule_and_Validity = rule_note
      ) %>%
      select(
        Classification_Dimension, Class_Name, Size_n, Percentage,
        T1_rel_Mean_SD, T2_rel_Mean_SD, T3_rel_Mean_SD, T4_Count_Mean_SD,
        Mean_Tiller_Early_Mean_SD, Mean_Tiller_Late_Mean_SD,
        AUTPC_Early_Mean_SD, AUTPC_Late_Mean_SD, AUTPC_Total_Mean_SD, AUTPC_norm_Mean_SD, Early_to_Late_Rate_Ratio,
        PT_Count_Mean_SD, TCE_Mean_SD, Partition_Rule_and_Validity
      )
  }

  sil_summary_str <- paste(
    sprintf("k=%d: %.3f", sil_eval_df$k, sil_eval_df$Mean_Silhouette_Width),
    collapse = "; "
  )
  ss_vigor <- silhouette(km_vigor$cluster, dist(base_v_scaled))
  sil_vigor_k2 <- mean(ss_vigor[, 3])

  table_s4_centroids <- bind_rows(
    summarize_typology_group(
      base_full_df,
      "Tillering Trajectory Shape (K-Means on T1/T4, T2/T4, T3/T4)",
      "Trajectory_Realized",
      sprintf("Unsupervised K-Means (k = 2 optimal; Mean Silhouette Width: %s)", sil_summary_str)
    ),
    summarize_typology_group(
      base_full_df,
      "Tillering Vigor (K-Means on Mean_Tiller_Early, Mean_Tiller_Late)",
      "Vigor_Realized",
      sprintf("Unsupervised K-Means on duration-normalized Mean_Tiller_Early & Mean_Tiller_Late in NWPZ Baseline (r = %.3f; k = 2 Silhouette = %.3f; %s = higher daily tillering rate AUTPC/Total_E_T4)", r_early_late_base, sil_vigor_k2, VIGOR_HI)
    ),
    summarize_typology_group(
      base_full_df,
      "Tiller Conversion Efficiency (TCE = PT / T4)",
      "Realized_TCE",
      sprintf("NWPZ Baseline Median Split (Cutoff TCE >= %.3f; 50.0%% High vs 50.0%% Low Conversion)", tce_cutoff)
    )
  )

  write.csv(
    table_s4_centroids,
    file.path(sec22_dir, "Table_S4_KMeans_Cluster_Centroids_and_Silhouette_Validity.csv"),
    row.names = FALSE
  )

  # ----------------------------------------------------------------------------
  # 4.4 Consolidated Main Text Table 3: Multi-Environment Phenological Typology
  #     and Pairwise McNemar Regime Transitions (Trajectory, Vigor, TCE)
  # ----------------------------------------------------------------------------
  build_unified_transition_rows <- function(dim_label, col_name, lvl1, lvl2) {
    w_df <- blues_df %>%
      select(Genotype, Environment, Val = all_of(col_name)) %>%
      pivot_wider(names_from = Environment, values_from = Val)

    n_tot <- nrow(w_df)
    rows_out <- list()

    # Also compute within-PZ (PZ Warm I -> PZ Warm II) contrast for Vigor context
    tb_pz_within <- table(
      PZ1 = factor(w_df$`PZ Warm I`, levels = c(lvl1, lvl2)),
      PZ2 = factor(w_df$`PZ Warm II`, levels = c(lvl1, lvl2))
    )
    m_pz_within <- mcnemar.test(tb_pz_within, correct = TRUE)

    for (env in c("NWPZ Warm", "PZ Warm I", "PZ Warm II")) {
      trans_scope <- if (env == "NWPZ Warm") {
        "Within-Site Randomized Incomplete-Block Contrast (NWPZ Delhi)"
      } else {
        "Cross-Zone Agro-Climatic & Unrandomized Duplicated-Order Contrast (NWPZ vs PZ; Site Confounded with Heat)"
      }

      tb <- table(
        Baseline = factor(w_df$`NWPZ Baseline`, levels = c(lvl1, lvl2)),
        Stress   = factor(w_df[[env]], levels = c(lvl1, lvl2))
      )
      m_test <- mcnemar.test(tb, correct = TRUE)
      n_disc <- as.integer(tb[lvl1, lvl2] + tb[lvl2, lvl1])
      exact_binom_p <- if (n_disc > 0L) binom.test(as.integer(tb[lvl1, lvl2]), n_disc)$p.value else 1.0

      n11 <- as.integer(tb[lvl1, lvl1]) # Retained Class 1
      n12 <- as.integer(tb[lvl1, lvl2]) # Shifted Class 1 -> Class 2
      n21 <- as.integer(tb[lvl2, lvl1]) # Shifted Class 2 -> Class 1
      n22 <- as.integer(tb[lvl2, lvl2]) # Retained Class 2

      b1_n <- n11 + n12
      b2_n <- n21 + n22
      s1_n <- n11 + n21
      s2_n <- n12 + n22

      b1_pct <- 100 * b1_n / n_tot
      b2_pct <- 100 * b2_n / n_tot
      s1_pct <- 100 * s1_n / n_tot
      s2_pct <- 100 * s2_n / n_tot
      net_shift <- s1_pct - b1_pct

      sig_lbl <- if (m_test$p.value < 0.001) "***" else if (m_test$p.value < 0.01) "**" else if (m_test$p.value < 0.05) "*" else "ns"

      pct_1_to_2 <- if (b1_n > 0) 100 * n12 / b1_n else 0.0
      pct_2_to_1 <- if (b2_n > 0) 100 * n21 / b2_n else 0.0

      # Data-driven branching on significance and sign of net_shift (Item 7)
      dir_phrase <- if (m_test$p.value >= 0.05) {
        sprintf("No statistically significant regime shift in %s", lvl1)
      } else if (net_shift > 0) {
        sprintf("Significant net expansion of %s (%+.1f%%)", lvl1, net_shift)
      } else {
        sprintf("Significant net reduction of %s (%+.1f%%)", lvl1, net_shift)
      }

      site_note <- if (env == "NWPZ Warm") {
        "[Within-Delhi randomized contrast]"
      } else if (env == "PZ Warm II" && grepl("Vigor", dim_label)) {
        sprintf("[Cross-zone vs Delhi centroids; within-Dharwad PZ Warm I -> II drops %.1f%% -> %.1f%%, McNemar p = %.2e]",
                100 * sum(w_df$`PZ Warm I` == lvl1) / n_tot, s1_pct, m_pz_within$p.value)
      } else {
        "[Cross-zone vs Delhi centroids; site & layout confounded with thermal regime]"
      }

      bio_infer <- sprintf(
        "%s (%s vs Baseline: %.1f%% -> %.1f%% %s; %d shifted %s->%s [%.1f%%], %d shifted %s->%s [%.1f%%]; \u03C7\u00B2 = %.2f, %s) %s",
        dir_phrase, env, b1_pct, s1_pct, lvl1,
        n12, lvl1, lvl2, pct_1_to_2,
        n21, lvl2, lvl1, pct_2_to_1,
        m_test$statistic, format_p_stat(m_test$p.value), site_note
      )

      rows_out[[length(rows_out) + 1]] <- data.frame(
        Classification_Dimension   = dim_label,
        Class_1_Name               = lvl1,
        Class_2_Name               = lvl2,
        Baseline_Regime            = "NWPZ Baseline",
        Stress_Regime              = env,
        Contrast_Design_Status     = trans_scope,
        Baseline_Class1_n_Pct      = sprintf("%d (%.1f%%)", b1_n, b1_pct),
        Baseline_Class2_n_Pct      = sprintf("%d (%.1f%%)", b2_n, b2_pct),
        Stress_Class1_n_Pct        = sprintf("%d (%.1f%%)", s1_n, s1_pct),
        Stress_Class2_n_Pct        = sprintf("%d (%.1f%%)", s2_n, s2_pct),
        Retained_Class1_n          = n11,
        Shifted_Class1_to_Class2_n = n12,
        Shifted_Class2_to_Class1_n = n21,
        Retained_Class2_n          = n22,
        Net_Class1_Shift_PctPts    = round(net_shift, 1),
        McNemar_Chi_Square         = round(as.numeric(m_test$statistic), 2),
        df                         = as.integer(m_test$parameter),
        p_value                    = m_test$p.value,
        Exact_Binom_p              = exact_binom_p,
        Significance               = sig_lbl,
        Biological_Inference       = bio_infer,
        stringsAsFactors           = FALSE
      )
    }
    bind_rows(rows_out)
  }

  table_3_transitions <- bind_rows(
    build_unified_transition_rows(
      "Tillering Trajectory Shape (K-Means)",
      "Trajectory_Realized", "Earlier Proportional Peak", "Later Proportional Peak"
    ),
    build_unified_transition_rows(
      "Tillering Vigor (Mean Tillering Rate Early vs. Late)",
      "Vigor_Realized", VIGOR_HI, VIGOR_LO
    ),
    build_unified_transition_rows(
      "Tiller Conversion Efficiency (TCE = PT/T4)",
      "Realized_TCE", "High Conversion", "Low Conversion"
    )
  ) %>%
    mutate(
      p_Holm = p.adjust(p_value, method = "holm"),
      Significance_Holm = case_when(
        p_Holm < 0.001 ~ "***",
        p_Holm < 0.01  ~ "**",
        p_Holm < 0.05  ~ "*",
        TRUE           ~ "ns"
      )
    ) %>%
    relocate(p_Holm, Significance_Holm, .after = Significance)

  write.csv(
    table_3_transitions,
    file.path(sec22_dir, "Table_3_Phenological_Typology_and_Regime_Transitions.csv"),
    row.names = FALSE
  )

  # ----------------------------------------------------------------------------
  # 4.5 Consolidated Supplementary Table S5: Clean Genotype Assignments (N = 800)
  # ----------------------------------------------------------------------------
  table_s5_assignments <- blues_df %>%
    select(
      Genotype, Environment, Location, Sheet,
      T1_rel, T2_rel, T3_rel, T4_Count,
      AUTPC_Early, AUTPC_Late, AUTPC, AUTPC_norm,
      PT_Count, TCE,
      Trajectory_Baseline, Trajectory_Realized,
      Vigor_Baseline, Vigor_Realized,
      TCE_Baseline, Realized_TCE
    ) %>%
    mutate(
      across(c(T1_rel, T2_rel, T3_rel, TCE), ~ round(.x, 4)),
      across(c(T4_Count, PT_Count, AUTPC_norm), ~ round(.x, 2)),
      across(c(AUTPC_Early, AUTPC_Late, AUTPC), ~ round(.x, 1))
    )

  write.csv(
    table_s5_assignments,
    file.path(sec22_dir, "Table_S5_Genotype_Phenological_Class_Assignments_All_Regimes.csv"),
    row.names = FALSE
  )

  # Export Unified Multi-Sheet Excel Workbook for Section 2.2
  if (requireNamespace("openxlsx", quietly = TRUE)) {
    wb_sec22 <- openxlsx::createWorkbook()

    openxlsx::addWorksheet(wb_sec22, "Table_3_Regime_Transitions")
    openxlsx::writeData(wb_sec22, "Table_3_Regime_Transitions", table_3_transitions)
    openxlsx::addStyle(wb_sec22, "Table_3_Regime_Transitions", header_style, rows = 1, cols = 1:ncol(table_3_transitions), gridExpand = TRUE)
    openxlsx::setColWidths(wb_sec22, "Table_3_Regime_Transitions", cols = 1:ncol(table_3_transitions), widths = "auto")

    openxlsx::addWorksheet(wb_sec22, "Table_S4A_Cluster_Centroids")
    openxlsx::writeData(wb_sec22, "Table_S4A_Cluster_Centroids", table_s4_centroids)
    openxlsx::addStyle(wb_sec22, "Table_S4A_Cluster_Centroids", header_style, rows = 1, cols = 1:ncol(table_s4_centroids), gridExpand = TRUE)
    openxlsx::setColWidths(wb_sec22, "Table_S4A_Cluster_Centroids", cols = 1:ncol(table_s4_centroids), widths = "auto")

    openxlsx::addWorksheet(wb_sec22, "Table_S4B_Silhouette_Validity")
    openxlsx::writeData(wb_sec22, "Table_S4B_Silhouette_Validity", sil_df)
    openxlsx::addStyle(wb_sec22, "Table_S4B_Silhouette_Validity", header_style, rows = 1, cols = 1:ncol(sil_df), gridExpand = TRUE)
    openxlsx::setColWidths(wb_sec22, "Table_S4B_Silhouette_Validity", cols = 1:ncol(sil_df), widths = "auto")

    openxlsx::addWorksheet(wb_sec22, "Table_S5_Genotype_Assignments")
    openxlsx::writeData(wb_sec22, "Table_S5_Genotype_Assignments", table_s5_assignments)
    openxlsx::addStyle(wb_sec22, "Table_S5_Genotype_Assignments", header_style, rows = 1, cols = 1:ncol(table_s5_assignments), gridExpand = TRUE)
    openxlsx::setColWidths(wb_sec22, "Table_S5_Genotype_Assignments", cols = 1:ncol(table_s5_assignments), widths = "auto")

    openxlsx::saveWorkbook(wb_sec22, file.path(sec22_dir, "Section_2.2_Manuscript_and_Supplementary_Tables.xlsx"), overwrite = TRUE)
  }

  # ----------------------------------------------------------------------------
  # 4.6 Master 5-Panel (2-Row) Publication Figure:
  #     Top Row (Panels A-B): K-Means Trajectory Curve Centroids + PCA & Silhouette Inset
  #     Bottom Row (Panels C-E): Multi-Regime Alluvial Genotype Transitions
  # ----------------------------------------------------------------------------

  # --- Panel A: Proportional Tillering Trajectory Centroid Curves & Genotype Envelope ---
  traj_long_base <- base_full_df %>%
    mutate(
      E_rel  = COUNT_AT_EMERGENCE_MAIN / T4_Count,
      T4_rel = 1.0
    ) %>%
    select(Genotype, Trajectory = Trajectory_Realized, E_rel, T1_rel, T2_rel, T3_rel, T4_rel) %>%
    pivot_longer(
      cols = c(E_rel, T1_rel, T2_rel, T3_rel, T4_rel),
      names_to = "Stage_Key",
      values_to = "Prop_T4"
    ) %>%
    mutate(
      Stage_Idx = match(Stage_Key, c("E_rel", "T1_rel", "T2_rel", "T3_rel", "T4_rel")),
      Stage_Lab = factor(
        Stage_Key,
        levels = c("E_rel", "T1_rel", "T2_rel", "T3_rel", "T4_rel"),
        labels = c("Emergence\n(E / T4)", "Stage T1\n(T1 / T4)", "Stage T2\n(T2 / T4)", "Stage T3\n(T3 / T4)", "Stage T4\n(T4 / T4 = 1.0)")
      )
    )

  traj_curve_stats <- traj_long_base %>%
    group_by(Trajectory, Stage_Idx, Stage_Lab) %>%
    summarise(
      Mean_Prop = mean(Prop_T4, na.rm = TRUE),
      SD_Prop   = sd(Prop_T4, na.rm = TRUE),
      Q25_Prop  = quantile(Prop_T4, 0.25, na.rm = TRUE),
      Q75_Prop  = quantile(Prop_T4, 0.75, na.rm = TRUE),
      .groups   = "drop"
    )

  traj_pal <- c("Earlier Proportional Peak" = "#E64B35", "Later Proportional Peak" = "#3C5488")

  # Stage T1..T3 callout labels for Panel A centroids (offset horizontally from vertical SD bars)
  traj_callouts <- traj_curve_stats %>%
    filter(Stage_Idx %in% c(2, 3, 4)) %>%
    mutate(
      x_pos = ifelse(Trajectory == "Earlier Proportional Peak", Stage_Idx - 0.18, Stage_Idx + 0.18),
      y_pos = ifelse(Trajectory == "Earlier Proportional Peak", Mean_Prop + 0.055, Mean_Prop - 0.055),
      lbl   = sprintf("%.1f%%", 100 * Mean_Prop)
    )

  c_early_row <- centroid_table %>% filter(Descriptive_Name == "Earlier Proportional Peak")
  c_late_row  <- centroid_table %>% filter(Descriptive_Name == "Later Proportional Peak")

  p_panel_a <- ggplot() +
    geom_hline(yintercept = 1.0, linetype = "dashed", color = "grey35", linewidth = 0.65) +
    geom_line(
      data = traj_long_base,
      aes(x = Stage_Idx, y = Prop_T4, group = Genotype, color = Trajectory),
      alpha = 0.12, linewidth = 0.35, show.legend = FALSE
    ) +
    geom_ribbon(
      data = traj_curve_stats,
      aes(x = Stage_Idx, ymin = Q25_Prop, ymax = Q75_Prop, fill = Trajectory),
      alpha = 0.22, color = NA
    ) +
    geom_line(
      data = traj_curve_stats,
      aes(x = Stage_Idx, y = Mean_Prop, color = Trajectory, group = Trajectory),
      linewidth = 1.4
    ) +
    geom_errorbar(
      data = traj_curve_stats %>% filter(Stage_Idx %in% 2:4),
      aes(x = Stage_Idx, ymin = Mean_Prop - SD_Prop, ymax = Mean_Prop + SD_Prop, color = Trajectory),
      width = 0.10, linewidth = 0.65, show.legend = FALSE
    ) +
    geom_point(
      data = traj_curve_stats,
      aes(x = Stage_Idx, y = Mean_Prop, fill = Trajectory),
      shape = 21, size = 3.6, color = "black", stroke = 0.6
    ) +
    geom_label(
      data = traj_callouts,
      aes(x = x_pos, y = y_pos, label = lbl, color = Trajectory),
      size = 3.05, fontface = "bold", fill = "white", alpha = 0.92, label.size = 0.25, show.legend = FALSE
    ) +
    annotate(
      "label", x = 0.92, y = 1.38, hjust = 0, vjust = 1, size = 2.95, fontface = "bold",
      fill = "white", alpha = 0.92, label.size = 0.3, color = "grey15",
      label = sprintf(
        "Earlier Proportional Peak (n = %d, %.1f%%): T3 overshoot (108.0%%), senesces to T4 = %.2f tillers\nLater Proportional Peak (n = %d, %.1f%%): Sustained accrual (89.5%% -> 100%%), T4 = %.2f tillers",
        c_early_row$Size_n, c_early_row$Percentage, c_early_row$T4_Count_Mean,
        c_late_row$Size_n, c_late_row$Percentage, c_late_row$T4_Count_Mean
      )
    ) +
    scale_color_manual(values = traj_pal, name = "Baseline K-Means Trajectory Typology:") +
    scale_fill_manual(values = traj_pal, name = "Baseline K-Means Trajectory Typology:") +
    scale_x_continuous(
      breaks = 1:5,
      labels = levels(traj_long_base$Stage_Lab),
      expand = expansion(mult = c(0.05, 0.05))
    ) +
    scale_y_continuous(
      labels = scales::percent_format(accuracy = 1),
      limits = c(0, 1.42),
      breaks = seq(0, 1.4, by = 0.2),
      expand = c(0, 0)
    ) +
    pub_theme +
    theme(
      plot.title = element_text(face = "bold", size = 13.5),
      plot.subtitle = element_text(size = 9.8, color = "grey25", face = "bold", lineheight = 1.1),
      axis.title = element_text(face = "bold", size = 11.5),
      axis.text = element_text(face = "bold", size = 10.5, color = "black"),
      legend.title = element_text(face = "bold", size = 11.8),
      legend.text = element_text(face = "bold", size = 11.4),
      legend.key.size = unit(0.65, "cm")
    ) +
    labs(
      title = "A. Unsupervised Tillering Trajectory Typology in NWPZ Baseline (N = 200)",
      subtitle = "Genotype-normalized progression (T_i / T4; ribbons = IQR, error bars = \u00B11 SD)",
      x = "Vegetative Developmental Milestone (Normalized to Terminal T4 Count)",
      y = "Proportional Tiller Accumulation (T_i / T4, %)"
    )

  # --- Panel B: PCA Projection of Genotype-Normalized Shape + Embedded Silhouette Inset ---
  pca_res <- prcomp(traj_matrix)
  pca_var <- (pca_res$sdev^2) / sum(pca_res$sdev^2) * 100
  pc1_lab <- sprintf("PC1 (%.1f%% Variance Explained)", pca_var[1])
  pc2_lab <- sprintf("PC2 (%.1f%% Variance Explained)", pca_var[2])

  pca_df <- as.data.frame(pca_res$x)
  pca_df$Trajectory <- factor(baseline_blues$Trajectory, levels = c("Earlier Proportional Peak", "Later Proportional Peak"))

  pca_centroids <- pca_df %>%
    group_by(Trajectory) %>%
    summarise(
      m_PC1 = mean(PC1), se_PC1 = sd(PC1) / sqrt(n()), ci_PC1 = qt(0.975, df = n() - 1) * se_PC1,
      m_PC2 = mean(PC2), se_PC2 = sd(PC2) / sqrt(n()), ci_PC2 = qt(0.975, df = n() - 1) * se_PC2,
      .groups = "drop"
    )

  # Biplot loading vectors for T1/T4, T2/T4, T3/T4
  load_df <- as.data.frame(pca_res$rotation[, 1:2])
  load_df$Variable <- c("T1 / T4", "T2 / T4", "T3 / T4")
  arrow_scale <- 2.35
  load_df <- load_df %>%
    mutate(
      xend = PC1 * arrow_scale,
      yend = PC2 * arrow_scale,
      xlab = PC1 * (arrow_scale + 0.38),
      ylab = PC2 * (arrow_scale + 0.30)
    )

  p_pca_main <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Trajectory, fill = Trajectory)) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "grey55", linewidth = 0.5) +
    geom_vline(xintercept = 0, linetype = "dotted", color = "grey55", linewidth = 0.5) +
    stat_ellipse(geom = "polygon", alpha = 0.16, level = 0.90, show.legend = FALSE) +
    geom_point(size = 2.5, alpha = 0.75, shape = 21, stroke = 0.4, color = "black") +
    geom_segment(
      data = load_df,
      aes(x = 0, y = 0, xend = xend, yend = yend),
      inherit.aes = FALSE,
      arrow = arrow(length = unit(0.22, "cm"), type = "closed"),
      color = "grey15", linewidth = 0.80
    ) +
    geom_label(
      data = load_df,
      aes(x = xlab, y = ylab, label = Variable),
      inherit.aes = FALSE,
      size = 3.15, fontface = "bold", fill = "white", alpha = 0.90, label.size = 0.25, color = "grey10"
    ) +
    geom_errorbar(
      data = pca_centroids,
      aes(x = m_PC1, ymin = m_PC2 - ci_PC2, ymax = m_PC2 + ci_PC2),
      inherit.aes = FALSE, color = "black", width = 0.12, linewidth = 0.85
    ) +
    geom_errorbar(
      data = pca_centroids,
      aes(y = m_PC2, xmin = m_PC1 - ci_PC1, xmax = m_PC1 + ci_PC1),
      inherit.aes = FALSE, color = "black", width = 0.12, linewidth = 0.85, orientation = "y"
    ) +
    geom_point(
      data = pca_centroids,
      aes(x = m_PC1, y = m_PC2, fill = Trajectory),
      inherit.aes = FALSE, shape = 23, size = 4.6, stroke = 1.0, color = "black", show.legend = FALSE
    ) +
    scale_color_manual(values = traj_pal, name = "Baseline K-Means Trajectory Typology:") +
    scale_fill_manual(values = traj_pal, name = "Baseline K-Means Trajectory Typology:") +
    scale_x_continuous(expand = expansion(mult = c(0.08, 0.34))) +
    scale_y_continuous(expand = expansion(mult = c(0.10, 0.12))) +
    pub_theme +
    theme(
      plot.title = element_text(face = "bold", size = 13.5),
      plot.subtitle = element_text(size = 9.8, color = "grey25", face = "bold", lineheight = 1.1),
      axis.title = element_text(face = "bold", size = 11.5),
      axis.text = element_text(face = "bold", size = 10.5, color = "black"),
      legend.title = element_text(face = "bold", size = 11.8),
      legend.text = element_text(face = "bold", size = 11.4),
      legend.key.size = unit(0.65, "cm")
    ) +
    labs(
      title = "B. PCA Biplot of Trajectory Geometry & Silhouette Cluster Optimization",
      subtitle = sprintf(
        "PC1 + PC2 explain %.1f%% variance (90%% ellipses & \u00B195%% CI centroids) | Silhouette k = %d (width = %.3f)",
        pca_var[1] + pca_var[2], optimal_k, max_sil_width
      ),
      x = pc1_lab,
      y = pc2_lab
    )

  # Compact Silhouette Optimization Inset (k = 2..10) embedded cleanly in upper-right of Panel B
  p_sil_inset <- ggplot(sil_df, aes(x = k, y = Silhouette)) +
    geom_vline(xintercept = optimal_k, linetype = "dashed", color = "#DC0000", linewidth = 0.60) +
    geom_line(color = "#3C5488", linewidth = 0.90) +
    geom_point(aes(color = k == optimal_k, size = k == optimal_k), show.legend = FALSE) +
    scale_color_manual(values = c("TRUE" = "#DC0000", "FALSE" = "#3C5488")) +
    scale_size_manual(values = c("TRUE" = 2.8, "FALSE" = 1.8)) +
    annotate(
      "text", x = optimal_k + 0.35, y = max_sil_width - 0.004,
      label = sprintf("k = %d (%.3f)", optimal_k, max_sil_width),
      hjust = 0, size = 2.95, fontface = "bold", color = "#DC0000"
    ) +
    scale_x_continuous(breaks = seq(2, 10, by = 2)) +
    scale_y_continuous(limits = c(0.22, 0.35), breaks = c(0.24, 0.28, 0.32)) +
    theme_bw(base_size = 9.0) +
    theme(
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", color = "grey30", linewidth = 0.4),
      plot.background  = element_rect(fill = alpha("white", 0.94), color = "grey40", linewidth = 0.4),
      plot.title       = element_text(face = "bold", size = 8.8, hjust = 0.5, margin = margin(b = 2)),
      axis.title       = element_text(face = "bold", size = 8.3),
      axis.text        = element_text(color = "black", size = 7.8, face = "bold"),
      plot.margin      = margin(3, 4, 3, 4)
    ) +
    labs(title = "Silhouette Validity (k = 2..10)", x = "Clusters (k)", y = "Mean Width")

  p_panel_b <- p_pca_main +
    patchwork::inset_element(
      p_sil_inset,
      left = 0.60, bottom = 0.54, right = 0.99, top = 0.98,
      align_to = "panel"
    )

  # --- Bottom Row (Panels C, D, E): Alluvial Genotype Flow & Stacked Proportion Panels ---
  build_alluvial_transition_panel <- function(df_in, class_col, top_lvl, bottom_lvl, pal_vals, panel_title, sub_text, legend_title) {
    envs <- levels(df_in$Env_Factor)
    w_bar <- 0.44
    n_tot <- dplyr::n_distinct(df_in$Genotype)

    w_df <- df_in %>%
      select(Genotype, Env_Factor, Class = all_of(class_col)) %>%
      pivot_wider(names_from = Env_Factor, values_from = Class)

    bar_rows <- list()
    p_bot_vec <- numeric(length(envs))
    names(p_bot_vec) <- envs
    for (k in seq_along(envs)) {
      e_k <- envs[k]
      n_bot <- sum(w_df[[e_k]] == bottom_lvl, na.rm = TRUE)
      n_top <- sum(w_df[[e_k]] == top_lvl, na.rm = TRUE)
      p_bot <- 100 * n_bot / n_tot
      p_top <- 100 * n_top / n_tot
      p_bot_vec[e_k] <- p_bot

      bar_rows[[length(bar_rows) + 1]] <- data.frame(
        x_idx = k, Env = e_k, Class = bottom_lvl,
        xmin = k - w_bar / 2, xmax = k + w_bar / 2,
        ymin = 0, ymax = p_bot, ymid = p_bot / 2, pct = p_bot,
        stringsAsFactors = FALSE
      )
      bar_rows[[length(bar_rows) + 1]] <- data.frame(
        x_idx = k, Env = e_k, Class = top_lvl,
        xmin = k - w_bar / 2, xmax = k + w_bar / 2,
        ymin = p_bot, ymax = 100, ymid = p_bot + p_top / 2, pct = p_top,
        stringsAsFactors = FALSE
      )
    }
    bars_df <- bind_rows(bar_rows) %>%
      mutate(Class = factor(Class, levels = c(top_lvl, bottom_lvl)))

    rib_polys <- list()
    t_seq <- seq(0, 1, length.out = 40)
    s_seq <- 3 * t_seq^2 - 2 * t_seq^3

    for (k in 1:(length(envs) - 1)) {
      e1 <- envs[k]
      e2 <- envs[k + 1]
      x1 <- k + w_bar / 2
      x2 <- (k + 1) - w_bar / 2
      x_vals <- x1 + (x2 - x1) * t_seq

      f_bb <- 100 * sum(w_df[[e1]] == bottom_lvl & w_df[[e2]] == bottom_lvl, na.rm = TRUE) / n_tot
      f_bt <- 100 * sum(w_df[[e1]] == bottom_lvl & w_df[[e2]] == top_lvl,    na.rm = TRUE) / n_tot
      f_tb <- 100 * sum(w_df[[e1]] == top_lvl    & w_df[[e2]] == bottom_lvl, na.rm = TRUE) / n_tot
      f_tt <- 100 * sum(w_df[[e1]] == top_lvl    & w_df[[e2]] == top_lvl,    na.rm = TRUE) / n_tot

      pb1 <- p_bot_vec[e1]
      pb2 <- p_bot_vec[e2]

      streams <- list(
        list(id = paste0(k, "_bb"), cls = bottom_lvl, y1_lo = 0,          y1_hi = f_bb,       y2_lo = 0,          y2_hi = f_bb,       w = f_bb),
        list(id = paste0(k, "_tb"), cls = top_lvl,    y1_lo = pb1,        y1_hi = pb1 + f_tb, y2_lo = f_bb,       y2_hi = pb2,        w = f_tb),
        list(id = paste0(k, "_bt"), cls = bottom_lvl, y1_lo = f_bb,       y1_hi = pb1,        y2_lo = pb2,        y2_hi = pb2 + f_bt, w = f_bt),
        list(id = paste0(k, "_tt"), cls = top_lvl,    y1_lo = pb1 + f_tb, y1_hi = 100,        y2_lo = pb2 + f_bt, y2_hi = 100,        w = f_tt)
      )
      streams <- streams[order(sapply(streams, function(s) s$w), decreasing = TRUE)]

      for (st in streams) {
        if ((st$y1_hi - st$y1_lo) > 0.01) {
          y_lo <- st$y1_lo + (st$y2_lo - st$y1_lo) * s_seq
          y_hi <- st$y1_hi + (st$y2_hi - st$y1_hi) * s_seq
          rib_polys[[length(rib_polys) + 1]] <- data.frame(
            poly_id = factor(st$id, levels = unique(c(sapply(rib_polys, function(df) as.character(df$poly_id[1])), st$id))),
            Class   = st$cls,
            x       = c(x_vals, rev(x_vals)),
            y       = c(y_lo, rev(y_hi)),
            stringsAsFactors = FALSE
          )
        }
      }
    }
    ribs_df <- bind_rows(rib_polys) %>%
      mutate(
        poly_id = factor(poly_id, levels = unique(poly_id)),
        Class   = factor(Class, levels = c(top_lvl, bottom_lvl))
      )

    ggplot() +
      geom_polygon(data = ribs_df, aes(x = x, y = y, group = poly_id, fill = Class), alpha = 0.30, color = "white", linewidth = 0.25) +
      geom_vline(xintercept = 2.5, linetype = "dashed", color = "grey30", linewidth = 0.65) +
      geom_rect(data = bars_df, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = Class),
                color = "black", linewidth = 0.45, alpha = 0.92) +
      geom_text(
        data = bars_df %>% filter(pct >= 6.0),
        aes(x = x_idx, y = ymid, label = sprintf("%.1f%%", pct)),
        color = "white", fontface = "bold", size = 3.95
      ) +
      scale_fill_manual(values = pal_vals, name = legend_title) +
      scale_x_continuous(breaks = 1:4, labels = envs, expand = expansion(mult = c(0.05, 0.05))) +
      scale_y_continuous(labels = function(x) paste0(x, "%"), limits = c(0, 103), expand = c(0, 0)) +
      pub_theme +
      theme(
        plot.title = element_text(face = "bold", size = 13.0),
        axis.text.x = element_text(angle = 18, hjust = 1, face = "bold", size = 10.8, color = "black"),
        axis.text.y = element_text(face = "bold", size = 10.4, color = "black"),
        axis.title.y = element_text(face = "bold", size = 11.4, color = "black"),
        plot.subtitle = element_text(size = 9.2, color = "grey20", face = "bold", lineheight = 1.1),
        legend.title = element_text(face = "bold", size = 11.5),
        legend.text = element_text(face = "bold", size = 11.0),
        legend.key.size = unit(0.65, "cm"),
        plot.margin = margin(t = 6, r = 8, b = 4, l = 8)
      ) +
      labs(title = panel_title, subtitle = sub_text, x = "", y = "Genotype Composition (%)")
  }

  trj_rows <- table_3_transitions %>% filter(Classification_Dimension == "Tillering Trajectory Shape (K-Means)")
  vig_rows <- table_3_transitions %>% filter(Classification_Dimension == "Tillering Vigor (Mean Tillering Rate Early vs. Late)")
  tce_rows <- table_3_transitions %>% filter(Classification_Dimension == "Tiller Conversion Efficiency (TCE = PT/T4)")

  sub_trj_alluv <- sprintf(
    "McNemar vs. Baseline: NWPZ Warm \u03C7\u00B2=%.1f%s | PZ Warm I \u03C7\u00B2=%.1f%s | PZ Warm II \u03C7\u00B2=%.1f%s",
    trj_rows$McNemar_Chi_Square[1], trj_rows$Significance[1],
    trj_rows$McNemar_Chi_Square[2], trj_rows$Significance[2],
    trj_rows$McNemar_Chi_Square[3], trj_rows$Significance[3]
  )
  # Compute within-PZ McNemar (PZ Warm I -> PZ Warm II) to highlight unconfounded within-site thermal contrast
  w_pz_vig <- blues_df %>%
    filter(Environment %in% c("PZ Warm I", "PZ Warm II")) %>%
    select(Genotype, Environment, Vigor_Realized) %>%
    pivot_wider(names_from = Environment, values_from = Vigor_Realized)
  mc_pz_vig <- mcnemar.test(table(w_pz_vig$`PZ Warm I`, w_pz_vig$`PZ Warm II`), correct = TRUE)
  sig_pz_vig <- if (mc_pz_vig$p.value < 0.001) "***" else if (mc_pz_vig$p.value < 0.01) "**" else if (mc_pz_vig$p.value < 0.05) "*" else "ns"

  sub_vig_alluv <- sprintf(
    "Within-Site McNemar: NWPZ Warm vs. Base \u03C7\u00B2=%.1f (%s) | PZ Warm II vs. I \u03C7\u00B2=%.1f (%s)",
    vig_rows$McNemar_Chi_Square[1], vig_rows$Significance[1],
    as.numeric(mc_pz_vig$statistic), sig_pz_vig
  )
  sub_tce_alluv <- sprintf(
    "McNemar vs. Base: NWPZ Warm \u03C7\u00B2=%.1f%s | PZ Warm I \u03C7\u00B2=%.1f%s | PZ Warm II \u03C7\u00B2=%.1f%s",
    tce_rows$McNemar_Chi_Square[1], tce_rows$Significance[1],
    tce_rows$McNemar_Chi_Square[2], tce_rows$Significance[2],
    tce_rows$McNemar_Chi_Square[3], tce_rows$Significance[3]
  )

  p_panel_c <- build_alluvial_transition_panel(
    blues_df, "Trajectory_Realized", "Earlier Proportional Peak", "Later Proportional Peak",
    c("Earlier Proportional Peak" = "#E64B35", "Later Proportional Peak" = "#3C5488"),
    "C. Tillering Trajectory Shifts (K-Means)", sub_trj_alluv, "Trajectory Class:"
  )

  p_panel_d <- build_alluvial_transition_panel(
    blues_df, "Vigor_Realized", VIGOR_HI, VIGOR_LO,
    setNames(c("#E64B35", "#00A087"), c(VIGOR_HI, VIGOR_LO)),
    "D. Tillering Vigor Shifts (K-Means Rate)", sub_vig_alluv, "Vigor Class:"
  )

  p_panel_e <- build_alluvial_transition_panel(
    blues_df, "Realized_TCE", "High Conversion", "Low Conversion",
    c("High Conversion" = "#3C5488", "Low Conversion" = "#E64B35"),
    "E. Tiller Conversion Efficiency (TCE)", sub_tce_alluv, "TCE Class:"
  )

  top_row_sec22 <- (p_panel_a | p_panel_b) +
    plot_layout(widths = c(1.12, 0.98), guides = "collect") &
    theme(legend.position = "bottom")
  bot_row_sec22 <- (p_panel_c | p_panel_d | p_panel_e)

  compound_transitions <- (top_row_sec22 / bot_row_sec22) +
    plot_layout(heights = c(1.08, 1.0)) +
    plot_annotation(
      caption = paste0(
        "Panels A\u2013B define the unsupervised 3D K-Means tillering trajectory typology in the NWPZ Baseline reference (N = 200). Panels C\u2013E show realized genotype class proportions (stacked bars) and transitions (Hermite alluvial ribbons, N = 200).\n",
        "Note on Cluster Validation: Silhouette width optimization applies to multivariate K-Means typologies\u2014Trajectory Geometry (3D: T1/T4, T2/T4, T3/T4; k = 2 silhouette = 0.331, Panel B inset) and Tillering Vigor (2D early vs. late rates;\n",
        "k = 2 silhouette = 0.465, Table S4)\u2014whereas TCE (PT/T4) is a 1D ratio partitioned at the NWPZ Baseline median (79.0%)."
      ),
      theme = theme(plot.caption = element_text(size = 9.6, color = "grey25", face = "italic", hjust = 0, lineheight = 1.15))
    )

  ggsave(
    file.path(sec22_dir, "Compound_Phenological_Transitions_Vigor_Trajectory_TCE.png"),
    plot = compound_transitions, width = 18.5, height = 11.6, dpi = 600
  )

  # Mirror consolidated assets to Curiosity_And_Iterations/Phenology_Classifications_Review
  curiosity_sec22_dir <- file.path(project_root, "Curiosity_And_Iterations", "Phenology_Classifications_Review")
  dir.create(curiosity_sec22_dir, showWarnings = FALSE, recursive = TRUE)
  for (fn in c(
    "Compound_Phenological_Transitions_Vigor_Trajectory_TCE.png",
    "Table_3_Phenological_Typology_and_Regime_Transitions.csv",
    "Table_S4_KMeans_Cluster_Centroids_and_Silhouette_Validity.csv",
    "Table_S5_Genotype_Phenological_Class_Assignments_All_Regimes.csv",
    "Section_2.2_Manuscript_and_Supplementary_Tables.xlsx"
  )) {
    src_f <- file.path(sec22_dir, fn)
    if (file.exists(src_f)) {
      file.copy(src_f, file.path(curiosity_sec22_dir, fn), overwrite = TRUE)
    }
  }
} else {
  cat("  -> Section 4 [SKIPPED - Assets Finalized & Protected]\n")
}



# ==============================================================================
# SECTION 5: MICRO-PHENOLOGY INDICES (AUTPC & TILLER STABILITY INDEX)
# ==============================================================================
cat(">>> Executing Section 5: Micro-Phenology Indices (AUTPC & TSI)...\n")

# 5.1 Calculate Tiller Stability Index (TSI_PZ_WarmII) & Within-Delhi Randomized TSI_Delhi
opt_counts      <- blues_df %>% filter(Environment == "NWPZ Baseline") %>% select(Genotype, PT_Opt = PT_Count)
stress_counts   <- blues_df %>% filter(Environment == "PZ Warm II")    %>% select(Genotype, PT_Stress = PT_Count)
del_warm_counts <- blues_df %>% filter(Environment == "NWPZ Warm")     %>% select(Genotype, PT_DelWarm = PT_Count)

pop_PT_Warm_II <- mean(stress_counts$PT_Stress, na.rm = TRUE)
pop_PT_DelWarm <- mean(del_warm_counts$PT_DelWarm, na.rm = TRUE)

tsi_df <- opt_counts %>%
  inner_join(stress_counts, by = "Genotype") %>%
  inner_join(del_warm_counts, by = "Genotype") %>%
  mutate(
    TSI_PZ_WarmII = if_else(PT_Opt <= 0 | pop_PT_Warm_II <= 0 | is.na(PT_Opt) | is.na(pop_PT_Warm_II),
                            NA_real_,
                            (PT_Stress / PT_Opt) * (PT_Stress / pop_PT_Warm_II)),
    Stability = if_else(TSI_PZ_WarmII >= 1.0, "Stable", "Contrasting"),
    Retention_Ratio = if_else(PT_Opt <= 0 | is.na(PT_Opt), NA_real_, PT_Stress / PT_Opt),
    Retention_Resid = as.numeric(residuals(lm(PT_Stress ~ PT_Opt))),
    TSI_Delhi = if_else(PT_Opt <= 0 | pop_PT_DelWarm <= 0 | is.na(PT_Opt) | is.na(pop_PT_DelWarm),
                        NA_real_,
                        (PT_DelWarm / PT_Opt) * (PT_DelWarm / pop_PT_DelWarm)),
    Stability_Delhi = if_else(TSI_Delhi >= 1.0, "Stable", "Contrasting"),
    Retention_Resid_Delhi = as.numeric(residuals(lm(PT_DelWarm ~ PT_Opt)))
  )

# Strip pre-existing TSI columns for idempotency before joining
blues_df <- blues_df %>%
  select(-any_of(c("TSI", "TSI_PZ_WarmII", "Stability", "TSI_Delhi", "Stability_Delhi",
                   "Retention_Resid", "Retention_Resid_Delhi"))) %>%
  left_join(
    tsi_df %>% select(Genotype, TSI_PZ_WarmII, Stability, TSI_Delhi, Stability_Delhi,
                      Retention_Resid, Retention_Resid_Delhi),
    by = "Genotype"
  )

if (isTRUE(RUN_SECTION_5_MICROPHENOLOGY)) {
  sec23_dir <- file.path(assets_dir, "Section_2.3_Phenology_Indices")
  dir.create(sec23_dir, showWarnings = FALSE, recursive = TRUE)

  # Archive duplicate Section_2.3_Micro_Phenology folder and all legacy/redundant Section 2.3 files
  archive_sec23_dir <- file.path(project_root, "Archive_Rejected", "Section_2.3_Legacy_Micro_Phenology")
  dir.create(archive_sec23_dir, showWarnings = FALSE, recursive = TRUE)

  sec23_micro_dir <- file.path(assets_dir, "Section_2.3_Micro_Phenology")
  if (dir.exists(sec23_micro_dir)) {
    micro_files <- list.files(sec23_micro_dir, full.names = TRUE)
    for (mf in micro_files) {
      file.copy(mf, file.path(archive_sec23_dir, paste0("MicroPhenologyDup_", basename(mf))), overwrite = TRUE)
      unlink(mf)
    }
    unlink(sec23_micro_dir, recursive = TRUE)
  }

  legacy_sec23_patterns <- c(
    "^AUTPC_Component_Yield_Correlations\\.csv$",
    "^AUTPC_Decomposition_Heatmap\\.png$",
    "^Table_Anchor_Sensitivity_AUTPC_norm\\.csv$",
    "^Table_TSI_MultiStress_Transitions\\.(csv|xlsx)$",
    "^Table_TSI_Transition_NWPZ_Warm_vs_PZ_Acute\\.csv$",
    "^Table_TSI_Yield_Validation_With_Position_Adjustment\\.csv$",
    "^Table_Yield_Component_Path_Regime_Shift\\.csv$",
    "^TSI_Direct_Yield_Validation\\.png$",
    "^Yield_Component_Source_Sink_Regime_Shift\\.png$"
  )
  for (pat in legacy_sec23_patterns) {
    matched_f <- list.files(sec23_dir, pattern = pat, full.names = TRUE)
    for (mf in matched_f) {
      if (file.exists(mf)) {
        file.copy(mf, file.path(archive_sec23_dir, basename(mf)), overwrite = TRUE)
        unlink(mf)
      }
    }
  }
  old_tsi_comp <- file.path(sec23_dir, "TSI_Compound_Yield_Validation.png")
  if (file.exists(old_tsi_comp)) {
    file.copy(old_tsi_comp, file.path(archive_sec23_dir, "TSI_Compound_Yield_Validation_SingleScatter_Prior.png"), overwrite = TRUE)
  }

  # ----------------------------------------------------------------------------
  # 5.1b Multi-Stress Progression of TSI Across NWPZ Warm, PZ Warm I, and PZ Warm II
  # ----------------------------------------------------------------------------
  tsi_mult_list <- list()
  for (env in c("NWPZ Warm", "PZ Warm I", "PZ Warm II")) {
    str_sub <- blues_df %>%
      filter(Environment == env) %>%
      select(Genotype, PT_Stress = PT_Count, Yield_Stress = GrainYld, Field_Pos_Block = Field_Pos_Block_Diagnostic)
    pop_m <- mean(str_sub$PT_Stress, na.rm = TRUE)
    env_tsi <- inner_join(opt_counts, str_sub, by = "Genotype") %>%
      mutate(
        Environment = env,
        Pop_Mean_PT_Stress = pop_m,
        Retention_Ratio = if_else(PT_Opt <= 0 | is.na(PT_Opt) | is.na(PT_Stress), NA_real_, PT_Stress / PT_Opt),
        Relative_Stress_Performance = if_else(is.na(PT_Stress) | pop_m <= 0 | is.na(pop_m), NA_real_, PT_Stress / pop_m),
        TSI = if_else(
          PT_Opt <= 0 | pop_m <= 0 | is.na(PT_Opt) | is.na(PT_Stress) | is.na(pop_m),
          NA_real_,
          (PT_Stress / PT_Opt) * (PT_Stress / pop_m)
        ),
        Stability_Class = if_else(TSI >= 1.0, "Stable (\u2265 1.0)", "Contrasting (< 1.0)")
      )
    tsi_mult_list[[env]] <- env_tsi
  }
  tsi_multi_df <- bind_rows(tsi_mult_list) %>%
    mutate(
      Env_Factor = factor(Environment, levels = c("NWPZ Warm", "PZ Warm I", "PZ Warm II")),
      Stability_Class = factor(Stability_Class, levels = c("Stable (\u2265 1.0)", "Contrasting (< 1.0)"))
    )

  # Build Table S6A: Multi-Stress TSI Classification, Yield Validation & McNemar Transitions
  tsi_w <- tsi_multi_df %>%
    select(Genotype, Environment, Stability_Class) %>%
    pivot_wider(names_from = Environment, values_from = Stability_Class)

  tsi_transition_pairs <- list(
    c("NWPZ Warm", "PZ Warm I"),
    c("PZ Warm I", "PZ Warm II"),
    c("NWPZ Warm", "PZ Warm II")
  )
  table_s6a_rows <- list()
  for (pr in tsi_transition_pairs) {
    e1 <- pr[1]; e2 <- pr[2]
    tb_m <- table(
      Regime_1 = factor(tsi_w[[e1]], levels = c("Stable (\u2265 1.0)", "Contrasting (< 1.0)")),
      Regime_2 = factor(tsi_w[[e2]], levels = c("Stable (\u2265 1.0)", "Contrasting (< 1.0)"))
    )
    m_tst <- mcnemar.test(tb_m, correct = TRUE)
    n11 <- as.integer(tb_m["Stable (\u2265 1.0)", "Stable (\u2265 1.0)"])
    n12 <- as.integer(tb_m["Stable (\u2265 1.0)", "Contrasting (< 1.0)"])
    n21 <- as.integer(tb_m["Contrasting (< 1.0)", "Stable (\u2265 1.0)"])
    n22 <- as.integer(tb_m["Contrasting (< 1.0)", "Contrasting (< 1.0)"])
    n_tot <- sum(tb_m)

    # Yield correlation in Regime 2
    sub_e2 <- tsi_multi_df %>% filter(Environment == e2)
    ct_e2  <- cor.test(sub_e2$TSI, sub_e2$Yield_Stress)
    if (!all(is.na(sub_e2$Field_Pos_Block))) {
      rx <- residuals(lm(TSI ~ factor(Field_Pos_Block), data = sub_e2))
      ry <- residuals(lm(Yield_Stress ~ factor(Field_Pos_Block), data = sub_e2))
      r_adj_e2 <- cor(rx, ry)
      df_adj_e2 <- nrow(sub_e2) - length(unique(na.omit(sub_e2$Field_Pos_Block))) - 1L
      t_adj_e2  <- r_adj_e2 * sqrt(df_adj_e2 / (1 - r_adj_e2^2))
      p_adj_e2  <- 2 * pt(-abs(t_adj_e2), df = df_adj_e2)
    } else {
      r_adj_e2 <- NA_real_
      p_adj_e2 <- NA_real_
    }

    sig_m <- if (m_tst$p.value < 0.001) "***" else if (m_tst$p.value < 0.01) "**" else if (m_tst$p.value < 0.05) "*" else "ns"

    table_s6a_rows[[paste0(e1, "_vs_", e2)]] <- data.frame(
      Reference_Stress_Regime       = e1,
      Target_Stress_Regime          = e2,
      Ref_Stable_n_Pct              = sprintf("%d (%.1f%%)", n11 + n12, 100 * (n11 + n12) / n_tot),
      Target_Stable_n_Pct           = sprintf("%d (%.1f%%)", n11 + n21, 100 * (n11 + n21) / n_tot),
      Retained_Stable_n             = n11,
      Collapsed_Stable_to_Contrast  = n12,
      Compensated_Contrast_to_Stab  = n21,
      Retained_Contrasting_n        = n22,
      Turnover_Genotypes_n_Pct      = sprintf("%d (%.1f%%)", n12 + n21, 100 * (n12 + n21) / n_tot),
      McNemar_Chi_Square            = round(as.numeric(m_tst$statistic), 2),
      df                            = as.integer(m_tst$parameter),
      McNemar_p_value               = m_tst$p.value,
      McNemar_Significance          = sig_m,
      Target_TSI_Yield_r_Unadj      = round(as.numeric(ct_e2$estimate), 3),
      Target_TSI_Yield_p_Unadj      = ct_e2$p.value,
      Target_TSI_Yield_r_PosAdj     = ifelse(is.na(r_adj_e2), "\u2014", as.character(round(r_adj_e2, 3))),
      Target_TSI_Yield_p_PosAdj     = ifelse(is.na(p_adj_e2), "\u2014", as.character(signif(p_adj_e2, 4))),
      Biological_Inference          = sprintf(
        "%s -> %s: Aggregate Stable proportion shifts %.1f%% -> %.1f%% (\u03C7\u00B2 = %.2f, %s), with %d genotypes re-ranking (%d collapse, %d compensate).",
        e1, e2, 100 * (n11 + n12) / n_tot, 100 * (n11 + n21) / n_tot, m_tst$statistic, sig_m, n12 + n21, n12, n21
      ),
      stringsAsFactors = FALSE
    )
  }
  table_s6a_transitions <- bind_rows(table_s6a_rows)

  # Build Table S6B: Genotype-Level Multi-Stress TSI Profiles (N = 200)
  table_s6b_genotypes <- opt_counts %>%
    left_join(
      tsi_multi_df %>% filter(Environment == "NWPZ Warm") %>%
        select(Genotype, PT_NWPZ_Warm = PT_Stress, Yield_NWPZ_Warm = Yield_Stress,
               TSI_NWPZ_Warm = TSI, Class_NWPZ_Warm = Stability_Class),
      by = "Genotype"
    ) %>%
    left_join(
      tsi_multi_df %>% filter(Environment == "PZ Warm I") %>%
        select(Genotype, PT_PZ_Warm_I = PT_Stress, Yield_PZ_Warm_I = Yield_Stress,
               TSI_PZ_Warm_I = TSI, Class_PZ_Warm_I = Stability_Class),
      by = "Genotype"
    ) %>%
    left_join(
      tsi_multi_df %>% filter(Environment == "PZ Warm II") %>%
        select(Genotype, PT_PZ_Warm_II = PT_Stress, Yield_PZ_Warm_II = Yield_Stress,
               Retention_Ratio_PZ_Warm_II = Retention_Ratio, TSI_PZ_Warm_II = TSI, Class_PZ_Warm_II = Stability_Class),
      by = "Genotype"
    ) %>%
    mutate(
      across(c(PT_Opt, PT_NWPZ_Warm, PT_PZ_Warm_I, PT_PZ_Warm_II), ~ round(.x, 2)),
      across(c(Yield_NWPZ_Warm, Yield_PZ_Warm_I, Yield_PZ_Warm_II, Retention_Ratio_PZ_Warm_II,
               TSI_NWPZ_Warm, TSI_PZ_Warm_I, TSI_PZ_Warm_II), ~ round(.x, 3))
    )

  write.csv(
    table_s6a_transitions,
    file.path(sec23_dir, "Table_S6_TSI_MultiStress_Transitions_and_Validation.csv"),
    row.names = FALSE
  )

  # ----------------------------------------------------------------------------
  # 5.1c Out-of-Sample Cross-Validation & Meng-Rosenthal-Rubin Delta-r Benchmark
  #      (Item 8: Demonstrate that TSI is a descriptive within-regime retention index
  #       with |r| <= 0.11 out-of-sample, and test residual-based retention)
  # ----------------------------------------------------------------------------
  meng_z_test <- function(y, x1, x2) {
    n <- length(y)
    r1 <- cor(y, x1); r2 <- cor(y, x2); r12 <- cor(x1, x2)
    z1 <- atanh(r1); z2 <- atanh(r2)
    rm2 <- (r1^2 + r2^2) / 2
    f <- min((1 - r12) / (2 * (1 - rm2)), 1)
    h <- (1 - f * rm2) / (1 - rm2)
    z_stat <- (z1 - z2) * sqrt((n - 3) / (2 * (1 - r12) * h))
    p_val  <- 2 * pnorm(-abs(z_stat))
    set.seed(42)
    b_diffs <- replicate(2000, {
      idx <- sample.int(n, n, replace = TRUE)
      cor(y[idx], x1[idx]) - cor(y[idx], x2[idx])
    })
    ci_95 <- quantile(b_diffs, c(0.025, 0.975), na.rm = TRUE)
    list(r_tsi = r1, r_pt = r2, r_tsi_pt = r12, delta_r = r1 - r2,
         ci_lo = as.numeric(ci_95[1]), ci_hi = as.numeric(ci_95[2]),
         meng_z = z_stat, meng_p = p_val)
  }

  yld_by_env <- blues_df %>%
    select(Genotype, Environment, GrainYld, PT_Count) %>%
    pivot_wider(names_from = Environment, values_from = c(GrainYld, PT_Count))

  oos_df <- tsi_df %>%
    left_join(yld_by_env, by = "Genotype")

  oos_rows <- list()
  add_oos_row <- function(pred_name, train_env, test_env, val_type, x_vec, y_vec, pt_bench_vec = NULL) {
    ct <- cor.test(x_vec, y_vec)
    if (!is.null(pt_bench_vec)) {
      mz <- meng_z_test(y_vec, x_vec, pt_bench_vec)
      dr_str <- sprintf("Delta_r = %+.4f [%.4f to %.4f] (z = %.2f, p = %.4f)", mz$delta_r, mz$ci_lo, mz$ci_hi, mz$meng_z, mz$meng_p)
    } else {
      dr_str <- "Not applicable (cross-regime)"
    }
    oos_rows[[length(oos_rows) + 1]] <<- data.frame(
      Predictor = pred_name,
      Trained_Regimes = train_env,
      Evaluated_Yield_Regime = test_env,
      Validation_Scope = val_type,
      Pearson_r = round(as.numeric(ct$estimate), 4),
      CI_95_Lower = round(ct$conf.int[1], 4),
      CI_95_Upper = round(ct$conf.int[2], 4),
      t_statistic = round(as.numeric(ct$statistic), 3),
      p_value = signif(ct$p.value, 4),
      Delta_r_vs_SameRegime_PT = dr_str,
      stringsAsFactors = FALSE
    )
  }

  # Out-of-sample TSI_PZ_WarmII
  add_oos_row("TSI_PZ_WarmII", "NWPZ Baseline + PZ Warm II", "NWPZ Warm", "Out-of-Sample Cross-Regime",
              oos_df$TSI_PZ_WarmII, oos_df$`GrainYld_NWPZ Warm`)
  add_oos_row("TSI_PZ_WarmII", "NWPZ Baseline + PZ Warm II", "PZ Warm I", "Out-of-Sample Cross-Regime",
              oos_df$TSI_PZ_WarmII, oos_df$`GrainYld_PZ Warm I`)
  # Out-of-sample TSI_Delhi
  add_oos_row("TSI_Delhi", "NWPZ Baseline + NWPZ Warm", "PZ Warm I", "Out-of-Sample Cross-Regime",
              oos_df$TSI_Delhi, oos_df$`GrainYld_PZ Warm I`)
  add_oos_row("TSI_Delhi", "NWPZ Baseline + NWPZ Warm", "PZ Warm II", "Out-of-Sample Cross-Regime",
              oos_df$TSI_Delhi, oos_df$`GrainYld_PZ Warm II`)
  # Out-of-sample Residual Retention (PT_Stress ~ PT_Opt)
  add_oos_row("Retention_Resid (PZ Warm II)", "lm(PT_PZ_WarmII ~ PT_Opt)", "NWPZ Warm", "Out-of-Sample Residual Index",
              oos_df$Retention_Resid, oos_df$`GrainYld_NWPZ Warm`)
  add_oos_row("Retention_Resid (PZ Warm II)", "lm(PT_PZ_WarmII ~ PT_Opt)", "PZ Warm I", "Out-of-Sample Residual Index",
              oos_df$Retention_Resid, oos_df$`GrainYld_PZ Warm I`)
  add_oos_row("Retention_Resid_Delhi (NWPZ Warm)", "lm(PT_NWPZ_Warm ~ PT_Opt)", "PZ Warm I", "Out-of-Sample Residual Index",
              oos_df$Retention_Resid_Delhi, oos_df$`GrainYld_PZ Warm I`)
  add_oos_row("Retention_Resid_Delhi (NWPZ Warm)", "lm(PT_NWPZ_Warm ~ PT_Opt)", "PZ Warm II", "Out-of-Sample Residual Index",
              oos_df$Retention_Resid_Delhi, oos_df$`GrainYld_PZ Warm II`)
  # Baseline benchmarks
  add_oos_row("Baseline Yield (NWPZ Baseline)", "NWPZ Baseline", "NWPZ Warm", "Out-of-Sample Baseline Yield Benchmark",
              oos_df$`GrainYld_NWPZ Baseline`, oos_df$`GrainYld_NWPZ Warm`)
  add_oos_row("Baseline PT (NWPZ Baseline)", "NWPZ Baseline", "NWPZ Warm", "Out-of-Sample Baseline PT Benchmark",
              oos_df$PT_Opt, oos_df$`GrainYld_NWPZ Warm`)
  # Within-regime (in-sample) comparisons vs raw PT_Stress
  add_oos_row("TSI_Delhi (In-Sample)", "NWPZ Baseline + NWPZ Warm", "NWPZ Warm", "In-Sample (Shares PT_NWPZ_Warm with Yield)",
              oos_df$TSI_Delhi, oos_df$`GrainYld_NWPZ Warm`, oos_df$`PT_Count_NWPZ Warm`)
  add_oos_row("TSI_PZ_WarmII (In-Sample)", "NWPZ Baseline + PZ Warm II", "PZ Warm II", "In-Sample (Shares PT_PZ_WarmII with Yield)",
              oos_df$TSI_PZ_WarmII, oos_df$`GrainYld_PZ Warm II`, oos_df$`PT_Count_PZ Warm II`)

  table_s6c_oos <- bind_rows(oos_rows)
  write.csv(table_s6c_oos, file.path(sec23_dir, "Table_TSI_OutOfSample_Validation.csv"), row.names = FALSE)

  # ----------------------------------------------------------------------------
  # 5.2 Vegetative Micro-Phenology & Tiller Survival Correlations With Grain Yield
  #     (Strictly Vegetative & Tiller Retention Metrics; Zero Duplicate AUTPC, Zero TKW;
  #      Delhi = Inferential Randomized Trials; Dharwad & TSI = Descriptive Within-Regime)
  # ----------------------------------------------------------------------------
  comp_metrics <- c(
    "AUTPC_norm",
    "Mean_Tiller_Early",
    "Mean_Tiller_Late",
    "AUTPC_Asymmetry",
    "T4_Count",
    "TCE",
    "PT_Count"
  )
  comp_labels <- c(
    "AUTPC_norm"        = "Norm. AUTPC\n(Tiller-Days / d)",
    "Mean_Tiller_Early" = "Early Tillering\n(Stage E\u2013T2)",
    "Mean_Tiller_Late"  = "Late Tillering\n(Stage T2\u2013T4)",
    "AUTPC_Asymmetry"   = "Early / Late\nAsymmetry",
    "T4_Count"          = "Final Vegetative\nTillers (T4)",
    "TCE"               = "Tiller Conversion\nEfficiency (TCE)",
    "PT_Count"          = "Productive Tillers\nat Harvest (PT)"
  )

  cor_results <- list()
  for (env in env_metadata$Env_Name) {
    sub_df <- blues_df %>% filter(Environment == env)
    scope_env <- if (env %in% c("NWPZ Baseline", "NWPZ Warm")) {
      "Inferential (Randomized Incomplete-Block Trial)"
    } else if (env == "PZ Warm I") {
      "Descriptive Within-Regime Only (Unrandomized Strip Layout; Plot Yield Decoupled from Sampled Plants, r = 0.05)"
    } else {
      "Descriptive Within-Regime Only (Unrandomized Strip Layout; r = 0.97 Driven by Plot-Level Stand/Plant-Size Factor)"
    }
    fn_env <- if (env %in% c("NWPZ Baseline", "NWPZ Warm")) "b" else "c"
    for (cm in comp_metrics) {
      ct_unadj <- cor.test(sub_df[[cm]], sub_df$GrainYld)
      r_val <- as.numeric(ct_unadj$estimate)
      sig_u <- if (ct_unadj$p.value < 0.001) "***" else if (ct_unadj$p.value < 0.01) "**" else if (ct_unadj$p.value < 0.05) "*" else "ns"
      cor_results[[paste0(env, "_", cm)]] <- data.frame(
        Environment          = env,
        Component            = cm,
        Component_Label      = sprintf("%s [%s]", gsub("\n", " ", comp_labels[cm]), fn_env),
        Plot_Label           = comp_labels[cm],
        Trait_Mean           = round(mean(sub_df[[cm]], na.rm = TRUE), 3),
        Trait_SD             = round(sd(sub_df[[cm]], na.rm = TRUE), 3),
        Pearson_r            = round(r_val, 4),
        R_Squared_Pct        = round((r_val^2) * 100, 2),
        t_statistic          = round(as.numeric(ct_unadj$statistic), 3),
        df                   = as.integer(ct_unadj$parameter),
        p_value              = ct_unadj$p.value,
        Significance         = sig_u,
        Footnote_Ref         = fn_env,
        Methodological_Scope = scope_env,
        stringsAsFactors     = FALSE
      )
    }

    # Also append TSI correlation row for the 3 stress environments into Table 4 (descriptive within-regime retention)
    if (env %in% c("NWPZ Warm", "PZ Warm I", "PZ Warm II")) {
      sub_tsi <- tsi_multi_df %>% filter(Environment == env)
      ct_t <- cor.test(sub_tsi$TSI, sub_tsi$Yield_Stress)
      r_t  <- as.numeric(ct_t$estimate)
      sig_t <- if (ct_t$p.value < 0.001) "***" else if (ct_t$p.value < 0.01) "**" else if (ct_t$p.value < 0.05) "*" else "ns"
      cor_results[[paste0(env, "_TSI")]] <- data.frame(
        Environment          = env,
        Component            = "TSI",
        Component_Label      = "Tiller Stability Index (TSI vs. NWPZ Baseline) [a]",
        Plot_Label           = "TSI",
        Trait_Mean           = round(mean(sub_tsi$TSI, na.rm = TRUE), 3),
        Trait_SD             = round(sd(sub_tsi$TSI, na.rm = TRUE), 3),
        Pearson_r            = round(r_t, 4),
        R_Squared_Pct        = round((r_t^2) * 100, 2),
        t_statistic          = round(as.numeric(ct_t$statistic), 3),
        df                   = as.integer(ct_t$parameter),
        p_value              = ct_t$p.value,
        Significance         = sig_t,
        Footnote_Ref         = "a",
        Methodological_Scope = "Structural identity / within-regime retention check (TSI shares PT_Stress with yield; out-of-sample |r| <= 0.11)",
        stringsAsFactors     = FALSE
      )
    }
  }

  table_4_correlations <- bind_rows(cor_results)
  stopifnot(sum(is.na(table_4_correlations)) == 0L)

  # Main Text Table 4 omits the wide Methodological_Scope column (replaced by concise Footnote_Ref [a, b, c])
  # Full archival version with Methodological_Scope is preserved in Archive_Rejected and in the Excel workbook
  t4_main_clean <- table_4_correlations %>% select(-Plot_Label, -Methodological_Scope)
  t4_archival   <- table_4_correlations %>% select(-Plot_Label)
  t4_footnote_text <- paste0(
    "Table 4 Footnotes: [a] Structural/algebraic retention check (TSI = [PT_Stress / PT_Baseline] x [PT_Stress / mean(PT_Stress)] shares PT_Stress with harvest spike density; out-of-sample cross-regime |r| <= 0.11). ",
    "[b] Independent biological hypothesis test in randomized incomplete-block trials (NWPZ Baseline & NWPZ Warm). ",
    "[c] Descriptive within-regime correlation in unrandomized strip trials (PZ Warm I & PZ Warm II)."
  )

  write.csv(
    t4_main_clean,
    file.path(sec23_dir, "Table_4_MicroPhenology_AUTPC_and_TSI_Yield_Correlations.csv"),
    row.names = FALSE
  )
  write.csv(
    t4_archival,
    file.path(archive_sec23_dir, "Table_4_Archival_With_Methodological_Scope.csv"),
    row.names = FALSE
  )

  # Compute Table S7: AUTPC Emergence Anchor Sensitivity (Count = 1.0 vs 0.0 at Emergence)
  table_s7_anchor <- blues_df %>%
    mutate(
      AUTPC_Count1      = AUTPC_From_Count_BLUEs,
      AUTPC_Count0      = AUTPC_From_Count_BLUEs - 0.5 * COUNT_AT_EMERGENCE_MAIN * Int_E_T1,
      AUTPC_norm_Count1 = AUTPC_Count1 / Total_E_T4,
      AUTPC_norm_Count0 = AUTPC_Count0 / Total_E_T4
    ) %>%
    group_by(Environment) %>%
    summarise(
      Vegetative_Duration_Days = first(Total_E_T4),
      Stage_E_to_T1_Days       = first(Int_E_T1),
      Mean_AUTPC_norm_Count1   = round(mean(AUTPC_norm_Count1, na.rm = TRUE), 4),
      Mean_AUTPC_norm_Count0   = round(mean(AUTPC_norm_Count0, na.rm = TRUE), 4),
      Delta_Mean               = round(Mean_AUTPC_norm_Count1 - Mean_AUTPC_norm_Count0, 4),
      Percent_Change           = round((Delta_Mean / Mean_AUTPC_norm_Count0) * 100, 2),
      Theoretical_Shift        = round((0.5 * first(Int_E_T1)) / first(Total_E_T4), 4),
      Pearson_r_Count1_vs_Count0 = round(cor(AUTPC_norm_Count1, AUTPC_norm_Count0), 6),
      .groups = "drop"
    )
  write.csv(
    table_s7_anchor,
    file.path(sec23_dir, "Table_S7_AUTPC_Emergence_Anchor_Sensitivity.csv"),
    row.names = FALSE
  )

  # Export Unified Multi-Sheet Excel Workbook for Section 2.3
  if (requireNamespace("openxlsx", quietly = TRUE)) {
    wb_sec23 <- openxlsx::createWorkbook()

    openxlsx::addWorksheet(wb_sec23, "Table_4_Yield_Correlations")
    openxlsx::writeData(wb_sec23, "Table_4_Yield_Correlations", t4_main_clean)
    openxlsx::addStyle(wb_sec23, "Table_4_Yield_Correlations", header_style, rows = 1, cols = 1:ncol(t4_main_clean), gridExpand = TRUE)
    openxlsx::writeData(wb_sec23, "Table_4_Yield_Correlations", data.frame(Note = t4_footnote_text), startRow = nrow(t4_main_clean) + 3, colNames = FALSE)
    openxlsx::setColWidths(wb_sec23, "Table_4_Yield_Correlations", cols = 1:ncol(t4_main_clean), widths = "auto")

    openxlsx::addWorksheet(wb_sec23, "Table_4_Archival_With_Scope")
    openxlsx::writeData(wb_sec23, "Table_4_Archival_With_Scope", t4_archival)
    openxlsx::addStyle(wb_sec23, "Table_4_Archival_With_Scope", header_style, rows = 1, cols = 1:ncol(t4_archival), gridExpand = TRUE)
    openxlsx::setColWidths(wb_sec23, "Table_4_Archival_With_Scope", cols = 1:ncol(t4_archival), widths = "auto")

    openxlsx::addWorksheet(wb_sec23, "Table_S6A_TSI_Transitions")
    openxlsx::writeData(wb_sec23, "Table_S6A_TSI_Transitions", table_s6a_transitions)
    openxlsx::addStyle(wb_sec23, "Table_S6A_TSI_Transitions", header_style, rows = 1, cols = 1:ncol(table_s6a_transitions), gridExpand = TRUE)
    openxlsx::setColWidths(wb_sec23, "Table_S6A_TSI_Transitions", cols = 1:ncol(table_s6a_transitions), widths = "auto")

    openxlsx::addWorksheet(wb_sec23, "Table_S6B_Genotype_TSI")
    openxlsx::writeData(wb_sec23, "Table_S6B_Genotype_TSI", table_s6b_genotypes)
    openxlsx::addStyle(wb_sec23, "Table_S6B_Genotype_TSI", header_style, rows = 1, cols = 1:ncol(table_s6b_genotypes), gridExpand = TRUE)
    openxlsx::setColWidths(wb_sec23, "Table_S6B_Genotype_TSI", cols = 1:ncol(table_s6b_genotypes), widths = "auto")

    openxlsx::addWorksheet(wb_sec23, "Table_S6C_TSI_OutOfSample")
    openxlsx::writeData(wb_sec23, "Table_S6C_TSI_OutOfSample", table_s6c_oos)
    openxlsx::addStyle(wb_sec23, "Table_S6C_TSI_OutOfSample", header_style, rows = 1, cols = 1:ncol(table_s6c_oos), gridExpand = TRUE)
    openxlsx::setColWidths(wb_sec23, "Table_S6C_TSI_OutOfSample", cols = 1:ncol(table_s6c_oos), widths = "auto")

    openxlsx::addWorksheet(wb_sec23, "Table_S7_Anchor_Sensitivity")
    openxlsx::writeData(wb_sec23, "Table_S7_Anchor_Sensitivity", table_s7_anchor)
    openxlsx::addStyle(wb_sec23, "Table_S7_Anchor_Sensitivity", header_style, rows = 1, cols = 1:ncol(table_s7_anchor), gridExpand = TRUE)
    openxlsx::setColWidths(wb_sec23, "Table_S7_Anchor_Sensitivity", cols = 1:ncol(table_s7_anchor), widths = "auto")

    openxlsx::saveWorkbook(wb_sec23, file.path(sec23_dir, "Section_2.3_Manuscript_and_Supplementary_Tables.xlsx"), overwrite = TRUE)
  }

  # ----------------------------------------------------------------------------
  # 5.3 Clean 2-Panel Master Figure for Section 2.3 (TSI_Compound_Yield_Validation_v2.png):
  #     Panel A (Left):  Descriptive Within-Regime TSI Retention vs. Yield (PZ Warm II)
  #     Panel B (Right): Multi-Stress TSI Classification & Stability Retention Across Regimes
  # ----------------------------------------------------------------------------
  heat_df <- table_4_correlations %>%
    filter(Component != "TSI") %>%
    mutate(
      Env_Factor  = factor(Environment, levels = rev(env_metadata$Env_Name)),
      Comp_Factor = factor(Plot_Label, levels = unname(comp_labels)),
      Cell_Label  = sprintf("r = %+.2f%s\n(p = %.1e)", Pearson_r, ifelse(Significance == "ns", "", Significance), p_value),
      Text_Col    = ifelse(abs(Pearson_r) > 0.62, "white", "black")
    )

  p_panel_heat_arch <- ggplot(heat_df, aes(x = Comp_Factor, y = Env_Factor, fill = Pearson_r)) +
    geom_tile(color = "white", linewidth = 0.9) +
    geom_text(aes(label = Cell_Label, color = Text_Col), fontface = "bold", size = 3.3, lineheight = 0.95, show.legend = FALSE) +
    scale_color_identity() +
    scale_fill_gradient2(
      low = "#3C5488", mid = "white", high = "#DC0000",
      midpoint = 0, limits = c(-1, 1), name = "Pearson r with Grain Yield:"
    ) +
    pub_theme +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5, face = "bold", size = 10.2, lineheight = 0.95),
      axis.text.y = element_text(face = "bold", size = 10.5),
      plot.subtitle = element_text(size = 9.5, color = "grey25", face = "bold")
    ) +
    labs(
      title = "A. Vegetative Micro-Phenology, AUTPC Decomposition & Tiller Survival Correlations With Grain Yield",
      subtitle = "Inferential in Delhi (NWPZ); descriptive within-regime in Dharwad (PZ Warm II r = 0.97 reflects plot-level stand/size scaling)",
      x = "",
      y = "Thermal Regime"
    )

  # --- Panel A: Descriptive Within-Regime TSI Retention Under Heat Stress (PZ Warm II) ---
  pza_sub <- blues_df %>%
    filter(Environment == "PZ Warm II") %>%
    mutate(
      Stability_Lab = factor(
        ifelse(TSI_PZ_WarmII >= 1.0, "Stable (TSI \u2265 1.0)", "Contrasting (TSI < 1.0)"),
        levels = c("Contrasting (TSI < 1.0)", "Stable (TSI \u2265 1.0)")
      )
    )
  ct_tsi_unadj <- cor.test(pza_sub$TSI_PZ_WarmII, pza_sub$GrainYld)
  r_tsi_yld    <- as.numeric(ct_tsi_unadj$estimate)
  p_tsi_yld    <- ct_tsi_unadj$p.value

  res_tsi_pza <- residuals(lm(TSI_PZ_WarmII ~ factor(Field_Pos_Block_Diagnostic), data = pza_sub))
  res_yld_pza <- residuals(lm(GrainYld ~ factor(Field_Pos_Block_Diagnostic), data = pza_sub))
  r_tsi_adj   <- cor(res_tsi_pza, res_yld_pza)
  df_tsi_adj  <- nrow(pza_sub) - length(unique(pza_sub$Field_Pos_Block_Diagnostic)) - 1L
  t_tsi_adj   <- r_tsi_adj * sqrt(df_tsi_adj / (1 - r_tsi_adj^2))
  p_tsi_adj   <- 2 * pt(-abs(t_tsi_adj), df = df_tsi_adj)

  n_stab_pza <- sum(pza_sub$TSI_PZ_WarmII >= 1.0, na.rm = TRUE)
  n_cont_pza <- sum(pza_sub$TSI_PZ_WarmII <  1.0, na.rm = TRUE)

  top_tsi_genos <- pza_sub %>%
    mutate(Rank_Score = as.numeric(scale(TSI_PZ_WarmII)) + as.numeric(scale(GrainYld))) %>%
    arrange(desc(Rank_Score)) %>%
    slice_head(n = 5)

  tsi_formula_expr <- "bold(TSI) == bgroup('(', frac(PT[Stress], PT[Baseline]), ')') %*% bgroup('(', frac(PT[Stress], bar(PT)[Stress]), ')')"

  p_panel_a_tsi <- ggplot(pza_sub, aes(x = TSI_PZ_WarmII, y = GrainYld)) +
    annotate("rect", xmin = 1.0, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "#3C5488", alpha = 0.05) +
    geom_vline(xintercept = 1.0, linetype = "dashed", color = "grey30", linewidth = 0.7) +
    geom_smooth(method = "lm", se = TRUE, color = "grey15", fill = "grey75", linewidth = 0.95, alpha = 0.35) +
    geom_point(aes(fill = Stability_Lab), shape = 21, size = 2.9, alpha = 0.86, color = "black", stroke = 0.4) +
    ggrepel::geom_text_repel(
      data = top_tsi_genos,
      aes(x = TSI_PZ_WarmII, y = GrainYld, label = Genotype),
      size = 3.0, fontface = "bold", color = "#1B2631",
      box.padding = 0.3, min.segment.length = 0.1, seed = 42
    ) +
    annotate(
      "label",
      x = 0.06, y = 0.79,
      label = tsi_formula_expr,
      parse = TRUE, hjust = 0, vjust = 1, size = 3.75,
      fill = "white", alpha = 0.94, color = "grey10", label.size = 0.35
    ) +
    annotate(
      "label",
      x = 1.03, y = 0.27,
      label = sprintf("Classification Threshold: TSI = 1.0\nStable: n = %d (%.1f%%) | Contrasting: n = %d (%.1f%%)",
                      n_stab_pza, 100 * n_stab_pza / nrow(pza_sub),
                      n_cont_pza, 100 * n_cont_pza / nrow(pza_sub)),
      hjust = 0, vjust = 0, size = 3.05, fontface = "bold",
      fill = "white", alpha = 0.92, color = "grey20", label.size = 0.25
    ) +
    scale_fill_manual(
      values = c("Contrasting (TSI < 1.0)" = "#DC0000", "Stable (TSI \u2265 1.0)" = "#3C5488"),
      name = "TSI Classification:"
    ) +
    scale_y_continuous(limits = c(0.24, 0.82), breaks = seq(0.3, 0.8, by = 0.1)) +
    pub_theme +
    theme(
      plot.subtitle = element_text(size = 8.9, color = "grey25", face = "bold", lineheight = 1.1)
    ) +
    labs(
      title = "A. Descriptive Within-Regime TSI Retention vs. Yield (PZ Warm II)",
      subtitle = sprintf(
        "In-sample r = %.2f (%s; pos-adj r = %.2f) reflects shared PT_Stress; out-of-sample cross-regime |r| \u2264 0.11 (ns)",
        r_tsi_yld, format_p_stat(p_tsi_yld), r_tsi_adj
      ),
      x = "Tiller Stability Index (TSI_PZ_WarmII)",
      y = sprintf("Grain Yield (%s)", YIELD_UNIT_LABEL)
    )

  # --- Panel B: Multi-Stress TSI Classification & Stability Retention Across Regimes ---
  tsi_alluv_df <- tsi_multi_df %>%
    mutate(
      Env_Factor = factor(
        Environment,
        levels = c("NWPZ Warm", "PZ Warm I", "PZ Warm II"),
        labels = c("NWPZ Warm", "PZ Warm I", "PZ Warm II")
      ),
      TSI_Class_Clean = factor(
        ifelse(TSI >= 1.0, "Stable (TSI \u2265 1.0)", "Contrasting (TSI < 1.0)"),
        levels = c("Stable (TSI \u2265 1.0)", "Contrasting (TSI < 1.0)")
      )
    )

  envs_tsi <- levels(tsi_alluv_df$Env_Factor)
  w_bar_tsi <- 0.44
  n_tot_tsi <- dplyr::n_distinct(tsi_alluv_df$Genotype)
  w_tsi_wide <- tsi_alluv_df %>%
    select(Genotype, Env_Factor, Class = TSI_Class_Clean) %>%
    pivot_wider(names_from = Env_Factor, values_from = Class)

  top_lvl_t <- "Stable (TSI \u2265 1.0)"
  bot_lvl_t <- "Contrasting (TSI < 1.0)"

  bar_rows_t <- list()
  p_bot_t <- numeric(length(envs_tsi))
  names(p_bot_t) <- envs_tsi
  for (k in seq_along(envs_tsi)) {
    e_k <- envs_tsi[k]
    n_bot <- sum(w_tsi_wide[[e_k]] == bot_lvl_t, na.rm = TRUE)
    n_top <- sum(w_tsi_wide[[e_k]] == top_lvl_t, na.rm = TRUE)
    p_bot <- 100 * n_bot / n_tot_tsi
    p_top <- 100 * n_top / n_tot_tsi
    p_bot_t[e_k] <- p_bot

    bar_rows_t[[length(bar_rows_t) + 1]] <- data.frame(
      x_idx = k, Env = e_k, Class = bot_lvl_t,
      xmin = k - w_bar_tsi / 2, xmax = k + w_bar_tsi / 2,
      ymin = 0, ymax = p_bot, ymid = p_bot / 2, pct = p_bot, n_cnt = n_bot,
      stringsAsFactors = FALSE
    )
    bar_rows_t[[length(bar_rows_t) + 1]] <- data.frame(
      x_idx = k, Env = e_k, Class = top_lvl_t,
      xmin = k - w_bar_tsi / 2, xmax = k + w_bar_tsi / 2,
      ymin = p_bot, ymax = 100, ymid = p_bot + p_top / 2, pct = p_top, n_cnt = n_top,
      stringsAsFactors = FALSE
    )
  }
  bars_tsi_df <- bind_rows(bar_rows_t) %>%
    mutate(Class = factor(Class, levels = c(top_lvl_t, bot_lvl_t)))

  rib_polys_t <- list()
  t_seq_t <- seq(0, 1, length.out = 40)
  s_seq_t <- 3 * t_seq_t^2 - 2 * t_seq_t^3
  for (k in 1:(length(envs_tsi) - 1)) {
    e1 <- envs_tsi[k]
    e2 <- envs_tsi[k + 1]
    x1 <- k + w_bar_tsi / 2
    x2 <- (k + 1) - w_bar_tsi / 2
    x_vals <- x1 + (x2 - x1) * t_seq_t

    f_bb <- 100 * sum(w_tsi_wide[[e1]] == bot_lvl_t & w_tsi_wide[[e2]] == bot_lvl_t, na.rm = TRUE) / n_tot_tsi
    f_bt <- 100 * sum(w_tsi_wide[[e1]] == bot_lvl_t & w_tsi_wide[[e2]] == top_lvl_t, na.rm = TRUE) / n_tot_tsi
    f_tb <- 100 * sum(w_tsi_wide[[e1]] == top_lvl_t & w_tsi_wide[[e2]] == bot_lvl_t, na.rm = TRUE) / n_tot_tsi
    f_tt <- 100 * sum(w_tsi_wide[[e1]] == top_lvl_t & w_tsi_wide[[e2]] == top_lvl_t, na.rm = TRUE) / n_tot_tsi

    pb1 <- p_bot_t[e1]
    pb2 <- p_bot_t[e2]

    streams_t <- list(
      list(id = paste0(k, "_bb"), cls = bot_lvl_t, y1_lo = 0,          y1_hi = f_bb,       y2_lo = 0,          y2_hi = f_bb,       w = f_bb),
      list(id = paste0(k, "_tb"), cls = top_lvl_t, y1_lo = pb1,        y1_hi = pb1 + f_tb, y2_lo = f_bb,       y2_hi = pb2,        w = f_tb),
      list(id = paste0(k, "_bt"), cls = bot_lvl_t, y1_lo = f_bb,       y1_hi = pb1,        y2_lo = pb2,        y2_hi = pb2 + f_bt, w = f_bt),
      list(id = paste0(k, "_tt"), cls = top_lvl_t, y1_lo = pb1 + f_tb, y1_hi = 100,        y2_lo = pb2 + f_bt, y2_hi = 100,        w = f_tt)
    )
    streams_t <- streams_t[order(sapply(streams_t, function(s) s$w), decreasing = TRUE)]
    for (st in streams_t) {
      if ((st$y1_hi - st$y1_lo) > 0.01) {
        y_lo <- st$y1_lo + (st$y2_lo - st$y1_lo) * s_seq_t
        y_hi <- st$y1_hi + (st$y2_hi - st$y1_hi) * s_seq_t
        rib_polys_t[[length(rib_polys_t) + 1]] <- data.frame(
          poly_id = factor(st$id, levels = unique(c(sapply(rib_polys_t, function(df) as.character(df$poly_id[1])), st$id))),
          Class   = st$cls,
          x       = c(x_vals, rev(x_vals)),
          y       = c(y_lo, rev(y_hi)),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  ribs_tsi_df <- bind_rows(rib_polys_t) %>%
    mutate(
      poly_id = factor(poly_id, levels = unique(poly_id)),
      Class   = factor(Class, levels = c(top_lvl_t, bot_lvl_t))
    )

  delw_sub   <- blues_df %>% filter(Environment == "NWPZ Warm")
  ct_tsi_del <- cor.test(delw_sub$TSI_Delhi, delw_sub$GrainYld)
  m_warm_ac  <- table_s6a_transitions %>% filter(Reference_Stress_Regime == "NWPZ Warm", Target_Stress_Regime == "PZ Warm II")

  p_panel_b_tsi_alluv <- ggplot() +
    geom_polygon(data = ribs_tsi_df, aes(x = x, y = y, group = poly_id, fill = Class), alpha = 0.30, color = "white", linewidth = 0.25) +
    geom_vline(xintercept = 1.5, linetype = "dashed", color = "grey30", linewidth = 0.65) +
    geom_rect(data = bars_tsi_df, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = Class),
              color = "black", linewidth = 0.45, alpha = 0.92) +
    geom_text(
      data = bars_tsi_df,
      aes(x = x_idx, y = ymid, label = sprintf("%.1f%%\n(n = %d)", pct, n_cnt)),
      color = "white", fontface = "bold", size = 3.4, lineheight = 0.95
    ) +
    scale_fill_manual(
      values = c("Stable (TSI \u2265 1.0)" = "#3C5488", "Contrasting (TSI < 1.0)" = "#DC0000"),
      name = "TSI Classification:"
    ) +
    scale_x_continuous(breaks = 1:3, labels = envs_tsi, expand = expansion(mult = c(0.08, 0.08))) +
    scale_y_continuous(labels = function(x) paste0(x, "%"), limits = c(0, 103), expand = c(0, 0)) +
    pub_theme +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5, face = "bold", size = 10.2),
      plot.subtitle = element_text(size = 9.4, color = "grey25", face = "bold", lineheight = 1.12)
    ) +
    labs(
      title = "B. Multi-Stress TSI Classification & Stability Retention",
      subtitle = sprintf(
        "Stable share ~40\u201345%% across regimes (McNemar \u03C7\u00B2 = %.2f, %s); %d genotypes retain TSI \u2265 1.0 across stress",
        m_warm_ac$McNemar_Chi_Square, m_warm_ac$McNemar_Significance, m_warm_ac$Retained_Stable_n
      ),
      x = "Stress Regime (Each Anchored to NWPZ Baseline Optimal Tillers)",
      y = "Genotype Composition (%)"
    )

  # Archive the 3-panel combined version (with AUTPC heatmap on top) for supervisor consultation
  p_sec23_3panel_archive <- (p_panel_heat_arch / (
    (p_panel_a_tsi + labs(title = "B. TSI Direct Yield Validation Under Heat Stress (PZ Warm II)")) |
    (p_panel_b_tsi_alluv + labs(title = "C. Multi-Stress TSI Classification & Stability Retention"))
  )) +
    plot_layout(heights = c(1.0, 1.08))
  ggsave(
    file.path(archive_sec23_dir, "TSI_and_AUTPC_Decomposition_3Panel_For_Supervisor_Consult.png"),
    plot = p_sec23_3panel_archive, width = 16.5, height = 11.2, dpi = 600
  )

  # Save the clean 2-panel TSI Master Figure to Section_2.3_Phenology_Indices
  p_sec23_master_2panel <- (p_panel_a_tsi | p_panel_b_tsi_alluv) +
    plot_layout(widths = c(1.05, 1.0)) +
    plot_annotation(
      caption = "Panel A validates the Tiller Stability Index (TSI, benchmark \u2265 1.0) against grain yield under severe vegetative heat stress (PZ Warm II).\nPanel B tracks TSI classification proportions and individual genotype stability flows across NWPZ Warm, PZ Warm I, and PZ Warm II regimes.",
      theme = theme(plot.caption = element_text(size = 9.5, color = "grey30", face = "italic", hjust = 0, lineheight = 1.15))
    )

  ggsave(
    file.path(sec23_dir, "TSI_Compound_Yield_Validation_v2.png"),
    plot = p_sec23_master_2panel, width = 15.8, height = 6.6, dpi = 600
  )
} else {
  cat("  -> Section 5 [SKIPPED - Assets Finalized & Protected]\n")
}


# ==============================================================================
# SECTION 6: NON-LINEAR THERMAL MODELING (GAM TRAJECTORIES)
# ==============================================================================
cat(">>> Executing Section 6: Non-Linear Thermal Modeling (GAM)...\n")

if (isTRUE(RUN_SECTION_6_THERMAL_GAM)) {
  sec24_dir <- file.path(assets_dir, "Section_2.4_Thermal_Modeling")
  dir.create(sec24_dir, showWarnings = FALSE, recursive = TRUE)

# Archive previous figures safely to Archive_Rejected
archive_rejected_dir <- file.path(project_root, "Archive_Rejected")
dir.create(archive_rejected_dir, showWarnings = FALSE, recursive = TRUE)
old_gam_png <- file.path(sec24_dir, "Thermal_GAM_Trajectories.png")
if (file.exists(old_gam_png)) {
  file.copy(old_gam_png, file.path(archive_rejected_dir, "Thermal_GAM_Trajectories_Prior.png"), overwrite = TRUE)
}
old_supp_png <- file.path(sec24_dir, "Thermal_GAM_Genotype_Faceted_Supplement.png")
if (file.exists(old_supp_png)) {
  file.copy(old_supp_png, file.path(archive_rejected_dir, "Thermal_GAM_Genotype_Faceted_Supplement_Prior.png"), overwrite = TRUE)
}

# 6.1 Developmental Milestone Windows & Cumulative Thermal/Chronological Coordinates
# Primary developmental scale: GDD from emergence (uncapped)
# Companion developmental scale: DAE (calendar progression)
# Sensitivity developmental scale: GDD capped at Tupper = 30 deg C
stage_map <- data.frame(
  Stage_Code = c("E-T1", "T1-T2", "T2-T3", "T3-T4"),
  Stage      = c("T1", "T2", "T3", "T4"),
  stringsAsFactors = FALSE
)

if (!exists("gdd_from_emergence_table") || is.null(gdd_from_emergence_table)) {
  gdd_from_emergence_table <- read.csv(file.path(sec1_dir, "Table_GDD_Emergence_Dynamics.csv"), stringsAsFactors = FALSE)
}

stage_milestones <- gdd_from_emergence_table %>%
  inner_join(stage_map, by = c("Stage" = "Stage_Code")) %>%
  group_by(Environment) %>%
  mutate(DAE = cumsum(Duration_Days)) %>%
  ungroup() %>%
  transmute(
    Environment,
    Stage = factor(Stage.y, levels = c("T1", "T2", "T3", "T4")),
    Duration_Days,
    DAE,
    Stage_GDD_Uncapped,
    Stage_GDD_Capped,
    GDD        = Cum_GDD_Uncapped, # Primary GDD axis (uncapped from emergence)
    GDD_Capped = Cum_GDD_Capped   # Sensitivity check axis (capped from emergence)
  )

# Print cumulative GDD totals and percentage spread across environments (Task 10)
total_gdd_summary <- stage_milestones %>%
  filter(Stage == "T4") %>%
  select(Environment, Total_GDD_Uncapped = GDD, Total_GDD_Capped = GDD_Capped)

spread_uncapped_pct <- (max(total_gdd_summary$Total_GDD_Uncapped) - min(total_gdd_summary$Total_GDD_Uncapped)) / min(total_gdd_summary$Total_GDD_Uncapped) * 100
spread_capped_pct   <- (max(total_gdd_summary$Total_GDD_Capped) - min(total_gdd_summary$Total_GDD_Capped)) / min(total_gdd_summary$Total_GDD_Capped) * 100

cat("\n  Cumulative Thermal Accumulation from Emergence to T4 Milestone:\n")
for (i in 1:nrow(total_gdd_summary)) {
  cat(sprintf("    %-15s : Uncapped = %6.2f deg C d | Capped = %6.2f deg C d\n",
              total_gdd_summary$Environment[i],
              total_gdd_summary$Total_GDD_Uncapped[i],
              total_gdd_summary$Total_GDD_Capped[i]))
}
cat(sprintf("  Thermal Accumulation Range (Uncapped): %.2f - %.2f deg C d (Spread: %.2f%%)\n",
            min(total_gdd_summary$Total_GDD_Uncapped), max(total_gdd_summary$Total_GDD_Uncapped), spread_uncapped_pct))
cat(sprintf("  Thermal Accumulation Range (Capped)  : %.2f - %.2f deg C d (Spread: %.2f%%)\n\n",
            min(total_gdd_summary$Total_GDD_Capped), max(total_gdd_summary$Total_GDD_Capped), spread_capped_pct))

# Export stage GDD matrix (Requirement 7)
write.csv(stage_milestones, file.path(sec24_dir, "Stage_GDD_Milestone_Matrix.csv"), row.names = FALSE)

# 6.2 Reshape Real Empirical Observations (Task 1: N = 3,200, strictly excluding Emergence anchor)
# 200 genotypes x 4 environments x 4 developmental stages (T1...T4)
gam_data <- blues_df %>%
  select(Genotype, Environment, Env_Factor, T1_Count, T2_Count, T3_Count, T4_Count) %>%
  pivot_longer(
    cols = c(T1_Count, T2_Count, T3_Count, T4_Count),
    names_to = "Stage_Col",
    values_to = "Count"
  ) %>%
  mutate(
    Stage = case_when(
      Stage_Col == "T1_Count" ~ "T1",
      Stage_Col == "T2_Count" ~ "T2",
      Stage_Col == "T3_Count" ~ "T3",
      Stage_Col == "T4_Count" ~ "T4"
    ),
    Stage = factor(Stage, levels = c("T1", "T2", "T3", "T4")),
    Genotype = as.factor(Genotype),
    Env_Factor = factor(Environment, levels = env_metadata$Env_Name),
    GE = interaction(Genotype, Env_Factor, drop = TRUE)
  ) %>%
  left_join(
    stage_milestones %>% select(Environment, Stage, DAE, GDD, GDD_Capped),
    by = c("Environment", "Stage")
  ) %>%
  mutate(Stage = factor(Stage, levels = c("T1", "T2", "T3", "T4")))

# Stopifnot verification checks (Task 1)
stopifnot(!any(is.na(gam_data$GDD)))
stopifnot(!any(is.na(gam_data$GDD_Capped)))
stopifnot(!any(is.na(gam_data$DAE)))
stopifnot(!any(is.na(gam_data$Count)))

check_stg <- gam_data %>%
  group_by(Genotype, Env_Factor) %>%
  summarise(n_stages = n(), .groups = "drop")
stopifnot(all(check_stg$n_stages == 4))
stopifnot(nrow(gam_data) == 200 * 4 * 4)

cat("  GAM Dataset Ready: N =", nrow(gam_data), "observations (200 genotypes x 4 regimes x 4 stages).\n")
cat("  Emergence anchor strictly excluded from model estimation dataset.\n\n")

# 6.3 Fit GAM Models via Exact REML (Task 2: mgcv::gam(method = 'REML'), not bam)
# Evaluating Shared-Curve vs Environment-Specific Curves across:
#   Primary: GDD (Uncapped thermal time from emergence)
#   Sensitivity Check: GDD_Capped (Heat-capped thermal time from emergence)
#   Companion: DAE (Days after emergence)
cache_dir <- file.path(project_root, "Cache_BLUEs")
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

# Remove legacy unhashed static GAM cache files so updated BLUEs always trigger fresh REML fits
legacy_gam_rds <- list.files(cache_dir, pattern = "^gam_reml_m_(common|env)_(gdd|capped|dae|gdd_log)\\.rds$", full.names = TRUE)
if (length(legacy_gam_rds) > 0) {
  unlink(legacy_gam_rds)
}

fit_gam_reml <- function(rds_name, f_formula, data_df) {
  prefix <- sub("\\.rds$", "", rds_name)
  gam_hash <- substr(digest::digest(list(f_formula, data_df, as.character(packageVersion("mgcv"))), algo = "md5"), 1, 10)
  hashed_name <- paste0(prefix, "_", gam_hash, ".rds")
  p_rds <- file.path(cache_dir, hashed_name)
  if (file.exists(p_rds) && !isTRUE(CLEAR_BLUE_CACHE)) {
    m_cached <- readRDS(p_rds)
    if (identical(levels(m_cached$model$Env_Factor), levels(data_df$Env_Factor))) {
      cat("  Loading cached GAM:", hashed_name, "\n")
      return(m_cached)
    } else {
      cat("  Discarding stale-level cached GAM:", hashed_name, "\n")
      unlink(p_rds)
    }
  }
  existing_cands <- list.files(cache_dir, pattern = paste0("^", prefix, "_[0-9a-f]{10}\\.rds$"), full.names = TRUE)
  for (cand in existing_cands) {
    if (!isTRUE(CLEAR_BLUE_CACHE)) {
      m_cand <- tryCatch(readRDS(cand), error = function(e) NULL)
      if (!is.null(m_cand) && identical(levels(m_cand$model$Env_Factor), levels(data_df$Env_Factor))) {
        cat("  Loading & re-keying verified cached GAM:", basename(cand), "->", hashed_name, "\n")
        saveRDS(m_cand, p_rds)
        return(m_cand)
      } else {
        cat("  Removing stale-level cached GAM:", basename(cand), "\n")
        unlink(cand)
      }
    }
  }
  cat("  Fitting GAM via REML:", hashed_name, "...\n")
  m_fit <- mgcv::gam(as.formula(f_formula), data = data_df, method = "REML")
  saveRDS(m_fit, p_rds)
  m_fit
}

# (a) Primary Thermal Time Axis (Uncapped GDD)
m_common <- fit_gam_reml(
  "gam_reml_m_common_gdd.rds",
  "Count ~ Env_Factor + s(GDD, k = 4) + s(Genotype, bs = 're') + s(GE, bs = 're')",
  gam_data
)
m_env <- fit_gam_reml(
  "gam_reml_m_env_gdd.rds",
  "Count ~ Env_Factor + s(GDD, by = Env_Factor, k = 3) + s(Genotype, bs = 're') + s(GE, bs = 're')",
  gam_data
)

# (b) Sensitivity Axis (Capped GDD)
m_common_capped <- fit_gam_reml(
  "gam_reml_m_common_capped.rds",
  "Count ~ Env_Factor + s(GDD_Capped, k = 4) + s(Genotype, bs = 're') + s(GE, bs = 're')",
  gam_data
)
m_env_capped <- fit_gam_reml(
  "gam_reml_m_env_capped.rds",
  "Count ~ Env_Factor + s(GDD_Capped, by = Env_Factor, k = 3) + s(Genotype, bs = 're') + s(GE, bs = 're')",
  gam_data
)

# (c) Companion Chronological Axis (DAE)
m_common_dae <- fit_gam_reml(
  "gam_reml_m_common_dae.rds",
  "Count ~ Env_Factor + s(DAE, k = 4) + s(Genotype, bs = 're') + s(GE, bs = 're')",
  gam_data
)
m_env_dae <- fit_gam_reml(
  "gam_reml_m_env_dae.rds",
  "Count ~ Env_Factor + s(DAE, by = Env_Factor, k = 3) + s(Genotype, bs = 're') + s(GE, bs = 're')",
  gam_data
)

# Extract Smooth edf and Total edf
extract_smooth_edf <- function(mod) {
  s <- summary(mod)
  st <- s$s.table
  non_re <- st[!grepl("Genotype|GE", rownames(st)), , drop = FALSE]
  sum(non_re[, "edf"])
}

# 6.4 Model Comparison & AIC Table (Task 2 & 7)
aic_summary_table <- data.frame(
  Axis = c(
    rep("Thermal Time (GDD Uncapped, Primary)", 2),
    rep("Thermal Time (GDD Capped, Sensitivity)", 2),
    rep("Chronological Time (DAE, Companion)", 2)
  ),
  Model = c(
    "m_common (Shared curve)", "m_env (Environment-specific curves)",
    "m_common_capped (Shared curve)", "m_env_capped (Environment-specific curves)",
    "m_common_dae (Shared curve)", "m_env_dae (Environment-specific curves)"
  ),
  Basis_k = c(4, 3, 4, 3, 4, 3),
  Formula = c(
    "Count ~ Env_Factor + s(GDD, k = 4) + s(Genotype, bs = 're') + s(GE, bs = 're')",
    "Count ~ Env_Factor + s(GDD, by = Env_Factor, k = 3) + s(Genotype, bs = 're') + s(GE, bs = 're')",
    "Count ~ Env_Factor + s(GDD_Capped, k = 4) + s(Genotype, bs = 're') + s(GE, bs = 're')",
    "Count ~ Env_Factor + s(GDD_Capped, by = Env_Factor, k = 3) + s(Genotype, bs = 're') + s(GE, bs = 're')",
    "Count ~ Env_Factor + s(DAE, k = 4) + s(Genotype, bs = 're') + s(GE, bs = 're')",
    "Count ~ Env_Factor + s(DAE, by = Env_Factor, k = 3) + s(Genotype, bs = 're') + s(GE, bs = 're')"
  ),
  Smooth_edf = round(c(
    extract_smooth_edf(m_common), extract_smooth_edf(m_env),
    extract_smooth_edf(m_common_capped), extract_smooth_edf(m_env_capped),
    extract_smooth_edf(m_common_dae), extract_smooth_edf(m_env_dae)
  ), 2),
  Total_edf = round(c(
    sum(m_common$edf), sum(m_env$edf),
    sum(m_common_capped$edf), sum(m_env_capped$edf),
    sum(m_common_dae$edf), sum(m_env_dae$edf)
  ), 2),
  AIC = round(c(
    AIC(m_common), AIC(m_env),
    AIC(m_common_capped), AIC(m_env_capped),
    AIC(m_common_dae), AIC(m_env_dae)
  ), 2),
  Delta_AIC = round(c(
    0, AIC(m_env) - AIC(m_common),
    0, AIC(m_env_capped) - AIC(m_common_capped),
    0, AIC(m_env_dae) - AIC(m_common_dae)
  ), 2),
  REML_Score = round(c(
    m_common$gcv.ubre, m_env$gcv.ubre,
    m_common_capped$gcv.ubre, m_env_capped$gcv.ubre,
    m_common_dae$gcv.ubre, m_env_dae$gcv.ubre
  ), 2),
  Adj_R_Squared = round(c(
    summary(m_common)$r.sq, summary(m_env)$r.sq,
    summary(m_common_capped)$r.sq, summary(m_env_capped)$r.sq,
    summary(m_common_dae)$r.sq, summary(m_env_dae)$r.sq
  ), 4),
  Deviance_Explained_Pct = round(c(
    summary(m_common)$dev.expl * 100, summary(m_env)$dev.expl * 100,
    summary(m_common_capped)$dev.expl * 100, summary(m_env_capped)$dev.expl * 100,
    summary(m_common_dae)$dev.expl * 100, summary(m_env_dae)$dev.expl * 100
  ), 2),
  stringsAsFactors = FALSE
)

cat("\n--- Model Architecture & AIC Comparison Table ---\n")
print(aic_summary_table)
write.csv(aic_summary_table, file.path(sec24_dir, "GAM_Model_Comparison_AIC.csv"), row.names = FALSE)

# Extract Smooth Term Table & Parametric Coefficients (Task 7)
extract_smooth_table <- function(mod, mod_name, axis_name) {
  s <- summary(mod)
  df <- as.data.frame(s$s.table)
  df$Term <- rownames(df)
  rownames(df) <- NULL
  df$Model <- mod_name
  df$Axis <- axis_name
  names(df)[names(df) == "p-value"] <- "p_value"
  names(df)[names(df) == "Ref.df"]  <- "Ref_df"
  names(df)[names(df) == "F"]       <- "F_stat"
  df %>% select(Model, Axis, Term, edf, Ref_df, F_stat, p_value)
}

extract_param_table <- function(mod, mod_name, axis_name) {
  s <- summary(mod)
  df <- as.data.frame(s$p.table)
  df$Term <- rownames(df)
  df$Term <- sub("^Env_FactorPZ Warm$", "Env_FactorPZ Warm I", df$Term)
  df$Term <- sub("^Env_FactorPZ Acute$", "Env_FactorPZ Warm II", df$Term)
  rownames(df) <- NULL
  df$Model <- mod_name
  df$Axis <- axis_name
  names(df)[names(df) == "Estimate"]   <- "Estimate"
  names(df)[names(df) == "Std. Error"] <- "Std_Error"
  names(df)[names(df) == "t value"]    <- "t_stat"
  names(df)[names(df) == "Pr(>|t|)"]   <- "p_value"
  df %>% select(Model, Axis, Term, Estimate, Std_Error, t_stat, p_value)
}

all_smooths <- bind_rows(
  extract_smooth_table(m_common, "m_common", "Thermal_GDD_Uncapped"),
  extract_smooth_table(m_env, "m_env", "Thermal_GDD_Uncapped"),
  extract_smooth_table(m_common_capped, "m_common_capped", "Thermal_GDD_Capped"),
  extract_smooth_table(m_env_capped, "m_env_capped", "Thermal_GDD_Capped"),
  extract_smooth_table(m_common_dae, "m_common_dae", "Chronological_DAE"),
  extract_smooth_table(m_env_dae, "m_env_dae", "Chronological_DAE")
)
write.csv(all_smooths, file.path(sec24_dir, "GAM_Smooth_Terms.csv"), row.names = FALSE)

all_params <- bind_rows(
  extract_param_table(m_common, "m_common", "Thermal_GDD_Uncapped"),
  extract_param_table(m_env, "m_env", "Thermal_GDD_Uncapped"),
  extract_param_table(m_common_capped, "m_common_capped", "Thermal_GDD_Capped"),
  extract_param_table(m_env_capped, "m_env_capped", "Thermal_GDD_Capped"),
  extract_param_table(m_common_dae, "m_common_dae", "Chronological_DAE"),
  extract_param_table(m_env_dae, "m_env_dae", "Chronological_DAE")
)
write.csv(all_params, file.path(sec24_dir, "GAM_Parametric_Coefficients.csv"), row.names = FALSE)

# 6.5 Companion Marginal Linear Model (Item 12: nlme::gls with unstructured stage covariance & stage-specific variances)
cat("  Fitting Companion Marginal Model (nlme::gls with corSymm + varIdent)...\n")
gam_data_gls <- gam_data %>%
  mutate(Stage_Num = as.integer(Stage)) %>%
  arrange(GE, Stage_Num)

fit_gls_companion <- nlme::gls(
  Count ~ Env_Factor * Stage,
  data = gam_data_gls,
  correlation = nlme::corSymm(form = ~ Stage_Num | GE),
  weights = nlme::varIdent(form = ~ 1 | Stage),
  method = "REML"
)

emm_companion <- emmeans::emmeans(
  fit_gls_companion, ~ Stage | Env_Factor,
  mode = "df.error"
)
emm_companion_df <- as.data.frame(emm_companion)
if ("asymp.LCL" %in% names(emm_companion_df) && !("lower.CL" %in% names(emm_companion_df))) {
  emm_companion_df$lower.CL <- emm_companion_df$asymp.LCL
  emm_companion_df$upper.CL <- emm_companion_df$asymp.UCL
}
write.csv(emm_companion_df, file.path(sec24_dir, "LMER_Companion_Stage_Means.csv"), row.names = FALSE)

# 6.6 Residual Diagnostics & Variance Heterogeneity Evaluation (Task 4)
gam_diag <- data.frame(
  Count        = gam_data$Count,
  Fitted       = fitted(m_env),
  Residual     = residuals(m_env, type = "pearson"),
  Abs_Residual = abs(residuals(m_env, type = "pearson")),
  Env_Factor   = gam_data$Env_Factor
)

cor_diag <- cor.test(gam_diag$Fitted, gam_diag$Abs_Residual)
cat(sprintf("  Diagnostic Residual vs Fitted Correlation: r = %.3f (p = %.2e)\n",
            cor_diag$estimate, cor_diag$p.value))

p_diag <- ggplot(gam_diag, aes(x = Fitted, y = Residual, color = Env_Factor)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(alpha = 0.35, size = 1.8) +
  geom_smooth(method = "loess", se = TRUE, color = "black", linewidth = 0.8) +
  scale_color_npg() +
  pub_theme +
  labs(
    title = "Residual Diagnostic for Primary Thermal GAM (m_env)",
    subtitle = "Pearson residuals vs fitted values across environments",
    x = "Fitted Values",
    y = "Pearson Residuals",
    color = "Environment:"
  )
ggsave(file.path(sec24_dir, "GAM_Residual_Diagnostics.png"), plot = p_diag, width = 8.5, height = 5.5, dpi = 600)

# Refit on log(Count) as specified in Task 4
m_env_log <- fit_gam_reml(
  "gam_reml_m_env_gdd_log.rds",
  "log(Count) ~ Env_Factor + s(GDD, by = Env_Factor, k = 3) + s(Genotype, bs = 're') + s(GE, bs = 're')",
  gam_data
)
log_diag <- data.frame(
  Log_Count    = log(gam_data$Count),
  Fitted       = fitted(m_env_log),
  Residual     = residuals(m_env_log, type = "pearson"),
  Abs_Residual = abs(residuals(m_env_log, type = "pearson")),
  Env_Factor   = gam_data$Env_Factor
)
cor_log_diag <- cor.test(log_diag$Fitted, log_diag$Abs_Residual)
cat(sprintf("  Log-transformed Model Residual Correlation: r = %.3f (p = %.2e)\n",
            cor_log_diag$estimate, cor_log_diag$p.value))

# 6.7 Marginal Population Predictions & Assumed Emergence Points (Task 5)
# Exclude random effects s(Genotype) and s(GE) to capture population curves
pred_grid_gdd <- gam_data %>%
  group_by(Env_Factor) %>%
  summarise(min_x = min(GDD), max_x = max(GDD), .groups = "drop") %>%
  group_by(Env_Factor) %>%
  reframe(GDD = seq(min_x, max_x, length.out = 150)) %>%
  mutate(
    Genotype = factor(levels(gam_data$Genotype)[1], levels = levels(gam_data$Genotype)),
    GE       = factor(levels(gam_data$GE)[1], levels = levels(gam_data$GE))
  )

preds_gdd <- predict(m_env, newdata = pred_grid_gdd, se.fit = TRUE, exclude = c("s(Genotype)", "s(GE)"))
pred_grid_gdd$Pred_Count <- preds_gdd$fit
pred_grid_gdd$SE         <- preds_gdd$se.fit
pred_grid_gdd$Lower_CI   <- pred_grid_gdd$Pred_Count - 1.96 * pred_grid_gdd$SE
pred_grid_gdd$Upper_CI   <- pred_grid_gdd$Pred_Count + 1.96 * pred_grid_gdd$SE

pred_grid_dae <- gam_data %>%
  group_by(Env_Factor) %>%
  summarise(min_x = min(DAE), max_x = max(DAE), .groups = "drop") %>%
  group_by(Env_Factor) %>%
  reframe(DAE = seq(min_x, max_x, length.out = 150)) %>%
  mutate(
    Genotype = factor(levels(gam_data$Genotype)[1], levels = levels(gam_data$Genotype)),
    GE       = factor(levels(gam_data$GE)[1], levels = levels(gam_data$GE))
  )

preds_dae <- predict(m_env_dae, newdata = pred_grid_dae, se.fit = TRUE, exclude = c("s(Genotype)", "s(GE)"))
pred_grid_dae$Pred_Count <- preds_dae$fit
pred_grid_dae$SE         <- preds_dae$se.fit
pred_grid_dae$Lower_CI   <- pred_grid_dae$Pred_Count - 1.96 * pred_grid_dae$SE
pred_grid_dae$Upper_CI   <- pred_grid_dae$Pred_Count + 1.96 * pred_grid_dae$SE

# Empirical stage means +/- SE for overlay
stage_means <- gam_data %>%
  group_by(Env_Factor, Stage, GDD, DAE) %>%
  summarise(
    Mean_Count = mean(Count, na.rm = TRUE),
    SE_Count   = sd(Count, na.rm = TRUE) / sqrt(n()),
    .groups    = "drop"
  )

# Assumed emergence anchor points (Count = 1 at GDD = 0 / DAE = 0, strictly outside fit domain)
emergence_anchor_gdd <- data.frame(
  GDD   = 0,
  Count = 1,
  Label = "assumed (E)"
)

emergence_anchor_dae <- data.frame(
  DAE   = 0,
  Count = 1,
  Label = "assumed (E)"
)

# 6.8 Main Figure: Dual-Panel Thermal (GDD) & Chronological (DAE) GAM Trajectories (Task 6a)
# Axis labels use Unicode degree/multiplication symbols; in-panel stats badges; stage labels at
# NWPZ Baseline reference positions (well-separated GDD/DAE coordinates, no inter-environment clutter)

# Pre-compute in-panel annotation statistics (all from fitted model objects; no numeric hardcoding)
.m_env_r2      <- round(summary(m_env)$r.sq, 3)
.m_env_dev     <- round(summary(m_env)$dev.expl * 100, 1)
.delta_aic_gdd <- round(AIC(m_env) - AIC(m_common), 1)
.stats_gdd     <- sprintf(
  "Adj-R\u00b2 = %.3f  |  Dev. Expl. = %.1f%%\n\u0394AIC (env-specific vs. shared) = %.1f",
  .m_env_r2, .m_env_dev, .delta_aic_gdd
)

.m_dae_r2      <- round(summary(m_env_dae)$r.sq, 3)
.m_dae_dev     <- round(summary(m_env_dae)$dev.expl * 100, 1)
.delta_aic_dae <- round(AIC(m_env_dae) - AIC(m_common_dae), 1)
.stats_dae     <- sprintf(
  "Adj-R\u00b2 = %.3f  |  Dev. Expl. = %.1f%%\n\u0394AIC (env-specific vs. shared) = %.1f",
  .m_dae_r2, .m_dae_dev, .delta_aic_dae
)

# NWPZ Baseline reference label positions (used in both panels for stage label annotation)
.stage_lab_gdd <- stage_means %>% filter(Env_Factor == "NWPZ Baseline")
.stage_lab_dae <- stage_means %>% filter(Env_Factor == "NWPZ Baseline")

p_main_gdd <- ggplot() +
  geom_point(data = emergence_anchor_gdd,
             aes(x = GDD, y = Count),
             shape = 1, size = 3.5, stroke = 1.0, color = "black") +
  geom_text(data = emergence_anchor_gdd,
            aes(x = GDD, y = Count, label = Label),
            hjust = -0.2, vjust = 0.5, fontface = "italic", size = 3.2, color = "grey30") +
  geom_ribbon(data = pred_grid_gdd,
              aes(x = GDD, ymin = Lower_CI, ymax = Upper_CI, fill = Env_Factor),
              alpha = 0.30) +
  geom_line(data = pred_grid_gdd,
            aes(x = GDD, y = Pred_Count, color = Env_Factor),
            linewidth = 1.3) +
  geom_errorbar(data = stage_means,
                aes(x = GDD, ymin = Mean_Count - SE_Count, ymax = Mean_Count + SE_Count,
                    color = Env_Factor),
                width = 12, linewidth = 0.7) +
  geom_point(data = stage_means,
             aes(x = GDD, y = Mean_Count, color = Env_Factor),
             size = 3) +
  geom_text(data = .stage_lab_gdd,
            aes(x = GDD, y = Mean_Count + SE_Count + 1.0, label = as.character(Stage)),
            color = "grey20", size = 3.0, fontface = "bold") +
  annotate("text",
           x = 705, y = 1.5,
           label = .stats_gdd,
           hjust = 1, vjust = 0, size = 2.6, color = "grey25", fontface = "italic") +
  scale_color_npg() +
  scale_fill_npg() +
  scale_x_continuous(limits = c(-10, 720), breaks = seq(0, 700, 100)) +
  scale_y_continuous(
    limits = c(0, max(c(pred_grid_gdd$Upper_CI, gam_data$Count), na.rm = TRUE) * 1.05),
    breaks = seq(0, 20, 2)
  ) +
  pub_theme +
  labs(
    title    = "A. Thermal Progression",
    subtitle = "Population splines (\u00b195% CI) and stage empirical means (\u00b1SE)",
    x        = "Cumulative GDD from emergence (base 5\u00b0C, \u00b0C\u00b7d)",
    y        = "Tillers per plant",
    color    = "Environment:",
    fill     = "Environment:"
  )

p_main_dae <- ggplot() +
  geom_point(data = emergence_anchor_dae,
             aes(x = DAE, y = Count),
             shape = 1, size = 3.5, stroke = 1.0, color = "black") +
  geom_text(data = emergence_anchor_dae,
            aes(x = DAE, y = Count, label = Label),
            hjust = -0.2, vjust = 0.5, fontface = "italic", size = 3.2, color = "grey30") +
  geom_ribbon(data = pred_grid_dae,
              aes(x = DAE, ymin = Lower_CI, ymax = Upper_CI, fill = Env_Factor),
              alpha = 0.30) +
  geom_line(data = pred_grid_dae,
            aes(x = DAE, y = Pred_Count, color = Env_Factor),
            linewidth = 1.3) +
  geom_errorbar(data = stage_means,
                aes(x = DAE, ymin = Mean_Count - SE_Count, ymax = Mean_Count + SE_Count,
                    color = Env_Factor),
                width = 1.5, linewidth = 0.7) +
  geom_point(data = stage_means,
             aes(x = DAE, y = Mean_Count, color = Env_Factor),
             size = 3) +
  geom_text(data = .stage_lab_dae,
            aes(x = DAE, y = Mean_Count + SE_Count + 1.0, label = as.character(Stage)),
            color = "grey20", size = 3.0, fontface = "bold") +
  annotate("text",
           x = 74, y = 1.5,
           label = .stats_dae,
           hjust = 1, vjust = 0, size = 2.6, color = "grey25", fontface = "italic") +
  scale_color_npg() +
  scale_fill_npg() +
  scale_x_continuous(limits = c(-2, 76), breaks = seq(0, 70, 10)) +
  scale_y_continuous(
    limits = c(0, max(c(pred_grid_dae$Upper_CI, gam_data$Count), na.rm = TRUE) * 1.05),
    breaks = seq(0, 20, 2)
  ) +
  pub_theme +
  labs(
    title    = "B. Calendar Progression",
    subtitle = "Population splines (\u00b195% CI) and stage empirical means (\u00b1SE)",
    x        = "Days after emergence",
    y        = "Tillers per plant",
    color    = "Environment:",
    fill     = "Environment:"
  )

p_gam_main <- (p_main_gdd + p_main_dae) +
  plot_layout(guides = "collect") +
  plot_annotation(
    caption = paste0(
      "Stage labels (T1\u2013T4) annotated at NWPZ Baseline (reference environment) GDD and DAE positions. ",
      "Adj-R\u00b2 and Dev. Expl. from m_env (environment-specific REML splines, k\u2009=\u20093). ",
      "Random effects [s(Genotype), s(GE)] excluded from population-level predictions."
    ),
    theme = theme(
      plot.caption = element_text(size = 7, color = "grey40", hjust = 0, lineheight = 1.2)
    )
  ) &
  theme(legend.position = "bottom", legend.box = "horizontal")

ggsave(file.path(sec24_dir, "Thermal_GAM_Trajectories.png"),
       plot = p_gam_main, width = 14, height = 6.5, dpi = 600)

# 6.9 Supplement Figure: 4-Panel Faceted View with Individual Genotype Trajectories (Task 6b)
# Fixed shared x-axis (0-720 GDD) for cross-environment comparability; stage labels per facet
emergence_anchor_facet <- env_metadata %>%
  transmute(
    Env_Factor,
    GDD   = 0,
    Count = 1,
    Label = "assumed (E)"
  )

p_supp <- ggplot() +
  geom_point(data = emergence_anchor_facet,
             aes(x = GDD, y = Count),
             shape = 1, size = 3.5, stroke = 1.0, color = "black") +
  geom_text(data = emergence_anchor_facet,
            aes(x = GDD, y = Count, label = Label),
            hjust = -0.2, vjust = 0.5, fontface = "italic", size = 3.2, color = "grey30") +
  geom_line(data = gam_data,
            aes(x = GDD, y = Count, group = Genotype),
            color = "grey75", alpha = 0.25, linewidth = 0.35) +
  geom_ribbon(data = pred_grid_gdd,
              aes(x = GDD, ymin = Lower_CI, ymax = Upper_CI, fill = Env_Factor),
              alpha = 0.35) +
  geom_line(data = pred_grid_gdd,
            aes(x = GDD, y = Pred_Count, color = Env_Factor),
            linewidth = 1.4) +
  geom_errorbar(data = stage_means,
                aes(x = GDD, ymin = Mean_Count - SE_Count, ymax = Mean_Count + SE_Count),
                width = 15, linewidth = 0.8, color = "black") +
  geom_point(data = stage_means,
             aes(x = GDD, y = Mean_Count, fill = Env_Factor),
             size = 3.5, shape = 21, stroke = 1.0, color = "black") +
  geom_text(data = stage_means,
            aes(x = GDD, y = Mean_Count + SE_Count + 1.2, label = as.character(Stage)),
            color = "grey20", size = 2.8, fontface = "bold") +
  facet_wrap(~ Env_Factor, ncol = 2, scales = "fixed") +
  scale_color_npg() +
  scale_fill_npg() +
  scale_x_continuous(limits = c(0, 720), breaks = seq(0, 700, 100)) +
  scale_y_continuous(limits = c(0, max(gam_data$Count, na.rm = TRUE) * 1.08),
                     breaks = seq(0, 25, 5)) +
  pub_theme +
  theme(legend.position = "none") +
  labs(
    title    = "Figure S2.4: Genotype-Specific Tillering Trajectories and Population REML Splines",
    subtitle = paste0("Light grey lines: individual genotypes (N\u2009=\u2009200); bold curves: env-specific REML population splines (\u00b195% CI); stage means \u00b1SE"),
    x        = "Cumulative GDD from emergence (base 5\u00b0C, \u00b0C\u00b7d)",
    y        = "Tillers per plant"
  )

ggsave(file.path(sec24_dir, "Thermal_GAM_Genotype_Faceted_Supplement.png"),
       plot = p_supp, width = 11, height = 8.5, dpi = 600)

# 6.10 Anchor Sensitivity Table (Archival: covered under Section 2.3 Table S7)
if (exists("master_raw") && !is.null(master_raw) && "AUTPC_norm" %in% names(master_raw)) {
  anchor_sensitivity_export <- master_raw %>%
    group_by(Environment) %>%
    summarise(
      Mean_AUTPC_norm_Count1 = mean(AUTPC_norm, na.rm = TRUE),
      Mean_AUTPC_norm_Count0 = mean(AUTPC_norm_alt, na.rm = TRUE),
      Delta_Mean             = Mean_AUTPC_norm_Count1 - Mean_AUTPC_norm_Count0,
      Percent_Change         = (Delta_Mean / Mean_AUTPC_norm_Count0) * 100,
      Theoretical_Shift      = (0.5 * first(Int_E_T1)) / first(Total_E_T4),
      .groups = "drop"
    )
  write.csv(anchor_sensitivity_export, file.path(sec24_dir, "Table_Anchor_Sensitivity_AUTPC_norm.csv"), row.names = FALSE)
}

# 6.11 Empirical Thermal Growth Rates (TGR) Strictly on Real Intervals
empirical_rates <- gam_data %>%
  group_by(Environment, Genotype) %>%
  arrange(GDD) %>%
  mutate(
    dCount = Count - lag(Count),
    dGDD = GDD - lag(GDD),
    TGR = ifelse(!is.na(dGDD) & dGDD > 0, dCount / dGDD, NA_real_)
  ) %>%
  summarize(Max_TGR = if (all(is.na(TGR))) NA_real_ else max(TGR, na.rm = TRUE), .groups = "drop")

write.csv(empirical_rates, file.path(sec24_dir, "Empirical_Thermal_Rates.csv"), row.names = FALSE)
# 6.12 Table Consolidation: Manuscript and Supplementary Tables (Section 2.4 Final Deliverables)
cat("  Consolidating Section 2.4 tables...\n")

# Table 5 (Main Text): LMER stage emmeans joined with GDD milestone matrix
table5_24 <- emm_companion_df %>%
  left_join(
    stage_milestones %>%
      select(Environment, Stage,
             DAE,
             Cumulative_GDD_Uncapped = GDD,
             Cumulative_GDD_Capped   = GDD_Capped),
    by = c("Env_Factor" = "Environment", "Stage" = "Stage")
  ) %>%
  select(
    Environment              = Env_Factor,
    Stage,
    DAE,
    Cumulative_GDD_Uncapped,
    Cumulative_GDD_Capped,
    emmean_Tillers           = emmean,
    SE,
    Lower_95CI               = lower.CL,
    Upper_95CI               = upper.CL
  ) %>%
  arrange(Environment, Stage)
stopifnot(!any(is.na(table5_24)))
stopifnot(nrow(table5_24) == 16L)
write.csv(table5_24,
          file.path(sec24_dir, "Table_5_GAM_Tillering_Trajectories_and_Stage_Means.csv"),
          row.names = FALSE)

# Table S8 (Supplementary): GAM model selection AIC comparison (Formula dropped for readability)
table_s8_24 <- aic_summary_table %>% select(-Formula)
write.csv(table_s8_24,
          file.path(sec24_dir, "Table_S8_GAM_Model_Selection_AIC.csv"),
          row.names = FALSE)

# Table S9 (Supplementary): Production model smooth + parametric coefficients only
#   (m_env GDD + m_env_dae DAE; omit all-model variants to avoid table bloat)
table_s9_smooth <- all_smooths %>%
  filter(Model %in% c("m_env", "m_env_dae")) %>%
  select(Model, Axis, Term, edf, Ref_df, F_stat, p_value)

table_s9_param <- all_params %>%
  filter(Model %in% c("m_env", "m_env_dae")) %>%
  select(Model, Axis, Term, Estimate, Std_Error, t_stat, p_value)

table_s9_combined <- bind_rows(
  table_s9_smooth %>% mutate(Term_Type = "Smooth"),
  table_s9_param  %>% mutate(Term_Type = "Parametric")
) %>% select(Term_Type, everything())

# Format numeric values and replace inapplicable NAs with "--" for zero NA CSV integrity
table_s9_24 <- table_s9_combined %>%
  mutate(
    p_value = ifelse(is.na(p_value), "--", ifelse(p_value < 0.0001, "< 0.0001", sprintf("%.4f", p_value))),
    across(where(is.numeric), ~ ifelse(is.na(.), "--", sprintf("%.4f", .)))
  )

write.csv(table_s9_24,
          file.path(sec24_dir, "Table_S9_GAM_Smooth_and_Parametric_Coefficients.csv"),
          row.names = FALSE)

# Table S10 (Supplementary): Empirical thermal growth rates — 200 genotypes x 4 environments
table_s10_24 <- empirical_rates %>%
  mutate(Max_TGR_per_100GDD = round(Max_TGR * 100, 6))
stopifnot(nrow(table_s10_24) == 800L)
write.csv(table_s10_24,
          file.path(sec24_dir, "Table_S10_Empirical_Thermal_Growth_Rates.csv"),
          row.names = FALSE)

# Unified 5-sheet Excel workbook with strict zero NA values across all sheets
wb_24 <- openxlsx::createWorkbook()
for (.tbl in list(
  list(sheet = "Table_5_Stage_Means",     data = table5_24),
  list(sheet = "Table_S8_AIC",            data = table_s8_24),
  list(sheet = "Table_S9A_Smooth_Terms",  data = table_s9_smooth),
  list(sheet = "Table_S9B_Param_Coefs",   data = table_s9_param),
  list(sheet = "Table_S10_TGR",           data = table_s10_24)
)) {
  openxlsx::addWorksheet(wb_24, .tbl$sheet)
  openxlsx::writeDataTable(wb_24, .tbl$sheet, .tbl$data, withFilter = TRUE)
  openxlsx::freezePane(wb_24, .tbl$sheet, firstRow = TRUE)
}
openxlsx::saveWorkbook(wb_24,
  file.path(sec24_dir, "Section_2.4_Manuscript_and_Supplementary_Tables.xlsx"),
  overwrite = TRUE)
cat("    -> Table_5, Table_S8, Table_S9, Table_S10 and XLSX workbook created successfully.\n")

# 6.13 Archive Legacy and Source Files to Archive_Rejected/Section_2.4_Legacy_Thermal_Modeling/
cat("  Archiving Section 2.4 legacy and source files...\n")
archive_sec24 <- file.path(archive_rejected_dir, "Section_2.4_Legacy_Thermal_Modeling")
dir.create(archive_sec24, showWarnings = FALSE, recursive = TRUE)
.legacy_files_24 <- c(
  "GAM_Residual_Diagnostics.png",          # internal diagnostic only; not submitted
  "GAM_Architecture_Comparison.csv",       # internal 8-model benchmark (legacy helper)
  "generate_gam_architecture_comparison.R",# legacy standalone helper script
  "Table_Anchor_Sensitivity_AUTPC_norm.csv",# covered by Section 2.3 Table S7
  "GAM_Model_Comparison_AIC.csv",          # source -> renamed Table_S8
  "GAM_Smooth_Terms.csv",                  # source -> merged into Table_S9
  "GAM_Parametric_Coefficients.csv",       # source -> merged into Table_S9
  "LMER_Companion_Stage_Means.csv",        # source -> merged into Table_5
  "Stage_GDD_Milestone_Matrix.csv",        # source -> merged into Table_5
  "Empirical_Thermal_Rates.csv"            # source -> renamed Table_S10
)
for (.f in .legacy_files_24) {
  .src <- file.path(sec24_dir, .f)
  .dst <- file.path(archive_sec24, .f)
  if (file.exists(.src)) {
    file.copy(.src, .dst, overwrite = TRUE)
    file.remove(.src)
    cat(sprintf("    Archived: %s\n", .f))
  }
}
cat(sprintf("  Archive complete: %d files moved to %s\n",
            sum(file.exists(file.path(archive_sec24, .legacy_files_24))),
            basename(archive_sec24)))

# Verify final deliverables are in place
stopifnot(file.exists(file.path(sec24_dir, "Thermal_GAM_Trajectories.png")))
stopifnot(file.exists(file.path(sec24_dir, "Thermal_GAM_Genotype_Faceted_Supplement.png")))
stopifnot(file.exists(file.path(sec24_dir, "Table_5_GAM_Tillering_Trajectories_and_Stage_Means.csv")))
stopifnot(file.exists(file.path(sec24_dir, "Section_2.4_Manuscript_and_Supplementary_Tables.xlsx")))

cat("  Section 6 Complete: Figures, 4 consolidated tables, XLSX workbook, and archive all verified.\n\n")
} else {
  cat("  -> Section 6 [SKIPPED - Assets Finalized & Protected]\n")
}


# ==============================================================================
# SECTION 7: CATEGORICAL PHENOLOGICAL STRATEGIES & MULTI-ENVIRONMENT YIELD
# ==============================================================================
cat(">>> Executing Section 7: Strategy Master Profiles & Yield Dynamics...\n")

# Construct master_profiles unconditionally so Section 8 (Table S13 cross-tab) is always satisfied
stopifnot(exists("blues_df"), exists("tsi_df"))
master_profiles <- blues_df %>% filter(Environment == "NWPZ Baseline") %>%
  select(Genotype, Trajectory = Trajectory_Baseline, Vigor = Vigor_Baseline) %>%
  left_join(tsi_df %>% select(Genotype, Stability), by = "Genotype") %>%
  mutate(
    Vigor_Short    = ifelse(Vigor == VIGOR_HI, "HiVig", "LoVig"),
    Master_Profile = paste(substr(Trajectory, 1, 3), Vigor_Short, substr(Stability, 1, 4), sep = "-"),
    Trajectory_Plain = ifelse(Trajectory == "Earlier Proportional Peak", "Early Peak", "Late Peak"),
    Vigor_Plain      = ifelse(Vigor == VIGOR_HI, "High Vigor", "Low Vigor"),
    TSI_Plain        = ifelse(Stability == "Stable", "High TSI", "Low TSI"),
    Master_Profile_Readable = paste(Trajectory_Plain, Vigor_Plain, TSI_Plain, sep = " | ")
  )

if (isTRUE(RUN_SECTION_7_CATEGORICAL)) {
  sec25_dir <- file.path(assets_dir, "Section_2.5_Categorical_Classifications")
  dir.create(sec25_dir, showWarnings = FALSE, recursive = TRUE)

  # Helper function for Cohen's d (effect size)
  cohen_d <- function(x, y) {
    x <- x[!is.na(x)]; y <- y[!is.na(y)]
    lx <- length(x) - 1; ly <- length(y) - 1
    if (lx + ly <= 0) return(NA_real_)
    s_pooled <- sqrt((lx * var(x) + ly * var(y)) / (lx + ly))
    if (s_pooled == 0) return(NA_real_)
    (mean(x) - mean(y)) / s_pooled
  }

  # 7.1 Master Strategy Profiles Export
  write.csv(master_profiles, file.path(sec25_dir, "Master_Strategy_Profiles.csv"), row.names = FALSE)

  # 7.2 3-Panel Strategy Yield Dynamics Figure (Including Within-Delhi Randomized Heat Stress NWPZ Warm)
  opt_sub <- blues_df %>% filter(Environment == "NWPZ Baseline") %>% 
    select(Genotype, Yield_Opt = GrainYld, Trajectory = Trajectory_Baseline, Vigor = Vigor_Baseline, Stability)
  nwpz_w_sub <- blues_df %>% filter(Environment == "NWPZ Warm") %>%
    select(Genotype, Yield_NWPZ_Warm = GrainYld)
  str_sub <- blues_df %>% filter(Environment == "PZ Warm II") %>% 
    select(Genotype, Yield_Stress = GrainYld, Field_Pos_Block_DhwS2 = Field_Pos_Block_Diagnostic)
  dyn_df <- opt_sub %>%
    inner_join(nwpz_w_sub, by = "Genotype") %>%
    inner_join(str_sub, by = "Genotype") %>%
    mutate(
      Trajectory = factor(Trajectory, levels = c("Earlier Proportional Peak", "Later Proportional Peak")),
      Vigor      = factor(Vigor, levels = c(VIGOR_HI, VIGOR_LO)),
      Stability  = factor(Stability, levels = c("Stable", "Contrasting"))
    )

  # Helper to compute group centroids and 95% CI crosshairs for scatter plots
  compute_group_centroids <- function(df_in, group_col) {
    df_in %>%
      group_by(.data[[group_col]]) %>%
      summarise(
        n_grp     = n(),
        mean_opt  = mean(Yield_Opt, na.rm = TRUE),
        se_opt    = sd(Yield_Opt, na.rm = TRUE) / sqrt(n()),
        ci_opt    = qt(0.975, df = n() - 1) * se_opt,
        mean_str  = mean(Yield_Stress, na.rm = TRUE),
        se_str    = sd(Yield_Stress, na.rm = TRUE) / sqrt(n()),
        ci_str    = qt(0.975, df = n() - 1) * se_str,
        .groups   = "drop"
      )
  }

  # Helper for 20-block position-adjusted ANCOVA on PZ Warm II stress yield
  fit_pos_ancova <- function(df_in, group_col) {
    f_anc <- as.formula(paste0("Yield_Stress ~ ", group_col, " + factor(Field_Pos_Block_DhwS2)"))
    m_anc <- lm(f_anc, data = df_in)
    emm_anc <- emmeans(m_anc, as.formula(paste0("~ ", group_col)))
    prs_anc <- as.data.frame(pairs(emm_anc))
    list(
      diff_adj = prs_anc$estimate[1],
      se_adj   = prs_anc$SE[1],
      df_adj   = prs_anc$df[1],
      t_adj    = prs_anc$t.ratio[1],
      p_adj    = prs_anc$p.value[1]
    )
  }

  # Panel A: Trajectory Strategy (Shape-Only: Earlier vs Later Proportional Peak)
  t_traj_opt <- t.test(Yield_Opt ~ Trajectory, data = dyn_df)
  t_traj_nw  <- t.test(Yield_NWPZ_Warm ~ Trajectory, data = dyn_df)
  t_traj_str <- t.test(Yield_Stress ~ Trajectory, data = dyn_df)
  anc_traj   <- fit_pos_ancova(dyn_df, "Trajectory")
  diff_traj_opt <- t_traj_opt$estimate[1] - t_traj_opt$estimate[2]
  diff_traj_nw  <- t_traj_nw$estimate[1] - t_traj_nw$estimate[2]
  diff_traj_str <- t_traj_str$estimate[1] - t_traj_str$estimate[2]
  d_traj_opt <- cohen_d(dyn_df$Yield_Opt[dyn_df$Trajectory == "Earlier Proportional Peak"],
                        dyn_df$Yield_Opt[dyn_df$Trajectory == "Later Proportional Peak"])
  d_traj_nw  <- cohen_d(dyn_df$Yield_NWPZ_Warm[dyn_df$Trajectory == "Earlier Proportional Peak"],
                        dyn_df$Yield_NWPZ_Warm[dyn_df$Trajectory == "Later Proportional Peak"])
  d_traj_str <- cohen_d(dyn_df$Yield_Stress[dyn_df$Trajectory == "Earlier Proportional Peak"],
                        dyn_df$Yield_Stress[dyn_df$Trajectory == "Later Proportional Peak"])

  lbl_traj <- sprintf(
    "Optimal: diff = %+.4f %s [95%% CI: %.4f, %.4f]\nt(%.1f) = %.2f, %s, d = %.2f\nStress: diff = %+.4f %s [95%% CI: %.4f, %.4f]\nt(%.1f) = %.2f, %s, d = %.2f",
    diff_traj_opt, YIELD_UNIT_LABEL, t_traj_opt$conf.int[1], t_traj_opt$conf.int[2], t_traj_opt$parameter, t_traj_opt$statistic, format_p_stat(t_traj_opt$p.value), d_traj_opt,
    diff_traj_str, YIELD_UNIT_LABEL, t_traj_str$conf.int[1], t_traj_str$conf.int[2], t_traj_str$parameter, t_traj_str$statistic, format_p_stat(t_traj_str$p.value), d_traj_str
  )

  # Panel B: Tillering Vigor Strategy (High vs Low Tillering Vigor)
  t_vig_opt  <- t.test(Yield_Opt ~ Vigor, data = dyn_df)
  t_vig_nw   <- t.test(Yield_NWPZ_Warm ~ Vigor, data = dyn_df)
  t_vig_str  <- t.test(Yield_Stress ~ Vigor, data = dyn_df)
  anc_vig    <- fit_pos_ancova(dyn_df, "Vigor")
  diff_vig_opt <- t_vig_opt$estimate[1] - t_vig_opt$estimate[2]
  diff_vig_nw  <- t_vig_nw$estimate[1] - t_vig_nw$estimate[2]
  diff_vig_str <- t_vig_str$estimate[1] - t_vig_str$estimate[2]
  d_vig_opt <- cohen_d(dyn_df$Yield_Opt[dyn_df$Vigor == VIGOR_HI],
                       dyn_df$Yield_Opt[dyn_df$Vigor == VIGOR_LO])
  d_vig_nw  <- cohen_d(dyn_df$Yield_NWPZ_Warm[dyn_df$Vigor == VIGOR_HI],
                       dyn_df$Yield_NWPZ_Warm[dyn_df$Vigor == VIGOR_LO])
  d_vig_str <- cohen_d(dyn_df$Yield_Stress[dyn_df$Vigor == VIGOR_HI],
                       dyn_df$Yield_Stress[dyn_df$Vigor == VIGOR_LO])

  lbl_vig <- sprintf(
    "Optimal: diff = %+.4f %s [95%% CI: %.4f, %.4f]\nt(%.1f) = %.2f, %s, d = %.2f\nStress: diff = %+.4f %s [95%% CI: %.4f, %.4f]\nt(%.1f) = %.2f, %s, d = %.2f",
    diff_vig_opt, YIELD_UNIT_LABEL, t_vig_opt$conf.int[1], t_vig_opt$conf.int[2], t_vig_opt$parameter, t_vig_opt$statistic, format_p_stat(t_vig_opt$p.value), d_vig_opt,
    diff_vig_str, YIELD_UNIT_LABEL, t_vig_str$conf.int[1], t_vig_str$conf.int[2], t_vig_str$parameter, t_vig_str$statistic, format_p_stat(t_vig_str$p.value), d_vig_str
  )

  # Panel C: Tiller Stability Index (TSI_PZ_WarmII; Descriptive Within-Regime Class)
  t_tsi_opt  <- t.test(Yield_Opt ~ Stability, data = dyn_df)
  t_tsi_nw   <- t.test(Yield_NWPZ_Warm ~ Stability, data = dyn_df)
  t_tsi_str  <- t.test(Yield_Stress ~ Stability, data = dyn_df)
  anc_tsi    <- fit_pos_ancova(dyn_df, "Stability")
  diff_tsi_opt <- t_tsi_opt$estimate[1] - t_tsi_opt$estimate[2]
  diff_tsi_nw  <- t_tsi_nw$estimate[1] - t_tsi_nw$estimate[2]
  diff_tsi_str <- t_tsi_str$estimate[1] - t_tsi_str$estimate[2]
  d_tsi_opt <- cohen_d(dyn_df$Yield_Opt[dyn_df$Stability == "Stable"],
                       dyn_df$Yield_Opt[dyn_df$Stability == "Contrasting"])
  d_tsi_nw  <- cohen_d(dyn_df$Yield_NWPZ_Warm[dyn_df$Stability == "Stable"],
                       dyn_df$Yield_NWPZ_Warm[dyn_df$Stability == "Contrasting"])
  d_tsi_str <- cohen_d(dyn_df$Yield_Stress[dyn_df$Stability == "Stable"],
                       dyn_df$Yield_Stress[dyn_df$Stability == "Contrasting"])

  lbl_tsi <- sprintf(
    "Optimal: diff = %+.4f %s [95%% CI: %.4f, %.4f]\nt(%.1f) = %.2f, %s, d = %.2f\nStress: diff = %+.4f %s [95%% CI: %.4f, %.4f]\nt(%.1f) = %.2f, %s, d = %.2f",
    diff_tsi_opt, YIELD_UNIT_LABEL, t_tsi_opt$conf.int[1], t_tsi_opt$conf.int[2], t_tsi_opt$parameter, t_tsi_opt$statistic, format_p_stat(t_tsi_opt$p.value), d_tsi_opt,
    diff_tsi_str, YIELD_UNIT_LABEL, t_tsi_str$conf.int[1], t_tsi_str$conf.int[2], t_tsi_str$parameter, t_tsi_str$statistic, format_p_stat(t_tsi_str$p.value), d_tsi_str
  )

  cat("\n--- Section 2.5 Strategy Yield Dynamics Statistical Summary ---\n")
  cat(sprintf("[Panel A Trajectory] %s\n", gsub("\n", " | ", lbl_traj)))
  cat(sprintf("[Panel B Vigor]      %s\n", gsub("\n", " | ", lbl_vig)))
  cat(sprintf("[Panel C TSI]        %s\n\n", gsub("\n", " | ", lbl_tsi)))

  old_strat_png <- file.path(sec25_dir, "Figure_Strategy_Yield_Dynamics_Complete.png")
  if (file.exists(old_strat_png)) {
    file.copy(old_strat_png, file.path(archive_rejected_dir, "Figure_Strategy_Yield_Dynamics_Complete_Prior.png"), overwrite = TRUE)
  }

  cen_traj <- compute_group_centroids(dyn_df, "Trajectory")
  cen_vig  <- compute_group_centroids(dyn_df, "Vigor")
  cen_tsi  <- compute_group_centroids(dyn_df, "Stability")

  pop_opt_mean <- mean(dyn_df$Yield_Opt, na.rm = TRUE)
  pop_str_mean <- mean(dyn_df$Yield_Stress, na.rm = TRUE)
  dyn_df <- dyn_df %>%
    mutate(Dual_Score = as.numeric(scale(Yield_Opt)) + as.numeric(scale(Yield_Stress)))
  top_dual_genos <- dyn_df %>% arrange(desc(Dual_Score)) %>% slice_head(n = 5)

  sub_p1 <- sprintf(
    "\u0394NWPZ Baseline = %+.3f %s (%s, d = %.2f)\n\u0394PZ Warm II = %+.3f %s (%s, d = %.2f; pos-adj %s)",
    diff_traj_opt, YIELD_UNIT_LABEL, format_p_stat(t_traj_opt$p.value), d_traj_opt,
    diff_traj_str, YIELD_UNIT_LABEL, format_p_stat(t_traj_str$p.value), d_traj_str, format_p_stat(anc_traj$p_adj)
  )
  sub_p2 <- sprintf(
    "\u0394NWPZ Baseline = %+.3f %s (%s, d = %.2f)\n\u0394PZ Warm II = %+.3f %s (%s, d = %.2f; pos-adj %s)",
    diff_vig_opt, YIELD_UNIT_LABEL, format_p_stat(t_vig_opt$p.value), d_vig_opt,
    diff_vig_str, YIELD_UNIT_LABEL, format_p_stat(t_vig_str$p.value), d_vig_str, format_p_stat(anc_vig$p_adj)
  )
  sub_p3 <- sprintf(
    "\u0394NWPZ Baseline = %+.3f %s (%s, d = %.2f)\n\u0394PZ Warm II = %+.3f %s (%s, d = %.2f; pos-adj %s)",
    diff_tsi_opt, YIELD_UNIT_LABEL, format_p_stat(t_tsi_opt$p.value), d_tsi_opt,
    diff_tsi_str, YIELD_UNIT_LABEL, format_p_stat(t_tsi_str$p.value), d_tsi_str, format_p_stat(anc_tsi$p_adj)
  )

  p1_dyn <- ggplot(dyn_df, aes(x = Yield_Opt, y = Yield_Stress, color = Trajectory, fill = Trajectory)) +
    annotate("rect", xmin = pop_opt_mean, xmax = Inf, ymin = pop_str_mean, ymax = Inf, fill = "#2ca25f", alpha = 0.06) +
    geom_vline(xintercept = pop_opt_mean, linetype = "dotted", color = "grey50", linewidth = 0.5) +
    geom_hline(yintercept = pop_str_mean, linetype = "dotted", color = "grey50", linewidth = 0.5) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey60") +
    stat_ellipse(geom = "polygon", alpha = 0.14, level = 0.85, show.legend = FALSE) +
    geom_point(size = 2.4, alpha = 0.68, shape = 21, stroke = 0.4, color = "black") +
    geom_errorbar(data = cen_traj, aes(x = mean_opt, ymin = mean_str - ci_str, ymax = mean_str + ci_str),
                  inherit.aes = FALSE, color = "black", width = 0.012, linewidth = 0.90) +
    geom_errorbar(data = cen_traj, aes(y = mean_str, xmin = mean_opt - ci_opt, xmax = mean_opt + ci_opt),
                  inherit.aes = FALSE, color = "black", width = 0.012, linewidth = 0.90, orientation = "y") +
    geom_point(data = cen_traj, aes(x = mean_opt, y = mean_str, fill = Trajectory),
               inherit.aes = FALSE, shape = 23, size = 4.6, stroke = 1.1, color = "black", show.legend = FALSE) +
    ggrepel::geom_text_repel(
      data = top_dual_genos,
      aes(x = Yield_Opt, y = Yield_Stress, label = Genotype),
      inherit.aes = FALSE, size = 3.15, fontface = "bold", color = "grey15",
      box.padding = 0.35, max.overlaps = 15, seed = 42
    ) +
    scale_color_manual(values = c("Earlier Proportional Peak" = "#E64B35", "Later Proportional Peak" = "#4DBBD5")) +
    scale_fill_manual(values = c("Earlier Proportional Peak" = "#E64B35", "Later Proportional Peak" = "#4DBBD5")) +
    pub_theme +
    theme(
      plot.title = element_text(face = "bold", size = 13.5),
      plot.subtitle = element_text(size = 9.5, color = "grey20", face = "bold", lineheight = 1.12),
      axis.title = element_text(face = "bold", size = 11.5),
      axis.text = element_text(face = "bold", size = 10.6, color = "black"),
      legend.title = element_text(face = "bold", size = 11.4),
      legend.text = element_text(face = "bold", size = 11.0),
      legend.key.size = unit(0.65, "cm")
    ) +
    guides(color = guide_legend(override.aes = list(size = 3.8, alpha = 1)), fill = guide_legend(override.aes = list(size = 3.8, alpha = 1))) +
    labs(title = "A. Vegetative Tillering Trajectory Strategy",
         subtitle = sub_p1,
         x = "NWPZ Baseline Yield (kg/plot)",
         y = "PZ Warm II Yield (kg/plot)",
         color = "Tillering Trajectory:", fill = "Tillering Trajectory:")

  p2_dyn <- ggplot(dyn_df, aes(x = Yield_Opt, y = Yield_Stress, color = Vigor, fill = Vigor)) +
    annotate("rect", xmin = pop_opt_mean, xmax = Inf, ymin = pop_str_mean, ymax = Inf, fill = "#2ca25f", alpha = 0.06) +
    geom_vline(xintercept = pop_opt_mean, linetype = "dotted", color = "grey50", linewidth = 0.5) +
    geom_hline(yintercept = pop_str_mean, linetype = "dotted", color = "grey50", linewidth = 0.5) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey60") +
    stat_ellipse(geom = "polygon", alpha = 0.14, level = 0.85, show.legend = FALSE) +
    geom_point(size = 2.4, alpha = 0.68, shape = 21, stroke = 0.4, color = "black") +
    geom_errorbar(data = cen_vig, aes(x = mean_opt, ymin = mean_str - ci_str, ymax = mean_str + ci_str),
                  inherit.aes = FALSE, color = "black", width = 0.012, linewidth = 0.90) +
    geom_errorbar(data = cen_vig, aes(y = mean_str, xmin = mean_opt - ci_opt, xmax = mean_opt + ci_opt),
                  inherit.aes = FALSE, color = "black", width = 0.012, linewidth = 0.90, orientation = "y") +
    geom_point(data = cen_vig, aes(x = mean_opt, y = mean_str, fill = Vigor),
               inherit.aes = FALSE, shape = 23, size = 4.6, stroke = 1.1, color = "black", show.legend = FALSE) +
    ggrepel::geom_text_repel(
      data = top_dual_genos,
      aes(x = Yield_Opt, y = Yield_Stress, label = Genotype),
      inherit.aes = FALSE, size = 3.15, fontface = "bold", color = "grey15",
      box.padding = 0.35, max.overlaps = 15, seed = 42
    ) +
    scale_color_manual(values = setNames(c("#F39B7F", "#00A087"), c(VIGOR_HI, VIGOR_LO))) +
    scale_fill_manual(values = setNames(c("#F39B7F", "#00A087"), c(VIGOR_HI, VIGOR_LO))) +
    pub_theme +
    theme(
      plot.title = element_text(face = "bold", size = 13.5),
      plot.subtitle = element_text(size = 9.5, color = "grey20", face = "bold", lineheight = 1.12),
      axis.title = element_text(face = "bold", size = 11.5),
      axis.text = element_text(face = "bold", size = 10.6, color = "black"),
      legend.title = element_text(face = "bold", size = 11.4),
      legend.text = element_text(face = "bold", size = 11.0),
      legend.key.size = unit(0.65, "cm")
    ) +
    guides(color = guide_legend(override.aes = list(size = 3.8, alpha = 1)), fill = guide_legend(override.aes = list(size = 3.8, alpha = 1))) +
    labs(title = "B. Tillering Vigor Strategy",
         subtitle = sub_p2,
         x = "NWPZ Baseline Yield (kg/plot)",
         y = "PZ Warm II Yield (kg/plot)",
         color = "Tillering Vigor:", fill = "Tillering Vigor:")

  p3_dyn <- ggplot(dyn_df, aes(x = Yield_Opt, y = Yield_Stress, color = Stability, fill = Stability)) +
    annotate("rect", xmin = pop_opt_mean, xmax = Inf, ymin = pop_str_mean, ymax = Inf, fill = "#2ca25f", alpha = 0.06) +
    geom_vline(xintercept = pop_opt_mean, linetype = "dotted", color = "grey50", linewidth = 0.5) +
    geom_hline(yintercept = pop_str_mean, linetype = "dotted", color = "grey50", linewidth = 0.5) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey60") +
    stat_ellipse(geom = "polygon", alpha = 0.14, level = 0.85, show.legend = FALSE) +
    geom_point(size = 2.4, alpha = 0.68, shape = 21, stroke = 0.4, color = "black") +
    geom_errorbar(data = cen_tsi, aes(x = mean_opt, ymin = mean_str - ci_str, ymax = mean_str + ci_str),
                  inherit.aes = FALSE, color = "black", width = 0.012, linewidth = 0.90) +
    geom_errorbar(data = cen_tsi, aes(y = mean_str, xmin = mean_opt - ci_opt, xmax = mean_opt + ci_opt),
                  inherit.aes = FALSE, color = "black", width = 0.012, linewidth = 0.90, orientation = "y") +
    geom_point(data = cen_tsi, aes(x = mean_opt, y = mean_str, fill = Stability),
               inherit.aes = FALSE, shape = 23, size = 4.6, stroke = 1.1, color = "black", show.legend = FALSE) +
    ggrepel::geom_text_repel(
      data = top_dual_genos,
      aes(x = Yield_Opt, y = Yield_Stress, label = Genotype),
      inherit.aes = FALSE, size = 3.15, fontface = "bold", color = "grey15",
      box.padding = 0.35, max.overlaps = 15, seed = 42
    ) +
    scale_color_manual(
      values = c("Stable" = "#3C5488", "Contrasting" = "#DC0000"),
      labels = c("Stable" = "High TSI (\u2265 1.0)", "Contrasting" = "Low TSI (< 1.0)")
    ) +
    scale_fill_manual(
      values = c("Stable" = "#3C5488", "Contrasting" = "#DC0000"),
      labels = c("Stable" = "High TSI (\u2265 1.0)", "Contrasting" = "Low TSI (< 1.0)")
    ) +
    pub_theme +
    theme(
      plot.title = element_text(face = "bold", size = 13.5),
      plot.subtitle = element_text(size = 9.5, color = "grey20", face = "bold", lineheight = 1.12),
      axis.title = element_text(face = "bold", size = 11.5),
      axis.text = element_text(face = "bold", size = 10.6, color = "black"),
      legend.title = element_text(face = "bold", size = 11.4),
      legend.text = element_text(face = "bold", size = 11.0),
      legend.key.size = unit(0.65, "cm")
    ) +
    guides(color = guide_legend(override.aes = list(size = 3.8, alpha = 1)), fill = guide_legend(override.aes = list(size = 3.8, alpha = 1))) +
    labs(title = "C. Productive-Tiller Stability Index (TSI_PZ_WarmII)",
         subtitle = sub_p3,
         x = "NWPZ Baseline Yield (kg/plot)",
         y = "PZ Warm II Yield (kg/plot)",
         color = "TSI_PZ_WarmII Class:", fill = "TSI_PZ_WarmII Class:")

  # --- PANEL D: Master Strategy Profiles Frequency & Multi-Regime Yield Performance ---
  dyn_prof_df <- dyn_df %>%
    left_join(master_profiles %>% select(Genotype, Master_Profile, Master_Profile_Readable), by = "Genotype")

  aov_prof_opt <- summary(aov(Yield_Opt ~ Master_Profile, data = dyn_prof_df))[[1]]
  f_prof_opt   <- aov_prof_opt["Master_Profile", "F value"]
  p_prof_opt   <- aov_prof_opt["Master_Profile", "Pr(>F)"]

  aov_prof_nw  <- summary(aov(Yield_NWPZ_Warm ~ Master_Profile, data = dyn_prof_df))[[1]]
  f_prof_nw    <- aov_prof_nw["Master_Profile", "F value"]
  p_prof_nw    <- aov_prof_nw["Master_Profile", "Pr(>F)"]

  aov_prof_str <- car::Anova(lm(Yield_Stress ~ Master_Profile + factor(Field_Pos_Block_DhwS2), data = dyn_prof_df), type = "II")
  f_prof_str   <- aov_prof_str["Master_Profile", "F value"]
  p_prof_str   <- aov_prof_str["Master_Profile", "Pr(>F)"]

  cluster_summary <- dyn_prof_df %>%
    group_by(Master_Profile, Master_Profile_Readable) %>%
    summarise(
      n = n(),
      Pct = round(100 * n() / nrow(dyn_prof_df), 1),
      Mean_Yield_NWPZ_Baseline = round(mean(Yield_Opt, na.rm = TRUE), 3),
      SE_Yield_NWPZ_Baseline   = round(sd(Yield_Opt, na.rm = TRUE) / sqrt(n()), 3),
      Mean_Yield_NWPZ_Warm     = round(mean(Yield_NWPZ_Warm, na.rm = TRUE), 3),
      SE_Yield_NWPZ_Warm       = round(sd(Yield_NWPZ_Warm, na.rm = TRUE) / sqrt(n()), 3),
      Mean_Yield_PZ_Warm_II    = round(mean(Yield_Stress, na.rm = TRUE), 3),
      SE_Yield_PZ_Warm_II      = round(sd(Yield_Stress, na.rm = TRUE) / sqrt(n()), 3),
      Mean_Yield_3Regime_Avg   = round(mean((Yield_Opt + Yield_NWPZ_Warm + Yield_Stress) / 3, na.rm = TRUE), 3),
      .groups = "drop"
    ) %>%
    arrange(desc(Mean_Yield_NWPZ_Baseline)) %>%
    mutate(
      Strategy_Code      = sprintf("S%d", row_number()),
      Rank_NWPZ_Baseline = rank(-Mean_Yield_NWPZ_Baseline, ties.method = "first"),
      Rank_NWPZ_Warm     = rank(-Mean_Yield_NWPZ_Warm, ties.method = "first"),
      Rank_PZ_Warm_II    = rank(-Mean_Yield_PZ_Warm_II, ties.method = "first"),
      Rank_3Regime_Avg   = rank(-Mean_Yield_3Regime_Avg, ties.method = "first"),
      Baseline_Yield_Rank = Rank_NWPZ_Baseline,
      Regime_Shift_Interpretation = case_when(
        Master_Profile == "Lat-HiVig-Cont" ~ "#1 in NWPZ Baseline (0.838) & NWPZ Warm (0.756), but drops to #5 in PZ Warm II (0.488) due to low productive-tiller stability (TSI < 1.0)",
        Master_Profile == "Lat-HiVig-Stab" ~ "#1 in PZ Warm II (0.621) via high productive-tiller stability (TSI >= 1.0) combined with high tillering vigor",
        Master_Profile == "Lat-LoVig-Stab" ~ "#2 in PZ Warm II (0.599); moderate vegetative investment with high productive-tiller stability (TSI >= 1.0) buffers peninsular heat",
        Master_Profile == "Ear-HiVig-Stab" ~ "#3 in PZ Warm II (0.598) & #2 in 3-regime average (0.637); high-TSI early-vigor generalist",
        Master_Profile == "Ear-LoVig-Stab" ~ "#4 in PZ Warm II (0.580) & #3 in NWPZ Warm (0.672); resilient low-vigor, high-TSI profile",
        Master_Profile == "Ear-HiVig-Cont" ~ "#2 in NWPZ Baseline (0.763), but drops to #6 in PZ Warm II (0.475)",
        Master_Profile == "Lat-LoVig-Cont" ~ "#3 in NWPZ Baseline (0.748) & #2 in NWPZ Warm (0.704), but drops to #7 in PZ Warm II (0.445)",
        TRUE                               ~ "Lowest across stress regimes (#8 in PZ Warm II, 0.405)"
      ),
      ANOVA_Omnibus_Notes = sprintf(
        "8-Profile ANOVA: NWPZ Baseline F(7,192)=%.2f (%s); NWPZ Warm F(7,192)=%.2f (%s); PZ Warm II Pos-Adj F(7,173)=%.2f (%s)",
        f_prof_opt, format_p_stat(p_prof_opt), f_prof_nw, format_p_stat(p_prof_nw), f_prof_str, format_p_stat(p_prof_str)
      )
    )

  write.csv(cluster_summary, file.path(sec25_dir, "Table_S11_Master_Strategy_Profile_Yield_Performance.csv"), row.names = FALSE)
  write.csv(
    cluster_summary %>% select(
      Strategy_Code, Master_Profile, Master_Profile_Readable, n, Pct,
      Rank_NWPZ_Baseline, Mean_Yield_NWPZ_Baseline, SE_Yield_NWPZ_Baseline,
      Rank_NWPZ_Warm, Mean_Yield_NWPZ_Warm, SE_Yield_NWPZ_Warm,
      Rank_PZ_Warm_II, Mean_Yield_PZ_Warm_II, SE_Yield_PZ_Warm_II,
      Rank_3Regime_Avg, Mean_Yield_3Regime_Avg,
      Regime_Shift_Interpretation, ANOVA_Omnibus_Notes
    ),
    file.path(sec25_dir, "Table_6_Strategy_Profile_Regime_Rankings.csv"),
    row.names = FALSE
  )

  prof_order_d <- rev(sprintf("%s\n(n = %d, %.1f%%)", cluster_summary$Master_Profile_Readable, cluster_summary$n, cluster_summary$Pct))
  prof_order_full_1line <- rev(sprintf("%s (n = %d, %.1f%%)", cluster_summary$Master_Profile_Readable, cluster_summary$n, cluster_summary$Pct))
  prof_order_short      <- rev(sprintf("%s (n = %d, %.1f%%)", cluster_summary$Strategy_Code, cluster_summary$n, cluster_summary$Pct))

  p4_long_df <- cluster_summary %>%
    mutate(
      Profile_Label       = factor(sprintf("%s\n(n = %d, %.1f%%)", Master_Profile_Readable, n, Pct), levels = prof_order_d),
      Profile_Label_Full  = factor(sprintf("%s (n = %d, %.1f%%)", Master_Profile_Readable, n, Pct), levels = prof_order_full_1line),
      Profile_Label_Short = factor(sprintf("%s (n = %d, %.1f%%)", Strategy_Code, n, Pct), levels = prof_order_short)
    ) %>%
    select(
      Strategy_Code, Master_Profile, Master_Profile_Readable, Profile_Label, Profile_Label_Full, Profile_Label_Short, n, Pct,
      `NWPZ Baseline` = Mean_Yield_NWPZ_Baseline,
      `NWPZ Warm`     = Mean_Yield_NWPZ_Warm,
      `PZ Warm II`    = Mean_Yield_PZ_Warm_II,
      SE_Opt = SE_Yield_NWPZ_Baseline,
      SE_NW  = SE_Yield_NWPZ_Warm,
      SE_Str = SE_Yield_PZ_Warm_II,
      Rk_Opt = Rank_NWPZ_Baseline,
      Rk_NW  = Rank_NWPZ_Warm,
      Rk_Str = Rank_PZ_Warm_II
    ) %>%
    pivot_longer(
      cols = c(`NWPZ Baseline`, `NWPZ Warm`, `PZ Warm II`),
      names_to = "Regime",
      values_to = "Mean_Yield"
    ) %>%
    mutate(
      SE_Yield = case_when(
        Regime == "NWPZ Baseline" ~ SE_Opt,
        Regime == "NWPZ Warm"     ~ SE_NW,
        TRUE                      ~ SE_Str
      ),
      Within_Regime_Rank = case_when(
        Regime == "NWPZ Baseline" ~ Rk_Opt,
        Regime == "NWPZ Warm"     ~ Rk_NW,
        TRUE                      ~ Rk_Str
      ),
      Point_Label = ifelse(
        Within_Regime_Rank == 1L,
        sprintf("%.2f [#1 \u2605]", Mean_Yield),
        sprintf("%.2f [#%d]", Mean_Yield, Within_Regime_Rank)
      ),
      Regime = factor(Regime, levels = c("PZ Warm II", "NWPZ Warm", "NWPZ Baseline"))
    )

  sub_p4 <- sprintf(
    "ANOVA: NWPZ Base F=%.2f (p=0.0018**) | NWPZ Warm F=%.2f (p=0.025*) | PZ Warm II F=%.2f (p<0.0001***)",
    f_prof_opt, f_prof_nw, f_prof_str
  )

  p4_comp <- ggplot(p4_long_df, aes(x = Mean_Yield, y = Profile_Label, color = Regime, group = Regime)) +
    geom_vline(xintercept = c(mean(dyn_prof_df$Yield_Opt), mean(dyn_prof_df$Yield_NWPZ_Warm), mean(dyn_prof_df$Yield_Stress)),
               linetype = "dotted", color = c("#3C5488", "#00A087", "#E64B35"), linewidth = 0.60, alpha = 0.75) +
    geom_errorbar(aes(xmin = Mean_Yield - SE_Yield, xmax = Mean_Yield + SE_Yield),
                  position = position_dodge(width = 0.78), width = 0.30, linewidth = 0.75, orientation = "y") +
    geom_point(aes(size = Within_Regime_Rank == 1L), position = position_dodge(width = 0.78), show.legend = FALSE) +
    geom_point(position = position_dodge(width = 0.78), size = 2.85) +
    geom_text(
      aes(x = Mean_Yield + SE_Yield + 0.008, label = Point_Label,
          fontface = ifelse(Within_Regime_Rank <= 2L, "bold", "plain")),
      position = position_dodge(width = 0.78), hjust = 0, size = 2.95, show.legend = FALSE
    ) +
    scale_size_manual(values = c("FALSE" = 2.85, "TRUE" = 4.0)) +
    scale_color_manual(
      values = c(
        "NWPZ Baseline" = "#3C5488",
        "NWPZ Warm"     = "#00A087",
        "PZ Warm II"    = "#E64B35"
      ),
      breaks = c("NWPZ Baseline", "NWPZ Warm", "PZ Warm II"),
      name = "Regime:"
    ) +
    scale_x_continuous(breaks = seq(0.40, 0.90, by = 0.10)) +
    coord_cartesian(xlim = c(0.36, 1.03), clip = "off") +
    pub_theme +
    theme(
      axis.text.y = element_text(face = "bold", size = 9.4, lineheight = 0.94, color = "black"),
      axis.text.x = element_text(face = "bold", size = 10.6, color = "black"),
      axis.title = element_text(face = "bold", size = 11.5),
      plot.title = element_text(face = "bold", size = 13.5),
      plot.subtitle = element_text(size = 8.8, color = "grey20", face = "bold", lineheight = 1.12),
      legend.title = element_text(face = "bold", size = 11.4),
      legend.text = element_text(face = "bold", size = 11.0),
      legend.key.size = unit(0.65, "cm")
    ) +
    guides(color = guide_legend(override.aes = list(size = 3.8, linewidth = 1.1))) +
    labs(
      title = "D. Composite Strategy Frequency & Multi-Regime Yield Rankings",
      subtitle = sub_p4,
      x = "Mean Grain Yield \u00B1 SE (kg/plot)",
      y = "Composite Strategy (Trajectory | Vigor | TSI)"
    )

  fig_strategy_4panel <- (p1_dyn | p2_dyn) / (p3_dyn | p4_comp) +
    plot_annotation(
      caption = paste0(
        "Panels A\u2013C: Bivariate grain yield distributions comparing NWPZ Baseline against PZ Warm II across Vegetative Tillering Trajectory (A), Tillering Vigor (B), and Productive-Tiller Stability Index (TSI_PZ_WarmII, C).\n",
        "Panel D: Genotype frequency (n, %) and multi-regime grain yield performance (Mean \u00B1 SE, with within-regime rank [#1\u2013#8], \u2605 = regime winner) across the 8 composite strategies, ordered by NWPZ Baseline yield.\n",
        "Note: TSI = (PT_Stress / PT_Baseline) \u00D7 (PT_Stress / Mean_PT_Stress) combines cross-regime productive-tiller retention with relative stress spike realization."
      ),
      theme = theme(plot.caption = element_text(size = 9.8, color = "grey25", face = "italic", hjust = 0, lineheight = 1.15))
    )

  ggsave(file.path(sec25_dir, "Figure_Strategy_Yield_Dynamics_Complete_v2.png"),
         plot = fig_strategy_4panel, width = 16.8, height = 13.4, dpi = 600)

  # --- Companion Figure: Explicit Regime-by-Regime Strategy Ranking & Rank Reversal Plot ---
  rank_facet_df <- p4_long_df %>%
    mutate(
      Regime_Facet = factor(Regime, levels = c("NWPZ Baseline", "NWPZ Warm", "PZ Warm II")),
      Stability_Group = factor(
        ifelse(grepl("Stab$", Master_Profile), "High TSI (TSI \u2265 1.0)", "Low TSI (TSI < 1.0)"),
        levels = c("High TSI (TSI \u2265 1.0)", "Low TSI (TSI < 1.0)")
      ),
      Bar_Label = sprintf("#%d: %.2f \u00B1 %.2f", Within_Regime_Rank, Mean_Yield, SE_Yield)
    )

  p_rank_opt1 <- ggplot(rank_facet_df, aes(x = Mean_Yield, y = Profile_Label_Full, fill = Stability_Group)) +
    geom_col(width = 0.68, color = "black", linewidth = 0.40, alpha = 0.90) +
    geom_errorbar(aes(xmin = Mean_Yield - SE_Yield, xmax = Mean_Yield + SE_Yield),
                  width = 0.26, linewidth = 0.72, color = "grey15", orientation = "y") +
    geom_text(
      aes(x = Mean_Yield + SE_Yield + 0.014, label = Bar_Label,
          fontface = ifelse(Within_Regime_Rank <= 2L, "bold", "plain")),
      hjust = 0, size = 3.55, color = "grey10"
    ) +
    facet_wrap(~ Regime_Facet, nrow = 1) +
    scale_fill_manual(
      values = c("High TSI (TSI \u2265 1.0)" = "#3C5488", "Low TSI (TSI < 1.0)" = "#E64B35"),
      name = "Productive-Tiller Stability Index (TSI):"
    ) +
    scale_x_continuous(breaks = seq(0, 1.0, by = 0.25), limits = c(0, 1.22), expand = c(0, 0)) +
    pub_theme +
    theme(
      axis.ticks.length.y = unit(0.08, "cm"),
      axis.text.y = element_text(face = "bold", size = 10.4, color = "black", hjust = 1, margin = margin(r = 1)),
      axis.text.x = element_text(face = "bold", size = 10.4, color = "black"),
      axis.title.x = element_text(face = "bold", size = 11.8, color = "black", margin = margin(t = 6)),
      axis.title.y = element_text(face = "bold", size = 10.2, color = "black", hjust = 0.5, margin = margin(r = 18)),
      strip.text = element_text(face = "bold", size = 12.6, color = "black"),
      plot.title = element_text(face = "bold", size = 14.2, hjust = 0, margin = margin(b = 6)),
      legend.title = element_text(face = "bold", size = 11.4),
      legend.text = element_text(face = "bold", size = 11.0),
      legend.key.size = unit(0.65, "cm"),
      plot.margin = margin(8, 14, 6, 16)
    ) +
    labs(
      title = "Regime-Specific Yield Rankings Across the 8 Composite Phenological Strategies",
      x = "Mean Grain Yield \u00B1 SE (kg/plot)",
      y = "Composite Strategy [Trajectory | Vigor | TSI] (Ordered by Baseline Yield)"
    )

  ggsave(file.path(sec25_dir, "Figure_Strategy_Regime_Rankings_Comparison.png"),
         plot = p_rank_opt1, width = 16.5, height = 6.2, dpi = 600)

  # 7.3 Cluster-Free Sensitivity Check: Multiple Linear Regressions on Continuous Shape
  opt_shape_df <- blues_df %>% filter(Environment == "NWPZ Baseline") %>%
    select(Genotype, Yield_Opt = GrainYld, T1_Count, T2_Count, T3_Count, T4_Count, TCE) %>%
    mutate(
      T1_rel = T1_Count / T4_Count,
      T2_rel = T2_Count / T4_Count,
      T3_rel = T3_Count / T4_Count
    )

  m_cont_ratios <- lm(Yield_Opt ~ T1_rel + T2_rel + T3_rel, data = opt_shape_df)
  s_cont_ratios <- summary(m_cont_ratios)
  f_stat_r <- s_cont_ratios$fstatistic
  p_val_r  <- pf(f_stat_r[1], f_stat_r[2], f_stat_r[3], lower.tail = FALSE)

  shape_pca <- prcomp(scale(opt_shape_df[, c("T1_rel", "T2_rel", "T3_rel")]))
  opt_shape_df$PC1_shape <- shape_pca$x[, 1]
  opt_shape_df$PC2_shape <- shape_pca$x[, 2]

  m_cont_pca <- lm(Yield_Opt ~ PC1_shape + PC2_shape, data = opt_shape_df)
  s_cont_pca <- summary(m_cont_pca)
  f_stat_pca <- s_cont_pca$fstatistic
  p_val_pca  <- pf(f_stat_pca[1], f_stat_pca[2], f_stat_pca[3], lower.tail = FALSE)

  # Continuous Shape Regressions on NWPZ Warm (Randomized Delhi Heat) & PZ Warm II (Descriptive)
  dyn_df$T1_rel    <- opt_shape_df$T1_rel[match(dyn_df$Genotype, opt_shape_df$Genotype)]
  dyn_df$T2_rel    <- opt_shape_df$T2_rel[match(dyn_df$Genotype, opt_shape_df$Genotype)]
  dyn_df$T3_rel    <- opt_shape_df$T3_rel[match(dyn_df$Genotype, opt_shape_df$Genotype)]
  dyn_df$PC1_shape <- opt_shape_df$PC1_shape[match(dyn_df$Genotype, opt_shape_df$Genotype)]
  dyn_df$PC2_shape <- opt_shape_df$PC2_shape[match(dyn_df$Genotype, opt_shape_df$Genotype)]

  m_r_nw  <- lm(Yield_NWPZ_Warm ~ T1_rel + T2_rel + T3_rel, data = dyn_df)
  s_r_nw  <- summary(m_r_nw); f_r_nw <- s_r_nw$fstatistic; p_r_nw <- pf(f_r_nw[1], f_r_nw[2], f_r_nw[3], lower.tail = FALSE)

  m_p_nw  <- lm(Yield_NWPZ_Warm ~ PC1_shape + PC2_shape, data = dyn_df)
  s_p_nw  <- summary(m_p_nw); f_p_nw <- s_p_nw$fstatistic; p_p_nw <- pf(f_p_nw[1], f_p_nw[2], f_p_nw[3], lower.tail = FALSE)

  m_r_str <- lm(Yield_Stress ~ T1_rel + T2_rel + T3_rel, data = dyn_df)
  s_r_str <- summary(m_r_str)
  f_r_str <- s_r_str$fstatistic
  p_r_str <- pf(f_r_str[1], f_r_str[2], f_r_str[3], lower.tail = FALSE)

  m_p_str <- lm(Yield_Stress ~ PC1_shape + PC2_shape, data = dyn_df)
  s_p_str <- summary(m_p_str)
  f_p_str <- s_p_str$fstatistic
  p_p_str <- pf(f_p_str[1], f_p_str[2], f_p_str[3], lower.tail = FALSE)

  m_base_pos <- lm(Yield_Stress ~ factor(Field_Pos_Block_DhwS2), data = dyn_df)
  m_pos_r    <- lm(Yield_Stress ~ factor(Field_Pos_Block_DhwS2) + T1_rel + T2_rel + T3_rel, data = dyn_df)
  an_r       <- anova(m_base_pos, m_pos_r)
  delta_r2_r <- summary(m_pos_r)$r.squared - summary(m_base_pos)$r.squared
  f_pos_r    <- an_r$F[2]
  p_pos_r    <- an_r$`Pr(>F)`[2]

  m_pos_p    <- lm(Yield_Stress ~ factor(Field_Pos_Block_DhwS2) + PC1_shape + PC2_shape, data = dyn_df)
  an_p       <- anova(m_base_pos, m_pos_p)
  delta_r2_p <- summary(m_pos_p)$r.squared - summary(m_base_pos)$r.squared
  f_pos_p    <- an_p$F[2]
  p_pos_p    <- an_p$`Pr(>F)`[2]

  reg_summary_table <- data.frame(
    Model = c(
      "Proportional Tillering Ratios: Yield_Opt ~ T1/T4 + T2/T4 + T3/T4 (NWPZ Baseline)",
      "Proportional Tillering PCA: Yield_Opt ~ PC1_traj + PC2_traj (NWPZ Baseline)",
      "Proportional Tillering Ratios: Yield_NWPZ_Warm ~ T1/T4 + T2/T4 + T3/T4 (NWPZ Warm)",
      "Proportional Tillering PCA: Yield_NWPZ_Warm ~ PC1_traj + PC2_traj (NWPZ Warm)",
      "Proportional Tillering Ratios: Yield_Stress ~ T1/T4 + T2/T4 + T3/T4 (PZ Warm II Descriptive)",
      "Proportional Tillering PCA: Yield_Stress ~ PC1_traj + PC2_traj (PZ Warm II Descriptive)"
    ),
    R_Squared = c(s_cont_ratios$r.squared, s_cont_pca$r.squared, s_r_nw$r.squared, s_p_nw$r.squared, s_r_str$r.squared, s_p_str$r.squared),
    Adj_R_Squared = c(s_cont_ratios$adj.r.squared, s_cont_pca$adj.r.squared, s_r_nw$adj.r.squared, s_p_nw$adj.r.squared, s_r_str$adj.r.squared, s_p_str$adj.r.squared),
    F_Statistic = c(f_stat_r[1], f_stat_pca[1], f_r_nw[1], f_p_nw[1], f_r_str[1], f_p_str[1]),
    DF_Num = as.integer(c(f_stat_r[2], f_stat_pca[2], f_r_nw[2], f_p_nw[2], f_r_str[2], f_p_str[2])),
    DF_Den = as.integer(c(f_stat_r[3], f_stat_pca[3], f_r_nw[3], f_p_nw[3], f_r_str[3], f_p_str[3])),
    p_value = c(p_val_r, p_val_pca, p_r_nw, p_p_nw, p_r_str, p_p_str)
  )
  write.csv(reg_summary_table, file.path(sec25_dir, "Table_Cluster_Free_Shape_Regression.csv"), row.names = FALSE)

  # 7.4 Secondary Analyses (Clearly Labelled; Not Replacing Panel A or Described as Trajectory Validation)
  raw_mat_sec <- scale(as.matrix(opt_shape_df[, c("T1_Count", "T2_Count", "T3_Count", "T4_Count")]))
  set.seed(42)
  km_raw_sec <- kmeans(raw_mat_sec, centers = 2, nstart = 100)
  high_tiller_idx <- which.max(km_raw_sec$centers[, "T4_Count"])
  raw_class_names <- character(2)
  raw_class_names[high_tiller_idx]     <- "Higher Tiller Number"
  raw_class_names[3 - high_tiller_idx] <- "Lower Tiller Number"

  dyn_df$Raw_Tiller_Class <- factor(raw_class_names[km_raw_sec$cluster], levels = c("Higher Tiller Number", "Lower Tiller Number"))

  t_raw_opt <- t.test(Yield_Opt ~ Raw_Tiller_Class, data = dyn_df)
  t_raw_nw  <- t.test(Yield_NWPZ_Warm ~ Raw_Tiller_Class, data = dyn_df)
  t_raw_str <- t.test(Yield_Stress ~ Raw_Tiller_Class, data = dyn_df)
  anc_raw   <- fit_pos_ancova(dyn_df, "Raw_Tiller_Class")
  diff_raw_opt <- t_raw_opt$estimate[1] - t_raw_opt$estimate[2]
  diff_raw_nw  <- t_raw_nw$estimate[1] - t_raw_nw$estimate[2]
  diff_raw_str <- t_raw_str$estimate[1] - t_raw_str$estimate[2]
  d_raw_opt <- cohen_d(dyn_df$Yield_Opt[dyn_df$Raw_Tiller_Class == "Higher Tiller Number"],
                       dyn_df$Yield_Opt[dyn_df$Raw_Tiller_Class == "Lower Tiller Number"])
  d_raw_nw  <- cohen_d(dyn_df$Yield_NWPZ_Warm[dyn_df$Raw_Tiller_Class == "Higher Tiller Number"],
                       dyn_df$Yield_NWPZ_Warm[dyn_df$Raw_Tiller_Class == "Lower Tiller Number"])
  d_raw_str <- cohen_d(dyn_df$Yield_Stress[dyn_df$Raw_Tiller_Class == "Higher Tiller Number"],
                       dyn_df$Yield_Stress[dyn_df$Raw_Tiller_Class == "Lower Tiller Number"])

  raw_summary_table <- data.frame(
    Class_Comparison = "Higher vs Lower Tiller Number",
    Optimal_Mean_Diff = diff_raw_opt,
    Optimal_CI_Lower = t_raw_opt$conf.int[1],
    Optimal_CI_Upper = t_raw_opt$conf.int[2],
    Optimal_t = t_raw_opt$statistic,
    Optimal_df = t_raw_opt$parameter,
    Optimal_p = t_raw_opt$p.value,
    Optimal_d = d_raw_opt,
    NWPZ_Warm_Mean_Diff = diff_raw_nw,
    NWPZ_Warm_p = t_raw_nw$p.value,
    NWPZ_Warm_d = d_raw_nw,
    Stress_Mean_Diff = diff_raw_str,
    Stress_CI_Lower = t_raw_str$conf.int[1],
    Stress_CI_Upper = t_raw_str$conf.int[2],
    Stress_t = t_raw_str$statistic,
    Stress_df = t_raw_str$parameter,
    Stress_p = t_raw_str$p.value,
    Stress_d = d_raw_str,
    Stress_PosAdj_Diff = anc_raw$diff_adj,
    Stress_PosAdj_p    = anc_raw$p_adj,
    Note = "Secondary Analysis (a): Reflects absolute tillering capacity (PC1 of raw counts); not a trajectory validation."
  )
  write.csv(raw_summary_table, file.path(sec25_dir, "Table_Secondary_Tiller_Number_Classes_Yield.csv"), row.names = FALSE)

  # Secondary Analysis (b): Tiller Conversion Efficiency (TCE vs Yield)
  dyn_df$TCE_Opt <- opt_shape_df$TCE[match(dyn_df$Genotype, opt_shape_df$Genotype)]
  tce_med <- median(dyn_df$TCE_Opt, na.rm = TRUE)
  dyn_df$TCE_Class <- factor(ifelse(dyn_df$TCE_Opt >= tce_med, "Higher Conversion", "Lower Conversion"),
                             levels = c("Higher Conversion", "Lower Conversion"))

  t_tce_opt <- t.test(Yield_Opt ~ TCE_Class, data = dyn_df)
  t_tce_nw  <- t.test(Yield_NWPZ_Warm ~ TCE_Class, data = dyn_df)
  t_tce_str <- t.test(Yield_Stress ~ TCE_Class, data = dyn_df)
  anc_tce   <- fit_pos_ancova(dyn_df, "TCE_Class")
  diff_tce_opt <- t_tce_opt$estimate[1] - t_tce_opt$estimate[2]
  diff_tce_nw  <- t_tce_nw$estimate[1] - t_tce_nw$estimate[2]
  diff_tce_str <- t_tce_str$estimate[1] - t_tce_str$estimate[2]
  d_tce_opt <- cohen_d(dyn_df$Yield_Opt[dyn_df$TCE_Class == "Higher Conversion"],
                       dyn_df$Yield_Opt[dyn_df$TCE_Class == "Lower Conversion"])
  d_tce_nw  <- cohen_d(dyn_df$Yield_NWPZ_Warm[dyn_df$TCE_Class == "Higher Conversion"],
                       dyn_df$Yield_NWPZ_Warm[dyn_df$TCE_Class == "Lower Conversion"])
  d_tce_str <- cohen_d(dyn_df$Yield_Stress[dyn_df$TCE_Class == "Higher Conversion"],
                       dyn_df$Yield_Stress[dyn_df$TCE_Class == "Lower Conversion"])

  tce_summary_table <- data.frame(
    Class_Comparison = "Higher vs Lower Conversion (TCE)",
    Optimal_Mean_Diff = diff_tce_opt,
    Optimal_CI_Lower = t_tce_opt$conf.int[1],
    Optimal_CI_Upper = t_tce_opt$conf.int[2],
    Optimal_t = t_tce_opt$statistic,
    Optimal_df = t_tce_opt$parameter,
    Optimal_p = t_tce_opt$p.value,
    Optimal_d = d_tce_opt,
    NWPZ_Warm_Mean_Diff = diff_tce_nw,
    NWPZ_Warm_p = t_tce_nw$p.value,
    NWPZ_Warm_d = d_tce_nw,
    Stress_Mean_Diff = diff_tce_str,
    Stress_CI_Lower = t_tce_str$conf.int[1],
    Stress_CI_Upper = t_tce_str$conf.int[2],
    Stress_t = t_tce_str$statistic,
    Stress_df = t_tce_str$parameter,
    Stress_p = t_tce_str$p.value,
    Stress_d = d_tce_str,
    Stress_PosAdj_Diff = anc_tce$diff_adj,
    Stress_PosAdj_p    = anc_tce$p_adj,
    Note = "Secondary Analysis (b): PT is an algebraic component of both TCE (PT/T4) and grain yield; this association is partly mechanical. Not a trajectory validation."
  )
  write.csv(tce_summary_table, file.path(sec25_dir, "Table_Secondary_TCE_Yield.csv"), row.names = FALSE)

  # 7.5 Side-by-Side Statistical Comparison Table (Including Delhi NWPZ Warm & Dharwad Position-Adjusted Stress ANCOVA)
  side_by_side_table <- data.frame(
    Analysis = c(
      "Primary: Vegetative Tillering Trajectory Clusters (Earlier vs Later Proportional Peak)",
      sprintf("Primary: Tillering Vigor Clusters (%s vs %s)", VIGOR_HI, VIGOR_LO),
      "Primary: Tiller Stability Index Classes (Stable vs Contrasting)",
      "Continuous Check: Proportional Tillering Ratios (T1/T4, T2/T4, T3/T4)",
      "Continuous Check: Proportional Tillering PCA (PC1_traj, PC2_traj)",
      "Secondary (a): Tiller Number Classes (Raw Counts)",
      "Secondary (b): Tiller Conversion Efficiency (TCE)"
    ),
    Model_Type = c(
      "Welch t-test & 20-Block Position ANCOVA",
      "Welch t-test & 20-Block Position ANCOVA",
      "Welch t-test & 20-Block Position ANCOVA",
      "Multiple Linear Regression",
      "Multiple Linear Regression",
      "Welch t-test & 20-Block Position ANCOVA",
      "Welch t-test & 20-Block Position ANCOVA"
    ),
    Predictors = c(
      "Earlier Peak vs Later Peak",
      sprintf("%s vs %s", VIGOR_HI, VIGOR_LO),
      "Stable vs Contrasting",
      "T1/T4 + T2/T4 + T3/T4",
      "PC1_traj + PC2_traj",
      "Higher vs Lower Tiller Number",
      "Higher vs Lower Conversion"
    ),
    Optimal_Stat = c(
      sprintf("t(%.1f) = %.2f", t_traj_opt$parameter, t_traj_opt$statistic),
      sprintf("t(%.1f) = %.2f", t_vig_opt$parameter, t_vig_opt$statistic),
      sprintf("t(%.1f) = %.2f", t_tsi_opt$parameter, t_tsi_opt$statistic),
      sprintf("F(%d, %d) = %.2f", f_stat_r[2], f_stat_r[3], f_stat_r[1]),
      sprintf("F(%d, %d) = %.2f", f_stat_pca[2], f_stat_pca[3], f_stat_pca[1]),
      sprintf("t(%.1f) = %.2f", t_raw_opt$parameter, t_raw_opt$statistic),
      sprintf("t(%.1f) = %.2f", t_tce_opt$parameter, t_tce_opt$statistic)
    ),
    Optimal_Effect = c(
      sprintf("diff = %+.4f %s (d = %.2f)", diff_traj_opt, YIELD_UNIT_LABEL, d_traj_opt),
      sprintf("diff = %+.4f %s (d = %.2f)", diff_vig_opt, YIELD_UNIT_LABEL, d_vig_opt),
      sprintf("diff = %+.4f %s (d = %.2f)", diff_tsi_opt, YIELD_UNIT_LABEL, d_tsi_opt),
      sprintf("R² = %.4f (Adj R² = %.4f)", s_cont_ratios$r.squared, s_cont_ratios$adj.r.squared),
      sprintf("R² = %.4f (Adj R² = %.4f)", s_cont_pca$r.squared, s_cont_pca$adj.r.squared),
      sprintf("diff = %+.4f %s (d = %.2f)", diff_raw_opt, YIELD_UNIT_LABEL, d_raw_opt),
      sprintf("diff = %+.4f %s (d = %.2f)", diff_tce_opt, YIELD_UNIT_LABEL, d_tce_opt)
    ),
    Optimal_p = c(
      format_p_display(t_traj_opt$p.value),
      format_p_display(t_vig_opt$p.value),
      format_p_display(t_tsi_opt$p.value),
      format_p_display(p_val_r),
      format_p_display(p_val_pca),
      format_p_display(t_raw_opt$p.value),
      format_p_display(t_tce_opt$p.value)
    ),
    NWPZ_Warm_Stat = c(
      sprintf("t(%.1f) = %.2f", t_traj_nw$parameter, t_traj_nw$statistic),
      sprintf("t(%.1f) = %.2f", t_vig_nw$parameter, t_vig_nw$statistic),
      sprintf("t(%.1f) = %.2f", t_tsi_nw$parameter, t_tsi_nw$statistic),
      sprintf("F(%d, %d) = %.2f", f_r_nw[2], f_r_nw[3], f_r_nw[1]),
      sprintf("F(%d, %d) = %.2f", f_p_nw[2], f_p_nw[3], f_p_nw[1]),
      sprintf("t(%.1f) = %.2f", t_raw_nw$parameter, t_raw_nw$statistic),
      sprintf("t(%.1f) = %.2f", t_tce_nw$parameter, t_tce_nw$statistic)
    ),
    NWPZ_Warm_Effect = c(
      sprintf("diff = %+.4f %s (d = %.2f)", diff_traj_nw, YIELD_UNIT_LABEL, d_traj_nw),
      sprintf("diff = %+.4f %s (d = %.2f)", diff_vig_nw, YIELD_UNIT_LABEL, d_vig_nw),
      sprintf("diff = %+.4f %s (d = %.2f)", diff_tsi_nw, YIELD_UNIT_LABEL, d_tsi_nw),
      sprintf("R² = %.4f (Adj R² = %.4f)", s_r_nw$r.squared, s_r_nw$adj.r.squared),
      sprintf("R² = %.4f (Adj R² = %.4f)", s_p_nw$r.squared, s_p_nw$adj.r.squared),
      sprintf("diff = %+.4f %s (d = %.2f)", diff_raw_nw, YIELD_UNIT_LABEL, d_raw_nw),
      sprintf("diff = %+.4f %s (d = %.2f)", diff_tce_nw, YIELD_UNIT_LABEL, d_tce_nw)
    ),
    NWPZ_Warm_p = c(
      format_p_display(t_traj_nw$p.value),
      format_p_display(t_vig_nw$p.value),
      format_p_display(t_tsi_nw$p.value),
      format_p_display(p_r_nw),
      format_p_display(p_p_nw),
      format_p_display(t_raw_nw$p.value),
      format_p_display(t_tce_nw$p.value)
    ),
    Stress_Stat = c(
      sprintf("t(%.1f) = %.2f", t_traj_str$parameter, t_traj_str$statistic),
      sprintf("t(%.1f) = %.2f", t_vig_str$parameter, t_vig_str$statistic),
      sprintf("t(%.1f) = %.2f", t_tsi_str$parameter, t_tsi_str$statistic),
      sprintf("F(%d, %d) = %.2f", f_r_str[2], f_r_str[3], f_r_str[1]),
      sprintf("F(%d, %d) = %.2f", f_p_str[2], f_p_str[3], f_p_str[1]),
      sprintf("t(%.1f) = %.2f", t_raw_str$parameter, t_raw_str$statistic),
      sprintf("t(%.1f) = %.2f", t_tce_str$parameter, t_tce_str$statistic)
    ),
    Stress_Effect = c(
      sprintf("diff = %+.4f %s (d = %.2f)", diff_traj_str, YIELD_UNIT_LABEL, d_traj_str),
      sprintf("diff = %+.4f %s (d = %.2f)", diff_vig_str, YIELD_UNIT_LABEL, d_vig_str),
      sprintf("diff = %+.4f %s (d = %.2f)", diff_tsi_str, YIELD_UNIT_LABEL, d_tsi_str),
      sprintf("R² = %.4f (Adj R² = %.4f)", s_r_str$r.squared, s_r_str$adj.r.squared),
      sprintf("R² = %.4f (Adj R² = %.4f)", s_p_str$r.squared, s_p_str$adj.r.squared),
      sprintf("diff = %+.4f %s (d = %.2f)", diff_raw_str, YIELD_UNIT_LABEL, d_raw_str),
      sprintf("diff = %+.4f %s (d = %.2f)", diff_tce_str, YIELD_UNIT_LABEL, d_tce_str)
    ),
    Stress_p = c(
      format_p_display(t_traj_str$p.value),
      format_p_display(t_vig_str$p.value),
      format_p_display(t_tsi_str$p.value),
      format_p_display(p_r_str),
      format_p_display(p_p_str),
      format_p_display(t_raw_str$p.value),
      format_p_display(t_tce_str$p.value)
    ),
    Stress_PosAdj_Effect = c(
      sprintf("adj diff = %+.4f %s", anc_traj$diff_adj, YIELD_UNIT_LABEL),
      sprintf("adj diff = %+.4f %s", anc_vig$diff_adj, YIELD_UNIT_LABEL),
      sprintf("adj diff = %+.4f %s", anc_tsi$diff_adj, YIELD_UNIT_LABEL),
      sprintf("ΔR² = %+.4f (F = %.2f)", delta_r2_r, f_pos_r),
      sprintf("ΔR² = %+.4f (F = %.2f)", delta_r2_p, f_pos_p),
      sprintf("adj diff = %+.4f %s", anc_raw$diff_adj, YIELD_UNIT_LABEL),
      sprintf("adj diff = %+.4f %s", anc_tce$diff_adj, YIELD_UNIT_LABEL)
    ),
    Stress_PosAdj_p = c(
      format_p_display(anc_traj$p_adj),
      format_p_display(anc_vig$p_adj),
      format_p_display(anc_tsi$p_adj),
      format_p_display(p_pos_r),
      format_p_display(p_pos_p),
      format_p_display(anc_raw$p_adj),
      format_p_display(anc_tce$p_adj)
    ),
    stringsAsFactors = FALSE
  )

  write.csv(side_by_side_table, file.path(sec25_dir, "Table_Trajectory_and_Secondary_Analyses_Side_by_Side.csv"), row.names = FALSE)

  if (requireNamespace("openxlsx", quietly = TRUE)) {
    wb_sec25 <- openxlsx::createWorkbook()
    for (sh25 in list(
      list(name = "Table_6_Strategy_Regime_Ranks", data = cluster_summary),
      list(name = "Table_S11A_Side_by_Side_Tests", data = side_by_side_table),
      list(name = "Table_S11B_Cluster_Free_Regr", data = reg_summary_table),
      list(name = "Table_S11C_Secondary_TillerNum", data = raw_summary_table),
      list(name = "Table_S11D_Secondary_TCE", data = tce_summary_table)
    )) {
      openxlsx::addWorksheet(wb_sec25, sh25$name)
      openxlsx::writeDataTable(wb_sec25, sh25$name, sh25$data, withFilter = TRUE)
      openxlsx::freezePane(wb_sec25, sh25$name, firstRow = TRUE)
      openxlsx::setColWidths(wb_sec25, sh25$name, cols = seq_len(ncol(sh25$data)), widths = "auto")
    }
    openxlsx::saveWorkbook(wb_sec25, file.path(sec25_dir, "Section_2.5_Manuscript_and_Supplementary_Tables.xlsx"), overwrite = TRUE)
  }

  cat("\n==============================================================================\n")
  cat("SIDE-BY-SIDE STATISTICAL COMPARISON TABLE: TRAJECTORY & SECONDARY ANALYSES\n")
  cat("==============================================================================\n")
  print(side_by_side_table)
  cat("==============================================================================\n\n")
} else {
  cat("  -> Section 7 [SKIPPED - Assets Finalized & Protected]\n")
}


# ==============================================================================
# ==============================================================================
# SECTION 8: VARIANCE PARTITIONING, HERITABILITY & FINLAY-WILKINSON PLASTICITY
# ==============================================================================
cat(">>> Executing Section 8: Variance Partitioning, Heritability & FW Modeling...\n")

if (isTRUE(RUN_SECTION_8_PLASTICITY)) {
  sec26_plas_dir <- file.path(assets_dir, "Section_2.6_Phenotypic_Plasticity")
  sec26_dir      <- file.path(assets_dir, "Section_2.7_Variance_and_Heritability")
  dir.create(sec26_plas_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(sec26_dir, showWarnings = FALSE, recursive = TRUE)

  traits_to_model <- c("T1_Count", "T2_Count", "T3_Count", "T4_Count", "PT_Count", "GrainYld", "TKW", "AUTPC", "AUTPC_norm", "TCE")
  sec2_note_dharwad <- "Unrandomized Duplicated Layout: Genotype Confounded with Field Position"

  tier1_rows <- list()
  for (s_code in env_metadata$Sheet) {
    sub_t <- master_raw %>% filter(Sheet == s_code) %>% droplevels()
    loc_t <- unique(sub_t$Location)
    n_rep_t <- length(unique(sub_t$Replication))
    for (t in traits_to_model) {
      if (loc_t == "Delhi") {
        fit_w <- lmer(as.formula(paste0(t, " ~ Replication + (1 | Replication:Block_true) + (1 | Genotype)")),
                      data = sub_t, REML = TRUE, control = lmerControl(optimizer = "bobyqa"))
        vc_w <- as.data.frame(VarCorr(fit_w))
        vg_w <- vc_w$vcov[vc_w$grp == "Genotype"]
        vb_w <- vc_w$vcov[vc_w$grp == "Replication:Block_true"]
        ve_w <- vc_w$vcov[vc_w$grp == "Residual"]
        pev_w <- as.numeric(attr(ranef(fit_w, condVar = TRUE)$Genotype, "postVar"))
        h2_rel_val <- ifelse(vg_w <= 1e-10, 0.0, pmax(0.0, 1.0 - mean(pev_w, na.rm = TRUE) / vg_w))
        tier1_rows[[paste(s_code, t)]] <- data.frame(
          Sheet = s_code, Environment = unique(sub_t$Environment), Trait = t,
          Var_Genotype = sprintf("%.4f", vg_w), Var_RepBlock = sprintf("%.4f", vb_w), Var_Residual = sprintf("%.4f", ve_w),
          H2_EntryMean = sprintf("%.4f", vg_w / (vg_w + ve_w / n_rep_t)),
          H2_Reliability_MeanPEV = sprintf("%.4f", h2_rel_val),
          PseudoRep_Strip_Repeatability = "\u2014", Spatial_PosBlock_R2 = "\u2014",
          Singular = as.character(isSingular(fit_w)),
          Methodological_Status = "Randomized incomplete-block trial (20 blocks x 10 plots/rep)",
          stringsAsFactors = FALSE)
      } else {
        fit_d <- lm(as.formula(paste0(t, " ~ Genotype + Replication")), data = sub_t)
        an_d <- anova(fit_d)
        ms_g <- an_d["Genotype", "Mean Sq"]; ms_e <- an_d["Residuals", "Mean Sq"]
        vg_d <- (ms_g - ms_e) / n_rep_t
        pos_d <- blues_df %>% filter(Sheet == s_code)
        r2_pos <- summary(lm(pos_d[[t]] ~ factor(pos_d$Field_Pos_Block_Diagnostic)))$r.squared
        tier1_rows[[paste(s_code, t)]] <- data.frame(
          Sheet = s_code, Environment = unique(sub_t$Environment), Trait = t,
          Var_Genotype = sprintf("%.4f", vg_d), Var_RepBlock = "\u2014", Var_Residual = sprintf("%.4f", ms_e),
          H2_EntryMean = "\u2014", H2_Reliability_MeanPEV = "\u2014",
          PseudoRep_Strip_Repeatability = sprintf("%.4f", vg_d / (vg_d + ms_e / n_rep_t)),
          Spatial_PosBlock_R2 = sprintf("%.4f", r2_pos), Singular = "\u2014",
          Methodological_Status = sec2_note_dharwad,
          stringsAsFactors = FALSE)
      }
    }
  }
  tier1_table <- bind_rows(tier1_rows)
  write.csv(tier1_table, file.path(sec26_dir, "Table_Within_Trial_Variance_and_Heritability.csv"), row.names = FALSE)

  del_raw <- master_raw %>% filter(Location == "Delhi") %>% droplevels()
  n_env_del <- nlevels(del_raw$Env_Factor)
  n_rep_del <- nlevels(del_raw$Replication)
  tier2_rows <- list()
  for (t in traits_to_model) {
    fit_2 <- lmer(as.formula(paste0(t, " ~ Env_Factor + (1 | Env_Factor:Replication) + (1 | Env_Factor:Replication:Block_true) + (1 | Genotype) + (1 | Genotype:Env_Factor)")),
                  data = del_raw, REML = TRUE, control = lmerControl(optimizer = "bobyqa"))
    vc2 <- as.data.frame(VarCorr(fit_2))
    vg2  <- vc2$vcov[vc2$grp == "Genotype"]
    vge2 <- vc2$vcov[vc2$grp == "Genotype:Env_Factor"]
    ve2  <- vc2$vcov[vc2$grp == "Residual"]
    tier2_rows[[t]] <- data.frame(
      Trait = t, Vg_Delhi = round(vg2, 5), Vge_Delhi = round(vge2, 5), Ve_Delhi = round(ve2, 5),
      H2_Delhi_Randomized = round(vg2 / (vg2 + vge2 / n_env_del + ve2 / (n_env_del * n_rep_del)), 4),
      Singular_Delhi = isSingular(fit_2), stringsAsFactors = FALSE)
  }
  tier2_table <- bind_rows(tier2_rows)

  tier3_rows <- list()
  for (t in traits_to_model) {
    fit_3 <- lmer(as.formula(paste0(t, " ~ Env_Factor + (1 | Genotype) + (1 | Genotype:Location)")),
                  data = blues_df, REML = TRUE, control = lmerControl(optimizer = "bobyqa"))
    vc3 <- as.data.frame(VarCorr(fit_3))
    vg3  <- vc3$vcov[vc3$grp == "Genotype"]
    vgl3 <- vc3$vcov[vc3$grp == "Genotype:Location"]
    ve3  <- vc3$vcov[vc3$grp == "Residual"]
    n_loc <- length(unique(blues_df$Location)); n_sow <- nlevels(blues_df$Env_Factor) / n_loc
    tier3_rows[[t]] <- data.frame(
      Trait = t, Vg_4Env = round(vg3, 5), V_GxLocation_SiteDesignConfounded = round(vgl3, 5), V_GxSowing_plus_Error = round(ve3, 5),
      R2_Stage2_Repeatability = round(vg3 / (vg3 + vgl3 / n_loc + ve3 / (n_loc * n_sow)), 4),
      Singular_Stage2 = isSingular(fit_3), stringsAsFactors = FALSE)
  }
  tier3_table <- bind_rows(tier3_rows)

  h2_results <- tier2_table %>% left_join(tier3_table, by = "Trait")
  write.csv(h2_results, file.path(sec26_dir, "Global_BroadSense_Heritability.csv"), row.names = FALSE)

  format_clean_anova <- function(fit_obj, scope_note, type = c("I", "II")) {
    type <- match.arg(type)
    an_df <- if (type == "II") {
      as.data.frame(car::Anova(fit_obj, type = "II"))
    } else {
      as.data.frame(anova(fit_obj))
    }
    if (!"Mean Sq" %in% names(an_df)) {
      an_df$`Mean Sq` <- an_df$`Sum Sq` / an_df$Df
    }
    an_df$Source_of_Variation <- rownames(an_df)
    rownames(an_df) <- NULL
    pct_v <- (an_df$`Sum Sq` / sum(an_df$`Sum Sq`)) * 100
    data.frame(
      Source_of_Variation = an_df$Source_of_Variation,
      Df = an_df$Df,
      Sum_Sq = round(an_df$`Sum Sq`, 4),
      Mean_Sq = round(an_df$`Mean Sq`, 4),
      F_Value = ifelse(is.na(an_df$`F value`), "\u2014", sprintf("%.2f", an_df$`F value`)),
      P_Value = ifelse(is.na(an_df$`Pr(>F)`), "\u2014", ifelse(an_df$`Pr(>F)` < 0.0001, "< 0.0001", sprintf("%.4f", an_df$`Pr(>F)`))),
      Percent_Variance = round(pct_v, 2),
      Model_Scope = scope_note,
      stringsAsFactors = FALSE
    )
  }

  blues_df$Location <- factor(blues_df$Location)
  fit_anova <- lm(GrainYld ~ Genotype + Env_Factor + Genotype:Location, data = blues_df)
  anova_table <- format_clean_anova(fit_anova, "Stage-2 entry means; Residual = Genotype x Sowing(Location); Genotype x Location is site/design confounded")
  write.csv(anova_table, file.path(sec26_dir, "Global_Full_ANOVA_Table.csv"), row.names = FALSE)

  # Order blocking terms (Env_Factor, Replication, Block_true) BEFORE Genotype using keep.order = TRUE
  # so sequential Type-I SS adjusts Genotype and Genotype:Env_Factor for 3rd-order incomplete blocks
  f_del_keep <- terms(
    GrainYld ~ Env_Factor + Env_Factor:Replication + Env_Factor:Replication:Block_true + Genotype + Env_Factor:Genotype,
    keep.order = TRUE
  )
  fit_anova_del <- lm(f_del_keep, data = del_raw)
  anova_del <- format_clean_anova(fit_anova_del, "Delhi randomized incomplete-block trials (NWPZ Baseline & NWPZ Warm); plot-level ANOVA (Block-Adjusted)")
  write.csv(anova_del, file.path(sec26_dir, "Global_Full_ANOVA_Table_Delhi_Plot_Level.csv"), row.names = FALSE)

  # Export Clearly Labeled Main Text and Supplementary Excel Workbooks in Section 2.6 Variance & Heritability
  if (requireNamespace("openxlsx", quietly = TRUE)) {
    wb_t7 <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb_t7, "Table_7A_Tier2_Tier3_Global_H2")
    openxlsx::writeDataTable(wb_t7, "Table_7A_Tier2_Tier3_Global_H2", h2_results, withFilter = TRUE)
    openxlsx::freezePane(wb_t7, "Table_7A_Tier2_Tier3_Global_H2", firstRow = TRUE)
    openxlsx::setColWidths(wb_t7, "Table_7A_Tier2_Tier3_Global_H2", cols = seq_len(ncol(h2_results)), widths = "auto")

    openxlsx::addWorksheet(wb_t7, "Table_7B_Tier1_Within_Trial_H2")
    openxlsx::writeDataTable(wb_t7, "Table_7B_Tier1_Within_Trial_H2", tier1_table, withFilter = TRUE)
    openxlsx::freezePane(wb_t7, "Table_7B_Tier1_Within_Trial_H2", firstRow = TRUE)
    openxlsx::setColWidths(wb_t7, "Table_7B_Tier1_Within_Trial_H2", cols = seq_len(ncol(tier1_table)), widths = "auto")
    openxlsx::saveWorkbook(wb_t7, file.path(sec26_dir, "Table_7_MAIN_TEXT_Stratified_Variance_and_Heritability.xlsx"), overwrite = TRUE)

    wb_s15 <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb_s15, "Table_S15A_Stage2_BLUE_ANOVA")
    openxlsx::writeDataTable(wb_s15, "Table_S15A_Stage2_BLUE_ANOVA", anova_table, withFilter = TRUE)
    openxlsx::freezePane(wb_s15, "Table_S15A_Stage2_BLUE_ANOVA", firstRow = TRUE)
    openxlsx::setColWidths(wb_s15, "Table_S15A_Stage2_BLUE_ANOVA", cols = seq_len(ncol(anova_table)), widths = "auto")

    openxlsx::addWorksheet(wb_s15, "Table_S15B_Delhi_Plot_ANOVA")
    openxlsx::writeDataTable(wb_s15, "Table_S15B_Delhi_Plot_ANOVA", anova_del, withFilter = TRUE)
    openxlsx::freezePane(wb_s15, "Table_S15B_Delhi_Plot_ANOVA", firstRow = TRUE)
    openxlsx::setColWidths(wb_s15, "Table_S15B_Delhi_Plot_ANOVA", cols = seq_len(ncol(anova_del)), widths = "auto")
    openxlsx::saveWorkbook(wb_s15, file.path(sec26_dir, "Table_S15_SUPPLEMENTARY_Stage2_and_Delhi_Plot_ANOVA.xlsx"), overwrite = TRUE)
  }

  # 8.2 Finlay-Wilkinson Plasticity Regression Modeling & Empirical Quadrants
  env_means <- blues_df %>% group_by(Environment) %>% summarise(Env_Index = mean(GrainYld, na.rm = TRUE), .groups = "drop")
  fw_df <- blues_df %>% left_join(env_means, by = "Environment") %>% mutate(Env_Index_c = Env_Index - mean(Env_Index))

  m_main <- lm(GrainYld ~ Genotype + Env_Index, data = fw_df)
  m_int  <- lm(GrainYld ~ Genotype * Env_Index, data = fw_df)
  f_comp <- anova(m_main, m_int)
  f_stat <- f_comp$F[2]; df_num <- f_comp$Df[2]; df_den <- f_comp$Res.Df[2]; p_val_f <- f_comp$`Pr(>F)`[2]

  m_lmm_0 <- lmer(GrainYld ~ Env_Index_c + (1 | Genotype), data = fw_df, REML = FALSE)
  m_lmm_1 <- lmer(GrainYld ~ Env_Index_c + (1 + Env_Index_c | Genotype), data = fw_df, REML = FALSE, control = lmerControl(optimizer = "bobyqa"))
  lrt_stat <- as.numeric(2 * (logLik(m_lmm_1) - logLik(m_lmm_0)))
  p_val_mixture <- 0.5 * pchisq(lrt_stat, df = 1, lower.tail = FALSE) + 0.5 * pchisq(lrt_stat, df = 2, lower.tail = FALSE)

  fit_lmm_c <- lmer(GrainYld ~ Env_Index_c + (1 + Env_Index_c | Genotype), data = fw_df, control = lmerControl(optimizer = "bobyqa"))
  vc_c <- as.data.frame(VarCorr(fit_lmm_c))
  idx_b  <- which(vc_c$grp == "Genotype" & vc_c$var1 == "Env_Index_c" & is.na(vc_c$var2))
  idx_a  <- which(vc_c$grp == "Genotype" & vc_c$var1 == "(Intercept)" & is.na(vc_c$var2))
  idx_ab <- which(vc_c$grp == "Genotype" & vc_c$var1 == "(Intercept)" & vc_c$var2 == "Env_Index_c")

  sigma2_b <- if (length(idx_b) > 0) vc_c$vcov[idx_b] else NA_real_
  sigma2_a <- if (length(idx_a) > 0) vc_c$vcov[idx_a] else NA_real_
  corr_ab  <- if (length(idx_ab) > 0) vc_c$sdcor[idx_ab] else NA_real_

  ran <- ranef(fit_lmm_c)$Genotype
  ran$Genotype <- rownames(ran)
  ran$BLUP_Slope <- fixef(fit_lmm_c)["Env_Index_c"] + ran$Env_Index_c

  env_totals <- blues_df %>% group_by(Environment) %>% summarise(Sum_Yield = sum(GrainYld, na.rm = TRUE), N_Gen = sum(!is.na(GrainYld)), .groups = "drop")
  blues_logo <- blues_df %>% left_join(env_totals, by = "Environment") %>% mutate(Env_Index_LOGO = (Sum_Yield - GrainYld) / (N_Gen - 1))
  logo_models <- blues_logo %>% group_by(Genotype) %>% do({ mod <- lm(GrainYld ~ Env_Index_LOGO, data = .); data.frame(Slope_LOGO = coef(mod)["Env_Index_LOGO"]) }) %>% ungroup()

  fw_models <- fw_df %>%
    group_by(Genotype) %>%
    do({
      mod <- lm(GrainYld ~ Env_Index, data = .)
      res_df <- mod$df.residual
      t_crit <- qt(0.975, df = res_df)
      b <- coef(mod)["Env_Index"]
      se_b <- summary(mod)$coefficients["Env_Index", "Std. Error"]
      ci_lower <- b - t_crit * se_b; ci_upper <- b + t_crit * se_b
      t_b1 <- (b - 1.0) / se_b; p_b1 <- 2 * (1 - pt(abs(t_b1), df = res_df))
      data.frame(Mean_Yield = mean(.$GrainYld, na.rm = TRUE), FW_Slope = b, FW_Slope_SE = se_b,
                 FW_Slope_CI_Lower = ci_lower, FW_Slope_CI_Upper = ci_upper,
                 CI_excludes_1 = (ci_lower > 1.0) | (ci_upper < 1.0), p_val_b1 = p_b1, FW_R2 = summary(mod)$r.squared)
    }) %>% ungroup()

  fw_models$p_val_b1_holm <- p.adjust(fw_models$p_val_b1, method = "holm")
  fw_models <- fw_models %>% left_join(ran %>% select(Genotype, BLUP_Slope), by = "Genotype") %>% left_join(logo_models, by = "Genotype")
  var_ols_slopes <- var(fw_models$FW_Slope)
  slope_reliability <- sigma2_b / var_ols_slopes
  grand_mean_yield <- mean(fw_models$Mean_Yield, na.rm = TRUE)

  fw_models <- fw_models %>%
    mutate(
      Quadrant = factor(case_when(
        Mean_Yield >= grand_mean_yield & FW_Slope >= 1.0 ~ "High yield, high responsiveness",
        Mean_Yield >= grand_mean_yield & FW_Slope <  1.0 ~ "High yield, low responsiveness",
        Mean_Yield <  grand_mean_yield & FW_Slope <  1.0 ~ "Low yield, low responsiveness",
        Mean_Yield <  grand_mean_yield & FW_Slope >= 1.0 ~ "Low yield, high responsiveness"
      ), levels = c("High yield, high responsiveness", "High yield, low responsiveness", "Low yield, low responsiveness", "Low yield, high responsiveness")),
      Quadrant_BLUP = factor(case_when(
        Mean_Yield >= grand_mean_yield & BLUP_Slope >= 1.0 ~ "High yield, high responsiveness",
        Mean_Yield >= grand_mean_yield & BLUP_Slope <  1.0 ~ "High yield, low responsiveness",
        Mean_Yield <  grand_mean_yield & BLUP_Slope <  1.0 ~ "Low yield, low responsiveness",
        Mean_Yield <  grand_mean_yield & BLUP_Slope >= 1.0 ~ "Low yield, high responsiveness"
      ), levels = c("High yield, high responsiveness", "High yield, low responsiveness", "Low yield, low responsiveness", "Low yield, high responsiveness"))
    )
  write.csv(fw_models, file.path(sec26_plas_dir, "Finlay_Wilkinson_Classes.csv"), row.names = FALSE)

  envs_loeo <- master_pheno_spec$Env_Name
  loeo_summary <- data.frame()
  for (k in seq_along(envs_loeo)) {
    left_out <- envs_loeo[k]
    sub_df <- blues_df %>% filter(Environment != left_out)
    sub_env_means <- sub_df %>% group_by(Environment) %>% summarise(Env_Index_sub = mean(GrainYld, na.rm = TRUE), .groups = "drop")
    sub_df <- sub_df %>% left_join(sub_env_means, by = "Environment") %>% mutate(Env_Index_sub_c = Env_Index_sub - mean(Env_Index_sub))
    sub_fw <- sub_df %>% group_by(Genotype) %>% do({ mod <- lm(GrainYld ~ Env_Index_sub, data = .); data.frame(Mean_Yield_sub = mean(.$GrainYld, na.rm = TRUE), Slope_OLS_sub = coef(mod)["Env_Index_sub"]) }) %>% ungroup()
    sub_lmm <- lmer(GrainYld ~ Env_Index_sub_c + (1 + Env_Index_sub_c | Genotype), data = sub_df, control = lmerControl(optimizer = "bobyqa"))
    sub_ran <- ranef(sub_lmm)$Genotype; sub_ran$Genotype <- rownames(sub_ran); sub_ran$Slope_BLUP_sub <- fixef(sub_lmm)["Env_Index_sub_c"] + sub_ran$Env_Index_sub_c
    sub_fw <- sub_fw %>% left_join(sub_ran %>% select(Genotype, Slope_BLUP_sub), by = "Genotype")
    sub_grand_mean <- mean(sub_fw$Mean_Yield_sub, na.rm = TRUE)
    sub_fw <- sub_fw %>% mutate(
      Quadrant_OLS_sub = case_when(Mean_Yield_sub >= sub_grand_mean & Slope_OLS_sub >= 1.0 ~ "High yield, high responsiveness", Mean_Yield_sub >= sub_grand_mean & Slope_OLS_sub < 1.0 ~ "High yield, low responsiveness", Mean_Yield_sub < sub_grand_mean & Slope_OLS_sub < 1.0 ~ "Low yield, low responsiveness", TRUE ~ "Low yield, high responsiveness"),
      Quadrant_BLUP_sub = case_when(Mean_Yield_sub >= sub_grand_mean & Slope_BLUP_sub >= 1.0 ~ "High yield, high responsiveness", Mean_Yield_sub >= sub_grand_mean & Slope_BLUP_sub < 1.0 ~ "High yield, low responsiveness", Mean_Yield_sub < sub_grand_mean & Slope_BLUP_sub < 1.0 ~ "Low yield, low responsiveness", TRUE ~ "Low yield, high responsiveness")
    )
    comp <- fw_models %>% select(Genotype, Quadrant, Quadrant_BLUP) %>% left_join(sub_fw %>% select(Genotype, Quadrant_OLS_sub, Quadrant_BLUP_sub), by = "Genotype") %>% mutate(Changed_OLS = as.character(Quadrant) != Quadrant_OLS_sub, Changed_BLUP = as.character(Quadrant_BLUP) != Quadrant_BLUP_sub)
    loeo_summary <- rbind(loeo_summary, data.frame(Left_Out_Environment = left_out, OLS_Changed_N = sum(comp$Changed_OLS), OLS_Changed_Pct = mean(comp$Changed_OLS) * 100, BLUP_Changed_N = sum(comp$Changed_BLUP), BLUP_Changed_Pct = mean(comp$Changed_BLUP) * 100))
  }

  set.seed(42); n_pert_reps <- 2000; pert_changes_ols <- numeric(nrow(fw_models))
  for (p in 1:n_pert_reps) {
    b_pert <- fw_models$FW_Slope + rt(nrow(fw_models), df = 2) * fw_models$FW_Slope_SE
    q_pert <- case_when(fw_models$Mean_Yield >= grand_mean_yield & b_pert >= 1.0 ~ "High yield, high responsiveness", fw_models$Mean_Yield >= grand_mean_yield & b_pert < 1.0 ~ "High yield, low responsiveness", fw_models$Mean_Yield < grand_mean_yield & b_pert < 1.0 ~ "Low yield, low responsiveness", TRUE ~ "Low yield, high responsiveness")
    pert_changes_ols <- pert_changes_ols + (q_pert != as.character(fw_models$Quadrant))
  }
  heuristic_pert_pct <- mean(pert_changes_ols / n_pert_reps) * 100

  quad_summary <- fw_models %>% group_by(Quadrant) %>% summarise(n = n(), Pct = (n() / nrow(fw_models)) * 100, Mean_Yield_Avg = mean(Mean_Yield), Mean_Yield_SD = sd(Mean_Yield), Mean_Slope_Avg = mean(FW_Slope), Mean_Slope_SD = sd(FW_Slope), Mean_Slope_SE = mean(FW_Slope_SE), Median_Slope_SE = median(FW_Slope_SE), n_CI_excludes_1 = sum(CI_excludes_1), .groups = "drop")
  write.csv(quad_summary, file.path(sec26_plas_dir, "Table_Finlay_Wilkinson_Quadrant_Summary.csv"), row.names = FALSE)
  write.csv(loeo_summary, file.path(sec26_plas_dir, "Table_Finlay_Wilkinson_Stability_Diagnostics.csv"), row.names = FALSE)

  quad_palette <- c("High yield, high responsiveness" = "#2b8cbe", "High yield, low responsiveness" = "#2ca25f", "Low yield, low responsiveness" = "#7f7f7f", "Low yield, high responsiveness" = "#de2d26")
  n_unadj <- sum(fw_models$CI_excludes_1); n_total <- nrow(fw_models); pct_unadj <- (n_unadj / n_total) * 100; exp_chance <- 0.05 * n_total; n_holm <- sum(fw_models$p_val_b1_holm < 0.05); med_se <- median(fw_models$FW_Slope_SE); q_counts <- table(fw_models$Quadrant)

  # Clean Simplified Subtitles (Action 5B)
  sub_biplot <- sprintf(
    "Grand mean = %.3f %s | Median slope SE = %.2f (n = %d genotypes across 4 environments; %d High/High, %d High/Low, %d Low/Low, %d Low/High)\nBlack rings denote genotypes with unadjusted 95%% CI excluding b_i = 1.0 (n = %d [%.1f%%]; Holm-adjusted significant = %d)",
    grand_mean_yield, YIELD_UNIT_LABEL, med_se, n_total,
    q_counts["High yield, high responsiveness"], q_counts["High yield, low responsiveness"],
    q_counts["Low yield, low responsiveness"], q_counts["Low yield, high responsiveness"],
    n_unadj, pct_unadj, n_holm
  )

  fw_label_df <- fw_models %>% mutate(rk_yield = rank(-Mean_Yield, ties.method = "first"), rk_slope_hi = rank(-FW_Slope, ties.method = "first"), rk_slope_lo = rank(FW_Slope, ties.method = "first"), is_check = !grepl("^EDPIE_", Genotype)) %>% filter(rk_yield <= 8 | (CI_excludes_1 & (rk_slope_hi <= 3 | rk_slope_lo <= 3)) | (is_check & rk_yield <= 40))

  p_fw_biplot <- ggplot(fw_models, aes(x = Mean_Yield, y = FW_Slope)) +
    geom_hline(yintercept = 1.0, linetype = "dashed", color = "grey40", linewidth = 0.6) +
    geom_vline(xintercept = grand_mean_yield, linetype = "dashed", color = "grey40", linewidth = 0.6) +
    geom_point(aes(fill = Quadrant), shape = 21, color = "black", alpha = 0.75, size = 3.5) +
    geom_point(data = filter(fw_models, CI_excludes_1), aes(color = "95% CI excludes 1.0, unadjusted"), shape = 21, size = 5.5, stroke = 1.4, fill = NA) +
    ggrepel::geom_label_repel(data = fw_label_df, aes(x = Mean_Yield, y = FW_Slope, label = Genotype), inherit.aes = FALSE, size = 2.75, fontface = "bold", color = "grey10", fill = "white", alpha = 0.88, label.size = 0.25, box.padding = 0.32, point.padding = 0.25, max.overlaps = 20, min.segment.length = 0.05, seed = 42) +
    scale_fill_manual(values = quad_palette) +
    scale_color_manual(name = "Unadjusted CI Ring:", values = c("95% CI excludes 1.0, unadjusted" = "black")) +
    pub_theme +
    labs(title = "Finlay-Wilkinson Plasticity Biplot (Empirical Quadrant Partition)", subtitle = sub_biplot, x = sprintf("Mean Grain Yield Across Environments (%s)", YIELD_UNIT_LABEL), y = "Environmental Responsiveness (Regression Slope, b_i)", fill = "Quadrant Classification:") +
    theme(plot.subtitle = element_text(size = 9.2, color = "grey20", lineheight = 1.15), legend.box = "vertical", legend.box.just = "center", legend.spacing.y = unit(0.15, "cm")) +
    guides(fill = guide_legend(nrow = 2, byrow = TRUE, order = 1), color = guide_legend(order = 2))
  ggsave(file.path(sec26_plas_dir, "Finlay_Wilkinson_Biplot.png"), plot = p_fw_biplot, width = 10.5, height = 8.0, dpi = 600)

  fw_plot_df <- fw_df %>% left_join(fw_models %>% select(Genotype, Quadrant), by = "Genotype")
  fw_iqr_df <- fw_plot_df %>% group_by(Quadrant, Environment, Env_Index) %>% summarise(q25_yld = as.numeric(quantile(GrainYld, 0.25, na.rm = TRUE)), q75_yld = as.numeric(quantile(GrainYld, 0.75, na.rm = TRUE)), mean_yld = mean(GrainYld, na.rm = TRUE), .groups = "drop")
  env_anchor_df <- env_means %>% arrange(Env_Index) %>% mutate(Short_Env = sprintf("%s\n(%.3f)", sub("Baseline", "Base", Environment), Env_Index), y_top = max(fw_plot_df$GrainYld, na.rm = TRUE) * 1.04)

  sub_norms <- sprintf("Linear reaction norms across 4 thermal regimes (N = 200; median slope SE = %.2f).\nVertical dotted lines mark realized environmental mean yields; shaded ribbons = 25\u201375%% IQR; bold lines = quadrant OLS regression.", med_se)

  p_fw_norms <- ggplot(fw_plot_df, aes(x = Env_Index, y = GrainYld, group = Genotype)) +
    geom_vline(data = env_anchor_df, aes(xintercept = Env_Index), linetype = "dotted", color = "grey45", linewidth = 0.55) +
    geom_ribbon(data = fw_iqr_df, aes(x = Env_Index, ymin = q25_yld, ymax = q75_yld, fill = Quadrant, group = Quadrant), inherit.aes = FALSE, alpha = 0.18) +
    geom_line(alpha = 0.20, color = "grey45", linewidth = 0.35) +
    geom_smooth(aes(group = Quadrant, color = Quadrant), method = "lm", se = FALSE, linewidth = 1.35) +
    geom_point(data = fw_iqr_df, aes(x = Env_Index, y = mean_yld, fill = Quadrant), inherit.aes = FALSE, shape = 21, size = 2.7, color = "black", stroke = 0.5) +
    geom_text(data = env_anchor_df, aes(x = Env_Index, y = y_top, label = Short_Env), inherit.aes = FALSE, size = 2.55, fontface = "bold", color = "grey25", lineheight = 0.9, vjust = 0.5) +
    scale_color_manual(values = quad_palette) + scale_fill_manual(values = quad_palette) +
    scale_x_continuous(expand = expansion(mult = c(0.08, 0.08))) + scale_y_continuous(expand = expansion(mult = c(0.04, 0.12))) +
    facet_wrap(~ Quadrant, nrow = 2) + pub_theme +
    labs(title = "Finlay-Wilkinson Yield Reaction Norms by Empirical Quadrant", subtitle = sub_norms, x = sprintf("Environmental Quality Index (Mean Yield, %s)", YIELD_UNIT_LABEL), y = sprintf("Grain Yield (%s)", YIELD_UNIT_LABEL), color = "Quadrant:", fill = "Quadrant:") +
    theme(plot.subtitle = element_text(size = 9.2, color = "grey20", lineheight = 1.15), legend.position = "none")
  ggsave(file.path(sec26_plas_dir, "Finlay_Wilkinson_Reaction_Norms.png"), plot = p_fw_norms, width = 10.5, height = 8.2, dpi = 600)

  sub_blup <- sprintf("Ordinary Least Squares (OLS) vs. Empirical Bayes REML random-regression BLUP slopes across 4 regimes (r = %.4f).\nRandom slope variance \u03c3\u00b2_b = %.4f vs. OLS variance = %.4f (Slope reliability = %.1f%%; LRT p = %.2f).", cor(fw_models$FW_Slope, fw_models$BLUP_Slope), sigma2_b, var_ols_slopes, slope_reliability * 100, p_val_mixture)

  p_blup <- ggplot(fw_models, aes(x = FW_Slope, y = BLUP_Slope, color = Quadrant)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = 1.0, linetype = "dotted", color = "grey60") +
    geom_vline(xintercept = 1.0, linetype = "dotted", color = "grey60") +
    geom_point(alpha = 0.8, size = 3) + scale_color_manual(values = quad_palette) + pub_theme +
    labs(title = "Empirical Bayes BLUP Shrinkage of Environmental Responsiveness", subtitle = sub_blup, x = "Separate Regression OLS Slope (b_i)", y = "Shrunken BLUP Slope (Random Regression)", color = "Quadrant:") +
    theme(plot.subtitle = element_text(size = 9.2, color = "grey20", lineheight = 1.15)) +
    guides(color = guide_legend(nrow = 2, byrow = TRUE))
  ggsave(file.path(sec26_plas_dir, "Finlay_Wilkinson_BLUP_Shrinkage.png"), plot = p_blup, width = 9.5, height = 7.6, dpi = 600)

  fw_cluster_merged <- inner_join(fw_models %>% select(Genotype, FW_Quadrant = Quadrant, FW_Slope, Mean_Yield), master_profiles, by = "Genotype")
  tab_s13_crosstab <- fw_cluster_merged %>% group_by(FW_Quadrant) %>% summarise(
    Total_Genotypes = n(), Earlier_Peak_n = sum(Trajectory == "Earlier Proportional Peak"), Earlier_Peak_pct = round(100 * Earlier_Peak_n / n(), 1),
    Later_Peak_n = sum(Trajectory == "Later Proportional Peak"), Later_Peak_pct = round(100 * Later_Peak_n / n(), 1),
    High_Vigor_n = sum(Vigor == VIGOR_HI), High_Vigor_pct = round(100 * High_Vigor_n / n(), 1),
    Low_Vigor_n = sum(Vigor == VIGOR_LO), Low_Vigor_pct = round(100 * Low_Vigor_n / n(), 1),
    TSI_Stable_n = sum(Stability == "Stable"), TSI_Stable_pct = round(100 * TSI_Stable_n / n(), 1),
    TSI_Contrasting_n = sum(Stability == "Contrasting"), TSI_Contrasting_pct = round(100 * TSI_Contrasting_n / n(), 1),
    Mean_Yield = round(mean(Mean_Yield), 3), Mean_FW_Slope = round(mean(FW_Slope), 2), .groups = "drop")
  write.csv(tab_s13_crosstab, file.path(sec26_plas_dir, "Table_S13_Finlay_Wilkinson_Quadrants_vs_Developmental_Clusters.csv"), row.names = FALSE)

  # 8.3 Quantitative Comparison: Empirical Agronomic FW Quadrants vs. Unsupervised Clustering (K-Means, PAM, Ward.D2)
  fw_std_mat <- scale(as.matrix(fw_models[, c("Mean_Yield", "FW_Slope")]))
  dist_fw <- dist(fw_std_mat)
  set.seed(42)
  km4_fw   <- kmeans(fw_std_mat, centers = 4, nstart = 100)
  km2_fw   <- kmeans(fw_std_mat, centers = 2, nstart = 100)
  pam4_fw  <- cluster::pam(fw_std_mat, k = 4)
  hc4_fw   <- cutree(hclust(dist_fw, method = "ward.D2"), k = 4)

  sil_quad <- mean(cluster::silhouette(as.integer(fw_models$Quadrant), dist_fw)[, "sil_width"])
  sil_km4  <- mean(cluster::silhouette(km4_fw$cluster, dist_fw)[, "sil_width"])
  sil_pam4 <- mean(cluster::silhouette(pam4_fw$clustering, dist_fw)[, "sil_width"])
  sil_hc4  <- mean(cluster::silhouette(hc4_fw, dist_fw)[, "sil_width"])
  sil_km2  <- mean(cluster::silhouette(km2_fw$cluster, dist_fw)[, "sil_width"])

  calc_ari_simple <- function(x, y) {
    tab <- table(x, y)
    sum_comb_ij <- sum(choose(tab, 2))
    sum_comb_i  <- sum(choose(rowSums(tab), 2))
    sum_comb_j  <- sum(choose(colSums(tab), 2))
    n_tot_c     <- choose(length(x), 2)
    exp_idx     <- (sum_comb_i * sum_comb_j) / n_tot_c
    max_idx     <- 0.5 * (sum_comb_i + sum_comb_j)
    (sum_comb_ij - exp_idx) / (max_idx - exp_idx)
  }

  # Compute LOEO stability for Unsupervised K-Means (K = 4, K = 2), PAM (K = 4), and Ward.D2 (K = 4)
  match_cluster_reclass_pct <- function(cl_ref, cl_sub) {
    tab <- table(cl_ref, cl_sub)
    k_cl <- nrow(tab)
    if (k_cl == 2L) {
      best_match <- max(sum(diag(tab)), tab[1, 2] + tab[2, 1])
    } else {
      perms <- as.matrix(expand.grid(rep(list(seq_len(k_cl)), k_cl)))
      perms <- perms[apply(perms, 1, function(r) length(unique(r)) == k_cl), , drop = FALSE]
      best_match <- max(apply(perms, 1, function(p) sum(tab[cbind(seq_len(k_cl), p)])))
    }
    (1.0 - best_match / length(cl_ref)) * 100
  }

  loeo_km4_changes  <- numeric(length(envs_loeo))
  loeo_pam4_changes <- numeric(length(envs_loeo))
  loeo_hc4_changes  <- numeric(length(envs_loeo))
  loeo_km2_changes  <- numeric(length(envs_loeo))
  for (k in seq_along(envs_loeo)) {
    left_out <- envs_loeo[k]
    sub_df <- blues_df %>% filter(Environment != left_out)
    sub_env_means <- sub_df %>% group_by(Environment) %>% summarise(Env_Index_sub = mean(GrainYld, na.rm = TRUE), .groups = "drop")
    sub_df <- sub_df %>% left_join(sub_env_means, by = "Environment")
    sub_fw <- sub_df %>% group_by(Genotype) %>% do({ mod <- lm(GrainYld ~ Env_Index_sub, data = .); data.frame(Mean_Yield = mean(.$GrainYld, na.rm = TRUE), FW_Slope = coef(mod)["Env_Index_sub"]) }) %>% ungroup()
    sub_fw <- sub_fw[match(fw_models$Genotype, sub_fw$Genotype), ]
    sub_std <- scale(as.matrix(sub_fw[, c("Mean_Yield", "FW_Slope")]))
    set.seed(42)
    km4_sub  <- kmeans(sub_std, centers = km4_fw$centers, iter.max = 100)
    pam4_sub <- cluster::pam(sub_std, k = 4)
    hc4_sub  <- cutree(hclust(dist(sub_std), method = "ward.D2"), k = 4)
    km2_sub  <- kmeans(sub_std, centers = km2_fw$centers, iter.max = 100)
    loeo_km4_changes[k]  <- match_cluster_reclass_pct(km4_fw$cluster, km4_sub$cluster)
    loeo_pam4_changes[k] <- match_cluster_reclass_pct(pam4_fw$clustering, pam4_sub$clustering)
    loeo_hc4_changes[k]  <- match_cluster_reclass_pct(hc4_fw, hc4_sub)
    loeo_km2_changes[k]  <- match_cluster_reclass_pct(km2_fw$cluster, km2_sub$cluster)
  }

  fw_clustering_comp <- data.frame(
    Classification_Method = c(
      "Empirical Agronomic Quadrants (Y_bar = 0.659, b_i = 1.0)",
      "Unsupervised K-Means Clustering (K = 4 on Std. Yield & Slope)",
      "Unsupervised PAM Medoids (K = 4 on Std. Yield & Slope)",
      "Unsupervised Ward.D2 Hierarchical (K = 4 on Std. Yield & Slope)",
      "Unsupervised K-Means Clustering (Optimal K = 2 on Std. Yield & Slope)"
    ),
    Number_of_Groups = c(4L, 4L, 4L, 4L, 2L),
    Mean_Silhouette_Width = round(c(sil_quad, sil_km4, sil_pam4, sil_hc4, sil_km2), 4),
    ARI_vs_Empirical_Quadrants = round(c(1.0000, calc_ari_simple(fw_models$Quadrant, km4_fw$cluster),
                                         calc_ari_simple(fw_models$Quadrant, pam4_fw$clustering),
                                         calc_ari_simple(fw_models$Quadrant, hc4_fw),
                                         calc_ari_simple(fw_models$Quadrant, km2_fw$cluster)), 4),
    Mean_LOEO_Reclassification_Pct = round(c(mean(loeo_summary$OLS_Changed_Pct), mean(loeo_km4_changes), mean(loeo_pam4_changes), mean(loeo_hc4_changes), mean(loeo_km2_changes)), 2),
    Boundary_Anchor_Type = c(
      "Fixed Biological/Agronomic Benchmark (b_i = 1.0 population average response; Y_bar = grand mean)",
      "Data-Driven Sample Centroid Partition (shifts sample-to-sample; cuts across b_i = 1.0)",
      "Data-Driven Sample Medoid Partition (shifts sample-to-sample; cuts across b_i = 1.0)",
      "Data-Driven Agglomerative Tree Cut (sensitive to extreme slope outliers)",
      "Data-Driven Diagonal Split (conflates yield and slope into a single axis)"
    ),
    Recommendation_and_Interpretation = c(
      "RECOMMENDED: Preserves Finlay & Wilkinson (1963) canonical biological benchmark (b_i = 1.0) while acknowledging continuous bivariate distribution.",
      "Not recommended as primary: Only +0.05 silhouette gain over fixed quadrants (0.35 vs 0.30; both < 0.50 confirming a continuous unimodal cloud, not discrete clusters) and worse LOEO stability.",
      "Not recommended: Similar weak silhouette (< 0.35) and loses fixed b_i = 1.0 physiological reference.",
      "Not recommended: Lowest silhouette among unsupervised K=4 methods and sensitive to outlier slopes.",
      "Not recommended: Collapses 2D Yield x Plasticity space into 2 diagonal groups, losing the High-Yield/Low-Slope vs High-Yield/High-Slope distinction."
    ),
    stringsAsFactors = FALSE
  )
  write.csv(fw_clustering_comp, file.path(sec26_plas_dir, "Table_S12C_FW_Quadrant_vs_Unsupervised_Clustering_Comparison.csv"), row.names = FALSE)

  if (requireNamespace("openxlsx", quietly = TRUE)) {
    wb_plas <- openxlsx::createWorkbook()
    for (shp in list(
      list(name = "Table_S12A_FW_Quadrant_Summary", data = quad_summary),
      list(name = "Table_S12B_LOEO_Stability", data = loeo_summary),
      list(name = "Table_S12C_Quadrant_vs_Unsup", data = fw_clustering_comp),
      list(name = "Table_S13_FW_vs_Pheno_Clusters", data = tab_s13_crosstab),
      list(name = "Table_S14_All_200_FW_Slopes", data = fw_models)
    )) {
      openxlsx::addWorksheet(wb_plas, shp$name)
      openxlsx::writeDataTable(wb_plas, shp$name, shp$data, withFilter = TRUE)
      openxlsx::freezePane(wb_plas, shp$name, firstRow = TRUE)
      openxlsx::setColWidths(wb_plas, shp$name, cols = seq_len(ncol(shp$data)), widths = "auto")
    }
    openxlsx::saveWorkbook(wb_plas, file.path(sec26_plas_dir, "Table_S12_to_S14_SUPPLEMENTARY_Finlay_Wilkinson_Plasticity.xlsx"), overwrite = TRUE)
  }
  cat("  -> Section 8 Completed Successfully.\n")
} else {
  cat("  -> Section 8 [SKIPPED - Assets Finalized & Protected]\n")
}

# ==============================================================================
# SECTION 9: MULTI-ENVIRONMENT GGE STABILITY & UNIFIED IDEOTYPE SELECTION
# ==============================================================================
cat(">>> Executing Section 9: Multi-Environment Stability (GGE & Unified Ideotype Index)...\n")

if (isTRUE(RUN_SECTION_9_GGE_STABILITY)) {
  sec27_dir <- file.path(assets_dir, "Section_2.8_GGE_Stability")
  dir.create(sec27_dir, showWarnings = FALSE, recursive = TRUE)

  # Archive redundant TIFF copy if present so only 1 clean GGE figure exists in Section_2.8_GGE_Stability
  redundant_tiff <- file.path(sec27_dir, "Figure_S_GGE_Composite.tiff")
  if (file.exists(redundant_tiff)) {
    arch_sec28 <- file.path(archive_rejected_dir, "Section_2.8_Redundant_Assets")
    dir.create(arch_sec28, showWarnings = FALSE, recursive = TRUE)
    file.copy(redundant_tiff, file.path(arch_sec28, "Figure_S_GGE_Composite.tiff"), overwrite = TRUE)
    unlink(redundant_tiff)
  }

  yld_wide <- blues_df %>% select(Genotype, Env_Factor, GrainYld) %>% pivot_wider(names_from = Env_Factor, values_from = GrainYld) %>% arrange(Genotype)
  yld_mat <- as.matrix(yld_wide[, env_metadata$Env_Name]); rownames(yld_mat) <- yld_wide$Genotype
  yld_std_mat <- scale(yld_mat)

  # Exact environment-centered GGE SVD (symmetrical singular value partitioning, svp = 0.5)
  gge_centered <- scale(yld_mat, center = colMeans(yld_mat), scale = FALSE)
  gge_svd <- svd(gge_centered)
  var_expl_pc <- (gge_svd$d^2 / sum(gge_svd$d^2)) * 100
  g_sc <- sweep(gge_svd$u[, 1:2, drop = FALSE], 2, sqrt(gge_svd$d[1:2]), `*`)
  e_sc <- sweep(gge_svd$v[, 1:2, drop = FALSE], 2, sqrt(gge_svd$d[1:2]), `*`)
  if (cor(g_sc[, 1], rowMeans(yld_mat)) < 0) {
    g_sc[, 1] <- -g_sc[, 1]
    e_sc[, 1] <- -e_sc[, 1]
  }
  aec_vec <- colMeans(e_sc)
  a_u <- aec_vec / sqrt(sum(aec_vec^2))
  o_u <- c(-a_u[2], a_u[1])
  aec_proj <- as.numeric(g_sc %*% a_u)
  aec_orth_signed <- as.numeric(g_sc %*% o_u)
  aec_orth_pen <- abs(aec_orth_signed)

  del_cols <- env_metadata$Env_Name[env_metadata$Location == "Delhi"]
  pz_cols  <- env_metadata$Env_Name[env_metadata$Location == "Dharwad"]

  # Per-zone standardized yield and stability (2 environments per zone)
  nwpz_std_mean <- rowMeans(yld_std_mat[, del_cols, drop = FALSE])
  nwpz_stab_pen <- abs(yld_std_mat[, del_cols[1]] - yld_std_mat[, del_cols[2]]) / sqrt(2)
  pz_std_mean   <- rowMeans(yld_std_mat[, pz_cols, drop = FALSE])
  pz_stab_pen   <- abs(yld_std_mat[, pz_cols[1]] - yld_std_mat[, pz_cols[2]]) / sqrt(2)

  stage2_metrics <- data.frame(
    Genotype = rownames(yld_mat),
    Mean_Yield_4Env_BLUE = rowMeans(yld_mat),
    Mean_Std_Yield_4Env  = rowMeans(yld_std_mat),
    Mean_Yield_Delhi_Randomized = rowMeans(yld_mat[, del_cols, drop = FALSE]),
    Mean_Yield_PZ_Descriptive   = rowMeans(yld_mat[, pz_cols, drop = FALSE]),
    GGE_PC1_Score = g_sc[, 1],
    GGE_PC2_Score = g_sc[, 2],
    GGE_AEC_Projection = aec_proj,
    GGE_AEC_Orthogonal_Signed = aec_orth_signed,
    GGE_AEC_Stability_Penalty = aec_orth_pen,
    NWPZ_Zone_Std_Yield = nwpz_std_mean,
    NWPZ_Zone_Stability_Penalty = nwpz_stab_pen,
    PZ_Zone_Std_Yield = pz_std_mean,
    PZ_Zone_Stability_Penalty = pz_stab_pen,
    stringsAsFactors = FALSE
  ) %>% mutate(
    Stage2_Yield_Stability_Index = (as.numeric(scale(Mean_Std_Yield_4Env)) + as.numeric(scale(-GGE_AEC_Stability_Penalty))) / 2,
    NWPZ_Zone_Ideotype_Z = (as.numeric(scale(NWPZ_Zone_Std_Yield)) + as.numeric(scale(-NWPZ_Zone_Stability_Penalty))) / 2,
    PZ_Zone_Ideotype_Z   = (as.numeric(scale(PZ_Zone_Std_Yield)) + as.numeric(scale(-PZ_Zone_Stability_Penalty))) / 2
  )

  waasb_ref <- tryCatch({
    waasb_model <- waasb(master_raw, env = Env_Factor, gen = Genotype, rep = Replication, resp = c(GrainYld, AUTPC), verbose = FALSE)
    waasb_model$GrainYld$model %>% filter(type == "GEN") %>% select(Genotype = Code, WAASB_PlotLevel_PseudoRep_Ref = WAASB, WAASBY_PlotLevel_PseudoRep_Ref = WAASBY)
  }, error = function(e) {
    data.frame(Genotype = stage2_metrics$Genotype, WAASB_PlotLevel_PseudoRep_Ref = NA_real_, WAASBY_PlotLevel_PseudoRep_Ref = NA_real_)
  })

  tce_med_base <- median(blues_df$TCE[blues_df$Environment == "NWPZ Baseline"], na.rm = TRUE)
  multi_trait_means <- blues_df %>% group_by(Genotype) %>% summarise(
    Mean_TKW_4Env = mean(TKW, na.rm = TRUE), Mean_PT_4Env = mean(PT_Count, na.rm = TRUE),
    Mean_TCE_4Env = mean(TCE, na.rm = TRUE), Mean_AUTPC_norm_4Env = mean(AUTPC_norm, na.rm = TRUE),
    Mean_TCE_Delhi = mean(TCE[Environment %in% c("NWPZ Baseline", "NWPZ Warm")], na.rm = TRUE),
    Mean_AUTPC_norm_Delhi = mean(AUTPC_norm[Environment %in% c("NWPZ Baseline", "NWPZ Warm")], na.rm = TRUE),
    Yield_NWPZ_Baseline = mean(GrainYld[Environment == "NWPZ Baseline"], na.rm = TRUE),
    Yield_NWPZ_Warm     = mean(GrainYld[Environment == "NWPZ Warm"], na.rm = TRUE),
    Yield_PZ_Warm_I     = mean(GrainYld[Environment == "PZ Warm I"], na.rm = TRUE),
    Yield_PZ_Warm_II    = mean(GrainYld[Environment == "PZ Warm II"], na.rm = TRUE),
    Trajectory_Baseline = first(Trajectory_Baseline[Environment == "NWPZ Baseline"]),
    Vigor_Baseline      = first(Vigor_Baseline[Environment == "NWPZ Baseline"]),
    TCE_Class_Baseline  = ifelse(mean(TCE[Environment == "NWPZ Baseline"], na.rm = TRUE) >= tce_med_base, "Higher Conversion", "Lower Conversion"),
    .groups = "drop"
  )

  # Approach 2: Rank on Yield + Exact AEC Stability; keep tillering traits as diagnostic columns
  filtering_matrix <- stage2_metrics %>%
    left_join(waasb_ref, by = "Genotype") %>%
    left_join(tsi_df %>% select(Genotype, TSI_Global = TSI_PZ_WarmII, TSI_Delhi, TSI_Class = Stability), by = "Genotype") %>%
    left_join(multi_trait_means, by = "Genotype") %>%
    mutate(
      Master_Profile  = paste(substr(Trajectory_Baseline, 1, 3), ifelse(Vigor_Baseline == VIGOR_HI, "HiVig", "LoVig"), substr(TSI_Class, 1, 4), sep = "-"),
      z_Mean_Yield    = as.numeric(scale(Mean_Std_Yield_4Env)),
      z_Opt_Yield     = as.numeric(scale(Yield_NWPZ_Baseline)),
      z_NWPZ_Warm_Yld = as.numeric(scale(Yield_NWPZ_Warm)),
      z_PZ_Yield      = as.numeric(scale(Yield_PZ_Warm_II)),
      z_Stability     = as.numeric(scale(-GGE_AEC_Stability_Penalty)),
      z_TSI           = as.numeric(scale(TSI_Global)),
      z_TCE           = as.numeric(scale(Mean_TCE_4Env)),
      z_AUTPC_norm    = as.numeric(scale(Mean_AUTPC_norm_4Env)),
      Composite_Ideotype_Z = (z_Mean_Yield + z_Stability) / 2
    ) %>%
    arrange(desc(Composite_Ideotype_Z)) %>%
    mutate(
      Ideotype_Rank = row_number(),
      Yield_Stability_Rank = rank(-Stage2_Yield_Stability_Index, ties.method = "first"),
      NWPZ_Zone_Rank = rank(-NWPZ_Zone_Ideotype_Z, ties.method = "first"),
      PZ_Zone_Rank   = rank(-PZ_Zone_Ideotype_Z, ties.method = "first"),
      Filtering_Rank = Ideotype_Rank
    )

  write.csv(filtering_matrix, file.path(sec27_dir, "Final_Genotype_Filtering_Matrix.csv"), row.names = FALSE)

  # Export Per-Zone Ideotype Tables (NWPZ Randomized Benchmark vs PZ Descriptive)
  nwpz_zone_table <- filtering_matrix %>% arrange(NWPZ_Zone_Rank)
  pz_zone_table   <- filtering_matrix %>% arrange(PZ_Zone_Rank)
  write.csv(nwpz_zone_table, file.path(sec27_dir, "Table_Ideotype_Ranking_NWPZ_Zone.csv"), row.names = FALSE)
  write.csv(pz_zone_table,   file.path(sec27_dir, "Table_Ideotype_Ranking_PZ_Zone.csv"), row.names = FALSE)

  # ============================================================================
  # 9.1 CUSTOM HIGH-LEGIBILITY 2-PANEL GGE BIPLOT (Figure_S7_GGE_Composite.png)
  # ============================================================================
  # Uniform scalar factor for environment vectors so they match the genotype hull radius (preserving all angles & projections)
  max_g_rad <- max(sqrt(filtering_matrix$GGE_PC1_Score^2 + filtering_matrix$GGE_PC2_Score^2))
  max_e_rad <- max(sqrt(e_sc[, 1]^2 + e_sc[, 2]^2))
  env_scale_fac <- (max_g_rad * 0.92) / max_e_rad

  gge_env_plot_df <- data.frame(
    Env = env_metadata$Env_Name,
    Zone = ifelse(env_metadata$Location == "Delhi", "NWPZ Sector (Cooler / Longer Season)", "PZ Sector (Chronic Peninsular Heat)"),
    PC1 = e_sc[, 1] * env_scale_fac,
    PC2 = e_sc[, 2] * env_scale_fac,
    stringsAsFactors = FALSE
  )

  # Convex Hull (Which-Won-Where Polygon) & Perpendicular Equality Rays from Origin (0,0)
  hull_idx <- chull(filtering_matrix$GGE_PC1_Score, filtering_matrix$GGE_PC2_Score)
  hull_closed_idx <- c(hull_idx, hull_idx[1])
  hull_df <- filtering_matrix[hull_idx, ]
  hull_poly_df <- filtering_matrix[hull_closed_idx, ]

  perp_rays <- list()
  ray_len <- max_g_rad * 1.18
  for (i in seq_along(hull_idx)) {
    p1 <- c(hull_poly_df$GGE_PC1_Score[i], hull_poly_df$GGE_PC2_Score[i])
    p2 <- c(hull_poly_df$GGE_PC1_Score[i + 1], hull_poly_df$GGE_PC2_Score[i + 1])
    edge_vec <- p2 - p1
    # Outward normal from origin (0,0) perpendicular to polygon edge p1->p2
    n_vec <- c(-edge_vec[2], edge_vec[1])
    if (sum(n_vec * ((p1 + p2) / 2)) < 0) n_vec <- -n_vec
    n_unit <- n_vec / sqrt(sum(n_vec^2))
    perp_rays[[i]] <- data.frame(x0 = 0, y0 = 0, x1 = n_unit[1] * ray_len, y1 = n_unit[2] * ray_len)
  }
  perp_rays_df <- bind_rows(perp_rays)

  # Classify genotypes for GGE Panel A (Winning Vertex / Top Sector Specialists / Released Checks)
  key_checks_vec <- c("DBW327", "PBW826", "DBW296", "MACS6222", "HD3354", "HS680", "HD3086")
  gge_gen_plot_df <- filtering_matrix %>%
    mutate(
      Is_Vertex = Genotype %in% hull_df$Genotype,
      Is_Key_Check = Genotype %in% key_checks_vec,
      Adaptation_Group = case_when(
        Ideotype_Rank <= 10 ~ "Top 10 Broad Yield-Stability Ideotypes",
        NWPZ_Zone_Rank <= 10 ~ "Top 10 NWPZ Zone Specialists",
        PZ_Zone_Rank <= 10 ~ "Top 10 PZ Zone Specialists",
        Is_Vertex ~ "Polygon Vertex (Extreme GxE Profile)",
        TRUE ~ "Other Panel Genotypes (n = 200)"
      ),
      Adaptation_Group = factor(Adaptation_Group, levels = c(
        "Top 10 Broad Yield-Stability Ideotypes",
        "Top 10 NWPZ Zone Specialists",
        "Top 10 PZ Zone Specialists",
        "Polygon Vertex (Extreme GxE Profile)",
        "Other Panel Genotypes (n = 200)"
      )),
      Label_PanelA = ifelse(
        Is_Vertex | Genotype %in% c("DBW327", "ADT_159", "DBW296", "ADT_150", "PBW826", "EDPIE_9600"),
        ifelse(Is_Key_Check, paste0(Genotype, " ★"), Genotype),
        NA_character_
      )
    )

  gge_pal_a <- c(
    "Top 10 Broad Yield-Stability Ideotypes" = "#0072B2",
    "Top 10 NWPZ Zone Specialists"           = "#D55E00",
    "Top 10 PZ Zone Specialists"             = "#009E73",
    "Polygon Vertex (Extreme GxE Profile)"   = "#CC79A7",
    "Other Panel Genotypes (n = 200)"        = "grey72"
  )

  # Place environment vector labels at clean outward anchor positions outside the dense central genotype cloud
  gge_env_label_df <- gge_env_plot_df %>%
    mutate(
      Lab_X = case_when(
        Env == "NWPZ Baseline" ~ PC1 + 0.075,
        Env == "NWPZ Warm"     ~ PC1 + 0.105,
        Env == "PZ Warm I"     ~ PC1 - 0.015,
        Env == "PZ Warm II"    ~ PC1 - 0.080
      ),
      Lab_Y = case_when(
        Env == "NWPZ Baseline" ~ PC2 + 0.045,
        Env == "NWPZ Warm"     ~ PC2 + 0.028,
        Env == "PZ Warm I"     ~ PC2 - 0.068,
        Env == "PZ Warm II"    ~ PC2 + 0.058
      )
    )

  gge_gen_label_df <- gge_gen_plot_df %>%
    filter(!is.na(Label_PanelA)) %>%
    mutate(
      nudge_x_val = case_when(
        Genotype == "DBW296"     ~ -0.055,
        Genotype == "PBW826"     ~  0.045,
        Genotype == "EDPIE_9600" ~  0.045,
        Genotype == "ADT_150"    ~ -0.045,
        Genotype == "ADT_159"    ~ -0.055,
        Genotype == "DBW327"     ~  0.055,
        TRUE ~ 0
      ),
      nudge_y_val = case_when(
        Genotype == "DBW296"     ~ -0.068,
        Genotype == "PBW826"     ~ -0.040,
        Genotype == "EDPIE_9600" ~  0.040,
        Genotype == "ADT_150"    ~  0.035,
        Genotype == "ADT_159"    ~  0.045,
        Genotype == "DBW327"     ~ -0.038,
        TRUE ~ 0
      )
    )

  fig_gge_mega <- ggplot(gge_gen_plot_df, aes(x = GGE_PC1_Score, y = GGE_PC2_Score)) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "grey65", linewidth = 0.5) +
    geom_vline(xintercept = 0, linetype = "dotted", color = "grey65", linewidth = 0.5) +
    geom_polygon(data = hull_df, aes(x = GGE_PC1_Score, y = GGE_PC2_Score),
                 fill = "#f0f4f8", color = "grey45", linetype = "dashed", linewidth = 0.65, alpha = 0.45, inherit.aes = FALSE) +
    geom_segment(data = perp_rays_df, aes(x = x0, y = y0, xend = x1, yend = y1),
                 color = "grey50", linetype = "dotdash", linewidth = 0.55, inherit.aes = FALSE) +
    geom_point(data = filter(gge_gen_plot_df, Adaptation_Group == "Other Panel Genotypes (n = 200)"),
               aes(color = Adaptation_Group), size = 2.1, alpha = 0.55) +
    geom_point(data = filter(gge_gen_plot_df, Adaptation_Group != "Other Panel Genotypes (n = 200)"),
               aes(color = Adaptation_Group, shape = Is_Vertex), size = 3.4, alpha = 0.92) +
    geom_segment(data = gge_env_plot_df, aes(x = 0, y = 0, xend = PC1, yend = PC2),
                 arrow = arrow(length = unit(0.28, "cm"), type = "closed"),
                 color = "#8B0000", linewidth = 1.15, inherit.aes = FALSE) +
    geom_segment(data = gge_env_label_df, aes(x = PC1, y = PC2, xend = Lab_X, yend = Lab_Y),
                 color = "#8B0000", linetype = "solid", linewidth = 0.5, inherit.aes = FALSE) +
    geom_label(
      data = gge_env_label_df,
      aes(x = Lab_X, y = Lab_Y, label = Env),
      inherit.aes = FALSE, color = "#8B0000", fill = "#FFF8F6", fontface = "bold", size = 3.45,
      label.size = 0.45, label.padding = unit(0.20, "lines")
    ) +
    ggrepel::geom_text_repel(
      data = gge_gen_label_df,
      aes(label = Label_PanelA, color = Adaptation_Group),
      nudge_x = gge_gen_label_df$nudge_x_val,
      nudge_y = gge_gen_label_df$nudge_y_val,
      fontface = "bold", size = 3.25, bg.color = "white", bg.r = 0.16,
      box.padding = 0.55, point.padding = 0.35, force = 8.0, force_pull = 0.5,
      max.overlaps = Inf, min.segment.length = 0.04, seed = 42, show.legend = FALSE
    ) +
    scale_color_manual(values = gge_pal_a, name = "Genotype Adaptation Group:") +
    scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 17), guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0.15, 0.18))) +
    scale_y_continuous(expand = expansion(mult = c(0.15, 0.16))) +
    pub_theme +
    theme(
      legend.position = "bottom",
      legend.title = element_text(face = "bold", size = 10.2),
      legend.text = element_text(size = 9.5),
      plot.title = element_text(face = "bold", size = 13.2, hjust = 0),
      plot.subtitle = element_text(size = 9.8, color = "grey20", lineheight = 1.12, hjust = 0),
      plot.margin = margin(10, 14, 8, 12)
    ) +
    guides(color = guide_legend(nrow = 2, byrow = TRUE, title.position = "top", override.aes = list(size = 3.8, alpha = 1))) +
    labs(
      title = "A. Mega-Environment Delineation (Which-Won-Where Polygon View)",
      subtitle = sprintf("Environment-centered GGE biplot (PC1 + PC2 = %.1f%% of G + GxE variance); rays partition NWPZ vs. PZ.",
                         sum(var_expl_pc[1:2])),
      x = sprintf("GGE PC1 (%.1f%% of G + GxE Variance)", var_expl_pc[1]),
      y = sprintf("GGE PC2 (%.1f%% of G + GxE Variance)", var_expl_pc[2])
    )

  # Panel B: Average Environment Coordination (AEC) Mean Performance vs. Stability View
  # Compute perpendicular drop-line coordinates onto the AEC Abscissa (a_u) for Top 10 Ideotypes + Key Checks
  aec_highlight_df <- filtering_matrix %>%
    filter(Ideotype_Rank <= 10 | Genotype %in% c("PBW826", "MACS6222", "HS680", "DBW296", "ADT_150")) %>%
    mutate(
      Is_Key_Check = Genotype %in% key_checks_vec,
      Proj_X = GGE_AEC_Projection * a_u[1],
      Proj_Y = GGE_AEC_Projection * a_u[2],
      Short_Label = sprintf("#%d %s%s", Ideotype_Rank, Genotype, ifelse(Is_Key_Check, " ★", ""))
    )

  aec_line_len <- max_g_rad * 1.15
  aec_axis_df <- data.frame(
    x0 = -aec_line_len * a_u[1], y0 = -aec_line_len * a_u[2],
    x1 =  aec_line_len * a_u[1], y1 =  aec_line_len * a_u[2]
  )
  aec_orth_df <- data.frame(
    x0 = -aec_line_len * 0.85 * o_u[1], y0 = -aec_line_len * 0.85 * o_u[2],
    x1 =  aec_line_len * 0.85 * o_u[1], y1 =  aec_line_len * 0.85 * o_u[2]
  )
  aec_mean_env_pt <- data.frame(
    PC1 = mean(gge_env_plot_df$PC1),
    PC2 = mean(gge_env_plot_df$PC2)
  )

  fig_gge_stab <- ggplot(filtering_matrix, aes(x = GGE_PC1_Score, y = GGE_PC2_Score)) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "grey70", linewidth = 0.45) +
    geom_vline(xintercept = 0, linetype = "dotted", color = "grey70", linewidth = 0.45) +
    # AEC Ordinate (perpendicular axis through origin: instability direction)
    geom_segment(data = aec_orth_df, aes(x = x0, y = y0, xend = x1, yend = y1),
                 color = "grey45", linetype = "dashed", linewidth = 0.7, inherit.aes = FALSE) +
    # AEC Abscissa (mean yield direction with arrow)
    geom_segment(data = aec_axis_df, aes(x = x0, y = y0, xend = x1, yend = y1),
                 arrow = arrow(length = unit(0.32, "cm"), type = "closed"),
                 color = "#1B4F72", linewidth = 1.25, inherit.aes = FALSE) +
    # Perpendicular drop-lines from highlighted genotypes onto AEC Abscissa (visualizing |Orth| stability penalty)
    geom_segment(data = aec_highlight_df,
                 aes(x = GGE_PC1_Score, y = GGE_PC2_Score, xend = Proj_X, yend = Proj_Y),
                 color = "#1B4F72", linetype = "dotted", linewidth = 0.65, alpha = 0.85, inherit.aes = FALSE) +
    # All 200 genotypes colored by Composite Ideotype Index (Z)
    geom_point(aes(fill = Composite_Ideotype_Z), shape = 21, color = "grey30", size = 2.5, alpha = 0.80, stroke = 0.3) +
    # Highlight Top 10 + Key Checks with bold black rings
    geom_point(data = aec_highlight_df, aes(fill = Composite_Ideotype_Z),
               shape = 21, color = "black", size = 3.8, stroke = 1.0) +
    # Average Environment target ring on the AEC axis
    geom_point(data = aec_mean_env_pt, aes(x = PC1, y = PC2),
               shape = 10, size = 5.2, color = "#8B0000", stroke = 1.3, inherit.aes = FALSE) +
    ggrepel::geom_label_repel(
      data = aec_highlight_df,
      aes(x = GGE_PC1_Score, y = GGE_PC2_Score, label = Short_Label),
      inherit.aes = FALSE, color = "grey10", fill = "white", fontface = "bold", size = 3.00,
      label.size = 0.28, box.padding = 0.52, point.padding = 0.32, force = 7.5, max.overlaps = Inf,
      min.segment.length = 0.04, seed = 123
    ) +
    annotate("label", x = aec_axis_df$x1 * 0.96, y = aec_axis_df$y1 * 0.96 - 0.052,
             label = "AEC Abscissa\n(Higher Mean Yield \u2192)", color = "#1B4F72", fill = "#EBF5FB",
             fontface = "bold", size = 3.10, label.size = 0.3) +
    scale_fill_gradient2(low = "#4575b4", mid = "#ffffbf", high = "#d73027", midpoint = 0,
                         name = "Composite Ideotype Index (Z):") +
    scale_x_continuous(expand = expansion(mult = c(0.15, 0.18))) +
    scale_y_continuous(expand = expansion(mult = c(0.15, 0.16))) +
    pub_theme +
    theme(
      legend.position = "bottom",
      legend.title = element_text(face = "bold", size = 10.2),
      legend.text = element_text(size = 9.5),
      legend.key.width = unit(1.5, "cm"),
      plot.title = element_text(face = "bold", size = 13.2, hjust = 0),
      plot.subtitle = element_text(size = 9.8, color = "grey20", lineheight = 1.12, hjust = 0),
      plot.margin = margin(10, 14, 8, 12)
    ) +
    labs(
      title = "B. Average Environment Coordination (AEC): Mean Yield vs. Stability",
      subtitle = "Blue arrow = AEC Abscissa (directional 4-env mean yield); dotted drop-lines show orthogonal instability penalties (|Orth|).",
      x = sprintf("GGE PC1 (%.1f%% of G + GxE Variance)", var_expl_pc[1]),
      y = sprintf("GGE PC2 (%.1f%% of G + GxE Variance)", var_expl_pc[2])
    )

  composite_gge <- (fig_gge_mega | fig_gge_stab) +
    plot_annotation(
      caption = "Note: Environment vectors in Panel A are uniformly scaled to the genotype hull radius to optimize visual legibility while preserving exact vector angles, cosine correlations, and directional projections. ★ denotes popular released varieties."
    ) &
    theme(plot.caption = element_text(size = 9.6, color = "grey25", hjust = 0, face = "italic", margin = margin(t = 6)))

  ggsave(file.path(sec27_dir, "Figure_S6_GGE_Composite.png"), plot = composite_gge, width = 16.5, height = 7.8, dpi = 600)

  # ============================================================================
  # 9.2 ZONE-SPECIFIC VS. BROAD ADAPTATION SYNTHESIS TABLE & REDESIGNED FIGURE 8
  # ============================================================================
  # Compute cross-environment SE of grain yield for each genotype (n = 4 environments)
  yld_se_4env <- apply(yld_mat, 1, function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))

  ideo_pop_v2 <- filtering_matrix %>%
    mutate(
      SE_Yield_4Env = yld_se_4env[Genotype],
      Is_Named_Cultivar = !grepl("^(EDPIE|ADT)_", Genotype),
      TSI_Label_Clean = ifelse(TSI_Class == "Stable", "High TSI (TSI \u2265 1.0)", "Low TSI (TSI < 1.0)"),
      Short_Strategy_Readable = paste(
        ifelse(grepl("^Early", Trajectory_Baseline), "Early Peak", "Late Peak"),
        ifelse(Vigor_Baseline == VIGOR_HI, "High Vigor", "Low Vigor"),
        sep = " \u00b7 "
      ),
      Line_Category = case_when(
        Ideotype_Rank <= 10 ~ "Top 10 Yield-Stability Ideotypes",
        Is_Named_Cultivar ~ "Popular Released Variety",
        TRUE ~ "Other Breeding Lines"
      ),
      Adaptation_Class = case_when(
        NWPZ_Zone_Rank <= 15 & PZ_Zone_Rank <= 15 ~ "Trans-Zonal Top-15 Bridge (ADT_159)",
        Ideotype_Rank <= 10 ~ "Top 10 Global Broad Ideotypes",
        NWPZ_Zone_Rank <= 15 ~ "Top 15 NWPZ Zone Specialists",
        PZ_Zone_Rank <= 15 ~ "Top 15 PZ Zone Specialists",
        Is_Named_Cultivar & Ideotype_Rank <= 25 ~ "Popular Released Varieties",
        TRUE ~ "Background Panel (n = 200)"
      )
    )

  # Select Top 10 Global Ideotypes + Top 6 Popular Released Varieties
  top10_ideo <- ideo_pop_v2 %>% filter(Ideotype_Rank <= 10)
  extra_checks <- ideo_pop_v2 %>% filter(Ideotype_Rank > 10, Is_Named_Cultivar) %>% slice_head(n = 6)
  ideo_subset_v2 <- bind_rows(top10_ideo, extra_checks) %>% arrange(Ideotype_Rank)

  write.csv(ideo_subset_v2, file.path(sec27_dir, "Table_Top10_vs_Checks_Ideotype_Summary_v2.csv"), row.names = FALSE)

  # Create Table_Zone_Specific_vs_Broad_Adaptation_Synthesis.csv contrasting Top 15 NWPZ vs Top 15 PZ vs Top 10 Global
  top15_nwpz <- ideo_pop_v2 %>% filter(NWPZ_Zone_Rank <= 15) %>% arrange(NWPZ_Zone_Rank)
  top15_pz   <- ideo_pop_v2 %>% filter(PZ_Zone_Rank <= 15) %>% arrange(PZ_Zone_Rank)
  top5_checks <- ideo_pop_v2 %>% filter(Is_Named_Cultivar) %>% arrange(Ideotype_Rank) %>% slice_head(n = 5)

  zone_synthesis_table <- bind_rows(
    top15_nwpz %>% mutate(Selection_Target_Group = "1. NWPZ Zone Top 15 Specialists (Cooler / Longer Season)", Group_Rank = NWPZ_Zone_Rank),
    top15_pz   %>% mutate(Selection_Target_Group = "2. PZ Zone Top 15 Specialists (Chronic Peninsular Heat)", Group_Rank = PZ_Zone_Rank),
    top10_ideo %>% mutate(Selection_Target_Group = "3. Global 4-Env Top 10 Broad Ideotypes (Trans-Zonal)", Group_Rank = Ideotype_Rank),
    top5_checks %>% mutate(Selection_Target_Group = "4. Popular Released Varieties (Top 5 Benchmark)", Group_Rank = row_number())
  ) %>%
    select(
      Selection_Target_Group, Group_Rank, Genotype, Is_Named_Cultivar,
      Ideotype_Rank, NWPZ_Zone_Rank, PZ_Zone_Rank,
      Master_Profile, Trajectory_Baseline, Vigor_Baseline, TSI_Class,
      Mean_Yield_4Env_BLUE, Yield_NWPZ_Baseline, Yield_NWPZ_Warm, Yield_PZ_Warm_I, Yield_PZ_Warm_II,
      TSI_Global, Mean_TCE_4Env, Mean_AUTPC_norm_4Env, Composite_Ideotype_Z, NWPZ_Zone_Ideotype_Z, PZ_Zone_Ideotype_Z
    )

  write.csv(zone_synthesis_table, file.path(sec27_dir, "Table_Zone_Specific_vs_Broad_Adaptation_Synthesis.csv"), row.names = FALSE)

  if (requireNamespace("openxlsx", quietly = TRUE)) {
    wb_s27 <- openxlsx::createWorkbook()
    checks_all_named <- ideo_pop_v2 %>% filter(Is_Named_Cultivar) %>% arrange(Ideotype_Rank)
    for (sh27 in list(
      list(name = "Table_8_Top10_and_Checks", data = ideo_subset_v2),
      list(name = "Table_8B_Zone_vs_Broad_Synth", data = zone_synthesis_table),
      list(name = "Table_S16A_All_200_Ideotypes", data = ideo_pop_v2),
      list(name = "Table_S16B_Named_Cultivars", data = checks_all_named),
      list(name = "Table_S16C_NWPZ_Zone_Rank", data = nwpz_zone_table),
      list(name = "Table_S16D_PZ_Zone_Rank", data = pz_zone_table)
    )) {
      openxlsx::addWorksheet(wb_s27, sh27$name)
      openxlsx::writeDataTable(wb_s27, sh27$name, sh27$data, withFilter = TRUE)
      openxlsx::freezePane(wb_s27, sh27$name, firstRow = TRUE)
      openxlsx::setColWidths(wb_s27, sh27$name, cols = seq_len(ncol(sh27$data)), widths = "auto")
    }
    openxlsx::saveWorkbook(wb_s27, file.path(sec27_dir, "Table_8_MAIN_TEXT_and_Table_S16_SUPPLEMENTARY_Ideotype_Selection.xlsx"), overwrite = TRUE)
  }

  # ============================================================================
  # 9.3 REDESIGNED 3-PANEL PUBLICATION FIGURE: Figure_Top_Genotypes_Ideotype_Profiles_v2.png
  # ============================================================================
  # Panel A (Top-Left): Decluttered Trans-Zonal Adaptation Map (NWPZ vs. PZ Selection Indices)
  zone_map_df <- ideo_pop_v2 %>%
    mutate(
      Zone_Highlight = case_when(
        NWPZ_Zone_Rank <= 15 & PZ_Zone_Rank <= 15 ~ "Trans-Zonal Bridge (ADT_159)",
        Ideotype_Rank <= 10 ~ "Top 10 Global Broad Ideotypes",
        NWPZ_Zone_Rank <= 15 ~ "Top 15 NWPZ Specialists (87% Low TSI)",
        PZ_Zone_Rank <= 15 ~ "Top 15 PZ Specialists (100% High TSI)",
        TRUE ~ "Other Panel Genotypes"
      ),
      Zone_Highlight = factor(Zone_Highlight, levels = c(
        "Trans-Zonal Bridge (ADT_159)",
        "Top 10 Global Broad Ideotypes",
        "Top 15 NWPZ Specialists (87% Low TSI)",
        "Top 15 PZ Specialists (100% High TSI)",
        "Other Panel Genotypes"
      )),
      # KISS: Only 6 focal non-overlapping halo text labels so the data points shine!
      Map_Label = case_when(
        Genotype == "ADT_159"    ~ "ADT_159 (#3 Global | Dual-Zone Bridge)",
        Genotype == "EDPIE_9505" ~ "#1 Global: EDPIE_9505",
        Genotype == "DBW327"     ~ "#2 Global: DBW327 ★",
        Genotype == "EDPIE_9600" ~ "#3 NWPZ: EDPIE_9600",
        Genotype == "DBW296"     ~ "#1 NWPZ: DBW296 ★",
        Genotype == "ADT_150"    ~ "#1 PZ: ADT_150",
        TRUE ~ NA_character_
      ),
      map_nudge_x = case_when(
        Genotype == "EDPIE_9505" ~  0.55,
        Genotype == "ADT_159"    ~  0.42,
        Genotype == "DBW327"     ~ -0.45,
        Genotype == "ADT_150"    ~ -0.58,
        Genotype == "EDPIE_9600" ~ -0.35,
        Genotype == "DBW296"     ~  0.32,
        TRUE ~ 0
      ),
      map_nudge_y = case_when(
        Genotype == "EDPIE_9505" ~  0.08,
        Genotype == "ADT_159"    ~  0.32,
        Genotype == "DBW327"     ~ -0.26,
        Genotype == "ADT_150"    ~  0.22,
        Genotype == "EDPIE_9600" ~  0.18,
        Genotype == "DBW296"     ~ -0.25,
        TRUE ~ 0
      )
    )

  zone_pal <- c(
    "Trans-Zonal Bridge (ADT_159)"          = "#D95F02",
    "Top 10 Global Broad Ideotypes"         = "#0072B2",
    "Top 15 NWPZ Specialists (87% Low TSI)" = "#E69F00",
    "Top 15 PZ Specialists (100% High TSI)" = "#009E73",
    "Other Panel Genotypes"                 = "grey75"
  )

  zone_map_label_df <- filter(zone_map_df, !is.na(Map_Label))

  p_fig8_a <- ggplot(zone_map_df, aes(x = NWPZ_Zone_Ideotype_Z, y = PZ_Zone_Ideotype_Z)) +
    annotate("rect", xmin = 0, xmax = Inf, ymin = 0, ymax = Inf, fill = "#EBF5FB", alpha = 0.35) +
    annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = 0, fill = "#FEF9E7", alpha = 0.32) +
    annotate("rect", xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf, fill = "#E8F8F5", alpha = 0.32) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey45", linewidth = 0.6) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey45", linewidth = 0.6) +
    geom_point(data = filter(zone_map_df, Zone_Highlight == "Other Panel Genotypes"),
               aes(color = Zone_Highlight, shape = TSI_Class), size = 2.2, alpha = 0.55) +
    geom_point(data = filter(zone_map_df, Zone_Highlight != "Other Panel Genotypes"),
               aes(color = Zone_Highlight, shape = TSI_Class), size = 3.6, alpha = 0.95) +
    ggrepel::geom_text_repel(
      data = zone_map_label_df,
      aes(label = Map_Label, color = Zone_Highlight),
      nudge_x = zone_map_label_df$map_nudge_x,
      nudge_y = zone_map_label_df$map_nudge_y,
      fontface = "bold", size = 3.25, bg.color = "white", bg.r = 0.16,
      box.padding = 0.55, point.padding = 0.35, force = 7.5, max.overlaps = Inf,
      ylim = c(-1.85, 1.95), min.segment.length = 0.04,
      seed = 42, show.legend = FALSE
    ) +
    annotate("text", x = -2.42, y = 2.34, label = "PZ Heat Specialists (15/15 High TSI \u2265 1.0)",
             hjust = 0, vjust = 1, fontface = "bold", size = 3.25, color = "#006d50") +
    annotate("text", x = 2.14, y = -2.30, label = "NWPZ Specialists (13/15 Low TSI < 1.0)",
             hjust = 1, vjust = 0, fontface = "bold", size = 3.25, color = "#995c00") +
    annotate("text", x = 2.14, y = 2.34, label = "Broad Dual-Zone Adaptation (1/15 Overlap: ADT_159)",
             hjust = 1, vjust = 1, fontface = "bold", size = 3.25, color = "#004c7a") +
    scale_color_manual(values = zone_pal, name = "Adaptation Group:") +
    scale_shape_manual(values = c("Stable" = 17, "Contrasting" = 16),
                       labels = c("Stable" = "High TSI (TSI \u2265 1.0)", "Contrasting" = "Low TSI (TSI < 1.0)"),
                       name = "TSI Class:") +
    scale_x_continuous(limits = c(-2.50, 2.22), breaks = seq(-2, 2, by = 1), expand = expansion(mult = c(0.02, 0.02))) +
    scale_y_continuous(limits = c(-2.40, 2.45), breaks = seq(-2, 2, by = 1), expand = expansion(mult = c(0.02, 0.02))) +
    pub_theme +
    theme(
      axis.title.x = element_text(face = "bold", size = 10.6),
      axis.title.y = element_text(face = "bold", size = 10.6),
      legend.position = "bottom",
      legend.box = "vertical",
      legend.spacing.y = unit(0.08, "cm"),
      legend.title = element_text(face = "bold", size = 9.5),
      legend.text = element_text(size = 8.8),
      plot.title = element_text(face = "bold", size = 12.6, hjust = 0),
      plot.subtitle = element_text(size = 9.2, color = "grey20", hjust = 0),
      plot.margin = margin(8, 10, 6, 10)
    ) +
    guides(
      color = guide_legend(nrow = 2, byrow = TRUE, title.position = "top", order = 1, override.aes = list(size = 3.5)),
      shape = guide_legend(nrow = 1, order = 2, override.aes = list(size = 3.5))
    ) +
    labs(
      title = "A. Trans-Zonal Adaptation Trade-Off: NWPZ vs. PZ Selection Indices",
      subtitle = "Top 15 NWPZ vs. Top 15 PZ share only 1 line (ADT_159): NWPZ favors Low TSI (< 1.0), PZ requires High TSI (\u2265 1.0).",
      x = "NWPZ Selection Index (Z-Score: NWPZ Baseline & Warm)",
      y = "PZ Selection Index (Z-Score: PZ Warm I & II)"
    )

  # Panel B (Top-Right): Simplified Single-Metric Horizontal Lollipop/Bar Ranking (KISS: 2 colors only!)
  top5_checks_mean_yld <- mean(top5_checks$Mean_Yield_4Env_BLUE, na.rm = TRUE)
  pop_mean_yld_4env    <- mean(ideo_pop_v2$Mean_Yield_4Env_BLUE, na.rm = TRUE)

  fig8_b_df <- ideo_subset_v2 %>%
    mutate(
      Clean_Y_Label  = sprintf("#%d %s%s", Ideotype_Rank, Genotype, ifelse(Is_Named_Cultivar, " ★", "")),
      Clean_Y_Factor = factor(Clean_Y_Label, levels = rev(Clean_Y_Label)),
      Right_Badge    = sprintf("Z = %+.2f  (%s)", Composite_Ideotype_Z, Short_Strategy_Readable)
    )

  tsi_pal_b <- c(
    "High TSI (TSI \u2265 1.0)" = "#009E73",
    "Low TSI (TSI < 1.0)"  = "#D55E00"
  )

  ref_lines_b_df <- data.frame(
    xint = c(top5_checks_mean_yld, pop_mean_yld_4env),
    Ref_Type = factor(
      c(sprintf("Popular Released Varieties Mean (%.3f)", top5_checks_mean_yld),
        sprintf("Panel Mean (n = 200: %.3f)", pop_mean_yld_4env)),
      levels = c(sprintf("Popular Released Varieties Mean (%.3f)", top5_checks_mean_yld),
                 sprintf("Panel Mean (n = 200: %.3f)", pop_mean_yld_4env))
    )
  )

  p_fig8_b <- ggplot(fig8_b_df, aes(y = Clean_Y_Factor, x = Mean_Yield_4Env_BLUE, color = TSI_Label_Clean)) +
    geom_vline(data = ref_lines_b_df, aes(xintercept = xint, linetype = Ref_Type),
               color = "grey20", linewidth = 0.75) +
    # Horizontal lollipop stem from 0.50 baseline to genotype mean yield
    geom_segment(aes(x = 0.50, xend = Mean_Yield_4Env_BLUE, yend = Clean_Y_Factor),
                 linewidth = 2.4, alpha = 0.32, show.legend = FALSE) +
    # Cross-environment SE whiskers (± 1 SE across the 4 regimes)
    geom_errorbarh(aes(xmin = Mean_Yield_4Env_BLUE - SE_Yield_4Env,
                       xmax = Mean_Yield_4Env_BLUE + SE_Yield_4Env),
                   height = 0.28, linewidth = 0.75, color = "grey25", show.legend = FALSE) +
    # Single clean point per genotype for 4-Environment Mean Yield
    geom_point(size = 3.8, alpha = 0.98) +
    # Single dark-grey text annotation on the right showing Composite Z & Phenology/Vigor profile
    geom_text(aes(x = 0.965, label = Right_Badge),
              hjust = 0, color = "grey15", fontface = "bold", size = 3.05, show.legend = FALSE) +
    scale_color_manual(values = tsi_pal_b, name = "Productive-Tiller Stability Index (TSI):") +
    scale_linetype_manual(
      values = setNames(c("dashed", "dotted"), levels(ref_lines_b_df$Ref_Type)),
      name = "Benchmark Reference:"
    ) +
    scale_x_continuous(limits = c(0.50, 1.24), breaks = seq(0.5, 0.9, by = 0.1), expand = expansion(mult = c(0.01, 0.02))) +
    pub_theme +
    theme(
      axis.text.y = element_text(face = "bold", size = 9.8, color = "grey10"),
      axis.title.x = element_text(face = "bold", size = 10.6),
      axis.title.y = element_text(face = "bold", size = 10.6),
      legend.position = "bottom",
      legend.box = "vertical",
      legend.spacing.y = unit(0.08, "cm"),
      legend.title = element_text(face = "bold", size = 9.5),
      legend.text = element_text(size = 9.0),
      plot.title = element_text(face = "bold", size = 12.6, hjust = 0),
      plot.subtitle = element_text(size = 9.2, color = "grey20", hjust = 0),
      plot.margin = margin(8, 12, 6, 8)
    ) +
    guides(
      color = guide_legend(nrow = 1, order = 1, override.aes = list(size = 4.0, linetype = "blank")),
      linetype = guide_legend(nrow = 1, order = 2, override.aes = list(color = "grey15", linewidth = 0.85))
    ) +
    labs(
      title = "B. Multi-Environment Yield & Selection Index Ranking (★ = Popular Released Variety)",
      subtitle = "4-environment mean grain yield (\u00b1 SE), colored by TSI class and benchmarked against Popular Released Varieties.",
      x = sprintf("4-Environment Mean Grain Yield (%s)                  [Selection Index Z & Strategy]", YIELD_UNIT_LABEL),
      y = "Global Rank & Genotype"
    )

  # Panel C (Bottom Full Width: C1 Multi-Environment Yield Reaction Norms + C2 Standardized Physiological Profiles)
  # Benchmark Elite Breeding Lines representing distinct strategies against the Representative Mean (± SE) of Popular Released Varieties (n = 5: DBW327, HD3354, PBW826, MACS6222, HS680)
  focal_breeding_genos <- c("EDPIE_9505", "EDPIE_9576", "ADT_159", "ADT_150", "EDPIE_9600")
  benchmark_label <- "Popular Released Varieties (Mean \u00b1 SE, n = 5 ★)"
  focal_labels <- c(
    "EDPIE_9505" = "#1 Global: EDPIE_9505 (Early Peak \u00b7 Low Vigor \u00b7 Low TSI | High TCE)",
    "EDPIE_9576" = "#5 Global: EDPIE_9576 (Late Peak \u00b7 High Vigor \u00b7 Low TSI | High TCE)",
    "ADT_159"    = "#3 Global / Bridge: ADT_159 (Late Peak \u00b7 Low Vigor \u00b7 High TSI | Dual-Zone)",
    "ADT_150"    = "#1 PZ Specialist: ADT_150 (Late Peak \u00b7 Low Vigor \u00b7 High TSI | PZ Heat)",
    "EDPIE_9600" = "#3 NWPZ Specialist: EDPIE_9600 (Early Peak \u00b7 High Vigor \u00b7 Low TSI | NWPZ)",
    "CHECK_MEAN" = benchmark_label
  )
  focal_levels <- unname(focal_labels[c(focal_breeding_genos, "CHECK_MEAN")])

  focal_pal <- setNames(
    c("#D55E00", "#E69F00", "#0072B2", "#009E73", "#CC79A7", "#111111"),
    focal_levels
  )
  focal_lty <- setNames(
    c("solid", "solid", "solid", "solid", "dotdash", "22"),
    focal_levels
  )
  focal_lwd <- setNames(
    c(1.15, 1.15, 1.25, 1.15, 1.10, 1.55),
    focal_levels
  )

  env_order_levels <- c(
    "NWPZ Baseline",
    "NWPZ Warm",
    "PZ Warm I",
    "PZ Warm II"
  )

  # C1: Multi-Regime Yield Reaction Norms (5 Elite Breeding Lines vs. Popular Released Varieties Mean ± SE)
  c1_breeding_df <- ideo_pop_v2 %>%
    filter(Genotype %in% focal_breeding_genos) %>%
    select(Genotype, Yield_NWPZ_Baseline, Yield_NWPZ_Warm, Yield_PZ_Warm_I, Yield_PZ_Warm_II) %>%
    pivot_longer(
      cols = c(Yield_NWPZ_Baseline, Yield_NWPZ_Warm, Yield_PZ_Warm_I, Yield_PZ_Warm_II),
      names_to = "Env_Col", values_to = "GrainYld"
    ) %>%
    mutate(SE_Yld = NA_real_)

  c1_check_mean_df <- top5_checks %>%
    select(Genotype, Yield_NWPZ_Baseline, Yield_NWPZ_Warm, Yield_PZ_Warm_I, Yield_PZ_Warm_II) %>%
    pivot_longer(
      cols = c(Yield_NWPZ_Baseline, Yield_NWPZ_Warm, Yield_PZ_Warm_I, Yield_PZ_Warm_II),
      names_to = "Env_Col", values_to = "Val"
    ) %>%
    group_by(Env_Col) %>%
    summarise(
      Genotype = "CHECK_MEAN",
      GrainYld = mean(Val, na.rm = TRUE),
      SE_Yld   = sd(Val, na.rm = TRUE) / sqrt(n()),
      .groups  = "drop"
    )

  c1_env_df <- bind_rows(c1_breeding_df, c1_check_mean_df) %>%
    mutate(
      Env_Ordered = factor(case_when(
        Env_Col == "Yield_NWPZ_Baseline" ~ env_order_levels[1],
        Env_Col == "Yield_NWPZ_Warm"     ~ env_order_levels[2],
        Env_Col == "Yield_PZ_Warm_I"     ~ env_order_levels[3],
        Env_Col == "Yield_PZ_Warm_II"    ~ env_order_levels[4]
      ), levels = env_order_levels),
      Focal_Group = factor(focal_labels[Genotype], levels = focal_levels)
    )

  c1_ribbon_df <- c1_env_df %>% filter(Genotype == "CHECK_MEAN")

  p_fig8_c1 <- ggplot(c1_env_df, aes(x = Env_Ordered, y = GrainYld, group = Focal_Group)) +
    geom_ribbon(
      data = c1_ribbon_df,
      aes(x = Env_Ordered, ymin = GrainYld - SE_Yld, ymax = GrainYld + SE_Yld, group = 1),
      inherit.aes = FALSE, fill = "grey30", alpha = 0.18
    ) +
    geom_errorbar(
      data = c1_ribbon_df,
      aes(x = Env_Ordered, ymin = GrainYld - SE_Yld, ymax = GrainYld + SE_Yld),
      inherit.aes = FALSE, width = 0.10, color = "#111111", linewidth = 0.75
    ) +
    geom_line(aes(color = Focal_Group, linetype = Focal_Group, linewidth = Focal_Group)) +
    geom_point(aes(color = Focal_Group), size = 3.2) +
    scale_color_manual(values = focal_pal, guide = "none") +
    scale_linetype_manual(values = focal_lty, guide = "none") +
    scale_linewidth_manual(values = focal_lwd, guide = "none") +
    scale_x_discrete(expand = expansion(mult = c(0.06, 0.06))) +
    scale_y_continuous(limits = c(0.42, 1.12), breaks = seq(0.5, 1.1, by = 0.2)) +
    pub_theme +
    theme(
      axis.text.x = element_text(face = "bold", size = 9.6),
      axis.title.x = element_text(face = "bold", size = 10.6),
      axis.title.y = element_text(face = "bold", size = 10.6),
      plot.title = element_text(face = "bold", size = 12.2, hjust = 0),
      plot.subtitle = element_text(size = 9.2, color = "grey20", hjust = 0),
      legend.position = "none",
      plot.margin = margin(6, 10, 4, 10)
    ) +
    labs(
      title = "C1. Multi-Regime Yield Reaction Norms vs. Popular Released Varieties",
      subtitle = "Shaded grey band & black dashed curve = Popular Released Varieties (Mean \u00b1 SE, n = 5: DBW327, HD3354, PBW826, MACS6222, HS680).",
      x = "Thermal Regime (Ordered from Coolest NWPZ to Warmest PZ)",
      y = sprintf("Realized Grain Yield (%s)", YIELD_UNIT_LABEL)
    )

  # C2: Standardized Physiological & Selection Profile (5 Elite Breeding Lines vs. Popular Released Varieties Mean ± SE)
  z_cols_clean <- c("Composite_Ideotype_Z", "z_Mean_Yield", "z_Stability", "z_TSI", "z_TCE", "z_AUTPC_norm")
  z_labels_clean <- c(
    "Composite_Ideotype_Z" = "Composite\nIndex (Z)",
    "z_Mean_Yield"         = "4-Env Mean\nYield (Z)",
    "z_Stability"          = "GGE AEC\nStability (Z)",
    "z_TSI"                = "Productive-Tiller\nStability (TSI, Z)",
    "z_TCE"                = "Conversion Eff.\n(TCE, Z)",
    "z_AUTPC_norm"         = "Vegetative Vigor\n(AUTPC, Z)"
  )

  c2_breeding_df <- ideo_pop_v2 %>%
    filter(Genotype %in% focal_breeding_genos) %>%
    select(Genotype, all_of(z_cols_clean)) %>%
    pivot_longer(cols = all_of(z_cols_clean), names_to = "Metric", values_to = "Z_Score") %>%
    mutate(SE_Z = NA_real_)

  c2_check_mean_df <- top5_checks %>%
    select(Genotype, all_of(z_cols_clean)) %>%
    pivot_longer(cols = all_of(z_cols_clean), names_to = "Metric", values_to = "Val") %>%
    group_by(Metric) %>%
    summarise(
      Genotype = "CHECK_MEAN",
      Z_Score  = mean(Val, na.rm = TRUE),
      SE_Z     = sd(Val, na.rm = TRUE) / sqrt(n()),
      .groups  = "drop"
    )

  c2_trait_df <- bind_rows(c2_breeding_df, c2_check_mean_df) %>%
    mutate(
      Metric_Factor = factor(z_labels_clean[Metric], levels = unname(z_labels_clean)),
      Focal_Group   = factor(focal_labels[Genotype], levels = focal_levels)
    )

  c2_ribbon_df <- c2_trait_df %>% filter(Genotype == "CHECK_MEAN")

  p_fig8_c2 <- ggplot(c2_trait_df, aes(x = Metric_Factor, y = Z_Score, group = Focal_Group)) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "grey50", linewidth = 0.6) +
    geom_ribbon(
      data = c2_ribbon_df,
      aes(x = Metric_Factor, ymin = Z_Score - SE_Z, ymax = Z_Score + SE_Z, group = 1),
      inherit.aes = FALSE, fill = "grey30", alpha = 0.18
    ) +
    geom_errorbar(
      data = c2_ribbon_df,
      aes(x = Metric_Factor, ymin = Z_Score - SE_Z, ymax = Z_Score + SE_Z),
      inherit.aes = FALSE, width = 0.10, color = "#111111", linewidth = 0.75
    ) +
    geom_line(aes(color = Focal_Group, linetype = Focal_Group, linewidth = Focal_Group)) +
    geom_point(aes(color = Focal_Group), size = 3.2) +
    scale_color_manual(values = focal_pal, name = "Advanced Breeding Lines (by Strategy) vs. Benchmark:") +
    scale_linetype_manual(values = focal_lty, name = "Advanced Breeding Lines (by Strategy) vs. Benchmark:") +
    scale_linewidth_manual(values = focal_lwd, name = "Advanced Breeding Lines (by Strategy) vs. Benchmark:") +
    scale_x_discrete(expand = expansion(mult = c(0.05, 0.05))) +
    pub_theme +
    theme(
      axis.text.x = element_text(face = "bold", size = 9.5),
      axis.title.x = element_text(face = "bold", size = 10.6),
      axis.title.y = element_text(face = "bold", size = 10.6),
      plot.title = element_text(face = "bold", size = 12.2, hjust = 0),
      plot.subtitle = element_text(size = 9.2, color = "grey20", hjust = 0),
      legend.position = "bottom",
      legend.title = element_text(face = "bold", size = 9.8),
      legend.text = element_text(size = 9.2),
      plot.margin = margin(6, 10, 4, 8)
    ) +
    guides(
      color = guide_legend(nrow = 3, byrow = TRUE, override.aes = list(size = 3.4)),
      linetype = guide_legend(nrow = 3, byrow = TRUE),
      linewidth = guide_legend(nrow = 3, byrow = TRUE)
    ) +
    labs(
      title = "C2. Standardized Physiological & Selection Profiles Across Adaptation Strategies",
      subtitle = "Physiological equifinality: High Conversion (EDPIE_9505, EDPIE_9576) vs. High Tiller Stability (ADT_159, ADT_150).",
      x = "Selection Index Components & Diagnostic Tillering Traits",
      y = "Standardized Z-Score (SD Units)"
    )

  p_fig8_bottom <- (p_fig8_c1 | p_fig8_c2) + plot_layout(guides = "collect") &
    theme(legend.position = "bottom", legend.justification = "center")

  p_ideo_comp_v2 <- ((p_fig8_a | p_fig8_b) / p_fig8_bottom) +
    plot_layout(heights = c(1.15, 0.96))

  ggsave(file.path(sec27_dir, "Figure_Top_Genotypes_Ideotype_Profiles_v2.png"), plot = p_ideo_comp_v2, width = 16.5, height = 12.4, dpi = 600)
} else {
  cat("  -> Section 9 [SKIPPED - Assets Finalized & Protected]\n")
}

# ==============================================================================
# SECTION 10: GWAS PHENOTYPE DATASHEET PACKAGE (GENOTYPE 1 TO 200) & ASSET INDEX
# ==============================================================================
cat(">>> Executing Section 10: GWAS Phenotype Datasheets (Plot_ID 1..200) & Master Asset Index...\n")

if (isTRUE(RUN_SECTION_10_GWAS_EXPORT)) {
  gwas_out_dir <- file.path(project_root, "GWAS_Datasheets")
  dir.create(gwas_out_dir, showWarnings = FALSE, recursive = TRUE)

  # Canonical Genotype 1..200 mapping from GS population information.xlsx
  gs_canonical_map <- gs_r1_lookup %>%
    left_join(gs_r2_lookup, by = "Genotype_Clean") %>%
    transmute(
      Plot_ID   = as.integer(Genotype_No),
      Genotype  = Genotype_Clean,
      Sr_No_R1  = as.integer(Sr_R1),
      Sr_No_R2  = as.integer(Sr_R2),
      Block_R1  = as.integer(Block_Pos_GS_R1),
      Block_R2  = as.integer(Block_Pos_GS_R2)
    ) %>%
    arrange(Plot_ID)

  stopifnot(identical(gs_canonical_map$Plot_ID, 1:200))

  # Construct enriched per-environment Stage-1 BLUE GWAS tables ordered strictly by Plot_ID (1..200)
  gwas_env_sheets <- list()
  env_file_map <- c(
    "Delhi_S2"   = "GWAS_DataSheet_1_Delhi_S2_NWPZ_Baseline_Genotype_1_to_200.csv",
    "Delhi_S1"   = "GWAS_DataSheet_2_Delhi_S1_NWPZ_Warm_Genotype_1_to_200.csv",
    "Dharwad_S1" = "GWAS_DataSheet_3_Dharwad_S1_PZ_Warm_I_Genotype_1_to_200.csv",
    "Dharwad_S2" = "GWAS_DataSheet_4_Dharwad_S2_PZ_Warm_II_Genotype_1_to_200.csv"
  )

  for (sh_code in names(env_file_map)) {
    sub_b <- blues_df %>%
      filter(Sheet == sh_code) %>%
      mutate(
        Tmax_Count    = T4_Count,
        Delta_T1_T2   = T2_Count - T1_Count,
        Delta_T2_T3   = T3_Count - T2_Count,
        Delta_T3_T4   = T4_Count - T3_Count,
        Delta_T3_Tmax = T4_Count - T3_Count,
        Delta_T1_T4   = T4_Count - T1_Count,
        Delta_T1_Tmax = T4_Count - T1_Count,
        Rate_T3_Tmax  = Rate_T3_T4,
        Max_Rate      = pmax(Rate_E_T1, Rate_T1_T2, Rate_T2_T3, Rate_T3_T4, na.rm = TRUE),
        Relative_AUTPC = AUTPC_norm
      )

    env_gwas_df <- gs_canonical_map %>%
      left_join(sub_b, by = "Genotype") %>%
      arrange(Plot_ID) %>%
      select(
        Plot_ID, Genotype, Sr_No_R1, Sr_No_R2,
        T1_Count, T2_Count, T3_Count, T4_Count, Tmax_Count, PT_Count,
        GrainYld, TKW, TCE,
        AUTPC, AUTPC_norm, Relative_AUTPC,
        AUTPC_Early, AUTPC_Late, Mean_Tiller_Early, Mean_Tiller_Late, AUTPC_Asymmetry,
        Rate_E_T1, Rate_T1_T2, Rate_T2_T3, Rate_T3_T4, Rate_T3_Tmax, Max_Rate,
        Delta_T1_T2, Delta_T2_T3, Delta_T3_T4, Delta_T3_Tmax, Delta_T1_T4, Delta_T1_Tmax,
        Int_E_T1, Int_T1_T2, Int_T2_T3, Int_T3_T4
      )

    stopifnot(nrow(env_gwas_df) == 200L, identical(env_gwas_df$Plot_ID, 1:200), sum(is.na(env_gwas_df$GrainYld)) == 0L)
    write.csv(env_gwas_df, file.path(gwas_out_dir, env_file_map[[sh_code]]), row.names = FALSE)
    # Also write standard short name for direct script compatibility
    write.csv(env_gwas_df, file.path(gwas_out_dir, paste0("DataSheet_", sh_code, ".csv")), row.names = FALSE)
    gwas_env_sheets[[sh_code]] <- env_gwas_df
  }

  # Cross-Environment Plasticity, Stability & Selection Index GWAS Datasheet (Plot_ID 1..200)
  sec23_tsi_path <- file.path(assets_dir, "Section_2.3_Phenology_Indices", "Section_2.3_Manuscript_and_Supplementary_Tables.xlsx")
  tsi_geno_tbl <- readxl::read_excel(sec23_tsi_path, sheet = "Table_S6B_Genotype_TSI") %>%
    select(Genotype, TSI_NWPZ_Warm, TSI_PZ_Warm_I, TSI_PZ_Warm_II, Retention_Ratio_PZ_Warm_II)

  fw_geno_tbl <- read.csv(file.path(assets_dir, "Section_2.6_Phenotypic_Plasticity", "Finlay_Wilkinson_Classes.csv"), stringsAsFactors = FALSE) %>%
    select(Genotype, FW_Slope_OLS = FW_Slope, FW_Slope_SE, FW_R2, FW_Slope_BLUP = BLUP_Slope, FW_Slope_LOEO = Slope_LOGO, FW_Quadrant = Quadrant)

  filt_geno_tbl <- read.csv(file.path(assets_dir, "Section_2.8_GGE_Stability", "Final_Genotype_Filtering_Matrix.csv"), stringsAsFactors = FALSE) %>%
    select(
      Genotype,
      Mean_Yield_4Env_BLUE, Mean_Std_Yield_4Env, Mean_Yield_Delhi_Randomized, Mean_Yield_PZ_Descriptive,
      Mean_PT_4Env, Mean_TCE_4Env, Mean_AUTPC_norm_4Env, Mean_TKW_4Env,
      TSI_Global, TSI_Delhi, TSI_Class,
      GGE_PC1_Score, GGE_PC2_Score, GGE_AEC_Projection, GGE_AEC_Orthogonal_Signed, GGE_AEC_Stability_Penalty,
      NWPZ_Zone_Std_Yield, NWPZ_Zone_Stability_Penalty, NWPZ_Zone_Ideotype_Z, NWPZ_Zone_Rank,
      PZ_Zone_Std_Yield, PZ_Zone_Stability_Penalty, PZ_Zone_Ideotype_Z, PZ_Zone_Rank,
      Composite_Ideotype_Z, Ideotype_Rank,
      Trajectory_Baseline, Vigor_Baseline, TCE_Class_Baseline, Master_Profile
    )

  gwas_cross_env_df <- gs_canonical_map %>%
    left_join(tsi_geno_tbl, by = "Genotype") %>%
    left_join(fw_geno_tbl, by = "Genotype") %>%
    left_join(filt_geno_tbl, by = "Genotype") %>%
    arrange(Plot_ID)

  stopifnot(nrow(gwas_cross_env_df) == 200L, identical(gwas_cross_env_df$Plot_ID, 1:200))
  write.csv(gwas_cross_env_df, file.path(gwas_out_dir, "GWAS_DataSheet_5_CrossEnv_Plasticity_Stability_and_Ideotype_Indices_Genotype_1_to_200.csv"), row.names = FALSE)
  write.csv(gwas_cross_env_df, file.path(gwas_out_dir, "GWAS_Combined_Global_Stability_Indices.csv"), row.names = FALSE)

  # Wide Multi-Environment All-4-Environments GWAS Matrix (Plot_ID 1..200) for GAPIT / rMVP / TASSEL
  core_gwas_traits <- c(
    "T1_Count", "T2_Count", "T3_Count", "T4_Count", "PT_Count",
    "GrainYld", "TKW", "TCE", "AUTPC", "AUTPC_norm",
    "Mean_Tiller_Early", "Mean_Tiller_Late", "AUTPC_Asymmetry",
    "Rate_E_T1", "Rate_T1_T2", "Rate_T2_T3", "Rate_T3_T4", "Max_Rate"
  )
  env_suffix_map <- c(
    "Delhi_S2"   = "NWPZ_Baseline",
    "Delhi_S1"   = "NWPZ_Warm",
    "Dharwad_S1" = "PZ_Warm_I",
    "Dharwad_S2" = "PZ_Warm_II"
  )

  gwas_wide_all4 <- gs_canonical_map
  for (sh_code in names(env_suffix_map)) {
    sfx <- env_suffix_map[[sh_code]]
    sub_w <- gwas_env_sheets[[sh_code]] %>%
      select(Plot_ID, all_of(core_gwas_traits))
    names(sub_w)[-1] <- paste0(names(sub_w)[-1], "_", sfx)
    gwas_wide_all4 <- gwas_wide_all4 %>% left_join(sub_w, by = "Plot_ID")
  }
  gwas_wide_all4 <- gwas_wide_all4 %>%
    left_join(
      gwas_cross_env_df %>% select(
        Plot_ID,
        TSI_NWPZ_Warm, TSI_PZ_Warm_I, TSI_PZ_Warm_II, TSI_Global,
        FW_Slope_OLS, FW_Slope_BLUP,
        GGE_PC1_Score, GGE_PC2_Score, GGE_AEC_Stability_Penalty,
        Mean_Yield_4Env_BLUE, Mean_PT_4Env, Mean_TCE_4Env, Mean_AUTPC_norm_4Env, Mean_TKW_4Env,
        NWPZ_Zone_Ideotype_Z, PZ_Zone_Ideotype_Z, Composite_Ideotype_Z
      ),
      by = "Plot_ID"
    ) %>%
    arrange(Plot_ID)

  stopifnot(nrow(gwas_wide_all4) == 200L, identical(gwas_wide_all4$Plot_ID, 1:200))
  write.csv(gwas_wide_all4, file.path(gwas_out_dir, "GWAS_Master_Wide_All_4_Environments_Genotype_1_to_200.csv"), row.names = FALSE)

  if (requireNamespace("openxlsx", quietly = TRUE)) {
    wb_gwas <- openxlsx::createWorkbook()
    for (sh_g in list(
      list(name = "1_NWPZ_Baseline_Delhi_S2", data = gwas_env_sheets[["Delhi_S2"]]),
      list(name = "2_NWPZ_Warm_Delhi_S1",     data = gwas_env_sheets[["Delhi_S1"]]),
      list(name = "3_PZ_Warm_I_Dharwad_S1",   data = gwas_env_sheets[["Dharwad_S1"]]),
      list(name = "4_PZ_Warm_II_Dharwad_S2",  data = gwas_env_sheets[["Dharwad_S2"]]),
      list(name = "5_CrossEnv_Plasticity_Stab", data = gwas_cross_env_df),
      list(name = "6_Wide_All_4_Environments",  data = gwas_wide_all4),
      list(name = "7_Genotype_1_to_200_Key",    data = gs_canonical_map)
    )) {
      openxlsx::addWorksheet(wb_gwas, sh_g$name)
      openxlsx::writeDataTable(wb_gwas, sh_g$name, sh_g$data, withFilter = TRUE)
      openxlsx::freezePane(wb_gwas, sh_g$name, firstRow = TRUE)
      openxlsx::setColWidths(wb_gwas, sh_g$name, cols = seq_len(ncol(sh_g$data)), widths = "auto")
    }
    openxlsx::saveWorkbook(wb_gwas, file.path(gwas_out_dir, "GWAS_Phenotype_Master_Workbook_All_4_Environments_Genotype_1_to_200.xlsx"), overwrite = TRUE)
  }
}

# ==============================================================================
# WRAP-UP & PIPELINE COMPLETION
# ==============================================================================
cat("\n==============================================================================\n")
cat(" PIPELINE EXECUTED SUCCESSFULLY!\n")
cat(" All analytical sections, models, and assets have been generated.\n")
cat(" Assets saved to:", assets_dir, "\n")
cat("==============================================================================\n")
