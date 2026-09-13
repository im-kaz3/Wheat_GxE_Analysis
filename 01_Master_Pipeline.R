# ==============================================================================
# MASTER PIPELINE: WHEAT GxE, MICRO-PHENOLOGY, AND CHRONIC VEGETATIVE STRESS
# Repository: Physiological and Molecular Basis of Tiller Number In Wheat Under Heat Stress
# Author: Apoorva Ashu, ICAR-IARI, New Delhi
# Version: 13-09-2026 (v5)
# Description: Evaluates 200 wheat genotypes under chronic vegetative thermal 
#              stress. Pre-calculates exact empirical tiller dynamics & AUTPC,
#              extracts BLUEs safely, computes TCE, TSI, and TCESI, dynamically
#              isolates relative heritability,
#              and executes multi-dimensional variance partitioning.
# ==============================================================================

# ==============================================================================
# TABLE OF CONTENTS / INDEX
# ==============================================================================
# 1. SETUP, LIBRARIES, PUBLICATION THEME & ENVIRONMENT METADATA
# 2. RAW DATA DERIVATIONS, DIAGNOSTICS & NORMALIZATION 
# 3. GLOBAL MERGE, CATEGORICAL TAGGING & COMPOSITES (5-Panel Rainclouds)
# 4. HIGH-IMPACT PHYSIOLOGICAL VISUALIZATIONS 
# 5. STEP-BY-STEP TCE, TSI VALIDATION & SELECTION MATRICES 
# 6. GLOBAL GxE VARIANCE PARTITIONING (ANOVA) - Dynamic Heritability Targets
# 7. SINGLE LOCATION ANALYSIS (LOCALIZED GxE) 
# 8. GGE MEAN VS STABILITY ANALYSIS 
# 9. METEOROLOGICAL ANALYSIS: THERMAL VELOCITY & GDD 
# 10. PHENOLOGICAL CLUSTER ANALYSIS (K-MEANS & PCA) + WINDOW-SPECIFIC GDD
# 11. SINK-CAPACITY RATIOS & CLUSTER SUPERIMPOSITION
# 12. THERMAL ESCAPE & DURATION COMPENSATION (ANOVA & EMMEANS)
# 13. STAGE-WISE PREDICTIVE POWER & INCREMENTAL R2 (Adjusted)
# 14. MASTER DATA EXPORT
# 15. STRUCTURAL EQUATION MODELING (SEM) / PATH ANALYSIS
# 16. NON-LINEAR THERMAL MODELING (LOGISTIC GDD CURVES)
# ==============================================================================

# ==========================================
# 0A. WORKING DIRECTORY & VERSIONED OUTPUT SETUP ####
# ==========================================
# Project root: always the GxE Analysis folder, regardless of how R was launched.
project_root <- "D:/PhD School/Research/Wheat_GxE_Analysis"

# Auto-generate a date-stamped output subdirectory (mirrors your existing naming convention).
# e.g. on 2026-09-13 this produces: Version_V5_13-09-26
output_version <- paste0("Version_V5_", format(Sys.Date(), "%d-%m-%y"))
output_dir     <- file.path(project_root, output_version)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Switch into the output directory so every relative path (Diagnostics/,
# SEM_Path_Analysis/, write.csv(), ggsave()) lands here automatically.
setwd(output_dir)
cat("=== Pipeline Output Directory ===\n", getwd(), "\n=================================\n")

# ==========================================
# 0. REPRODUCIBILITY: PACKAGE CHECK/INSTALL ####
# ==========================================
required_packages <- c(
  "lme4", "car", "emmeans", "readxl", "dplyr", "tidyr", "tibble", "purrr",
  "stringr", "lubridate", "metan", "broom", "scales",
  "ggplot2", "ggrepel", "ggpubr", "patchwork", "cowplot", "ggdist",
  "ggalluvial", "ggExtra", "ggsci", "GGally"
)
missing_packages <- setdiff(required_packages, rownames(installed.packages()))
if (length(missing_packages) > 0) install.packages(missing_packages)
invisible(lapply(required_packages, library, character.only = TRUE))

pub_theme <- theme_classic(base_size = 12, base_family = "Helvetica") +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 12),
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

# ==========================================
# 1B. CENTRAL ENVIRONMENT METADATA ####
# ==========================================
env_metadata <- data.frame(
  Sheet          = c("Delhi_S1", "Delhi_S2", "Dharwad_S1", "Dharwad_S2"),
  Location       = c("Delhi", "Delhi", "Dharwad", "Dharwad"),
  Sowing         = c("S1", "S2", "S1", "S2"),
  Int_E_T1       = c(26, 27, 10, 12),
  Int_T1_T2      = c(15, 15, 10, 10),
  Int_T2_T3      = c(15, 15, 10, 10),
  Int_T3_Tmax    = c(15, 15, 10, 10),
  Total_T1_Tmax  = c(45, 45, 30, 30),
  Start_Date     = as.Date(c("2025-11-12", "2025-12-04", "2025-11-15", "2025-11-27")),
  End_Date       = as.Date(c("2026-02-05", "2026-02-20", "2026-01-14", "2026-01-29")),
  stringsAsFactors = FALSE
)
env_metadata$Environment <- paste0(env_metadata$Location, "_", env_metadata$Sowing)


# ==========================================
# 2. RAW DATA DERIVATIONS, DIAGNOSTICS & NORMALIZATION ####
# ==========================================
# (Refactored to loop over all environments automatically and cache REML models)

cache_dir <- file.path(project_root, "Cache_BLUEs")
dir.create(cache_dir, showWarnings = FALSE)
dir.create("Diagnostics", showWarnings = FALSE)

all_blues <- list()

for (current_sheet in env_metadata$Sheet) {
  cache_file <- file.path(cache_dir, paste0("BLUEs_", current_sheet, ".rds"))
  
  if (file.exists(cache_file)) {
    message(paste("\n--- Loading cached BLUEs for", current_sheet, "---"))
    blues_matrix <- readRDS(cache_file)
  } else {
    message(paste("\n--- Computing BLUEs for", current_sheet, "---"))
    
    trial_data <- read_excel("D:/PhD School/Research/FinalDataYear1.xlsx", sheet = current_sheet)
    trial_data <- trial_data %>% filter(!is.na(Genotype))
    
    trial_data$Genotype <- as.factor(trial_data$Genotype)
    trial_data$Replication <- as.factor(trial_data$Replication)
    trial_data$Block <- as.factor(trial_data$Block)
    
    meta <- env_metadata %>% filter(Sheet == current_sheet)
    current_location <- meta$Location
    current_sowing   <- meta$Sowing
    
    trial_data <- trial_data %>%
      mutate(
        Int_E_T1      = meta$Int_E_T1,
        Int_T1_T2     = meta$Int_T1_T2,
        Int_T2_T3     = meta$Int_T2_T3,
        Int_T3_Tmax   = meta$Int_T3_Tmax,
        Total_T1_Tmax = meta$Total_T1_Tmax,
        
        Delta_T1_T2   = T2_Count - T1_Count,
        Delta_T2_T3   = T3_Count - T2_Count,
        Delta_T3_Tmax = Tmax_Count - T3_Count,
        Delta_T1_Tmax = Tmax_Count - T1_Count,
        
        AUTPC = (0.5 * (0 + T1_Count) * Int_E_T1) +
          (0.5 * (T1_Count + T2_Count) * Int_T1_T2) +
          (0.5 * (T2_Count + T3_Count) * Int_T2_T3) +
          (0.5 * (T3_Count + Tmax_Count) * Int_T3_Tmax),
        
        Rate_E_T1    = T1_Count / Int_E_T1,
        Rate_T1_T2   = Delta_T1_T2 / Int_T1_T2,
        Rate_T2_T3   = Delta_T2_T3 / Int_T2_T3,
        Rate_T3_Tmax = Delta_T3_Tmax / Int_T3_Tmax
      )
    
    design_cols <- c("Sr_No", "Genotype", "Replication", "Block", "Int_E_T1", "Int_T1_T2", "Int_T2_T3", "Int_T3_Tmax", "Total_T1_Tmax")
    trait_cols <- setdiff(colnames(trial_data), design_cols)
    trait_cols <- trait_cols[sapply(trial_data[trait_cols], is.numeric)]
    
    blues_matrix <- data.frame(Genotype = unique(trial_data$Genotype))
    
    normality_results <- data.frame()
    for (trait in trait_cols) {
      shapiro_p <- tryCatch(shapiro.test(trial_data[[trait]])$p.value, error = function(e) NA)
      levene_p <- tryCatch(leveneTest(trial_data[[trait]], trial_data$Replication)$`Pr(>F)`[1], error = function(e) NA)
      normality_results <- rbind(normality_results, data.frame(Trait = trait, Shapiro_P = shapiro_p, Levene_P = levene_p))
    }
    write.csv(normality_results, paste0("Diagnostics/Normality_Tests_", current_location, "_", current_sowing, ".csv"), row.names = FALSE)
    
    model_tier_log <- data.frame()
    for (trait in trait_cols) {
      formula_alpha <- paste(trait, "~ Genotype + (1 | Replication/Block)")
      formula_rcbd  <- paste(trait, "~ Genotype + (1 | Replication)")
      formula_lm    <- paste(trait, "~ Genotype + Replication")
      
      model_alpha <- try(lmer(as.formula(formula_alpha), data = trial_data, na.action = na.exclude,
                              control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))), silent = TRUE)
      
      if (!inherits(model_alpha, "try-error") && !isSingular(model_alpha)) {
        model <- model_alpha
        tier_used <- "Alpha-Lattice"
      } else {
        model_rcbd <- try(lmer(as.formula(formula_rcbd), data = trial_data, na.action = na.exclude,
                               control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))), silent = TRUE)
        if (!inherits(model_rcbd, "try-error") && !isSingular(model_rcbd)) {
          model <- model_rcbd
          tier_used <- "RCBD"
        } else {
          model <- try(lm(as.formula(formula_lm), data = trial_data, na.action = na.exclude), silent = TRUE)
          tier_used <- "Fixed-effects LM (last resort)"
          if (inherits(model, "try-error")) {
            warning(paste("Trait", trait, "could not be modeled by any tier - skipping."))
            model_tier_log <- rbind(model_tier_log, data.frame(Trait = trait, Model_Used = "FAILED - skipped"))
            next
          }
        }
      }
      
      model_tier_log <- rbind(model_tier_log, data.frame(Trait = trait, Model_Used = tier_used))
      
      res <- residuals(model)
      p_qq <- ggplot(data.frame(res), aes(sample = res)) +
        stat_qq() + stat_qq_line(color = "red") +
        theme_bw() + labs(title = paste("Q-Q Plot:", trait))
      p_hist <- ggplot(data.frame(res), aes(x = res)) +
        geom_histogram(bins = 30, fill = "steelblue", color = "black", alpha = 0.7) +
        theme_bw() + labs(title = paste("Residuals:", trait))
      
      diag_plot <- ggarrange(p_qq, p_hist, ncol = 2)
      ggsave(paste0("Diagnostics/", trait, "_Normalization_", current_location, "_", current_sowing, ".png"), plot = diag_plot, width = 8, height = 4)
      
      em_out <- suppressMessages(as.data.frame(emmeans(model, "Genotype")))
      trait_blues <- em_out[, c("Genotype", "emmean")]
      colnames(trait_blues)[2] <- trait
      blues_matrix <- merge(blues_matrix, trait_blues, by = "Genotype", all.x = TRUE)
    }
    
    write.csv(model_tier_log, paste0("Diagnostics/Model_Tier_Used_", current_location, "_", current_sowing, ".csv"), row.names = FALSE)
    
    blues_matrix$Location <- current_location
    blues_matrix$Sowing <- current_sowing
    blues_matrix$Environment <- paste0(current_location, "_", current_sowing)
    
    # Save to cache
    saveRDS(blues_matrix, cache_file)
  }
  
  # Append to master list
  all_blues[[current_sheet]] <- blues_matrix
}


# ==========================================
# 3. GLOBAL MERGE, CATEGORICAL TAGGING & COMPOSITES ####
# ==========================================
master_gxe <- do.call(rbind, all_blues)
master_gxe$Genotype <- as.factor(master_gxe$Genotype)
master_gxe$Location <- as.factor(master_gxe$Location)
master_gxe$Sowing <- as.factor(master_gxe$Sowing)
master_gxe$Environment <- as.factor(master_gxe$Environment)

# ---------------------------------------------------------
# 3A. PEAK-PHASE CLASSIFICATION 
# ---------------------------------------------------------
master_gxe <- master_gxe %>%
  mutate(
    Max_Rate = pmax(Rate_E_T1, Rate_T1_T2, Rate_T2_T3, Rate_T3_Tmax, na.rm = TRUE),
    Peak_Phase = case_when(
      Max_Rate <= 0 ~ "Continuous Decline", # INTERCEPT: Genotypes actively dying/aborting across all intervals
      Rate_E_T1 == Max_Rate ~ "Very Early (E-T1)",
      Rate_T1_T2 == Max_Rate ~ "Early Vigor (T1-T2)",
      Rate_T2_T3 == Max_Rate ~ "Late Vigor (T2-T3)",
      Rate_T3_Tmax == Max_Rate ~ "Terminal Vigor (T3-Tmax)",
      TRUE ~ "Uniform"   
    ),
    Phase_T1_T2   = case_when(Delta_T1_T2 > 0.05 ~ "Active Emergence", Delta_T1_T2 < -0.05 ~ "Active Abortion", TRUE ~ "Stagnation"),
    Phase_T2_T3   = case_when(Delta_T2_T3 > 0.05 ~ "Active Emergence", Delta_T2_T3 < -0.05 ~ "Active Abortion", TRUE ~ "Stagnation"),
    Phase_T3_Tmax = case_when(Delta_T3_Tmax > 0.05 ~ "Active Emergence", Delta_T3_Tmax < -0.05 ~ "Active Abortion", TRUE ~ "Stagnation"),
    Phase_T1_Tmax = case_when(Delta_T1_Tmax > 0.05 ~ "Active Emergence", Delta_T1_Tmax < -0.05 ~ "Active Abortion", TRUE ~ "Stagnation"),
    TCE = PT_Count / Tmax_Count,
    Relative_AUTPC = AUTPC / Tmax_Count
  )

n_all_declining <- sum(master_gxe$Max_Rate <= 0, na.rm = TRUE)
message(paste(n_all_declining, "of", nrow(master_gxe),
              "genotype-environment rows have Max_Rate <= 0 (net decline across all intervals)."))

phase_definitions <- data.frame(
  Phase_Classification = c("Continuous Decline", "Very Early (E-T1)", "Early Vigor (T1-T2)", "Late Vigor (T2-T3)", "Terminal Vigor (T3-Tmax)"),
  Biological_Definition = c(
    "Tillering velocity is negative or zero across all measured intervals (net abortion).",
    "Maximum tillering velocity occurs during the initial emergence phase (E to T1).",
    "Maximum tillering velocity occurs between T1 and T2.",
    "Maximum tillering velocity is delayed until the T2 to T3 interval.",
    "Maximum tillering velocity occurs during the final pre-anthesis interval (T3 to Tmax)."
  )
)
write.csv(phase_definitions, "Peak_Phase_Methodological_Definitions.csv", row.names = FALSE)
write.csv(master_gxe %>% count(Environment, Peak_Phase, name = "N_Genotypes"), "Peak_Phase_Counts_By_Environment.csv", row.names = FALSE)

# ---------------------------------------------------------
# 3B. COMPREHENSIVE PHENOTYPIC SUMMARY TABLE
# ---------------------------------------------------------
summary_table <- data.frame()
numeric_traits <- c("PT_Count", "Tmax_Count", "GrainYld", "TKW", "AUTPC", "TCE")

for (trait in numeric_traits) {
  formula_anova <- paste(trait, "~ Environment")
  model <- tryCatch(lm(as.formula(formula_anova), data = master_gxe), error = function(e) NULL)
  p_val <- ifelse(!is.null(model), anova(model)$"Pr(>F)"[1], NA)
  
  trait_summary <- master_gxe %>%
    group_by(Environment) %>%
    summarise(
      Trait = trait, Mean = mean(get(trait), na.rm = TRUE),
      SD = sd(get(trait), na.rm = TRUE), SE = SD / sqrt(n()),
      Range = paste(round(min(get(trait), na.rm = TRUE), 2), "-", round(max(get(trait), na.rm = TRUE), 2)),
      `CV(%)` = (SD / Mean) * 100, ANOVA_p = p_val, .groups = 'drop'
    )
  summary_table <- rbind(summary_table, trait_summary)
}
write.csv(summary_table, "Phenotypic_Summary_Table.csv", row.names = FALSE)

master_gxe$Environment <- factor(master_gxe$Environment, levels = c("Delhi_S1", "Delhi_S2", "Dharwad_S1", "Dharwad_S2"))

# ---------------------------------------------------------
# 3C. FIGURE 1: MACROSCOPIC PHENOTYPIC RESPONSE (5-PANEL RAINCLOUDS)
# ---------------------------------------------------------
dir.create("Raincloud_Plots_All_Traits", showWarnings = FALSE)
for (trait in numeric_traits) {
  p_rain <- ggplot(master_gxe, aes(x = Environment, y = .data[[trait]], fill = Environment)) +
    stat_halfeye(adjust = 0.5, width = 0.6, .width = 0, justification = -0.3, point_colour = NA) +
    geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.5) +
    geom_point(aes(color = Environment), size = 1.3, alpha = 0.3, position = position_jitter(width = 0.1, seed = 42)) +
    scale_fill_npg() + scale_color_npg() + pub_theme + theme(legend.position = "none") +
    labs(title = paste("Phenotypic Distribution:", trait), x = "Environment", y = trait)
  ggsave(paste0("Raincloud_Plots_All_Traits/", trait, "_Distribution.png"), plot = p_rain, width = 8, height = 6, dpi = 300)
}

fig1a <- ggplot(master_gxe, aes(x = Environment, y = GrainYld, fill = Environment)) +
  stat_halfeye(adjust = 0.5, width = 0.6, .width = 0, justification = -0.3, point_colour = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.5) +
  geom_point(aes(color = Environment), size = 1.3, alpha = 0.3, position = position_jitter(width = 0.1, seed = 42)) +
  scale_fill_npg() + scale_color_npg() + pub_theme + theme(legend.position = "none") +
  labs(title = "A. Grain Yield", y = "Grain Yield (kg/plot)", x = "")

fig1b <- ggplot(master_gxe, aes(x = Environment, y = Tmax_Count, fill = Environment)) +
  stat_halfeye(adjust = 0.5, width = 0.6, .width = 0, justification = -0.3, point_colour = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.5) +
  geom_point(aes(color = Environment), size = 1.3, alpha = 0.3, position = position_jitter(width = 0.1, seed = 42)) +
  scale_fill_npg() + scale_color_npg() + pub_theme + theme(legend.position = "none") +
  labs(title = expression(bold("B. Maximum Tillers ("*T[max]*")")), y = "Tiller Count", x = "")

fig1c <- ggplot(master_gxe, aes(x = Environment, y = PT_Count, fill = Environment)) +
  stat_halfeye(adjust = 0.5, width = 0.6, .width = 0, justification = -0.3, point_colour = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.5) +
  geom_point(aes(color = Environment), size = 1.3, alpha = 0.3, position = position_jitter(width = 0.1, seed = 42)) +
  scale_fill_npg() + scale_color_npg() + pub_theme + theme(legend.position = "none") +
  labs(title = "C. Productive Tillers (PT)", y = "PT Count", x = "")

fig1d <- ggplot(master_gxe, aes(x = Environment, y = TCE, fill = Environment)) +
  stat_halfeye(adjust = 0.5, width = 0.6, .width = 0, justification = -0.3, point_colour = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.5) +
  geom_point(aes(color = Environment), size = 1.3, alpha = 0.3, position = position_jitter(width = 0.1, seed = 42)) +
  scale_fill_npg() + scale_color_npg() + pub_theme + theme(legend.position = "none") +
  labs(title = "D. Tiller Conversion Efficiency", y = "Ratio (PT/Tmax)", x = "")

fig1e <- ggplot(master_gxe, aes(x = Environment, y = TKW, fill = Environment)) +
  stat_halfeye(adjust = 0.5, width = 0.6, .width = 0, justification = -0.3, point_colour = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.5) +
  geom_point(aes(color = Environment), size = 1.3, alpha = 0.3, position = position_jitter(width = 0.1, seed = 42)) +
  scale_fill_npg() + scale_color_npg() + pub_theme + theme(legend.position = "none") +
  labs(title = "E. Thousand Kernel Weight", y = "TKW (g)", x = "")

composite_fig1 <- (fig1a | fig1b) / (fig1c | fig1d) / (fig1e | plot_spacer())
ggsave("Figure_1_Raincloud_Distributions_Composite.png", plot = composite_fig1, width = 12, height = 15, dpi = 600)


# ==========================================
# 4. HIGH-IMPACT PHYSIOLOGICAL VISUALIZATIONS ####
# ==========================================

# ---------------------------------------------------------
# 4A. FIGURE 2A: MICRO-PHENOLOGY (ALLUVIAL) + FIGURE 2C: TRANSITION HEATMAP
# ---------------------------------------------------------
alluvial_data <- master_gxe %>%
  mutate(across(starts_with("Phase_"), ~ case_when(
    . == "Active Emergence" ~ "Emergence",
    . == "Active Abortion"  ~ "Abortion",
    . == "Stagnation"       ~ "Stagnation",
    TRUE ~ .
  ))) %>%
  group_by(Environment, Phase_T1_T2, Phase_T2_T3, Phase_T3_Tmax) %>%
  tally()

fig2a <- ggplot(alluvial_data, aes(y = n, axis1 = Phase_T1_T2, axis2 = Phase_T2_T3, axis3 = Phase_T3_Tmax)) +
  geom_alluvium(aes(fill = Phase_T1_T2), width = 1/12, alpha = 0.7) +
  geom_stratum(width = 1/12, fill = "grey90", color = "black") +
  geom_text(stat = "stratum",
            aes(label = str_wrap(after_stat(stratum), width = 10),
                hjust = after_stat(ifelse(x == 1, 0, ifelse(x == 3, 1, 0.5)))),
            size = 3.2, fontface = "bold", lineheight = 0.85) +
  facet_wrap(~ Environment, scales = "free_y") +
  scale_fill_npg() + pub_theme +
  theme(legend.position = "none", axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.line.y = element_blank(), panel.spacing.x = unit(2, "lines"), panel.spacing.y = unit(2.5, "lines")) +
  scale_x_discrete(limits = c("T1-T2", "T2-T3", "T3-Tmax"), expand = expansion(mult = c(0.22, 0.22))) +
  labs(title = "Developmental Phase Transitions (Crossover Flow)",
       x = "Phenological Interval", y = "Number of Genotypes (n)")
ggsave("Figure_2A_Alluvial_Transitions.png", plot = fig2a, width = 13, height = 8, dpi = 600)

transition_1 <- master_gxe %>%
  count(Environment, From = Phase_T1_T2, To = Phase_T2_T3) %>%
  group_by(Environment, From) %>%
  mutate(Pct = n / sum(n) * 100) %>%
  ungroup() %>%
  mutate(Transition = "T1-T2 -> T2-T3")

transition_2 <- master_gxe %>%
  count(Environment, From = Phase_T2_T3, To = Phase_T3_Tmax) %>%
  group_by(Environment, From) %>%
  mutate(Pct = n / sum(n) * 100) %>%
  ungroup() %>%
  mutate(Transition = "T2-T3 -> T3-Tmax")

transition_data <- bind_rows(transition_1, transition_2) %>%
  mutate(
    From = factor(From, levels = c("Active Emergence", "Stagnation", "Active Abortion")),
    To = factor(To, levels = c("Active Emergence", "Stagnation", "Active Abortion")),
    Transition = factor(Transition, levels = c("T1-T2 -> T2-T3", "T2-T3 -> T3-Tmax"))
  )
write.csv(transition_data, "Phase_Transition_Probability_Table.csv", row.names = FALSE)

transition_heatmap <- ggplot(transition_data, aes(x = To, y = From, fill = Pct)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = paste0(round(Pct, 0), "%\n(n=", n, ")")), size = 3, fontface = "bold", lineheight = 0.85) +
  facet_grid(Transition ~ Environment) +
  scale_fill_gradientn(colors = c("white", "#fdae61", "#d73027"), name = "% of genotypes\nfrom this state") +
  pub_theme + theme(axis.text.x = element_text(angle = 30, hjust = 1), panel.grid = element_blank()) +
  labs(title = "Developmental Phase Transition Probabilities",
       subtitle = "Row = starting state, Column = ending state",
       x = "State at end of interval", y = "State at start of interval")
ggsave("Figure_2C_Transition_Heatmap.png", plot = transition_heatmap, width = 14, height = 7, dpi = 600)

# ---------------------------------------------------------
# 4B. FIGURE 3: YIELD DETERMINANTS WITH MARGINAL DENSITY
# ---------------------------------------------------------
env_list <- levels(master_gxe$Environment)
plot_list <- list()

for (env in env_list) {
  p <- ggplot(master_gxe %>% filter(Environment == env), aes(x = Relative_AUTPC, y = GrainYld)) +
    geom_point(aes(color = Peak_Phase), alpha = 0.6, size = 1.5) +
    geom_smooth(aes(group = 1), method = "lm", color = "black", linewidth = 1) +
    stat_cor(aes(group = 1, label = paste(after_stat(r.label), after_stat(p.label), sep = "~`,`~")),
             label.x.npc = "left", label.y.npc = "top", size = 4) +
    scale_color_npg() + pub_theme + theme(legend.position = "none") +
    labs(title = env, x = "Relative AUTPC", y = "Grain Yield (kg/plot)")
  plot_list[[env]] <- ggMarginal(p, type = "density", groupColour = TRUE, groupFill = TRUE, alpha = 0.4)
}

shared_legend <- get_legend(
  ggplot(master_gxe, aes(x = Relative_AUTPC, y = GrainYld, color = Peak_Phase)) +
    geom_point() + scale_color_npg() + pub_theme
)

fig3_grid <- plot_grid(plot_list[[1]], plot_list[[2]], plot_list[[3]], plot_list[[4]], ncol = 2, labels = "AUTO")
fig3_final <- plot_grid(fig3_grid, shared_legend, ncol = 1, rel_heights = c(1, 0.1))
save_plot("Figure_3_Marginal_Yield_Determinants.png", fig3_final, base_height = 9, base_width = 12, dpi = 600)

# ---------------------------------------------------------
# 4C. AUTPC DECOMPOSITION - SHAPE, NOT JUST AREA
# ---------------------------------------------------------
master_gxe <- master_gxe %>%
  select(-any_of(c("Int_E_T1", "Int_T1_T2", "Int_T2_T3", "Int_T3_Tmax", "Total_T1_Tmax"))) %>%
  left_join(env_metadata %>% select(Environment, Int_E_T1, Int_T1_T2, Int_T2_T3, Int_T3_Tmax), by = "Environment")

master_gxe$Environment <- factor(master_gxe$Environment, levels = c("Delhi_S1", "Delhi_S2", "Dharwad_S1", "Dharwad_S2"))
master_gxe <- master_gxe %>%
  rowwise() %>%
  mutate(
    AUTPC_Early = 0.5 * (0 + T1_Count) * Int_E_T1 + 0.5 * (T1_Count + T2_Count) * Int_T1_T2,
    AUTPC_Late  = 0.5 * (T2_Count + T3_Count) * Int_T2_T3 + 0.5 * (T3_Count + Tmax_Count) * Int_T3_Tmax,
    AUTPC_Asymmetry = (AUTPC_Early - AUTPC_Late) / AUTPC
  ) %>%
  ungroup()

autpc_components <- c("AUTPC", "AUTPC_Early", "AUTPC_Late", "AUTPC_Asymmetry", "Relative_AUTPC")
autpc_yield_corr <- map_dfr(levels(master_gxe$Environment), function(env) {
  d <- master_gxe %>% filter(Environment == env)
  map_dfr(autpc_components, function(v) {
    ct <- suppressWarnings(cor.test(d[[v]], d$GrainYld))
    tibble(Environment = env, Component = v, R = unname(ct$estimate), p_value = ct$p.value)
  })
})
write.csv(autpc_yield_corr, "AUTPC_Component_Yield_Correlations.csv", row.names = FALSE)

autpc_component_plot <- ggplot(autpc_yield_corr, aes(x = Component, y = R, fill = Environment)) +
  geom_col(position = "dodge", color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_npg() + pub_theme + theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  labs(title = "Decomposing AUTPC: Does Shape Predict Yield Better Than Total Area?",
       subtitle = "AUTPC_Asymmetry isolates front-loaded vs. back-loaded tillering, which the whole-curve integral cannot",
       x = "", y = "Pearson R vs. Grain Yield")
ggsave("AUTPC_Decomposition_Plot.png", plot = autpc_component_plot, width = 10, height = 6, dpi = 300)


# ==========================================
# 5. TCE/TSI/TCESI VALIDATION & CLASSICAL CLOSURE ####
# ==========================================
efficiency_base <- master_gxe %>%
  mutate(Stress_Regime = case_when(
    Environment %in% c("Delhi_S2", "Dharwad_S1") ~ "Optimal_Combined",
    Environment %in% c("Delhi_S1", "Dharwad_S2") ~ "Stress_Combined"
  ))

dir.create("Selection_Indices_Exports", showWarnings = FALSE)

# ---------------------------------------------------------
# 5A. REGIONAL BASELINE ANALYSIS (DELHI ONLY)
# ---------------------------------------------------------
delhi_metrics <- efficiency_base %>% filter(Location == "Delhi")
delhi_summary <- delhi_metrics %>%
  group_by(Genotype, Sowing) %>%
  summarise(Mean_PT = mean(PT_Count, na.rm = TRUE), Mean_Yield = mean(GrainYld, na.rm = TRUE), Mean_TCE = mean(TCE, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = Sowing, values_from = c(Mean_PT, Mean_Yield, Mean_TCE))

pop_PT_Delhi_S1 <- mean(delhi_summary$Mean_PT_S1, na.rm = TRUE)
pop_TCE_Delhi_S2 <- mean(delhi_summary$Mean_TCE_S2, na.rm = TRUE)
pop_TCE_Delhi_S1 <- mean(delhi_summary$Mean_TCE_S1, na.rm = TRUE)
mean_Yp_delhi <- mean(delhi_summary$Mean_Yield_S2, na.rm = TRUE)
mean_Ys_delhi <- mean(delhi_summary$Mean_Yield_S1, na.rm = TRUE)
D_delhi <- 1 - (mean_Ys_delhi / mean_Yp_delhi)

delhi_index <- delhi_summary %>%
  mutate(
    Actual_Stress_Yield = Mean_Yield_S1,
    STI_Delhi = (Mean_Yield_S2 * Mean_Yield_S1) / (mean_Yp_delhi^2),
    SSI_Delhi = (1 - (Mean_Yield_S1 / Mean_Yield_S2)) / D_delhi,
    TSI_Delhi = (Mean_PT_S1 / Mean_PT_S2) * (Mean_PT_S1 / pop_PT_Delhi_S1),
    TCESI_Delhi = (Mean_TCE_S1 / Mean_TCE_S2) * (Mean_TCE_S1 / pop_TCE_Delhi_S1),
    Category_Delhi = case_when(
      Mean_TCE_S2 >= pop_TCE_Delhi_S2 & Mean_TCE_S1 >= pop_TCE_Delhi_S1 ~ "Stable High Efficiency",
      Mean_TCE_S2 >= pop_TCE_Delhi_S2 & Mean_TCE_S1 < pop_TCE_Delhi_S1 ~ "Contrasting Susceptible",
      Mean_TCE_S2 < pop_TCE_Delhi_S2 & Mean_TCE_S1 >= pop_TCE_Delhi_S1 ~ "Contrasting Tolerant",
      TRUE ~ "Stable Low Efficiency"
    )
  )

tsi_delhi_plot <- ggplot(delhi_index, aes(x = TSI_Delhi, y = Actual_Stress_Yield)) +
  geom_point(aes(color = Category_Delhi), alpha = 0.7, size = 2.5) +
  geom_smooth(method = "lm", color = "black", fill = "lightgrey") +
  stat_cor(aes(label = paste(after_stat(r.label), after_stat(p.label), sep = "~`,`~")), size = 5, fontface = "bold", label.x.npc = "left", label.y.npc = "top") +
  scale_color_npg() + pub_theme +
  labs(title = "TSI Predicting Actual Stress Yield (Delhi)", x = "TSI (Capacity Stability)", y = "Actual Grain Yield (kg/plot)", color = "TCE Category")
ggsave("Selection_Indices_Exports/TSI_Direct_Yield_Validation_Delhi.png", plot = tsi_delhi_plot, width = 9, height = 6, dpi = 600)

tcesi_delhi_plot <- ggplot(delhi_index, aes(x = TCESI_Delhi, y = Actual_Stress_Yield)) +
  geom_point(aes(color = Category_Delhi), alpha = 0.7, size = 2.5) +
  geom_smooth(method = "lm", color = "black", fill = "lightgrey") +
  stat_cor(aes(label = paste(after_stat(r.label), after_stat(p.label), sep = "~`,`~")), size = 5, fontface = "bold", label.x.npc = "left", label.y.npc = "top") +
  scale_color_npg() + pub_theme +
  labs(title = "TCESI Predicting Actual Stress Yield (Delhi)", x = "TCESI (Efficiency Stability)", y = "Actual Grain Yield (kg/plot)", color = "TCE Category")
ggsave("Selection_Indices_Exports/TCESI_Direct_Yield_Validation_Delhi.png", plot = tcesi_delhi_plot, width = 9, height = 6, dpi = 600)
write.csv(delhi_index, "Selection_Indices_Exports/Selection_Matrix_Delhi.csv", row.names = FALSE)

# ---------------------------------------------------------
# 5B. REGIONAL BASELINE ANALYSIS (DHARWAD ONLY) 
# ---------------------------------------------------------
dharwad_metrics <- efficiency_base %>% filter(Location == "Dharwad")
dharwad_summary <- dharwad_metrics %>%
  group_by(Genotype, Sowing) %>%
  summarise(Mean_PT = mean(PT_Count, na.rm = TRUE), Mean_Yield = mean(GrainYld, na.rm = TRUE), Mean_TCE = mean(TCE, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = Sowing, values_from = c(Mean_PT, Mean_Yield, Mean_TCE))

pop_PT_Dharwad_S2 <- mean(dharwad_summary$Mean_PT_S2, na.rm = TRUE)
pop_TCE_Dharwad_S1 <- mean(dharwad_summary$Mean_TCE_S1, na.rm = TRUE)
pop_TCE_Dharwad_S2 <- mean(dharwad_summary$Mean_TCE_S2, na.rm = TRUE)
mean_Yp_dharwad <- mean(dharwad_summary$Mean_Yield_S1, na.rm = TRUE)
mean_Ys_dharwad <- mean(dharwad_summary$Mean_Yield_S2, na.rm = TRUE)
D_dharwad <- 1 - (mean_Ys_dharwad / mean_Yp_dharwad)

dharwad_index <- dharwad_summary %>%
  mutate(
    Actual_Stress_Yield = Mean_Yield_S2,
    STI_Dharwad = (Mean_Yield_S1 * Mean_Yield_S2) / (mean_Yp_dharwad^2),
    SSI_Dharwad = (1 - (Mean_Yield_S2 / Mean_Yield_S1)) / D_dharwad,
    TSI_Dharwad = (Mean_PT_S2 / Mean_PT_S1) * (Mean_PT_S2 / pop_PT_Dharwad_S2),
    TCESI_Dharwad = (Mean_TCE_S2 / Mean_TCE_S1) * (Mean_TCE_S2 / pop_TCE_Dharwad_S2),
    Category_Dharwad = case_when(
      Mean_TCE_S1 >= pop_TCE_Dharwad_S1 & Mean_TCE_S2 >= pop_TCE_Dharwad_S2 ~ "Stable High Efficiency",
      Mean_TCE_S1 >= pop_TCE_Dharwad_S1 & Mean_TCE_S2 < pop_TCE_Dharwad_S2 ~ "Contrasting Susceptible",
      Mean_TCE_S1 < pop_TCE_Dharwad_S1 & Mean_TCE_S2 >= pop_TCE_Dharwad_S2 ~ "Contrasting Tolerant",
      TRUE ~ "Stable Low Efficiency"
    )
  )

tsi_dharwad_plot <- ggplot(dharwad_index, aes(x = TSI_Dharwad, y = Actual_Stress_Yield)) +
  geom_point(aes(color = Category_Dharwad), alpha = 0.7, size = 2.5) +
  geom_smooth(method = "lm", color = "black", fill = "lightgrey") +
  stat_cor(aes(label = paste(after_stat(r.label), after_stat(p.label), sep = "~`,`~")), size = 5, fontface = "bold", label.x.npc = "left", label.y.npc = "top") +
  scale_color_npg() + pub_theme +
  labs(title = "TSI Predicting Actual Stress Yield (Dharwad)", x = "TSI (Capacity Stability)", y = "Actual Grain Yield (kg/plot)", color = "TCE Category")
ggsave("Selection_Indices_Exports/TSI_Direct_Yield_Validation_Dharwad.png", plot = tsi_dharwad_plot, width = 9, height = 6, dpi = 600)

tcesi_dharwad_plot <- ggplot(dharwad_index, aes(x = TCESI_Dharwad, y = Actual_Stress_Yield)) +
  geom_point(aes(color = Category_Dharwad), alpha = 0.7, size = 2.5) +
  geom_smooth(method = "lm", color = "black", fill = "lightgrey") +
  stat_cor(aes(label = paste(after_stat(r.label), after_stat(p.label), sep = "~`,`~")), size = 5, fontface = "bold", label.x.npc = "left", label.y.npc = "top") +
  scale_color_npg() + pub_theme +
  labs(title = "TCESI Predicting Actual Stress Yield (Dharwad)", x = "TCESI (Efficiency Stability)", y = "Actual Grain Yield (kg/plot)", color = "TCE Category")
ggsave("Selection_Indices_Exports/TCESI_Direct_Yield_Validation_Dharwad.png", plot = tcesi_dharwad_plot, width = 9, height = 6, dpi = 600)
write.csv(dharwad_index, "Selection_Indices_Exports/Selection_Matrix_Dharwad.csv", row.names = FALSE)

# ---------------------------------------------------------
# 5C. GLOBAL CROSS-ZONAL ANALYSIS (COMBINED REGIMES)
# ---------------------------------------------------------
combined_summary <- efficiency_base %>%
  group_by(Genotype, Stress_Regime) %>%
  summarise(Mean_PT = mean(PT_Count, na.rm = TRUE), Mean_Yield = mean(GrainYld, na.rm = TRUE), Mean_TCE = mean(TCE, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = Stress_Regime, values_from = c(Mean_PT, Mean_Yield, Mean_TCE))

pop_PT_Stress_Combined <- mean(combined_summary$Mean_PT_Stress_Combined, na.rm = TRUE)
pop_TCE_Opt_Combined <- mean(combined_summary$Mean_TCE_Optimal_Combined, na.rm = TRUE)
pop_TCE_Stress_Combined <- mean(combined_summary$Mean_TCE_Stress_Combined, na.rm = TRUE)
mean_Yp_global <- mean(combined_summary$Mean_Yield_Optimal_Combined, na.rm = TRUE)
mean_Ys_global <- mean(combined_summary$Mean_Yield_Stress_Combined, na.rm = TRUE)
D_global <- 1 - (mean_Ys_global / mean_Yp_global)

global_index <- combined_summary %>%
  mutate(
    Actual_Stress_Yield_Global = Mean_Yield_Stress_Combined,
    STI_Global = (Mean_Yield_Optimal_Combined * Mean_Yield_Stress_Combined) / (mean_Yp_global^2),
    SSI_Global = (1 - (Mean_Yield_Stress_Combined / Mean_Yield_Optimal_Combined)) / D_global,
    TSI_Global = (Mean_PT_Stress_Combined / Mean_PT_Optimal_Combined) * (Mean_PT_Stress_Combined / pop_PT_Stress_Combined),
    TCESI_Global = (Mean_TCE_Stress_Combined / Mean_TCE_Optimal_Combined) * (Mean_TCE_Stress_Combined / pop_TCE_Stress_Combined),
    Category_Global = case_when(
      Mean_TCE_Optimal_Combined >= pop_TCE_Opt_Combined & Mean_TCE_Stress_Combined >= pop_TCE_Stress_Combined ~ "Stable High Efficiency",
      Mean_TCE_Optimal_Combined >= pop_TCE_Opt_Combined & Mean_TCE_Stress_Combined < pop_TCE_Stress_Combined ~ "Contrasting Susceptible",
      Mean_TCE_Optimal_Combined < pop_TCE_Opt_Combined & Mean_TCE_Stress_Combined >= pop_TCE_Stress_Combined ~ "Contrasting Tolerant",
      TRUE ~ "Stable Low Efficiency"
    )
  )

tsi_global_plot <- ggplot(global_index, aes(x = TSI_Global, y = Actual_Stress_Yield_Global)) +
  geom_point(aes(color = Category_Global), alpha = 0.7, size = 2.5) +
  geom_smooth(method = "lm", color = "black", fill = "lightgrey") +
  stat_cor(aes(label = paste(after_stat(r.label), after_stat(p.label), sep = "~`,`~")), size = 5, fontface = "bold", label.x.npc = "left", label.y.npc = "top") +
  scale_color_npg() + pub_theme +
  labs(title = "Global Validation: TSI Predicting Actual Stress Yield", x = "Global TSI", y = "Actual Grain Yield (kg/plot)", color = "TCE Category")
ggsave("Selection_Indices_Exports/TSI_Direct_Yield_Validation_Global.png", plot = tsi_global_plot, width = 9, height = 6, dpi = 600)

tcesi_global_plot <- ggplot(global_index, aes(x = TCESI_Global, y = Actual_Stress_Yield_Global)) +
  geom_point(aes(color = Category_Global), alpha = 0.7, size = 2.5) +
  geom_smooth(method = "lm", color = "black", fill = "lightgrey") +
  stat_cor(aes(label = paste(after_stat(r.label), after_stat(p.label), sep = "~`,`~")), size = 5, fontface = "bold", label.x.npc = "left", label.y.npc = "top") +
  scale_color_npg() + pub_theme +
  labs(title = "Global Validation: TCESI Predicting Actual Stress Yield", x = "Global TCESI", y = "Actual Grain Yield (kg/plot)", color = "TCE Category")
ggsave("Selection_Indices_Exports/TCESI_Direct_Yield_Validation_Global.png", plot = tcesi_global_plot, width = 9, height = 6, dpi = 600)
write.csv(global_index, "Selection_Indices_Exports/Selection_Matrix_Combined_Global.csv", row.names = FALSE)

# ---------------------------------------------------------
# 5D. TCE vs YIELD ACROSS ALL ENVIRONMENTS
# ---------------------------------------------------------
tce_yield_list <- list()
for (env in levels(master_gxe$Environment)) {
  p <- ggplot(master_gxe %>% filter(Environment == env), aes(x = TCE, y = GrainYld)) +
    geom_point(alpha = 0.6, size = 1.5, color = "#4575b4") +
    geom_smooth(method = "lm", color = "black", linewidth = 1) +
    stat_cor(aes(label = paste(after_stat(r.label), after_stat(p.label), sep = "~`,`~")), label.x.npc = "left", label.y.npc = "top", size = 4) +
    pub_theme +
    labs(title = env, x = "Tiller Conversion Efficiency (TCE)", y = "Grain Yield (kg/plot)")
  tce_yield_list[[env]] <- p
}
tce_yield_grid <- wrap_plots(tce_yield_list, ncol = 2) + plot_annotation(title = "Tiller Conversion Efficiency vs. Grain Yield, All Environments")
ggsave("Selection_Indices_Exports/TCE_vs_Yield_All_Environments.png", plot = tce_yield_grid, width = 12, height = 10, dpi = 600)

# ---------------------------------------------------------
# 5E. THE BAKE-OFF: TSI vs TCESI, WHICH ACTUALLY PREDICTS YIELD BETTER?
# ---------------------------------------------------------
hotelling_dependent_corr_test <- function(r_jk, r_jh, r_kh, n) {
  R_det <- 1 + 2 * r_jk * r_jh * r_kh - r_jk^2 - r_jh^2 - r_kh^2
  t_stat <- (r_jk - r_jh) * sqrt((n - 3) * (1 + r_kh)) / sqrt(2 * R_det)
  df <- n - 3
  p_value <- 2 * pt(-abs(t_stat), df)
  list(t = t_stat, df = df, p_value = p_value)
}

run_index_bakeoff <- function(data, tsi_col, tcesi_col, yield_col, label) {
  d <- data %>% select(all_of(c(tsi_col, tcesi_col, yield_col))) %>% drop_na()
  r_jk <- cor(d[[tsi_col]], d[[yield_col]])     
  r_jh <- cor(d[[tcesi_col]], d[[yield_col]])   
  r_kh <- cor(d[[tsi_col]], d[[tcesi_col]])     
  n <- nrow(d)
  
  test_result <- hotelling_dependent_corr_test(r_jk, r_jh, r_kh, n)
  data.frame(
    Comparison = label, N = n,
    r_TSI_Yield = round(r_jk, 3), r_TCESI_Yield = round(r_jh, 3), r_TSI_TCESI = round(r_kh, 3),
    Hotelling_t = round(test_result$t, 3), df = test_result$df, p_value = round(test_result$p_value, 4),
    Winner = case_when(
      test_result$p_value >= 0.05 ~ "No significant difference",
      abs(r_jk) > abs(r_jh) ~ "TSI (significantly better)",
      TRUE ~ "TCESI (significantly better)"
    )
  )
}

bakeoff_results <- bind_rows(
  run_index_bakeoff(delhi_index, "TSI_Delhi", "TCESI_Delhi", "Actual_Stress_Yield", "Delhi"),
  run_index_bakeoff(dharwad_index, "TSI_Dharwad", "TCESI_Dharwad", "Actual_Stress_Yield", "Dharwad"),
  run_index_bakeoff(global_index, "TSI_Global", "TCESI_Global", "Actual_Stress_Yield_Global", "Global (Combined)")
)
write.csv(bakeoff_results, "Selection_Indices_Exports/TSI_vs_TCESI_Bakeoff_HotellingsT.csv", row.names = FALSE)

bakeoff_plot_data <- bakeoff_results %>%
  select(Comparison, r_TSI_Yield, r_TCESI_Yield) %>%
  pivot_longer(cols = c(r_TSI_Yield, r_TCESI_Yield), names_to = "Index", values_to = "r") %>%
  mutate(Index = ifelse(Index == "r_TSI_Yield", "TSI", "TCESI"))

bakeoff_plot <- ggplot(bakeoff_plot_data, aes(x = Comparison, y = r, fill = Index)) +
  geom_col(position = "dodge", color = "black", alpha = 0.85) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("TSI" = "#4575b4", "TCESI" = "#d73027")) +
  pub_theme +
  labs(title = "TSI vs. TCESI: Correlation with Actual Stress Yield",
       subtitle = "See TSI_vs_TCESI_Bakeoff_HotellingsT.csv for the significance test behind each gap",
       x = "", y = "Pearson r vs. Actual Stress Yield", fill = "Index")
ggsave("Selection_Indices_Exports/TSI_vs_TCESI_Bakeoff_Plot.png", plot = bakeoff_plot, width = 8, height = 6, dpi = 300)

# ---------------------------------------------------------
# 5F. INDEX ROBUSTNESS CHECKS 
# ---------------------------------------------------------
full_indices <- global_index %>%
  mutate(
    MP  = (Mean_Yield_Optimal_Combined + Mean_Yield_Stress_Combined) / 2,
    GMP = sqrt(Mean_Yield_Optimal_Combined * Mean_Yield_Stress_Combined),
    YI  = Mean_Yield_Stress_Combined / mean_Ys_global,
    YSI = Mean_Yield_Stress_Combined / Mean_Yield_Optimal_Combined
  )
index_matrix <- full_indices %>% select(TSI_Global, TCESI_Global, STI_Global, SSI_Global, MP, GMP, YI, YSI)
write.csv(cor(index_matrix, use = "pairwise.complete.obs"), "TSI_Index_Correlation_Matrix.csv")

index_pairs_plot <- ggpairs(index_matrix, progress = FALSE) + pub_theme
ggsave("TSI_Index_Correlation_Pairs.png", plot = index_pairs_plot, width = 11, height = 11, dpi = 300)

denominator_check <- full_indices %>%
  mutate(Pctile_PT_Optimal = percent_rank(Mean_PT_Optimal_Combined)) %>%
  filter(Pctile_PT_Optimal < 0.05) %>%
  select(Genotype, Mean_PT_Optimal_Combined, Mean_PT_Stress_Combined, TSI_Global) %>%
  arrange(desc(TSI_Global))
write.csv(denominator_check, "TSI_Denominator_Risk_Genotypes.csv", row.names = FALSE)

denominator_check_tcesi <- full_indices %>%
  mutate(Pctile_TCE_Optimal = percent_rank(Mean_TCE_Optimal_Combined)) %>%
  filter(Pctile_TCE_Optimal < 0.05) %>%
  select(Genotype, Mean_TCE_Optimal_Combined, Mean_TCE_Stress_Combined, TCESI_Global) %>%
  arrange(desc(TCESI_Global))
write.csv(denominator_check_tcesi, "TCESI_Denominator_Risk_Genotypes.csv", row.names = FALSE)

genotypes_all <- unique(combined_summary$Genotype)
jackknife_topN_tsi <- function(leave_out) {
  d <- combined_summary %>% filter(Genotype != leave_out)
  pop_pt <- mean(d$Mean_PT_Stress_Combined, na.rm = TRUE)
  d %>% mutate(TSI_jk = (Mean_PT_Stress_Combined / Mean_PT_Optimal_Combined) * (Mean_PT_Stress_Combined / pop_pt)) %>%
    slice_max(TSI_jk, n = 10) %>% pull(Genotype)
}
jackknife_topN_tcesi <- function(leave_out) {
  d <- combined_summary %>% filter(Genotype != leave_out)
  pop_tce <- mean(d$Mean_TCE_Stress_Combined, na.rm = TRUE)
  d %>% mutate(TCESI_jk = (Mean_TCE_Stress_Combined / Mean_TCE_Optimal_Combined) * (Mean_TCE_Stress_Combined / pop_tce)) %>%
    slice_max(TCESI_jk, n = 10) %>% pull(Genotype)
}

original_top10_tsi <- global_index %>% slice_max(TSI_Global, n = 10) %>% pull(Genotype)
original_top10_tcesi <- global_index %>% slice_max(TCESI_Global, n = 10) %>% pull(Genotype)

set.seed(42)
jk_sample <- sample(genotypes_all, min(30, length(genotypes_all)))
jk_overlap_tsi <- map_dbl(map(jk_sample, jackknife_topN_tsi), ~length(intersect(.x, original_top10_tsi)) / 10)
jk_overlap_tcesi <- map_dbl(map(jk_sample, jackknife_topN_tcesi), ~length(intersect(.x, original_top10_tcesi)) / 10)

jk_summary <- data.frame(
  Metric = c("TSI_Top10_Jackknife_Stability", "TCESI_Top10_Jackknife_Stability"),
  Mean_Overlap_Fraction = round(c(mean(jk_overlap_tsi), mean(jk_overlap_tcesi)), 3),
  N_Genotypes_Jackknifed = length(jk_sample),
  Note = "Closer to 1.0 = more stable top-10 under leave-one-out resampling"
)
write.csv(jk_summary, "TSI_Jackknife_Stability.csv", row.names = FALSE)


# ==========================================
# 6. GLOBAL GxE VARIANCE PARTITIONING (ANOVA) ####
# ==========================================
variance_results <- data.frame()
global_anova_full <- data.frame()
core_traits <- c("PT_Count", "Tmax_Count", "GrainYld", "TKW", "AUTPC", "Relative_AUTPC")

for (trait in core_traits) {
  formula_anova <- paste(trait, "~ Genotype + Environment")
  model_anova <- tryCatch(lm(as.formula(formula_anova), data = master_gxe), error = function(e) return(NULL))
  if (!is.null(model_anova)) {
    anova_table <- car::Anova(model_anova, type = 2)
    global_anova_full <- rbind(global_anova_full, tidy(anova_table) %>% mutate(Trait = trait))
    
    ss_genotype <- anova_table["Genotype", "Sum Sq"]
    ss_environment <- anova_table["Environment", "Sum Sq"]
    ss_residual <- anova_table["Residuals", "Sum Sq"]
    ss_total <- sum(ss_genotype, ss_environment, ss_residual)
    
    variance_results <- rbind(variance_results, data.frame(
      Trait = trait,
      Pct_Genotype = (ss_genotype / ss_total) * 100,
      Pct_Environment = (ss_environment / ss_total) * 100,
      Pct_GxE_and_Residual = (ss_residual / ss_total) * 100
    ))
  }
}
write.csv(global_anova_full, "Global_Full_ANOVA_Table.csv", row.names = FALSE)

# DYNAMIC HERITABILITY: relative ranking within the panel (median split).
variance_results <- variance_results %>%
  mutate(Breeder_Target = ifelse(Pct_Genotype >= quantile(Pct_Genotype, 0.5, na.rm = TRUE),
                                 "Primary Heritable Target", "Environmentally Driven"))
write.csv(variance_results, "Global_GxE_Variance_Partitioning.csv", row.names = FALSE)

# v5 FIX: a pure median split always calls exactly half your traits "Primary
# Heritable Target" even if none of them are meaningfully heritable in
# absolute terms - this check makes that visible instead of letting it pass
# silently.
if (max(variance_results$Pct_Genotype, na.rm = TRUE) < 20) {
  warning("Even the single best trait has only ", round(max(variance_results$Pct_Genotype), 1),
          "% genotypic variance - 'Primary Heritable Target' below reflects RELATIVE ranking",
          " within a genuinely weak-heritability panel, not strong absolute heritability.",
          " Interpret Elite_Breeding_Targets_Heritability_Based.csv accordingly.")
}

elite_genotypes <- global_index %>%
  filter(Category_Global == "Stable High Efficiency", Actual_Stress_Yield_Global > mean_Ys_global) %>%
  arrange(desc(Actual_Stress_Yield_Global)) %>%
  select(Genotype, Actual_Stress_Yield_Global, TSI_Global, TCESI_Global, STI_Global, Category_Global)
write.csv(elite_genotypes, "Elite_Genotypes_Master_List.csv", row.names = FALSE)

# v5 FIX: elite_traits now reads the Breeder_Target classification computed
# just above, instead of an unconditional top-3 that ignored it entirely.
elite_traits <- variance_results %>% filter(Breeder_Target == "Primary Heritable Target") %>% pull(Trait)

if (length(elite_traits) == 0) {
  warning("No trait was classified as 'Primary Heritable Target' - Elite_Breeding_Targets_Heritability_Based.csv was not generated.")
} else {
  genotype_means <- master_gxe %>% group_by(Genotype) %>%
    summarise(across(all_of(elite_traits), \(x) mean(x, na.rm = TRUE)), .groups = "drop")
  z_scored <- genotype_means %>% mutate(across(all_of(elite_traits), \(x) as.numeric(scale(x)), .names = "Z_{.col}"))
  z_cols <- paste0("Z_", elite_traits)
  z_scored$Composite_Elite_Score <- rowMeans(z_scored[, z_cols, drop = FALSE], na.rm = TRUE)
  
  elite_breeding_targets <- z_scored %>% arrange(desc(Composite_Elite_Score)) %>%
    mutate(Rank = row_number()) %>% relocate(Rank, Genotype, Composite_Elite_Score)
  write.csv(elite_breeding_targets, "Elite_Breeding_Targets_Heritability_Based.csv", row.names = FALSE)
  
  top_heritable <- elite_breeding_targets %>% slice_max(Composite_Elite_Score, n = round(nrow(elite_breeding_targets) * 0.2), with_ties = FALSE)
  cross_validated_elite <- elite_genotypes %>% filter(Genotype %in% top_heritable$Genotype)
  write.csv(cross_validated_elite, "Elite_Genotypes_CrossValidated.csv", row.names = FALSE)
}

var_plot_data <- variance_results %>%
  pivot_longer(cols = starts_with("Pct_"), names_to = "Variance_Component", values_to = "Percentage") %>%
  mutate(Variance_Component = factor(Variance_Component, levels = c("Pct_GxE_and_Residual", "Pct_Genotype", "Pct_Environment")))

var_plot <- ggplot(var_plot_data, aes(x = reorder(Trait, -Percentage), y = Percentage, fill = Variance_Component)) +
  geom_bar(stat = "identity", position = "stack", color = "black") +
  scale_fill_manual(values = c("Pct_Environment" = "#d73027", "Pct_Genotype" = "#4575b4", "Pct_GxE_and_Residual" = "#fdae61"),
                    labels = c("G×E + Residual", "Genotype", "Environment")) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Genotypic Repeatability & Variance Partitioning", x = "Phenotypic Trait", y = "Proportion of Variance (%)", fill = "Component")
ggsave("Genotypic_Repeatability_Partitioning.png", plot = var_plot, width = 12, height = 6, dpi = 600)


# ==========================================
# 7. SINGLE LOCATION ANALYSIS (LOCALIZED GxE & ANOVAs) ####
# ==========================================
run_local_analysis <- function(loc_name) {
  loc_data <- subset(master_gxe, Location == loc_name)
  local_variance <- data.frame()
  anova_exports <- data.frame()
  dir.create(paste0("Local_Analysis_", loc_name), showWarnings = FALSE)
  
  for (trait in core_traits) {
    formula_local <- paste(trait, "~ Genotype + Sowing")
    model_local <- tryCatch(lm(as.formula(formula_local), data = loc_data), error = function(e) return(NULL))
    if (!is.null(model_local)) {
      anova_table <- car::Anova(model_local, type = 2)
      anova_exports <- rbind(anova_exports, tidy(anova_table) %>% mutate(Trait = trait, Location = loc_name))
      
      ss_genotype <- anova_table["Genotype", "Sum Sq"]
      ss_sowing <- anova_table["Sowing", "Sum Sq"]
      ss_residual <- anova_table["Residuals", "Sum Sq"]
      ss_total <- sum(ss_genotype, ss_sowing, ss_residual)
      
      local_variance <- rbind(local_variance, data.frame(
        Location = loc_name, Trait = trait,
        Pct_Genotype = (ss_genotype / ss_total) * 100, Pct_Sowing = (ss_sowing / ss_total) * 100, Pct_GxS_and_Residual = (ss_residual / ss_total) * 100
      ))
      
      p_local <- ggplot(loc_data, aes(x = Sowing, y = .data[[trait]], fill = Sowing)) +
        geom_boxplot(alpha = 0.8, outlier.shape = 21) +
        scale_fill_npg() + pub_theme +
        labs(title = paste(loc_name, "Local Variance:", trait), x = "Sowing Regime", y = trait) +
        stat_compare_means(method = "t.test", label = "p.signif")
      ggsave(paste0("Local_Analysis_", loc_name, "/", trait, "_Sowing_Contrast.png"), plot = p_local, width = 6, height = 5, dpi = 300)
    }
  }
  write.csv(anova_exports, paste0("Local_Analysis_", loc_name, "/ANOVA_Tables_Complete.csv"), row.names = FALSE)
  write.csv(local_variance, paste0("Local_Analysis_", loc_name, "/Variance_Partitioning.csv"), row.names = FALSE)
  
  var_plot_data_local <- local_variance %>%
    pivot_longer(cols = starts_with("Pct_"), names_to = "Component", values_to = "Percentage") %>%
    mutate(Component = factor(Component, levels = c("Pct_GxS_and_Residual", "Pct_Genotype", "Pct_Sowing")))
  
  p_var <- ggplot(var_plot_data_local, aes(x = reorder(Trait, -Percentage), y = Percentage, fill = Component)) +
    geom_bar(stat = "identity", position = "stack", color = "black") +
    scale_fill_manual(values = c("Pct_Sowing" = "#d73027", "Pct_Genotype" = "#4575b4", "Pct_GxS_and_Residual" = "#fdae61"),
                      labels = c("G×S + Residual", "Genotype", "Sowing Regime")) +
    theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = paste("Local Genotypic Repeatability:", loc_name),
         subtitle = "Evaluating structural stability exclusively across local sowing treatments",
         x = "Phenotypic Trait", y = "Proportion of Variance (%)", fill = "Component")
  ggsave(paste0("Local_Analysis_", loc_name, "/_Local_Variance_Barplot.png"), plot = p_var, width = 10, height = 6, dpi = 600)
}

run_local_analysis("Delhi")
run_local_analysis("Dharwad")


# ==========================================
# 8. GGE MEAN VS STABILITY ANALYSIS ####
# ==========================================
gge_yield <- gge(master_gxe, env = Environment, gen = Genotype, resp = GrainYld)
png("GGE_Mean_vs_Stability_Yield.png", width = 3000, height = 2400, res = 300)
plot(gge_yield, type = 2, title = "GGE Biplot: Mean Yield vs. Environmental Stability")
dev.off()

gge_pt <- gge(master_gxe, env = Environment, gen = Genotype, resp = PT_Count)
png("GGE_Mean_vs_Stability_PT.png", width = 3000, height = 2400, res = 300)
plot(gge_pt, type = 2, title = "GGE Biplot: Mean Productive Tillers vs. Environmental Stability")
dev.off()


# ==========================================
# 9. METEOROLOGICAL ANALYSIS: THERMAL VELOCITY & GDD ####
# ==========================================
weather_data <- read_excel("D:/PhD School/Research/Weather_Master.xlsx", sheet = "Sheet1")
weather_data$Date <- as.Date(weather_data$Date)
base_temp <- 5

weather_data <- weather_data %>%
  mutate(Tmean = (Tmax + Tmin) / 2, Daily_GDD = ifelse(Tmean > base_temp, Tmean - base_temp, 0))

phase_intervals <- env_metadata %>% select(Environment, Location, Sowing, Start_Date, End_Date)

thermal_summary <- data.frame()
for (i in 1:nrow(phase_intervals)) {
  env_weather <- weather_data %>% filter(Location == phase_intervals$Location[i], Date >= phase_intervals$Start_Date[i] & Date <= phase_intervals$End_Date[i])
  thermal_summary <- rbind(thermal_summary, data.frame(
    Environment = phase_intervals$Environment[i], Total_Days = nrow(env_weather),
    Cumulative_GDD = sum(env_weather$Daily_GDD, na.rm = TRUE),
    GDD_Velocity = sum(env_weather$Daily_GDD, na.rm = TRUE) / nrow(env_weather),
    Mean_Tmax = mean(env_weather$Tmax, na.rm = TRUE), Mean_Tmin = mean(env_weather$Tmin, na.rm = TRUE)
  ))
}
write.csv(thermal_summary, "Meteorological_Thermal_Summary.csv", row.names = FALSE)

weather_plot <- ggplot() +
  geom_rect(data = phase_intervals, aes(xmin = Start_Date, xmax = End_Date, ymin = -Inf, ymax = Inf, fill = Sowing), alpha = 0.15) +
  geom_line(data = weather_data, aes(x = Date, y = Tmax, color = "Tmax"), linewidth = 1) +
  geom_line(data = weather_data, aes(x = Date, y = Tmin, color = "Tmin"), linewidth = 1) +
  geom_hline(yintercept = base_temp, linetype = "dashed", color = "black") +
  facet_wrap(~ Location, scales = "free_x") +
  scale_color_manual(values = c("Tmax" = "#d73027", "Tmin" = "#4575b4")) +
  scale_fill_npg() + pub_theme +
  labs(title = "Meteorological Trajectories with Crop Windows",
       subtitle = "Shaded regions represent in-field durations (Emergence to PI)",
       x = "Date", y = "Temperature (°C)", color = "Parameter", fill = "Sowing Regime")
ggsave("Meteorological_Thermal_Trends.png", plot = weather_plot, width = 12, height = 5, dpi = 600)

gdd_plot <- ggplot(thermal_summary, aes(x = Environment, y = Cumulative_GDD, fill = Environment)) +
  geom_col(color = "black", alpha = 0.8) +
  geom_text(aes(label = round(Cumulative_GDD, 0)), vjust = -0.5, fontface = "bold") +
  scale_fill_npg() + pub_theme + theme(legend.position = "none") +
  labs(title = "Cumulative Growing Degree Days (GDD)", x = "Testing Regime", y = "Total GDD (°C·d)")
ggsave("Meteorological_Cumulative_GDD.png", plot = gdd_plot, width = 8, height = 5, dpi = 600)


# ==========================================
# 10. PHENOLOGICAL CLUSTER ANALYSIS (K-MEANS & PCA) ####
# ==========================================
cluster_vars <- c("Rate_E_T1", "Rate_T1_T2", "Rate_T2_T3", "Rate_T3_Tmax")

pheno_df <- master_gxe %>%
  filter(Environment == "Delhi_S2") %>%
  select(Genotype, all_of(cluster_vars)) %>%
  drop_na() %>%
  column_to_rownames("Genotype")

pheno_scaled <- scale(pheno_df)

set.seed(42)
kmeans_result <- kmeans(pheno_scaled, centers = 3, nstart = 25)
pheno_df$Raw_Cluster <- as.factor(kmeans_result$cluster)

cluster_profiles <- pheno_df %>% group_by(Raw_Cluster) %>% summarise(Speed = mean(Rate_E_T1)) %>% arrange(desc(Speed))
cluster_map <- data.frame(Raw_Cluster = cluster_profiles$Raw_Cluster,
                          Cluster = c("Early/Fast Establishers", "Intermediate Establishers", "Late/Slow Establishers"))

pheno_df <- pheno_df %>%
  rownames_to_column("Genotype") %>%
  left_join(cluster_map, by = "Raw_Cluster") %>%
  select(-Raw_Cluster) %>%
  mutate(Cluster = as.factor(Cluster)) %>%
  column_to_rownames("Genotype")

pca_result <- prcomp(pheno_scaled, center = FALSE, scale. = FALSE)
pheno_df$PC1 <- pca_result$x[, 1]
pheno_df$PC2 <- pca_result$x[, 2]

cluster_export <- pheno_df %>% rownames_to_column("Genotype")
write.csv(cluster_export, "Phenological_Cluster_Assignments.csv", row.names = FALSE)

var_explained <- round(100 * pca_result$sdev^2 / sum(pca_result$sdev^2), 1)

pheno_df_labeled <- pheno_df %>% rownames_to_column("Genotype")
centroids <- pheno_df_labeled %>% group_by(Cluster) %>% summarise(PC1_c = mean(PC1), PC2_c = mean(PC2))
pheno_df_labeled <- pheno_df_labeled %>% left_join(centroids, by = "Cluster") %>%
  mutate(Dist = sqrt((PC1 - PC1_c)^2 + (PC2 - PC2_c)^2))
top_extremes <- pheno_df_labeled %>% group_by(Cluster) %>% slice_max(Dist, n = 5)

cluster_plot <- ggplot(pheno_df_labeled, aes(x = PC1, y = PC2, fill = Cluster)) +
  geom_point(shape = 21, size = 4, alpha = 0.8, color = "black") +
  stat_ellipse(aes(color = Cluster), type = "norm", linetype = "dashed", linewidth = 1) +
  geom_text_repel(data = top_extremes, aes(label = Genotype), size = 3.5, fontface = "bold", box.padding = 0.5, max.overlaps = Inf) +
  scale_fill_npg() + scale_color_npg() + pub_theme +
  labs(title = "K-Means Clustering of Genotypic Phenological Strategies (Delhi S2 Baseline)",
       subtitle = "Partitioning multidimensional velocities across constitutive development (Top Extremes Labeled)",
       x = paste0("PC1 (", var_explained[1], "%)"), y = paste0("PC2 (", var_explained[2], "%)"))
ggsave("Phenological_Clusters_PCA.png", plot = cluster_plot, width = 10, height = 8, dpi = 600)

cluster_mapping <- pheno_df %>% rownames_to_column("Genotype") %>% select(Genotype, Cluster)

if ("Cluster" %in% colnames(master_gxe)) master_gxe$Cluster <- NULL
master_gxe <- master_gxe %>% left_join(cluster_mapping, by = "Genotype")

trajectory_data <- master_gxe %>%
  filter(!is.na(Cluster)) %>%
  group_by(Cluster, Environment) %>%
  summarise(Mean_Yield = mean(GrainYld, na.rm = TRUE), SE_Yield = sd(GrainYld, na.rm = TRUE) / sqrt(n()), .groups = 'drop') %>%
  mutate(Environment = factor(Environment, levels = c("Delhi_S2", "Delhi_S1", "Dharwad_S1", "Dharwad_S2")))

trajectory_plot <- ggplot(trajectory_data, aes(x = Environment, y = Mean_Yield, color = Cluster, group = Cluster)) +
  geom_line(linewidth = 1.5, alpha = 0.8) +
  geom_point(size = 4, shape = 21, fill = "white", stroke = 1.5) +
  geom_errorbar(aes(ymin = Mean_Yield - SE_Yield, ymax = Mean_Yield + SE_Yield), width = 0.1, linewidth = 0.8) +
  scale_color_npg() + pub_theme +
  labs(title = "Reaction Norms: Phenological Clusters Across Zonal Stress Regimes",
       subtitle = "Tracking yield stability based on inherent (Delhi S2 Optimal) baseline strategies",
       x = "Environmental Regime", y = "Mean Grain Yield (kg/plot) ± SE", color = "Constitutive Baseline Strategy")
ggsave("Reaction_Norms_Cluster_Overlay.png", plot = trajectory_plot, width = 10, height = 6, dpi = 600)

# ---------------------------------------------------------
# 10C. WINDOW-SPECIFIC GDD [v5 RESTORED - Section 16 depends on window_gdd]
# ---------------------------------------------------------
compute_window_gdd <- function(env_row, weather) {
  loc <- env_row$Location; start <- env_row$Start_Date
  w <- weather %>% filter(Location == loc, Date >= start) %>% arrange(Date)
  bounds <- cumsum(c(env_row$Int_E_T1, env_row$Int_T1_T2, env_row$Int_T2_T3, env_row$Int_T3_Tmax))
  if (nrow(w) < bounds[4]) {
    warning(paste("Insufficient weather records for", env_row$Environment, "- window GDD may be incomplete."))
  }
  tibble(
    Environment = env_row$Environment,
    Window = c("E_T1", "T1_T2", "T2_T3", "T3_Tmax"),
    Window_GDD = c(
      sum(w$Daily_GDD[1:bounds[1]], na.rm = TRUE),
      sum(w$Daily_GDD[(bounds[1] + 1):bounds[2]], na.rm = TRUE),
      sum(w$Daily_GDD[(bounds[2] + 1):bounds[3]], na.rm = TRUE),
      sum(w$Daily_GDD[(bounds[3] + 1):bounds[4]], na.rm = TRUE)
    )
  )
}
window_gdd <- map_dfr(1:nrow(env_metadata), ~compute_window_gdd(env_metadata[.x, ], weather_data))
write.csv(window_gdd, "Window_Specific_GDD_By_Environment.csv", row.names = FALSE)

window_gdd_plot <- ggplot(window_gdd, aes(x = factor(Window, levels = c("E_T1", "T1_T2", "T2_T3", "T3_Tmax")), y = Window_GDD, fill = Environment)) +
  geom_col(position = "dodge", color = "black") +
  scale_fill_npg() + pub_theme +
  labs(title = "Thermal Load by Developmental Window", x = "Phenological Window", y = "GDD (°C·day)")
ggsave("Window_Specific_GDD_Plot.png", plot = window_gdd_plot, width = 9, height = 6, dpi = 300)

trajectory_thermal <- trajectory_data %>%
  left_join(thermal_summary %>% select(Environment, Cumulative_GDD, Mean_Tmax), by = "Environment")

thermal_reaction_norm <- ggplot(trajectory_thermal, aes(x = Cumulative_GDD, y = Mean_Yield, color = Cluster)) +
  geom_line(linewidth = 1.2, alpha = 0.8) +
  geom_point(size = 4, shape = 21, fill = "white", stroke = 1.5) +
  geom_errorbar(aes(ymin = Mean_Yield - SE_Yield, ymax = Mean_Yield + SE_Yield), width = 15) +
  geom_text_repel(aes(label = Environment), size = 3, color = "black", show.legend = FALSE) +
  scale_color_npg() + pub_theme +
  labs(title = "Thermal Reaction Norms: Cluster Yield vs. Cumulative Thermal Load",
       subtitle = "Continuous GDD axis (vs. categorical Environment) shows the actual dose-response",
       x = "Cumulative GDD (°C·day)", y = "Mean Grain Yield ± SE", color = "Baseline Strategy")
ggsave("Thermal_Reaction_Norm_By_Cluster.png", plot = thermal_reaction_norm, width = 9, height = 6, dpi = 300)


# ==========================================
# 11. SINK-CAPACITY RATIOS (DELHI BASELINE PCA & DHARWAD STRESS OVERLAY) ####
# ==========================================
ratio_df <- master_gxe %>%
  mutate(Ratio_TCE = PT_Count / Tmax_Count, Ratio_Total_Yield = Tmax_Count / GrainYld, Ratio_PT_Yield = PT_Count / GrainYld) %>%
  drop_na(Ratio_TCE, Ratio_Total_Yield, Ratio_PT_Yield)

ratio_baseline <- ratio_df %>% filter(Environment == "Delhi_S2") %>% column_to_rownames("Genotype") %>% select(starts_with("Ratio_"))
ratio_pca <- prcomp(ratio_baseline, center = TRUE, scale. = TRUE)

pca_baseline_coords <- as.data.frame(ratio_pca$x) %>% rownames_to_column("Genotype") %>% mutate(Regime = "Delhi S2 (Optimal Baseline)")
ratio_stress <- ratio_df %>% filter(Environment == "Dharwad_S2") %>% column_to_rownames("Genotype") %>% select(starts_with("Ratio_"))
pca_stress_coords <- as.data.frame(predict(ratio_pca, newdata = ratio_stress)) %>% rownames_to_column("Genotype") %>% mutate(Regime = "Dharwad S2 (Extreme Stress)")

ratio_plot_data <- bind_rows(pca_baseline_coords, pca_stress_coords) %>%
  left_join(pheno_df %>% rownames_to_column("Genotype") %>% select(Genotype, Cluster), by = "Genotype")

shift_data <- pca_baseline_coords %>% select(Genotype, PC1_base = PC1, PC2_base = PC2) %>%
  inner_join(pca_stress_coords %>% select(Genotype, PC1_stress = PC1, PC2_stress = PC2), by = "Genotype") %>%
  mutate(Shift_Dist = sqrt((PC1_stress - PC1_base)^2 + (PC2_stress - PC2_base)^2)) %>%
  arrange(desc(Shift_Dist)) %>% head(10)
write.csv(shift_data, "Sink_Capacity_Shift_Magnitude.csv", row.names = FALSE)

label_data_overlay <- ratio_plot_data %>% filter(Regime == "Dharwad S2 (Extreme Stress)", Genotype %in% shift_data$Genotype)

ratio_overlay_plot <- ggplot(ratio_plot_data, aes(x = PC1, y = PC2, color = Regime, shape = Regime)) +
  geom_point(size = 3, alpha = 0.7) +
  stat_ellipse(type = "norm", linetype = "dashed", linewidth = 1) +
  geom_text_repel(data = label_data_overlay, aes(label = Genotype), size = 3.5, fontface = "bold", box.padding = 0.5, show.legend = FALSE, max.overlaps = Inf) +
  scale_color_npg() + pub_theme +
  labs(title = "Sink-Capacity Efficiency Shift (Delhi Baseline vs Dharwad Overlay)",
       subtitle = "Projecting extreme terminal heat response into the optimal baseline PCA space (Largest Shifters Labeled)",
       x = paste0("PC1: Vegetative Burden (", round(summary(ratio_pca)$importance[2, 1] * 100, 1), "%)"),
       y = paste0("PC2: Conversion Efficiency (", round(summary(ratio_pca)$importance[2, 2] * 100, 1), "%)"))
ggsave("Ratio_Shift_Delhi_Dharwad_Overlay.png", plot = ratio_overlay_plot, width = 10, height = 7, dpi = 600)


# ==========================================
# 12. THERMAL ESCAPE & DURATION COMPENSATION (ANOVA & EMMEANS) ####
# ==========================================
duration_analysis <- master_gxe %>%
  filter(!Peak_Phase %in% c("Uniform", "Continuous Decline")) %>% 
  mutate(
    Duration_Class = case_when(
      Peak_Phase %in% c("Very Early (E-T1)", "Early Vigor (T1-T2)") ~ "Early (Short Duration)",
      Peak_Phase %in% c("Late Vigor (T2-T3)", "Terminal Vigor (T3-Tmax)") ~ "Late (Long Duration)"
    ),
    Environment = factor(Environment, levels = c("Delhi_S2", "Delhi_S1", "Dharwad_S1", "Dharwad_S2"))
  )

yield_compensation_plot <- ggplot(duration_analysis, aes(x = Environment, y = GrainYld, fill = Duration_Class)) +
  geom_boxplot(alpha = 0.85, outlier.shape = 21, color = "black") +
  scale_fill_npg() + pub_theme +
  labs(title = "Yield Compensation Across Temporal Stress Gradients",
       subtitle = "Full 4-environment comparison; see emmeans post-hoc contrasts below for all pairwise tests",
       x = "Testing Regime", y = "Grain Yield (kg/plot)", fill = "Genotypic Strategy") +
  stat_compare_means(aes(group = Duration_Class), label = "p.signif", method = "t.test")
ggsave("Yield_Compensation_Temporal_Stress.png", plot = yield_compensation_plot, width = 10, height = 7, dpi = 600)

abortion_penalty_plot <- ggplot(duration_analysis, aes(x = Duration_Class, y = Delta_T3_Tmax, fill = Duration_Class)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.2, fill = "white", outlier.shape = NA) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 1) +
  facet_wrap(~ Environment, nrow = 1) +
  scale_fill_npg() + pub_theme +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  labs(title = "Late-Stage Tiller Abortion Dynamics",
       subtitle = "Negative Delta indicates active sink abortion driven by regime-specific thermal loads",
       x = "", y = "Tiller Change (Delta T3 to Tmax)")
ggsave("Late_Stage_Abortion_Dynamics.png", plot = abortion_penalty_plot, width = 12, height = 6, dpi = 600)

compensation_model <- lm(GrainYld ~ Duration_Class * Environment, data = duration_analysis)

sink("Thermal_Escape_ANOVA_Results.txt")
cat("--- Thermal Escape Omnibus ANOVA (Type II), all 4 environments ---\n")
print(car::Anova(compensation_model, type = 2))

cat("\n\n--- Post-Hoc Contrasts (Tukey Adjusted) ---\n")
cat("1. Comparing Environments within each Duration Strategy:\n")
env_contrasts <- emmeans(compensation_model, ~ Environment | Duration_Class)
print(contrast(env_contrasts, "pairwise"))

cat("\n2. Comparing Duration Strategies within each Environment:\n")
dur_contrasts <- emmeans(compensation_model, ~ Duration_Class | Environment)
print(contrast(dur_contrasts, "pairwise"))
sink()

# ==========================================
# 13. STAGE-WISE PREDICTIVE POWER OF YIELD ####
# ==========================================
stage_traits <- c("T1_Count", "T2_Count", "T3_Count", "Tmax_Count", "PT_Count")

stage_predictive_power <- map_dfr(levels(master_gxe$Environment), function(env) {
  d <- master_gxe %>% filter(Environment == env)
  map_dfr(stage_traits, function(st) {
    m <- lm(as.formula(paste("GrainYld ~", st)), data = d)
    tibble(Environment = env, Stage = st, R2 = summary(m)$r.squared,
           Slope = coef(m)[2], p_value = summary(m)$coefficients[2, 4])
  })
})
write.csv(stage_predictive_power, "Stage_Predictive_Power_By_Environment.csv", row.names = FALSE)

stage_power_plot <- ggplot(stage_predictive_power, aes(x = factor(Stage, levels = stage_traits), y = R2, fill = Environment)) +
  geom_col(position = "dodge", color = "black") +
  scale_fill_npg() + pub_theme +
  labs(title = "Predictive Power of Each Developmental Stage for Grain Yield",
       subtitle = "Highest bar per environment = the stage most worth measuring in that thermal regime",
       x = "Developmental Stage", y = expression(R^2 ~ "(single-stage linear model)"))
ggsave("Stage_Predictive_Power_Plot.png", plot = stage_power_plot, width = 10, height = 6, dpi = 300)

# ==========================================
# 13B. INCREMENTAL PREDICTIVE VALUE OF EACH STAGE (ADJUSTED R-SQUARED) ####
# ==========================================
incremental_stages <- c("T1_Count", "T2_Count", "T3_Count", "Tmax_Count")

incremental_r2 <- map_dfr(levels(master_gxe$Environment), function(env) {
  d <- master_gxe %>% filter(Environment == env)
  formula_strs <- accumulate(incremental_stages, ~paste(.x, "+", .y))
  formula_strs <- paste("GrainYld ~", formula_strs)
  
  r2_vals <- map_dbl(formula_strs, ~summary(lm(as.formula(.x), data = d))$adj.r.squared)
  
  tibble(
    Environment = env,
    Stage_Added = incremental_stages,
    Cumulative_Adj_R2 = r2_vals,
    Incremental_Adj_R2 = c(r2_vals[1], diff(r2_vals))
  )
})
write.csv(incremental_r2, "Incremental_Stage_R2_By_Environment.csv", row.names = FALSE)

incremental_plot <- ggplot(incremental_r2, aes(x = factor(Stage_Added, levels = incremental_stages), y = Incremental_Adj_R2, fill = Environment)) +
  geom_col(position = "dodge", color = "black") +
  scale_fill_npg() + pub_theme +
  labs(title = "Unique Predictive Contribution of Each Stage (Beyond Earlier Stages)",
       subtitle = "How much does each successive stage add to Adjusted R^2, once earlier stages are already in the model?",
       x = "Stage Added to the Model", y = expression(Delta ~ R[adj]^2))
ggsave("Incremental_Stage_R2_Plot.png", plot = incremental_plot, width = 10, height = 6, dpi = 300)


# ==========================================
# 14. MASTER DATA EXPORT & COMPLETION ####
# ==========================================
write.csv(master_gxe, "Master_GxE_Matrix_Complete.csv", row.names = FALSE)
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
print("--- Master Pipeline Execution Complete ---")

# ==========================================
# 15. STRUCTURAL EQUATION MODELING (SEM) / PATH ANALYSIS ####
# ==========================================
required_sem_pkgs <- c("lavaan", "semPlot")
missing_sem <- setdiff(required_sem_pkgs, rownames(installed.packages()))
if(length(missing_sem) > 0) install.packages(missing_sem)
library(lavaan)
library(semPlot)

print("--- Executing Structural Equation Modeling (Path Analysis) ---")

# ------------------------------------------------------------------
# Guard: recompute AUTPC_Asymmetry if missing
# ------------------------------------------------------------------
if (!"AUTPC_Asymmetry" %in% colnames(master_gxe)) {
  message("AUTPC_Asymmetry not found — recomputing.")
  master_gxe <- master_gxe %>%
    select(-any_of(c("Int_E_T1","Int_T1_T2","Int_T2_T3","Int_T3_Tmax","Total_T1_Tmax"))) %>%
    left_join(
      env_metadata %>% select(Environment, Int_E_T1, Int_T1_T2, Int_T2_T3, Int_T3_Tmax),
      by = "Environment"
    ) %>%
    rowwise() %>%
    mutate(
      AUTPC_Early     = 0.5*(0 + T1_Count)*Int_E_T1 +
                        0.5*(T1_Count + T2_Count)*Int_T1_T2,
      AUTPC_Late      = 0.5*(T2_Count + T3_Count)*Int_T2_T3 +
                        0.5*(T3_Count + Tmax_Count)*Int_T3_Tmax,
      AUTPC_Asymmetry = ifelse(is.na(AUTPC) | AUTPC == 0, NA_real_,
                               (AUTPC_Early - AUTPC_Late) / AUTPC)
    ) %>%
    ungroup()
}

master_gxe <- master_gxe %>%
  mutate(AUTPC_Asymmetry = ifelse(
    is.infinite(AUTPC_Asymmetry) | is.nan(AUTPC_Asymmetry),
    NA_real_, AUTPC_Asymmetry
  ))

# ------------------------------------------------------------------
# 15A. Extended Biophysical Path Model
# ------------------------------------------------------------------
path_model <- '
  # --- Sequential development chain ---
  T2_Count   ~ T1_Count
  Tmax_Count ~ T2_Count + T1_Count

  # --- Sink survival ---
  PT_Count   ~ Tmax_Count + T1_Count + T2_Count + AUTPC_Asymmetry

  # --- Final yield ---
  GrainYld   ~ PT_Count + TKW + T1_Count + AUTPC_Asymmetry

  # --- Methodological / Structural Covariances ---
  # We only covary AUTPC_Asymmetry with exogenous T1 and the peak Tmax.
  T1_Count   ~~ AUTPC_Asymmetry
  Tmax_Count ~~ AUTPC_Asymmetry
  
  # Contemporaneous and MI-recommended
  PT_Count   ~~ TKW
  T2_Count   ~~ PT_Count
  T1_Count   ~~ Tmax_Count
'

dir.create("SEM_Path_Analysis", showWarnings = FALSE)
sem_results_list     <- list()
sem_fit_indices_list <- list()

env_list_sem <- levels(master_gxe$Environment)

# ------------------------------------------------------------------
# 15B. Fit across all four environments
# ------------------------------------------------------------------
for (env in env_list_sem) {

  env_data <- master_gxe %>%
    filter(Environment == env) %>%
    select(T1_Count, T2_Count, T3_Count, Tmax_Count,
           PT_Count, AUTPC_Asymmetry, GrainYld, TKW) %>%
    drop_na()

  if (nrow(env_data) < 10) {
    warning(paste("SEM skipped for", env, "— fewer than 10 complete cases."))
    next
  }

  fit <- tryCatch(
    suppressWarnings(sem(path_model, data = env_data, estimator = "MLR")),
    error = function(e) {
      warning(paste("SEM hard error for", env, ":", conditionMessage(e)))
      NULL
    }
  )

  if (is.null(fit)) next

  converged_ok    <- isTRUE(lavInspect(fit, "converged"))
  post_raw        <- tryCatch(lavInspect(fit, "post.check"), error = function(e) NA)
  proper_solution <- isTRUE(post_raw)

  if (!converged_ok || !proper_solution)
    warning(paste0("SEM [", env, "]: converged=", converged_ok,
                   " | proper_solution=", proper_solution))

  fi <- fitMeasures(fit, c("chisq","df","pvalue","cfi","tli","rmsea","srmr"))

  sem_fit_indices_list[[env]] <- data.frame(
    Environment     = env,
    chisq           = unname(fi["chisq"]),
    df              = unname(fi["df"]),
    pvalue          = unname(fi["pvalue"]),
    cfi             = unname(fi["cfi"]),
    tli             = unname(fi["tli"]),
    rmsea           = unname(fi["rmsea"]),
    srmr            = unname(fi["srmr"]),
    Converged       = converged_ok,
    Proper_Solution = proper_solution,
    Fit_Adequate    = unname(fi["cfi"]) >= 0.90 &
                      unname(fi["rmsea"]) < 0.08 &
                      unname(fi["srmr"]) < 0.08
  )

  mi <- tryCatch(modindices(fit, sort. = TRUE, maximum.number = 5),
                 error = function(e) NULL)
  if (!is.null(mi) && nrow(mi) > 0) {
    cat("\n--- Top modification indices for", env, "---\n")
    print(mi[, c("lhs","op","rhs","mi","epc")])
  }

  param_ests <- tryCatch(standardizedsolution(fit), error = function(e) NULL)
  if (!is.null(param_ests)) {
    sem_results_list[[env]] <- param_ests %>% filter(op %in% c("~", "~~")) %>% mutate(Environment = env)
  }

  tryCatch({
    png(paste0("SEM_Path_Analysis/Path_Diagram_", env, ".png"), width = 2400, height = 1800, res = 300)
    semPaths(fit, what = "std", layout = "tree2", edge.label.cex = 1)
    dev.off()
  }, error = function(e) {
    if (dev.cur() > 1) dev.off()
    warning(paste("semPaths diagram failed for", env, ":", conditionMessage(e)))
  })
}

# ------------------------------------------------------------------
# 15C. Export
# ------------------------------------------------------------------
if (length(sem_results_list) == 0) {
  warning("SEM: no environment produced usable results.")
  write.csv(data.frame(),
            "SEM_Path_Analysis/SEM_Comparative_Coefficients_All_Environments.csv",
            row.names = FALSE)
  write.csv(data.frame(),
            "SEM_Path_Analysis/SEM_Fit_Indices_By_Environment.csv",
            row.names = FALSE)
} else {
  sem_comparison <- bind_rows(sem_results_list) %>%
    select(Environment,
           LHS             = lhs,
           Operator        = op,
           RHS             = rhs,
           Std_Coefficient = est.std,
           SE              = se,
           z_value         = z,
           p_value         = pvalue) %>%
    mutate(Significance = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE            ~ "ns"
    )) %>%
    arrange(LHS, RHS, Environment)

  write.csv(sem_comparison,
            "SEM_Path_Analysis/SEM_Comparative_Coefficients_All_Environments.csv",
            row.names = FALSE)

  sem_fit_summary <- bind_rows(sem_fit_indices_list)
  write.csv(sem_fit_summary,
            "SEM_Path_Analysis/SEM_Fit_Indices_By_Environment.csv",
            row.names = FALSE)

  cat("\n========== SEM Fit Indices Summary ==========\n")
  print(sem_fit_summary)
  cat("=============================================\n")
  cat("Adequate fit: CFI >= 0.90, RMSEA < 0.08, SRMR < 0.08\n")

  n_adequate <- sum(sem_fit_summary$Fit_Adequate, na.rm = TRUE)
  message(n_adequate, " of ", nrow(sem_fit_summary),
          " environments show adequate model fit.")
  message("Section 15 complete — results in SEM_Path_Analysis/")
}

# ------------------------------------------------------------------
# 15D. Location-Specific SEM Models (Pooled working resource per location)
# ------------------------------------------------------------------
print("--- Executing Location-Specific Pooled SEM ---")
dir.create("SEM_Path_Analysis/Location_Pooled", showWarnings = FALSE)

location_fit_indices <- list()
location_results <- list()

for (loc in unique(master_gxe$Location)) {
  loc_data <- master_gxe %>% filter(Location == loc) %>%
    select(T1_Count, T2_Count, T3_Count, Tmax_Count,
           PT_Count, AUTPC_Asymmetry, GrainYld, TKW) %>%
    drop_na()

  if (nrow(loc_data) > 10) {
    loc_fit <- try(sem(path_model, data = loc_data, estimator = "MLR"), silent = TRUE)
    if (!inherits(loc_fit, "try-error")) {
      
      fi <- fitMeasures(loc_fit, c("chisq","df","pvalue","cfi","tli","rmsea","srmr"))
      location_fit_indices[[loc]] <- data.frame(
        Location = loc, chisq = unname(fi["chisq"]), df = unname(fi["df"]),
        pvalue = unname(fi["pvalue"]), cfi = unname(fi["cfi"]),
        tli = unname(fi["tli"]), rmsea = unname(fi["rmsea"]), srmr = unname(fi["srmr"])
      )
      
      param_ests <- tryCatch(standardizedsolution(loc_fit), error = function(e) NULL)
      if (!is.null(param_ests)) {
        location_results[[loc]] <- param_ests %>% filter(op %in% c("~", "~~")) %>% mutate(Location = loc)
      }
      
      tryCatch({
        png(paste0("SEM_Path_Analysis/Location_Pooled/Path_Diagram_", loc, "_Pooled.png"), width = 2400, height = 1800, res = 300)
        p <- semPaths(loc_fit, what = "std", layout = "tree2", edge.label.cex = 1)
        plot(p)
        dev.off()
      }, error = function(e) { if (dev.cur() > 1) dev.off() })
    }
  }
}

if(length(location_fit_indices) > 0) {
  write.csv(bind_rows(location_fit_indices), "SEM_Path_Analysis/Location_Pooled/Pooled_Fit_Indices.csv", row.names = FALSE)
  write.csv(bind_rows(location_results), "SEM_Path_Analysis/Location_Pooled/Pooled_Coefficients.csv", row.names = FALSE)
  message("Section 15D complete — Location-specific resources in Location_Pooled/")
}

# ==========================================
# 16. NON-LINEAR THERMAL MODELING (LOGISTIC GDD CURVES) ####
# ==========================================
print("--- Executing Non-Linear Thermal Modeling ---")
dir.create("Thermal_Logistic_Models", showWarnings = FALSE)

if (!exists("window_gdd")) {
  message("window_gdd not found — recomputing from env_metadata + weather_data.")
  if (!exists("weather_data")) {
    weather_data <- readxl::read_excel("D:/PhD School/Research/Weather_Master.xlsx", sheet = "Sheet1")
    weather_data$Date <- as.Date(weather_data$Date)
    weather_data <- weather_data %>% mutate(Tmean = (Tmax + Tmin) / 2, Daily_GDD = ifelse(Tmean > 5, Tmean - 5, 0))
  }
  compute_window_gdd <- function(env_row, weather) {
    loc   <- env_row$Location
    start <- env_row$Start_Date
    w     <- weather %>% filter(Location == loc, Date >= start) %>% arrange(Date)
    bounds <- cumsum(c(env_row$Int_E_T1, env_row$Int_T1_T2, env_row$Int_T2_T3, env_row$Int_T3_Tmax))
    tibble(
      Environment = env_row$Environment,
      Window      = c("E_T1","T1_T2","T2_T3","T3_Tmax"),
      Window_GDD  = c(
        sum(w$Daily_GDD[seq_len(bounds[1])], na.rm = TRUE),
        sum(w$Daily_GDD[(bounds[1]+1):bounds[2]], na.rm = TRUE),
        sum(w$Daily_GDD[(bounds[2]+1):bounds[3]], na.rm = TRUE),
        sum(w$Daily_GDD[(bounds[3]+1):min(bounds[4], nrow(w))], na.rm = TRUE)
      )
    )
  }
  window_gdd <- map_dfr(seq_len(nrow(env_metadata)), ~compute_window_gdd(env_metadata[.x, ], weather_data))
  write.csv(window_gdd, "Window_Specific_GDD_By_Environment.csv", row.names = FALSE)
}

gdd_timeline <- window_gdd %>%
  mutate(Window = factor(Window, levels = c("E_T1","T1_T2","T2_T3","T3_Tmax"))) %>%
  arrange(Environment, Window) %>%
  group_by(Environment) %>%
  mutate(Cum_GDD = cumsum(Window_GDD)) %>%
  ungroup()

gdd_lookup <- gdd_timeline %>%
  select(Environment, Window, Cum_GDD) %>%
  pivot_wider(names_from = Window, values_from = Cum_GDD) %>%
  rename(GDD_T1 = E_T1, GDD_T2 = T1_T2, GDD_T3 = T2_T3, GDD_Tmax = T3_Tmax) %>%
  mutate(GDD_E = 0)

long_tillers <- master_gxe %>%
  select(Environment, Genotype, T1_Count, T2_Count, T3_Count, Tmax_Count) %>%
  left_join(gdd_lookup, by = "Environment") %>%
  pivot_longer(
    cols      = c(T1_Count, T2_Count, T3_Count, Tmax_Count),
    names_to  = "Stage",
    values_to = "Count"
  ) %>%
  mutate(
    GDD = case_when(
      Stage == "T1_Count"   ~ GDD_T1,
      Stage == "T2_Count"   ~ GDD_T2,
      Stage == "T3_Count"   ~ GDD_T3,
      Stage == "Tmax_Count" ~ GDD_Tmax
    )
  ) %>%
  select(Environment, Genotype, Stage, Count, GDD)

origin_points <- master_gxe %>%
  select(Environment, Genotype) %>%
  mutate(Stage = "E", Count = 0, GDD = 0)

long_tillers <- bind_rows(long_tillers, origin_points) %>%
  arrange(Environment, Genotype, GDD)

# ------------------------------------------------------------------
# NEW APPROACH: Fit Logistic by Yield Group instead of Genotype
# (Resolves the 200 implausible/fragile 5-point fits)
# ------------------------------------------------------------------
extremes_global <- master_gxe %>%
  group_by(Environment) %>%
  arrange(desc(GrainYld)) %>%
  mutate(Yield_Rank = row_number()) %>%
  filter(Yield_Rank <= 5 | Yield_Rank >= (n() - 4)) %>%
  mutate(Group = ifelse(Yield_Rank <= 5, "Top 5 Yielding (Elite)", "Bottom 5 Yielding (Susceptible)")) %>%
  select(Environment, Genotype, Group) %>%
  ungroup()

write.csv(extremes_global, "Thermal_Logistic_Models/Top_Bottom_Yielders_By_Environment.csv", row.names = FALSE)

# Pool data by Environment + Group
pooled_tillers <- long_tillers %>%
  inner_join(extremes_global, by = c("Environment","Genotype")) %>%
  filter(!is.na(Count), !is.na(GDD))

extract_pooled_logistic <- function(df_sub) {
  # Now fitting on 25 points per curve (5 genotypes * 5 points) -> robust fit
  fit <- try(nls(Count ~ SSlogis(GDD, Asym, xmid, scal), data = df_sub), silent = TRUE)
  if (!inherits(fit, "try-error")) {
    p <- coef(fit)
    tibble(Thermal_Ceiling_Asym = p["Asym"], Thermal_Inflection_xmid = p["xmid"], Growth_Window_scal = p["scal"], Curve_Fit = "Success")
  } else {
    tibble(Thermal_Ceiling_Asym = NA_real_, Thermal_Inflection_xmid = NA_real_, Growth_Window_scal = NA_real_, Curve_Fit = "Failed")
  }
}

logistic_params_grouped <- pooled_tillers %>%
  group_by(Environment, Group) %>%
  nest() %>%
  mutate(Params = map(data, extract_pooled_logistic)) %>%
  unnest(Params) %>%
  select(-data)

write.csv(logistic_params_grouped, "Thermal_Logistic_Models/Group_Logistic_Parameters.csv", row.names = FALSE)

# Generate smooth predicted curves for plotting
if (!exists("thermal_summary")) {
  thermal_summary <- window_gdd %>% group_by(Environment) %>% summarise(Cumulative_GDD = sum(Window_GDD), .groups="drop")
}

max_gdd_by_env <- thermal_summary %>% select(Environment, Cumulative_GDD)

predicted_curves <- logistic_params_grouped %>%
  filter(Curve_Fit == "Success") %>%
  inner_join(max_gdd_by_env, by = "Environment") %>%
  group_by(Environment, Group, Thermal_Ceiling_Asym, Thermal_Inflection_xmid, Growth_Window_scal, Cumulative_GDD) %>%
  reframe(
    GDD = seq(0, Cumulative_GDD, length.out = 120),
    Count_Pred = Thermal_Ceiling_Asym / (1 + exp(-(GDD - Thermal_Inflection_xmid) / Growth_Window_scal))
  ) %>%
  ungroup()

logistic_plot <- ggplot() +
  geom_point(data = pooled_tillers, aes(x = GDD, y = Count, color = Group), alpha = 0.35, size = 1.8, position = position_jitter(width=10, height=0)) +
  geom_line(data = predicted_curves, aes(x = GDD, y = Count_Pred, color = Group), linewidth = 1.2, alpha = 0.95) +
  scale_color_manual(values = c("Top 5 Yielding (Elite)" = "#4575b4", "Bottom 5 Yielding (Susceptible)" = "#d73027")) +
  facet_wrap(~ Environment, scales = "free_x") +
  pub_theme +
  labs(
    title = "Non-Linear Thermal Growth Trajectories Across Stress Gradients",
    subtitle = "Top 5 vs Bottom 5 yielders pooled directly for robust SSlogis parameterisation",
    x = "Cumulative Thermal Load (GDD, °C·day)", y = "Tiller Count", color = "Yield Group"
  )

ggsave("Thermal_Logistic_Models/Elite_vs_Susceptible_Thermal_Curves_All_Envs.png", plot = logistic_plot, width = 12, height = 8, dpi = 600)

print("--- Advanced Analysis Module (Sections 15 & 16) Complete ---")
