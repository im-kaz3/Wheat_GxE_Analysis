# ==============================================================================
# MASTER PIPELINE: WHEAT GxE, MICRO-PHENOLOGY, AND CHRONIC VEGETATIVE STRESS
# Repository: Physiological and Molecular Basis of Tiller Number In Wheat Under Heat Stress
# Author: Apoorva Ashu, ICAR-IARI, New Delhi
# Date: September 2026
# Description: Evaluates 217 wheat genotypes under chronic vegetative thermal 
#              stress. Pre-calculates exact empirical tiller dynamics & AUTPC,
#              extracts BLUEs safely (Alpha-Lattice -> RCBD -> LM), computes TCE 
#              and TSI, and executes global and localized variance partitioning.
# ==============================================================================
# 
# ==============================================================================
# TABLE OF CONTENTS / INDEX
# ==============================================================================
# 1. SETUP, LIBRARIES & PUBLICATION THEME 
#    - Loads required statistical and advanced visualization packages.
# 2. RAW DATA DERIVATIONS, DIAGNOSTICS & NORMALIZATION 
#    - Pre-calculates empirical dynamics, AUTPC, and extracts BLUEs.
# 3. GLOBAL MERGE, CATEGORICAL TAGGING & COMPOSITES
#    - Merges environments, creates phenotypic summary table, plots Figure 1.
# 4. HIGH-IMPACT PHYSIOLOGICAL VISUALIZATIONS 
#    - Plots AUTPC vs Yield and Phase Transition Alluvial diagrams.
# 5. STEP-BY-STEP TCE, TSI VALIDATION & SELECTION MATRICES 
#    - Calculates local/global Tiller Stability Indices and categorized matrices.
# 6. GLOBAL GxE VARIANCE PARTITIONING (ANOVA) 
#    - Partitions Genotypic Repeatability, isolates Elite Targets (>30% heritability).
# 7. SINGLE LOCATION ANALYSIS (LOCALIZED GxE) 
#    - Executes specific structural stability tests locally across sowing treatments.
# 8. GGE MEAN VS STABILITY ANALYSIS 
#    - Generates GGE biplots for grain yield and productive tillers.
# 9. METEOROLOGICAL ANALYSIS: THERMAL VELOCITY & GDD 
#    - Quantifies cumulative thermal load and dark respiration across macro-phases.
# 10. PHENOLOGICAL CLUSTER ANALYSIS (K-MEANS & PCA)
#    - Partitions multidimensional physiological velocities across vegetative development.
# 11. SINK-CAPACITY RATIOS & CLUSTER SUPERIMPOSITION
#    - Multivariate mapping of reproductive efficiency vs. vegetative burden.
# 12. THERMAL ESCAPE & DURATION COMPENSATION (ANOVA)
#    - Proves short-duration genotypes protect yield against late-stage thermal abortion.
# ==============================================================================

# ==========================================
# 1. SETUP, LIBRARIES & PUBLICATION THEME ####
# ==========================================
# Core Data & Stats
library(lme4)
library(car)
library(emmeans)
library(readxl)
library(dplyr)
library(tidyr)
library(tibble)
library(lubridate)
library(metan)

# Advanced Visualization Suite
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(patchwork)
library(cowplot)
library(ggdist)      # For Raincloud plots
library(ggbump)      # For Crossover slope graphs
library(ggalluvial)  # For Phase Transition flow diagrams
library(ggExtra)     # For Marginal Density plots
library(ggsci)       # For Nature Publishing Group (NPG) color palettes

# Define Universal High-Impact Publication Theme
pub_theme <- theme_classic(base_size = 12, base_family = "Helvetica") + 
  theme(
    strip.background = element_blank(), 
    strip.text = element_text(face = "bold", size = 12),
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

# ==========================================
# 2. RAW DATA DERIVATIONS, DIAGNOSTICS & NORMALIZATION ####
# ==========================================
# INSTRUCTION: Execute this section 4 times independently. 
# Update the target sheet and output object name for each run.

trial_data <- read_excel("D:/PhD School/Research/FinalDataYear1.xlsx", sheet = "Delhi_S1") # UPDATE PER RUN

# A. Enforce Factor Structure
trial_data$Genotype <- as.factor(trial_data$Genotype)
trial_data$Replication <- as.factor(trial_data$Replication)
trial_data$Block <- as.factor(trial_data$Block)

# Set Location and Sowing manually for this specific run
current_location <- "Delhi"  # UPDATE PER RUN
current_sowing <- "S1"       # UPDATE PER RUN

# B. PRE-CALCULATE DERIVED TRAITS (To ensure they undergo spatial normalization)
trial_data <- trial_data %>%
  mutate(
    # 1. Define exact calendar durations based on empirical DAS timeline
    Int_E_T1 = case_when(
      current_location == "Delhi" & current_sowing == "S1" ~ 26, 
      current_location == "Delhi" & current_sowing == "S2" ~ 27, 
      current_location == "Dharwad" & current_sowing == "S1" ~ 10, 
      current_location == "Dharwad" & current_sowing == "S2" ~ 12  
    ),
    
    # Subsequent intervals are uniform per location
    Int_T1_T2 = ifelse(current_location == "Delhi", 15, 10),
    Int_T2_T3 = ifelse(current_location == "Delhi", 15, 10),
    Int_T3_Tmax = ifelse(current_location == "Delhi", 15, 10),
    Total_T1_Tmax = ifelse(current_location == "Delhi", 45, 30), 
    
    # 2. Absolute Change (Delta)
    Delta_T1_T2 = T2_Count - T1_Count,
    Delta_T2_T3 = T3_Count - T2_Count,
    Delta_T3_Tmax = Tmax_Count - T3_Count,
    Delta_T1_Tmax = Tmax_Count - T1_Count,
    
    # 3. Area Under the Tillering Progress Curve (AUTPC) via Trapezoidal Rule
    AUTPC = (0.5 * (0 + T1_Count) * Int_E_T1) +         
      (0.5 * (T1_Count + T2_Count) * Int_T1_T2) +        
      (0.5 * (T2_Count + T3_Count) * Int_T2_T3) +        
      (0.5 * (T3_Count + Tmax_Count) * Int_T3_Tmax), 
    
    Relative_AUTPC = AUTPC / Tmax_Count,
    
    # 4. Fractional Rates (Exact physiological velocity)
    Rate_E_T1 = T1_Count / Int_E_T1,
    Rate_T1_T2 = Delta_T1_T2 / Int_T1_T2,
    Rate_T2_T3 = Delta_T2_T3 / Int_T2_T3
  )

# C. Define Traits to Model (Isolate strictly numeric variables)
design_cols <- c("Sr_No", "Genotype", "Replication", "Block", "Int_E_T1", "Int_T1_T2", "Int_T2_T3", "Int_T3_Tmax", "Total_T1_Tmax")
trait_cols <- setdiff(colnames(trial_data), design_cols)
trait_cols <- trait_cols[sapply(trial_data[trait_cols], is.numeric)] 

blues_matrix <- data.frame(Genotype = unique(trial_data$Genotype))
dir.create("Diagnostics", showWarnings = FALSE)

# D. Generate Shapiro-Wilk and Levene's Test Diagnostics Table
normality_results <- data.frame()
for (trait in trait_cols) {
  shapiro_p <- tryCatch(shapiro.test(trial_data[[trait]])$p.value, error = function(e) NA)
  levene_p <- tryCatch(leveneTest(trial_data[[trait]], trial_data$Replication)$`Pr(>F)`[1], error = function(e) NA)
  normality_results <- rbind(normality_results, data.frame(Trait = trait, Shapiro_P = shapiro_p, Levene_P = levene_p))
}
write.csv(normality_results, paste0("Diagnostics/Normality_Tests_", current_location, "_", current_sowing, ".csv"), row.names = FALSE)

# E. Master Loop: Spatial Normalization (BLUEs) & Q-Q Plots
for (trait in trait_cols) {
  # Step-down spatial modeling to prevent singularity crashes
  formula_alpha <- paste(trait, "~ Genotype + (1 | Replication/Block)")
  formula_rcbd  <- paste(trait, "~ Genotype + (1 | Replication)")
  formula_lm    <- paste(trait, "~ Genotype + Replication") 
  
  model_alpha <- try(lmer(as.formula(formula_alpha), data = trial_data, na.action = na.exclude, 
                          control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))), silent = TRUE)
  
  if (!inherits(model_alpha, "try-error") && !isSingular(model_alpha)) {
    model <- model_alpha
  } else {
    model_rcbd <- try(lmer(as.formula(formula_rcbd), data = trial_data, na.action = na.exclude,
                           control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))), silent = TRUE)
    if (!inherits(model_rcbd, "try-error") && !isSingular(model_rcbd)) {
      model <- model_rcbd
    } else {
      model <- lm(as.formula(formula_lm), data = trial_data, na.action = na.exclude)
    }
  }
  
  # Diagnostics Export
  res <- residuals(model)
  p_qq <- ggplot(data.frame(res), aes(sample = res)) + 
    stat_qq() + stat_qq_line(color = "red") + 
    theme_bw() + labs(title = paste("Q-Q Plot:", trait))
  p_hist <- ggplot(data.frame(res), aes(x = res)) + 
    geom_histogram(bins = 30, fill = "steelblue", color = "black", alpha = 0.7) + 
    theme_bw() + labs(title = paste("Residuals:", trait))
  
  diag_plot <- ggarrange(p_qq, p_hist, ncol = 2)
  ggsave(paste0("Diagnostics/", trait, "_Normalization_", current_location, "_", current_sowing, ".png"), plot = diag_plot, width = 8, height = 4)
  
  # Extract BLUEs safely
  em_out <- suppressMessages(as.data.frame(emmeans(model, "Genotype")))
  trait_blues <- em_out[, c("Genotype", "emmean")]
  colnames(trait_blues)[2] <- trait
  blues_matrix <- merge(blues_matrix, trait_blues, by = "Genotype", all.x = TRUE)
}

blues_matrix$Location <- current_location
blues_matrix$Sowing <- current_sowing
blues_matrix$Environment <- paste0(current_location, "_", current_sowing)

# SAVE OBJECT 
blues_Delhi_S1 <- blues_matrix 

# [STOP & REPEAT SECTION 2 UNTIL ALL 4 ENVIRONMENTS ARE PROCESSED]


# ==========================================
# 3. GLOBAL MERGE, CATEGORICAL TAGGING & COMPOSITES ####
# ==========================================
master_gxe <- rbind(blues_Delhi_S1, blues_Delhi_S2, blues_Dharwad_S1, blues_Dharwad_S2)
master_gxe$Genotype <- as.factor(master_gxe$Genotype)
master_gxe$Location <- as.factor(master_gxe$Location)
master_gxe$Sowing <- as.factor(master_gxe$Sowing)
master_gxe$Environment <- as.factor(master_gxe$Environment)

# Apply Biological Categorical Tags
master_gxe <- master_gxe %>%
  mutate(
    Peak_Phase = case_when(
      Rate_E_T1 > Rate_T1_T2 & Rate_E_T1 > Rate_T2_T3 ~ "Very Early (E-T1)",
      Rate_T1_T2 >= Rate_E_T1 & Rate_T1_T2 > Rate_T2_T3 ~ "Early Vigor (T1-T2)",
      Rate_T2_T3 >= Rate_T1_T2 & Rate_T2_T3 >= Rate_E_T1 ~ "Late Vigor (T2-T3)",
      TRUE ~ "Uniform"
    ),
    Phase_T1_T2 = case_when(Delta_T1_T2 > 0.05 ~ "Active Emergence", Delta_T1_T2 < -0.05 ~ "Active Abortion", TRUE ~ "Stagnation"),
    Phase_T2_T3 = case_when(Delta_T2_T3 > 0.05 ~ "Active Emergence", Delta_T2_T3 < -0.05 ~ "Active Abortion", TRUE ~ "Stagnation"),
    Phase_T3_Tmax = case_when(Delta_T3_Tmax > 0.05 ~ "Active Emergence", Delta_T3_Tmax < -0.05 ~ "Active Abortion", TRUE ~ "Stagnation"),
    Phase_T1_Tmax = case_when(Delta_T1_Tmax > 0.05 ~ "Active Emergence", Delta_T1_Tmax < -0.05 ~ "Active Abortion", TRUE ~ "Stagnation")
  )

# ---------------------------------------------------------
# COMPREHENSIVE PHENOTYPIC SUMMARY TABLE
# ---------------------------------------------------------
summary_table <- data.frame()
numeric_traits <- c("PT_Count", "Tmax_Count", "GrainYld", "TKW", "AUTPC", "TCE")
master_gxe$TCE <- master_gxe$PT_Count / master_gxe$Tmax_Count

for (trait in numeric_traits) {
  # Calculate ANOVA p-value for Environmental differences
  formula_anova <- paste(trait, "~ Environment")
  model <- tryCatch(lm(as.formula(formula_anova), data = master_gxe), error = function(e) NULL)
  p_val <- ifelse(!is.null(model), anova(model)$"Pr(>F)"[1], NA)
  
  # Calculate Descriptive Statistics per Environment
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

# ---------------------------------------------------------
# FIGURE 1: MACROSCOPIC PHENOTYPIC RESPONSE (RAINCLOUD PLOTS)
# ---------------------------------------------------------
# Ensure environment factor is ordered chronologically
master_gxe$Environment <- factor(master_gxe$Environment, levels = c("Delhi_S1", "Delhi_S2", "Dharwad_S1", "Dharwad_S2"))

fig1a <- ggplot(master_gxe, aes(x = Environment, y = GrainYld, fill = Environment)) +
  stat_halfeye(adjust = 0.5, width = 0.6, .width = 0, justification = -0.3, point_colour = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.5) +
  geom_point(aes(color = Environment), size = 1.3, alpha = 0.3, position = position_jitter(width = 0.1, seed = 42)) +
  scale_fill_npg() + scale_color_npg() + pub_theme + theme(legend.position = "none") +
  labs(title = "A. Grain Yield Distribution Across Regimes", y = "Grain Yield (t/ha)", x = "")

fig1b <- ggplot(master_gxe, aes(x = Environment, y = Tmax_Count, fill = Environment)) +
  stat_halfeye(adjust = 0.5, width = 0.6, .width = 0, justification = -0.3, point_colour = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.5) +
  geom_point(aes(color = Environment), size = 1.3, alpha = 0.3, position = position_jitter(width = 0.1, seed = 42)) +
  scale_fill_npg() + scale_color_npg() + pub_theme + theme(legend.position = "none") +
  labs(title = expression(bold("B. Maximum Tillers ("*T[max]*") Distribution")), y = "Tiller Count", x = "")

composite_fig1 <- fig1a / fig1b
ggsave("Figure_1_Raincloud_Distributions.png", plot = composite_fig1, width = 10, height = 12, dpi = 600)


# ==========================================
# 4. HIGH-IMPACT PHYSIOLOGICAL VISUALIZATIONS ####
# ==========================================
# ---------------------------------------------------------
# FIGURE 2: MICRO-PHENOLOGY (ALLUVIAL PHASE TRANSITIONS)
# ---------------------------------------------------------
alluvial_data <- master_gxe %>%
  group_by(Environment, Phase_T1_T2, Phase_T2_T3, Phase_T3_Tmax) %>%
  tally()

fig2a <- ggplot(alluvial_data, aes(y = n, axis1 = Phase_T1_T2, axis2 = Phase_T2_T3, axis3 = Phase_T3_Tmax)) +
  geom_alluvium(aes(fill = Phase_T1_T2), width = 1/12, alpha = 0.7) +
  geom_stratum(width = 1/12, fill = "grey90", color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3, fontface = "bold") +
  facet_wrap(~ Environment, scales = "free_y") +
  scale_fill_npg() + pub_theme + 
  theme(legend.position = "none", axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank()) +
  scale_x_discrete(limits = c("T1-T2", "T2-T3", "T3-Tmax"), expand = c(.05, .05)) +
  labs(title = "Developmental Phase Transitions (Crossover Flow)", x = "Phenological Interval", y = "Genotype Density")

ggsave("Figure_2_Alluvial_Transitions.png", plot = fig2a, width = 12, height = 8, dpi = 600)

# ---------------------------------------------------------
# FIGURE 3: YIELD DETERMINANTS WITH MARGINAL DENSITY
# ---------------------------------------------------------
p_scatter <- ggplot(master_gxe, aes(x = Relative_AUTPC, y = GrainYld, color = Peak_Phase)) +
  geom_point(alpha = 0.6, size = 2) + 
  geom_smooth(method = "lm", color = "black", linewidth = 1) +
  facet_grid(Location ~ Sowing) +
  scale_color_npg() + pub_theme +
  labs(title = "Dependence of Grain Yield on Early Tiller Establishment",
       x = "Relative Area Under Tillering Progress Curve (AUTPC)", y = "Grain Yield (t/ha)", color = "Peak Phase")

# Apply ggMarginal to add distribution density curves to the axes
fig3_marginal <- ggMarginal(p_scatter, type = "density", groupColour = TRUE, groupFill = TRUE, alpha = 0.4)
save_plot("Figure_3_Marginal_Yield_Determinants.png", fig3_marginal, base_height = 8, base_width = 11, dpi = 600)


# ==========================================
# 5. STEP-BY-STEP TCE, TSI VALIDATION & SELECTION MATRICES ####
# ==========================================
# Calculate base Tiller Conversion Efficiency (TCE) and set biological regimes
efficiency_base <- master_gxe %>%
  mutate(
    TCE = PT_Count / Tmax_Count,
    Stress_Regime = case_when(
      Environment %in% c("Delhi_S2", "Dharwad_S1") ~ "Optimal_Combined",
      Environment %in% c("Delhi_S1", "Dharwad_S2") ~ "Stress_Combined"
    )
  )

# ---------------------------------------------------------
# STEP 5A: REGIONAL BASELINE ANALYSIS (DELHI ONLY)
# ---------------------------------------------------------
delhi_metrics <- efficiency_base %>% filter(Location == "Delhi")

delhi_summary <- delhi_metrics %>%
  group_by(Genotype, Sowing) %>%
  summarise(Mean_PT = mean(PT_Count, na.rm = TRUE), 
            Mean_Yield = mean(GrainYld, na.rm = TRUE), 
            Mean_TCE = mean(TCE, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = Sowing, values_from = c(Mean_PT, Mean_Yield, Mean_TCE))

pop_PT_Delhi_S1 <- mean(delhi_summary$Mean_PT_S1, na.rm = TRUE)       
pop_Yield_Delhi_S2 <- mean(delhi_summary$Mean_Yield_S2, na.rm = TRUE) 
pop_TCE_Delhi_S2 <- mean(delhi_summary$Mean_TCE_S2, na.rm = TRUE)     
pop_TCE_Delhi_S1 <- mean(delhi_summary$Mean_TCE_S1, na.rm = TRUE)     

delhi_index <- delhi_summary %>%
  mutate(
    TSI_Delhi = (Mean_PT_S1 / Mean_PT_S2) * (Mean_PT_S1 / pop_PT_Delhi_S1),
    Yield_STI_Delhi = (Mean_Yield_S2 * Mean_Yield_S1) / (pop_Yield_Delhi_S2^2),
    Category_Delhi = case_when(
      Mean_TCE_S2 >= pop_TCE_Delhi_S2 & Mean_TCE_S1 >= pop_TCE_Delhi_S1 ~ "Stable High Efficiency",
      Mean_TCE_S2 >= pop_TCE_Delhi_S2 & Mean_TCE_S1 < pop_TCE_Delhi_S1 ~ "Contrasting Susceptible",
      Mean_TCE_S2 < pop_TCE_Delhi_S2 & Mean_TCE_S1 >= pop_TCE_Delhi_S1 ~ "Contrasting Tolerant",
      TRUE ~ "Stable Low Efficiency"
    )
  )

tsi_delhi_plot <- ggplot(delhi_index, aes(x = TSI_Delhi, y = Yield_STI_Delhi)) +
  geom_point(aes(color = Category_Delhi), alpha = 0.7, size = 2.5) +
  geom_smooth(method = "lm", color = "black", fill = "lightgrey") +
  stat_cor(size = 5, fontface = "bold", label.x.npc = "left", label.y.npc = "top") +
  scale_color_npg() + pub_theme +
  labs(title = "Regional Validation: Tiller Stability Index (Delhi)",
       subtitle = "Predicting terminal yield stability (Optimal S2 vs Stress S1)",
       x = "Tiller Stability Index (Delhi)", y = "Yield STI (Delhi)", color = "TCE Category")
ggsave("TSI_Validation_Delhi.png", plot = tsi_delhi_plot, width = 9, height = 6, dpi = 600)
write.csv(delhi_index, "Molecular_Selection_Matrix_Delhi.csv", row.names = FALSE)

# ---------------------------------------------------------
# STEP 5B: GLOBAL CROSS-ZONAL ANALYSIS (COMBINED REGIMES)
# ---------------------------------------------------------
combined_summary <- efficiency_base %>%
  group_by(Genotype, Stress_Regime) %>%
  summarise(Mean_PT = mean(PT_Count, na.rm = TRUE), 
            Mean_Yield = mean(GrainYld, na.rm = TRUE), 
            Mean_TCE = mean(TCE, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = Stress_Regime, values_from = c(Mean_PT, Mean_Yield, Mean_TCE))

pop_PT_Stress_Combined <- mean(combined_summary$Mean_PT_Stress_Combined, na.rm = TRUE)
pop_Yield_Opt_Combined <- mean(combined_summary$Mean_Yield_Optimal_Combined, na.rm = TRUE)
pop_TCE_Opt_Combined <- mean(combined_summary$Mean_TCE_Optimal_Combined, na.rm = TRUE)
pop_TCE_Stress_Combined <- mean(combined_summary$Mean_TCE_Stress_Combined, na.rm = TRUE)

global_index <- combined_summary %>%
  mutate(
    TSI_Global = (Mean_PT_Stress_Combined / Mean_PT_Optimal_Combined) * (Mean_PT_Stress_Combined / pop_PT_Stress_Combined),
    Yield_STI_Global = (Mean_Yield_Optimal_Combined * Mean_Yield_Stress_Combined) / (pop_Yield_Opt_Combined^2),
    Category_Global = case_when(
      Mean_TCE_Optimal_Combined >= pop_TCE_Opt_Combined & Mean_TCE_Stress_Combined >= pop_TCE_Stress_Combined ~ "Stable High Efficiency",
      Mean_TCE_Optimal_Combined >= pop_TCE_Opt_Combined & Mean_TCE_Stress_Combined < pop_TCE_Stress_Combined ~ "Contrasting Susceptible",
      Mean_TCE_Optimal_Combined < pop_TCE_Opt_Combined & Mean_TCE_Stress_Combined >= pop_TCE_Stress_Combined ~ "Contrasting Tolerant",
      TRUE ~ "Stable Low Efficiency"
    )
  )

tsi_global_plot <- ggplot(global_index, aes(x = TSI_Global, y = Yield_STI_Global)) +
  geom_point(aes(color = Category_Global), alpha = 0.7, size = 2.5) +
  geom_smooth(method = "lm", color = "black", fill = "lightgrey") +
  stat_cor(size = 5, fontface = "bold", label.x.npc = "left", label.y.npc = "top") +
  scale_color_npg() + pub_theme +
  labs(title = "Global Validation: Cross-Zonal Tiller Stability Index",
       subtitle = "Predicting macro-environmental yield stability across extreme thermal gradients",
       x = "Global Tiller Stability Index (Optimal vs Stress)", y = "Global Yield STI", color = "TCE Category")
ggsave("TSI_Validation_Combined_Global.png", plot = tsi_global_plot, width = 9, height = 6, dpi = 600)
write.csv(global_index, "Molecular_Selection_Matrix_Combined_Global.csv", row.names = FALSE)


# ==========================================
# 6. GLOBAL GxE VARIANCE PARTITIONING (ANOVA) ####
# ==========================================
variance_results <- data.frame()
core_traits <- c("PT_Count", "Tmax_Count", "GrainYld", "TKW", "AUTPC", "Relative_AUTPC")

for (trait in core_traits) {
  # Standard linear model for BLUEs. The residual captures GxE.
  formula_anova <- paste(trait, "~ Genotype + Environment")
  model_anova <- tryCatch(lm(as.formula(formula_anova), data = master_gxe), error = function(e) return(NULL))
  
  if(!is.null(model_anova)) {
    anova_table <- car::Anova(model_anova, type = 2)
    ss_genotype <- anova_table["Genotype", "Sum Sq"]
    ss_environment <- anova_table["Environment", "Sum Sq"]
    ss_residual <- anova_table["Residuals", "Sum Sq"] 
    ss_total <- sum(ss_genotype, ss_environment, ss_residual)
    
    trait_variance <- data.frame(
      Trait = trait,
      Pct_Genotype = (ss_genotype / ss_total) * 100,
      Pct_Environment = (ss_environment / ss_total) * 100,
      Pct_GxE_and_Residual = (ss_residual / ss_total) * 100
    )
    variance_results <- rbind(variance_results, trait_variance)
  }
}

# Isolate Elite Breeding Targets (>30% Genotypic Variance)
variance_results <- variance_results %>%
  mutate(Breeder_Target = ifelse(Pct_Genotype > 30, "Highly Heritable (>30%)", "Environmentally Driven"))
write.csv(variance_results, "Global_GxE_Variance_Partitioning.csv", row.names = FALSE)

elite_targets <- variance_results %>% filter(Pct_Genotype > 30)
write.csv(elite_targets, "Elite_Breeding_Targets.csv", row.names = FALSE)

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
# 7. SINGLE LOCATION ANALYSIS (LOCALIZED GxE) ####
# ==========================================
run_local_analysis <- function(loc_name) {
  loc_data <- subset(master_gxe, Location == loc_name)
  local_variance <- data.frame()
  
  for (trait in core_traits) {
    formula_local <- paste(trait, "~ Genotype + Sowing")
    model_local <- tryCatch(lm(as.formula(formula_local), data = loc_data), error = function(e) return(NULL))
    
    if(!is.null(model_local)) {
      anova_table <- car::Anova(model_local, type = 2)
      ss_genotype <- anova_table["Genotype", "Sum Sq"]
      ss_sowing <- anova_table["Sowing", "Sum Sq"]
      ss_residual <- anova_table["Residuals", "Sum Sq"] 
      ss_total <- sum(ss_genotype, ss_sowing, ss_residual)
      
      trait_var <- data.frame(
        Location = loc_name, Trait = trait,
        Pct_Genotype = (ss_genotype / ss_total) * 100,
        Pct_Sowing = (ss_sowing / ss_total) * 100,
        Pct_GxS_and_Residual = (ss_residual / ss_total) * 100
      )
      local_variance <- rbind(local_variance, trait_var)
    }
  }
  write.csv(local_variance, paste0(loc_name, "_Local_Variance_Partitioning.csv"), row.names = FALSE)
  
  var_plot_data_local <- local_variance %>%
    pivot_longer(cols = starts_with("Pct_"), names_to = "Component", values_to = "Percentage") %>%
    mutate(Component = factor(Component, levels = c("Pct_GxS_and_Residual", "Pct_Genotype", "Pct_Sowing")))
  
  p <- ggplot(var_plot_data_local, aes(x = reorder(Trait, -Percentage), y = Percentage, fill = Component)) +
    geom_bar(stat = "identity", position = "stack", color = "black") +
    scale_fill_manual(values = c("Pct_Sowing" = "#d73027", "Pct_Genotype" = "#4575b4", "Pct_GxS_and_Residual" = "#fdae61"),
                      labels = c("G×S + Residual", "Genotype", "Sowing Regime")) +
    theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = paste("Local Genotypic Repeatability:", loc_name), 
         subtitle = "Evaluating structural stability exclusively across local sowing treatments",
         x = "Phenotypic Trait", y = "Proportion of Variance (%)", fill = "Component")
  
  ggsave(paste0(loc_name, "_Local_Variance.png"), plot = p, width = 10, height = 6, dpi = 600)
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
weather_data <- read_excel("D:/PhD School/Research/Weather_Master.xlsx")
weather_data$Date <- as.Date(weather_data$Date)
base_temp <- 5 

weather_data <- weather_data %>%
  mutate(Tmean = (Tmax + Tmin) / 2, Daily_GDD = ifelse(Tmean > base_temp, Tmean - base_temp, 0))

phase_intervals <- data.frame(
  Environment = c("Delhi_S1", "Delhi_S2", "Dharwad_S1", "Dharwad_S2"),
  Location = c("Delhi", "Delhi", "Dharwad", "Dharwad"),
  Start_Date = as.Date(c("2025-11-12", "2025-12-04", "2025-11-15", "2025-11-27")), 
  End_Date = as.Date(c("2026-02-05", "2026-02-20", "2026-01-14", "2026-01-29"))    
)

thermal_summary <- data.frame()
for(i in 1:nrow(phase_intervals)) {
  env_weather <- weather_data %>% filter(Location == phase_intervals$Location[i], Date >= phase_intervals$Start_Date[i] & Date <= phase_intervals$End_Date[i])
  summary_row <- data.frame(
    Environment = phase_intervals$Environment[i], Total_Days = nrow(env_weather),
    Cumulative_GDD = sum(env_weather$Daily_GDD, na.rm = TRUE),
    GDD_Velocity = sum(env_weather$Daily_GDD, na.rm = TRUE) / nrow(env_weather), 
    Mean_Tmax = mean(env_weather$Tmax, na.rm = TRUE), Mean_Tmin = mean(env_weather$Tmin, na.rm = TRUE) 
  )
  thermal_summary <- rbind(thermal_summary, summary_row)
}
write.csv(thermal_summary, "Meteorological_Thermal_Summary.csv", row.names = FALSE)


# ==========================================
# 10. PHENOLOGICAL CLUSTER ANALYSIS (K-MEANS & PCA) ####
# ==========================================
cluster_vars <- c("Rate_E_T1", "Rate_T1_T2", "Rate_T2_T3")

pheno_df <- master_gxe %>%
  group_by(Genotype) %>%
  summarise(across(all_of(cluster_vars), mean, na.rm = TRUE)) %>%
  drop_na() %>%
  column_to_rownames("Genotype")

pheno_scaled <- scale(pheno_df)
set.seed(42) 
kmeans_result <- kmeans(pheno_scaled, centers = 3, nstart = 25)
pheno_df$Cluster <- as.factor(kmeans_result$cluster)

pca_result <- prcomp(pheno_scaled, center = FALSE, scale. = FALSE)
pheno_df$PC1 <- pca_result$x[, 1]
pheno_df$PC2 <- pca_result$x[, 2]

cluster_export <- pheno_df %>% rownames_to_column("Genotype")
write.csv(cluster_export, "Phenological_Cluster_Assignments.csv", row.names = FALSE)

var_explained <- round(100 * pca_result$sdev^2 / sum(pca_result$sdev^2), 1)
cluster_plot <- ggplot(pheno_df, aes(x = PC1, y = PC2, fill = Cluster)) +
  geom_point(shape = 21, size = 4, alpha = 0.8, color = "black") +
  stat_ellipse(aes(color = Cluster), type = "norm", linetype = "dashed", linewidth = 1) +
  scale_fill_npg() + scale_color_npg() + pub_theme +
  labs(title = "K-Means Clustering of Genotypic Phenological Strategies",
       subtitle = "Partitioning multidimensional physiological velocities across vegetative development",
       x = paste0("PC1 (", var_explained[1], "%)"), y = paste0("PC2 (", var_explained[2], "%)"))

ggsave("Phenological_Clusters_PCA.png", plot = cluster_plot, width = 9, height = 7, dpi = 600)


# ==========================================
# 11. SINK-CAPACITY RATIOS & CLUSTER SUPERIMPOSITION ####
# ==========================================
ratio_df <- master_gxe %>%
  mutate(
    Ratio_TCE = PT_Count / Tmax_Count,           
    Ratio_Total_Yield = Tmax_Count / GrainYld,   
    Ratio_PT_Yield = PT_Count / GrainYld         
  ) %>%
  drop_na(Ratio_TCE, Ratio_Total_Yield, Ratio_PT_Yield) 

ratio_summary <- ratio_df %>%
  group_by(Genotype) %>%
  summarise(across(starts_with("Ratio_"), mean)) %>%
  column_to_rownames("Genotype")

ratio_pca <- prcomp(ratio_summary, center = TRUE, scale. = TRUE)

ratio_plot_data <- as.data.frame(ratio_pca$x) %>%
  rownames_to_column("Genotype") %>%
  left_join(pheno_df %>% rownames_to_column("Genotype") %>% select(Genotype, Cluster), by = "Genotype")

ratio_biplot <- ggplot(ratio_plot_data, aes(x = PC1, y = PC2, fill = Cluster)) +
  geom_point(shape = 21, size = 4, alpha = 0.8, color = "black") +
  stat_ellipse(aes(color = Cluster), type = "norm", linetype = "dashed", linewidth = 1) +
  scale_fill_npg() + scale_color_npg() + pub_theme +
  labs(title = "Superimposition of Sink-Capacity Ratios",
       subtitle = "Multivariate mapping of reproductive efficiency vs. vegetative burden",
       x = "PC1: Vegetative Burden (Total Tillers / Yield)", y = "PC2: Conversion Efficiency (PT / Total)", fill = "Pheno-Cluster")

ggsave("Ratio_Cluster_Superimposition.png", plot = ratio_biplot, width = 9, height = 7, dpi = 600)


# ==========================================
# 12. THERMAL ESCAPE & DURATION COMPENSATION (ANOVA) ####
# ==========================================
duration_analysis <- master_gxe %>%
  filter(Peak_Phase != "Uniform") %>% 
  mutate(
    Duration_Class = case_when(
      Peak_Phase %in% c("Very Early (E-T1)", "Early Vigor (T1-T2)") ~ "Early (Short Duration)",
      Peak_Phase == "Late Vigor (T2-T3)" ~ "Late (Long Duration)"
    ),
    Environment = factor(Environment, levels = c("Delhi_S2", "Delhi_S1", "Dharwad_S1", "Dharwad_S2"))
  )

yield_compensation_plot <- ggplot(duration_analysis, aes(x = Environment, y = GrainYld, fill = Duration_Class)) +
  geom_boxplot(alpha = 0.85, outlier.shape = 21, color = "black") +
  scale_fill_npg() + pub_theme +
  labs(title = "Yield Compensation via Thermal Escape",
       subtitle = "Short-duration genotypes protecting yield under severe vegetative stress",
       x = "Testing Regime", y = "Grain Yield (t/ha)", fill = "Genotypic Strategy") +
  stat_compare_means(aes(group = Duration_Class), label = "p.signif", method = "t.test")

ggsave("Yield_Compensation_Thermal_Escape.png", plot = yield_compensation_plot, width = 10, height = 7, dpi = 600)

abortion_penalty_plot <- ggplot(duration_analysis, aes(x = Duration_Class, y = Delta_T3_Tmax, fill = Duration_Class)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.2, fill = "white", outlier.shape = NA) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 1) + 
  facet_wrap(~ Environment, nrow = 1) +
  scale_fill_npg() + pub_theme +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  labs(title = "The Temperature Trap: Late-Stage Tiller Abortion",
       subtitle = "Negative Delta indicates active abortion as fixed-interval heat hits late-developing lines",
       x = "", y = "Tiller Change (Delta T3 to Tmax)")

ggsave("Late_Stage_Abortion_Penalty.png", plot = abortion_penalty_plot, width = 12, height = 6, dpi = 600)

# Statistical Proof (ANOVA for Interaction)
sink("Thermal_Escape_ANOVA_Results.txt")
cat("--- Thermal Escape ANOVA (Duration x Environment) ---\n")
compensation_model <- aov(GrainYld ~ Duration_Class * Environment, data = duration_analysis)
print(summary(compensation_model))
sink()


# ==========================================
# MASTER DATA EXPORT & COMPLETION
# ==========================================
write.csv(master_gxe, "Master_GxE_Matrix_Complete.csv", row.names = FALSE)
print("--- Master Pipeline Execution Complete ---")
