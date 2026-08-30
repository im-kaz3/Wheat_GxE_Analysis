# ==============================================================================
# MASTER PIPELINE: WHEAT GxE, ABIOTIC STRESS, AND PHYSIOLOGICAL INDICES
# ==============================================================================

# ==========================================
# 1. SETUP & LIBRARIES ####
# ==========================================
# Install missing packages if necessary: 
# install.packages(c("lme4", "car", "emmeans", "readxl", "ggplot2", "tidyr", "dplyr", "lubridate", "ggrepel", "metan", "ggpubr"))

library(lme4)
library(car)
library(emmeans)
library(readxl)
library(ggplot2)
library(tidyr)
library(dplyr)
library(lubridate)
library(ggrepel)
library(metan)
library(ggpubr)

# ==========================================
# 2. SINGLE ENVIRONMENT EXTRACTION & DIAGNOSTICS ####
# ==========================================
# INSTRUCTION: You must run this entire section 4 separate times (once for each dataset).
# Update the file path/sheet and the metadata tags at the end of this section each time.

# A. Import Raw Data
trial_data <- read_excel("Data_Dharwad.xlsx", sheet = "Normal_Final") # UPDATE THIS LINE FOR EACH RUN

# B. Enforce Factor Structure
trial_data$Genotype <- as.factor(trial_data$Genotype)
trial_data$Replication <- as.factor(trial_data$Replication)
trial_data$Block <- as.factor(trial_data$Block)

# C. Initialize Matrices and Define Exclusions (Note: Sr_No spelling corrected)
design_cols <- c("Sr_No", "Genotype", "Replication", "Block")
trait_cols <- setdiff(colnames(trial_data), design_cols)

blues_matrix <- data.frame(Genotype = unique(trial_data$Genotype))
diagnostics_matrix <- data.frame(Trait = character(), Shapiro_P = numeric(), Levene_P = numeric())

# D. Master Loop: Alpha Lattice, Diagnostics, Visuals, and BLUEs Extraction
for (trait in trait_cols) {
  
  # Define models
  formula_alpha <- paste(trait, "~ Genotype + (1 | Replication/Block)")
  formula_rcbd  <- paste(trait, "~ Genotype + (1 | Replication)")
  
  # Fit Model: Try Alpha Lattice, fallback to RCBD if matrix collapses (zero spatial variance)
  model <- tryCatch({
    lmer(as.formula(formula_alpha), data = trial_data, na.action = na.exclude)
  }, error = function(e) {
    message(paste("Alpha Lattice failed for", trait, "- Falling back to RCBD."))
    lmer(as.formula(formula_rcbd), data = trial_data, na.action = na.exclude)
  })
  
  # Extract Diagnostics
  resids <- residuals(model)
  fitted_vals <- fitted(model)
  shapiro_p <- tryCatch(shapiro.test(resids)$p.value, error = function(e) NA)
  levene_p <- tryCatch(leveneTest(resids ~ trial_data$Replication)$`Pr(>F)`[1], error = function(e) NA)
  
  diagnostics_matrix <- rbind(diagnostics_matrix, 
                              data.frame(Trait = trait, Shapiro_P = shapiro_p, Levene_P = levene_p))
  
  # Export Diagnostic Plots
  png_name <- paste0("Diagnostic_", trait, ".png")
  png(filename = png_name, width = 1200, height = 500, res = 120)
  par(mfrow = c(1, 2))
  
  qqnorm(resids, main = sprintf("Q-Q Plot: %s\n(Shapiro p = %.4f)", trait, shapiro_p), pch = 16, col = "darkblue")
  qqline(resids, col = "red", lwd = 2)
  plot(fitted_vals, resids, main = sprintf("Residuals vs Fitted: %s\n(Levene p = %.4f)", trait, levene_p),
       xlab = "Fitted Values", ylab = "Residuals", pch = 16, col = "darkcyan")
  abline(h = 0, col = "red", lwd = 2)
  
  dev.off()
  
  # Extract Raw BLUEs
  em_out <- suppressMessages(as.data.frame(emmeans(model, "Genotype")))
  trait_blues <- em_out[, c("Genotype", "emmean")]
  colnames(trait_blues)[2] <- trait
  blues_matrix <- merge(blues_matrix, trait_blues, by = "Genotype", all.x = TRUE)
}

# E. Apply Universal Transformations (Based on committee review of diagnostics)
traits_to_log <- c("T1", "T2", "T_Max") # UPDATE THESE ARRAYS BASED ON YOUR FINAL DECISION
traits_to_sqrt <- c("PT", "TKW")        

for (trait in traits_to_log) {
  trial_data[[paste0(trait, "_Log")]] <- log(trial_data[[trait]] + 1)
  model_log <- lmer(as.formula(paste0(trait, "_Log ~ Genotype + (1 | Replication/Block)")), data = trial_data, na.action = na.exclude)
  blues_matrix[[trait]] <- suppressMessages(as.data.frame(emmeans(model_log, "Genotype")))$emmean 
}

for (trait in traits_to_sqrt) {
  trial_data[[paste0(trait, "_Sqrt")]] <- sqrt(trial_data[[trait]])
  model_sqrt <- lmer(as.formula(paste0(trait, "_Sqrt ~ Genotype + (1 | Replication/Block)")), data = trial_data, na.action = na.exclude)
  blues_matrix[[trait]] <- suppressMessages(as.data.frame(emmeans(model_sqrt, "Genotype")))$emmean 
}

# F. Append Metadata & Save Matrix
blues_matrix$Location <- "Dharwad"        # UPDATE THIS LINE FOR EACH RUN
blues_matrix$Sowing <- "Normal"           # UPDATE THIS LINE FOR EACH RUN
blues_matrix$Environment <- "Dharwad_Normal" # UPDATE THIS LINE FOR EACH RUN

# SAVE OBJECT - UPDATE NAME FOR EACH RUN (blues_D_E, blues_D_N, blues_DH_E, blues_DH_N)
blues_DH_N <- blues_matrix 

# [RETURN TO SECTION 2 AND REPEAT UNTIL ALL 4 ENVIRONMENTS ARE COMPLETED]


# ==========================================
# 3. GxE VARIANCE & HERITABILITY ####
# ==========================================
# Run this section ONLY after all 4 matrices exist in the environment

# Stack the environments
master_gxe <- rbind(blues_D_E, blues_D_N, blues_DH_E, blues_DH_N)
master_gxe$Genotype <- as.factor(master_gxe$Genotype)
master_gxe$Location <- as.factor(master_gxe$Location)
master_gxe$Sowing <- as.factor(master_gxe$Sowing)
master_gxe$Environment <- as.factor(master_gxe$Environment)

meta_cols <- c("Genotype", "Location", "Sowing", "Environment")
master_trait_cols <- setdiff(colnames(master_gxe), meta_cols)
variance_results <- data.frame()

for (trait in master_trait_cols) {
  formula_gxe <- paste(trait, "~ (1|Genotype) + (1|Location) + (1|Sowing) + 
                                  (1|Genotype:Location) + (1|Genotype:Sowing) + 
                                  (1|Location:Sowing)")
  model_gxe <- lmer(as.formula(formula_gxe), data = master_gxe, na.action = na.exclude)
  var_comps <- as.data.frame(VarCorr(model_gxe))
  
  trait_variance <- data.frame(
    Trait = trait,
    Var_Genotype = var_comps$vcov[var_comps$grp == "Genotype"],
    Var_Location = var_comps$vcov[var_comps$grp == "Location"],
    Var_Sowing = var_comps$vcov[var_comps$grp == "Sowing"],
    Var_G_L = var_comps$vcov[var_comps$grp == "Genotype:Location"],
    Var_G_S = var_comps$vcov[var_comps$grp == "Genotype:Sowing"],
    Var_L_S = var_comps$vcov[var_comps$grp == "Location:Sowing"],
    Var_Residual = var_comps$vcov[var_comps$grp == "Residual"]
  )
  variance_results <- rbind(variance_results, trait_variance)
}

# Calculate Broad-Sense Heritability (H^2)
l <- 2 # 2 Locations
s <- 2 # 2 Sowing Dates
variance_results$Heritability <- with(variance_results, 
                                      Var_Genotype / (Var_Genotype + (Var_G_L / l) + (Var_G_S / s) + (Var_Residual / (l * s)))
)
variance_results$Heritability_Percent <- round(variance_results$Heritability * 100, 2)


# ==========================================
# 4. GxE VISUALIZATIONS ####
# ==========================================

# A. Environmental Phenotypic Distributions (Violin Plots)
for (trait in master_trait_cols) {
  plot_name <- paste0("Distribution_", trait, ".png")
  p <- ggplot(master_gxe, aes(x = Environment, y = .data[[trait]], fill = Sowing)) +
    geom_violin(alpha = 0.5, trim = FALSE) +
    geom_boxplot(width = 0.2, fill = "white", color = "black", alpha = 0.8) +
    scale_fill_manual(values = c("Normal" = "#2c7bb6", "Early" = "#d7191c")) +
    theme_minimal() +
    labs(title = paste("Phenotypic Distribution of", trait), x = "Environment", y = trait) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"))
  ggsave(plot_name, plot = p, width = 8, height = 6, dpi = 300)
}

# B. Variance Partitioning Stacked Bar Chart
var_percent <- variance_results
var_cols <- c("Var_Genotype", "Var_Location", "Var_Sowing", "Var_G_L", "Var_G_S", "Var_L_S", "Var_Residual")
var_percent$Total_Var <- rowSums(var_percent[, var_cols])
var_percent[, var_cols] <- (var_percent[, var_cols] / var_percent$Total_Var) * 100

var_long <- pivot_longer(var_percent, cols = all_of(var_cols), names_to = "Variance_Component", values_to = "Percentage")
var_long$Variance_Component <- gsub("Var_", "", var_long$Variance_Component)

var_plot <- ggplot(var_long, aes(x = Trait, y = Percentage, fill = Variance_Component)) +
  geom_bar(stat = "identity", position = "stack", color = "black", linewidth = 0.2) +
  scale_fill_brewer(palette = "Set3") +
  theme_classic() +
  labs(title = "Phenotypic Variance Partitioning (GxE)", x = "Traits", y = "Variance (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"))
ggsave("Variance_Partitioning_Summary.png", plot = var_plot, width = 10, height = 6, dpi = 300)


# ==========================================
# 5. MOLECULAR PANEL SELECTION & STRESS INDICES ####
# ==========================================

# Calculate Means and Indices (Using PT as the primary classification trait)
pt_selection <- master_gxe %>%
  group_by(Genotype, Sowing) %>%
  summarise(Mean_PT = mean(PT, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = Sowing, values_from = Mean_PT)

pop_mean_normal <- mean(pt_selection$Normal, na.rm = TRUE)
pop_mean_early <- mean(pt_selection$Early, na.rm = TRUE)
D <- 1 - (pop_mean_early / pop_mean_normal)

pt_selection <- pt_selection %>%
  mutate(
    Category = case_when(
      Normal >= pop_mean_normal & Early >= pop_mean_early ~ "Stable High",
      Normal >= pop_mean_normal & Early < pop_mean_early ~ "Contrasting Susceptible",
      Normal < pop_mean_normal & Early >= pop_mean_early ~ "Contrasting Tolerant",
      Normal < pop_mean_normal & Early < pop_mean_early ~ "Stable Low"
    ),
    STI = (Normal * Early) / (pop_mean_normal^2),
    SSI = (1 - (Early / Normal)) / D
  ) %>%
  arrange(desc(STI))

# Export Panel Data
write.csv(pt_selection, "Molecular_Panel_with_Indices_PT.csv", row.names = FALSE)

# Generate Four-Quadrant Selection Plot
top_candidates <- pt_selection %>% group_by(Category) %>% top_n(3, wt = abs(Normal - Early))

pt_plot <- ggplot(pt_selection, aes(x = Normal, y = Early, color = Category)) +
  geom_point(alpha = 0.6, size = 3) +
  geom_vline(xintercept = pop_mean_normal, linetype = "dashed", color = "black") +
  geom_hline(yintercept = pop_mean_early, linetype = "dashed", color = "black") +
  geom_text_repel(data = top_candidates, aes(label = Genotype), size = 4, fontface = "bold", show.legend = FALSE) +
  scale_color_manual(values = c("Stable High" = "forestgreen", "Contrasting Susceptible" = "firebrick", 
                                "Contrasting Tolerant" = "darkorange", "Stable Low" = "grey50")) +
  theme_bw() +
  labs(title = "PT: Genotype Selection Matrix", x = "Normal PT", y = "Early Sowing PT")
ggsave("PT_Molecular_Selection_Matrix.png", plot = pt_plot, width = 9, height = 7, dpi = 300)

# Generate GGE Biplot for Stability Analysis
gge_model_pt <- gge(master_gxe, env = Environment, gen = Genotype, resp = PT)
png("GGE_Mean_vs_Stability_PT.png", width = 2400, height = 2000, res = 300)
plot(gge_model_pt, type = 2, title = "GGE Biplot: Mean Performance vs. Stability (PT)")
dev.off()


# ==========================================
# 6. METEOROLOGICAL DATA & BIOPHYSICS ####
# ==========================================

# A. Import Master Weather (Format: Date, Location, Tmax, Tmin, Rainfall, Sunshine_hrs, RH_max, RH_min)
# weather_raw <- read_excel("Weather_Master.xlsx")
# weather_raw$Date <- as.Date(weather_raw$Date)
# weather_raw$Location <- as.factor(weather_raw$Location)

# Define Sowing Milestones (UPDATE DATES)
sowing_dates <- data.frame(
  Location = c("Delhi", "Delhi", "Dharwad", "Dharwad"),
  Sowing = c("Early", "Normal", "Early", "Normal"),
  Sowing_Date = as.Date(c("2025-10-25", "2025-11-15", "2025-10-20", "2025-11-10"))
)

# B. Absolute Calendar Milestone Plot
timeline_plot <- ggplot(weather_raw, aes(x = Date, y = Tmax, color = Location)) +
  geom_line(linewidth = 1, alpha = 0.8) +
  geom_hline(yintercept = 32, linetype = "dashed", color = "darkred", linewidth = 0.8) +
  geom_vline(data = sowing_dates, aes(xintercept = as.numeric(Sowing_Date), color = Location, linetype = Sowing), linewidth = 1) +
  scale_color_manual(values = c("Delhi" = "#2c7bb6", "Dharwad" = "#d7191c")) +
  theme_classic() +
  labs(title = "Thermal Regimes and Trial Milestones", x = "Calendar Date", y = "Maximum Temperature (°C)")
ggsave("Calendar_Milestone_Weather.png", plot = timeline_plot, width = 10, height = 5, dpi = 300)

# C. Calculate Biophysical Indices (VPD, GDD, HTU)
weather_trials <- weather_raw %>%
  inner_join(sowing_dates, by = "Location") %>%
  mutate(DAS = as.numeric(Date - Sowing_Date)) %>%
  filter(DAS >= 0 & DAS <= 140) %>%
  mutate(
    # Vapor Pressure Deficit (VPD)
    es = (0.6108 * exp((17.27 * Tmax) / (Tmax + 237.3)) + 0.6108 * exp((17.27 * Tmin) / (Tmin + 237.3))) / 2,
    ea = (0.6108 * exp((17.27 * Tmin) / (Tmin + 237.3)) * (RH_max / 100) + 
            0.6108 * exp((17.27 * Tmax) / (Tmax + 237.3)) * (RH_min / 100)) / 2,
    VPD = es - ea,
    
    # Heliothermal Units
    Mean_Temp = (Tmax + Tmin) / 2,
    GDD = ifelse(Mean_Temp > 5, Mean_Temp - 5, 0),
    HTU = GDD * Sunshine_hrs
  )


# ==========================================
# 7. NOVEL TILLER STABILITY INDEX (TSI) ####
# ==========================================

trait_summary <- master_gxe %>%
  group_by(Genotype, Sowing) %>%
  summarise(Mean_PT = mean(PT, na.rm = TRUE), Mean_Yield = mean(GrainYld, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = Sowing, values_from = c(Mean_PT, Mean_Yield))

pop_PT_early <- mean(trait_summary$Mean_PT_Early, na.rm = TRUE)
pop_Yield_normal <- mean(trait_summary$Mean_Yield_Normal, na.rm = TRUE)

# Calculate Indices
index_validation <- trait_summary %>%
  mutate(
    TSI = (Mean_PT_Early / Mean_PT_Normal) * (Mean_PT_Early / pop_PT_early),
    Yield_STI = (Mean_Yield_Normal * Mean_Yield_Early) / (pop_Yield_normal^2)
  )

# Validation Regression Plot
validation_plot <- ggplot(index_validation, aes(x = TSI, y = Yield_STI)) +
  geom_point(alpha = 0.6, color = "darkblue", size = 2) +
  geom_smooth(method = "lm", color = "red", fill = "lightpink", alpha = 0.5) +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), 
           label.x = min(index_validation$TSI, na.rm = TRUE), 
           label.y = max(index_validation$Yield_STI, na.rm = TRUE), size = 5, fontface = "bold") +
  theme_bw() +
  labs(title = "Validation of the Novel Tiller Stability Index (TSI)", x = "TSI", y = "Conventional Yield STI")
ggsave("TSI_Validation_Regression.png", plot = validation_plot, width = 7, height = 6, dpi = 300)