# ==============================================================================
# MASTER PIPELINE: WHEAT GxE, MICRO-PHENOLOGY, AND CHRONIC VEGETATIVE STRESS
# Repository: [Insert GitHub Repo Name]
# Author: Apoorva Ashu, ICAR-IARI
# Date: September 2026
# Description: Evaluates 217 wheat genotypes under chronic vegetative thermal 
#              stress across varying sowing dates (S1, S2) and locations 
#              (Delhi, Dharwad). Extracts BLUEs using REML Alpha-Lattice models, 
#              calculates AUTPC, TCE, and localized variance components.
# ==============================================================================

# ==========================================
# 1. SETUP & LIBRARIES ####
# ==========================================
# Ensure renv is deactivated if encountering vctrs/dplyr conflicts
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
# INSTRUCTION: Execute this section 4 times independently. 
# Update the target sheet and output object name for each run:
# Delhi_S1, Delhi_S2, Dharwad_S1, Dharwad_S2

# A. Import Raw Data
trial_data <- read_excel("Data_Dharwad.xlsx", sheet = "S2_Final") # UPDATE PER RUN

# B. Enforce Factor Structure (CRITICAL for REML)
trial_data$Genotype <- as.factor(trial_data$Genotype)
trial_data$Replication <- as.factor(trial_data$Replication)
trial_data$Block <- as.factor(trial_data$Block)

# C. Initialize Matrices
design_cols <- c("Sr_No", "Genotype", "Replication", "Block")
trait_cols <- setdiff(colnames(trial_data), design_cols)
blues_matrix <- data.frame(Genotype = unique(trial_data$Genotype))

# D. Master Loop: Spatial Normalization & BLUEs Extraction
for (trait in trait_cols) {
  # Define Spatial Models
  formula_alpha <- paste(trait, "~ Genotype + (1 | Replication/Block)")
  formula_rcbd  <- paste(trait, "~ Genotype + (1 | Replication)")
  
  # Try Alpha Lattice; Fallback to RCBD if spatial variance collapses (zero variance)
  model <- tryCatch({
    lmer(as.formula(formula_alpha), data = trial_data, na.action = na.exclude)
  }, error = function(e) {
    lmer(as.formula(formula_rcbd), data = trial_data, na.action = na.exclude)
  })
  
  # Extract Best Linear Unbiased Estimators (BLUEs)
  em_out <- suppressMessages(as.data.frame(emmeans(model, "Genotype")))
  trait_blues <- em_out[, c("Genotype", "emmean")]
  colnames(trait_blues)[2] <- trait
  
  # Append to specific environment matrix
  blues_matrix <- merge(blues_matrix, trait_blues, by = "Genotype", all.x = TRUE)
}

# E. Append Metadata & Save Environment Object
blues_matrix$Location <- "Dharwad"        # UPDATE PER RUN
blues_matrix$Sowing <- "S2"               # UPDATE PER RUN
blues_matrix$Environment <- "Dharwad_S2"  # UPDATE PER RUN

# SAVE OBJECT - UPDATE NAME FOR EACH RUN (blues_Delhi_S1, blues_Delhi_S2, blues_Dharwad_S1, blues_Dharwad_S2)
blues_Dharwad_S2 <- blues_matrix 

# [STOP & REPEAT SECTION 2 UNTIL ALL 4 ENVIRONMENTS ARE PROCESSED]


# ==========================================
# 3. GLOBAL GxE MERGE & LOCALIZED VARIANCE ####
# ==========================================
# Combine the 4 independent environments into the Master Matrix
master_gxe <- rbind(blues_Delhi_S1, blues_Delhi_S2, blues_Dharwad_S1, blues_Dharwad_S2)
master_gxe$Genotype <- as.factor(master_gxe$Genotype)
master_gxe$Location <- as.factor(master_gxe$Location)
master_gxe$Sowing <- as.factor(master_gxe$Sowing)

meta_cols <- c("Genotype", "Location", "Sowing", "Environment")
master_trait_cols <- setdiff(colnames(master_gxe), meta_cols)

# Function: Localized Variance Partitioning (Two-Way Mixed Model)
run_local_variance <- function(local_data, location_name) {
  local_variance_results <- data.frame()
  for (trait in master_trait_cols) {
    formula_local <- paste(trait, "~ (1|Genotype) + (1|Sowing) + (1|Genotype:Sowing)")
    model_local <- tryCatch(lmer(as.formula(formula_local), data = local_data, na.action = na.exclude), 
                            error = function(e) return(NULL))
    
    if(!is.null(model_local)) {
      var_comps <- as.data.frame(VarCorr(model_local))
      trait_var <- data.frame(
        Trait = trait,
        Var_Genotype = var_comps$vcov[var_comps$grp == "Genotype"],
        Var_Sowing = var_comps$vcov[var_comps$grp == "Sowing"],
        Var_G_S = var_comps$vcov[var_comps$grp == "Genotype:Sowing"],
        Var_Residual = var_comps$vcov[var_comps$grp == "Residual"]
      )
      local_variance_results <- rbind(local_variance_results, trait_var)
    }
  }
  # Calculate Genotypic Repeatability for Single-Year Trials
  local_variance_results$Genotypic_Repeatability <- with(local_variance_results, 
    Var_Genotype / (Var_Genotype + (Var_G_S / 2) + (Var_Residual / 2)))
  
  write.csv(local_variance_results, paste0(location_name, "_Local_Variance_Analysis.csv"), row.names = FALSE)
  return(local_variance_results)
}

# Execute Local Models to establish regional baseline variance
delhi_variance <- run_local_variance(subset(master_gxe, Location == "Delhi"), "Delhi")
dharwad_variance <- run_local_variance(subset(master_gxe, Location == "Dharwad"), "Dharwad")


# ==========================================
# 4. MICRO-PHENOLOGY: AUTPC & MATURITY CATEGORIZATION ####
# ==========================================
# Evaluates early vigor uncoupled from absolute maximum capacity
# NOTE: Replace T1, T2, T3, Tmax_trait with your exact column names
master_gxe <- master_gxe %>%
  mutate(
    # Interval lengths based on empirical meteorological observations
    Int_Days = ifelse(Location == "Dharwad", 10, 15),
    
    # Area Under the Tillering Progress Curve (AUTPC) via Trapezoidal Rule
    AUTPC = (0.5 * (0 + T1) * Int_Days) +         
            (0.5 * (T1 + T2) * Int_Days) +        
            (0.5 * (T2 + T3) * Int_Days) +        
            (0.5 * (T3 + Tmax_trait) * Int_Days), 
    
    # Standardize to compare across treatments with varying final capacities
    Relative_AUTPC = AUTPC / Tmax_trait,
    
    # Identify the specific interval of Peak Emergence Velocity
    Rate_E_T1 = T1 / Int_Days,
    Rate_T1_T2 = (T2 - T1) / Int_Days,
    Rate_T2_T3 = (T3 - T2) / Int_Days,
    
    Peak_Phase = case_when(
      Rate_E_T1 > Rate_T1_T2 & Rate_E_T1 > Rate_T2_T3 ~ "Very Early (E-T1)",
      Rate_T1_T2 >= Rate_E_T1 & Rate_T1_T2 > Rate_T2_T3 ~ "Early Vigor (T1-T2)",
      Rate_T2_T3 >= Rate_T1_T2 & Rate_T2_T3 >= Rate_E_T1 ~ "Late Vigor (T2-T3)",
      TRUE ~ "Uniform"
    )
  )

autpc_yield_plot <- ggplot(master_gxe, aes(x = Relative_AUTPC, y = GrainYld)) +
  geom_point(aes(color = Peak_Phase), alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", color = "black", linewidth = 1) +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x.npc = "left", label.y.npc = "bottom") +
  facet_grid(Location ~ Sowing) +
  scale_color_brewer(palette = "Set1") +
  theme_bw() +
  labs(title = "Dependence of Grain Yield on Early Tiller Establishment (AUTPC)",
       x = "Relative Area Under Tillering Progress Curve (AUTPC)", y = "Grain Yield")
ggsave("AUTPC_vs_Yield.png", plot = autpc_yield_plot, width = 11, height = 8, dpi = 300)


# ==========================================
# 5. TILLER DYNAMICS: ASYMPTOTES & PHASE TRANSITIONS ####
# ==========================================
# Safely calculates the Inverse Emergence Metric (Days/Tiller) 
tiller_dynamics <- master_gxe %>%
  select(Genotype, Location, Sowing, T1, T2) %>% # Evaluates specific micro-fractions (e.g., T1 to T2)
  mutate(
    Interval_Days = ifelse(Location == "Dharwad", 10, 15), 
    Delta_Tillers = T2 - T1,
    
    # Tag Biological State to handle mathematical asymptotes
    Tiller_Phase = case_when(
      Delta_Tillers > 0.05 ~ "Active Emergence",
      Delta_Tillers < -0.05 ~ "Active Abortion",
      TRUE ~ "Stagnation"
    ),
    
    # Calculate Conditional Metric
    Days_Per_Event = case_when(
      Tiller_Phase == "Active Emergence" ~ Interval_Days / Delta_Tillers,
      Tiller_Phase == "Active Abortion" ~ Interval_Days / abs(Delta_Tillers),
      Tiller_Phase == "Stagnation" ~ NA_real_ # Prevents Divide-by-Zero
    )
  )
write.csv(tiller_dynamics, "Tiller_Dynamics_Phase_Transitions.csv", row.names = FALSE)


# ==========================================
# 6. METEOROLOGY & PHENOLOGICAL COMPRESSION ####
# ==========================================
weather_raw <- read_excel("Weather_Master.xlsx")
weather_raw$Date <- as.Date(weather_raw$Date)

# Macro-Boundaries (TI and PI)
phase_dates <- data.frame(
  Location = c("Delhi", "Delhi", "Dharwad", "Dharwad"),
  Sowing = c("S1", "S2", "S1", "S2"),
  Sowing_Date = as.Date(c("2025-11-03", "2025-11-21", "2025-11-10", "2025-11-20")),
  TI_Date = as.Date(c("2025-12-10", "2025-12-30", "2025-12-05", "2025-12-25")), 
  PI_Date = as.Date(c("2026-02-05", "2026-02-20", "2026-01-14", "2026-01-29"))  
) %>%
  mutate(Emergence_Date = case_when(
    Location == "Delhi" & Sowing == "S1" ~ Sowing_Date + 9,
    Location == "Delhi" & Sowing == "S2" ~ Sowing_Date + 13,
    Location == "Dharwad" & Sowing == "S1" ~ Sowing_Date + 5,
    Location == "Dharwad" & Sowing == "S2" ~ Sowing_Date + 7
  ))

# Micro-Milestones to extract specific Tmax values for plotting
micro_milestones_das <- data.frame(
  Location = rep(c("Delhi", "Delhi", "Dharwad", "Dharwad"), each = 4),
  Sowing = rep(c("S1", "S2", "S1", "S2"), each = 4),
  Stage = rep(c("T1", "T2", "T3", "Tmax"), 4),
  DAS = c(35, 50, 65, 80, 40, 55, 70, 85, 15, 25, 35, 45, 19, 29, 39, 49)
)

micro_milestones_dates <- micro_milestones_das %>%
  left_join(phase_dates %>% select(Location, Sowing, Sowing_Date), by = c("Location", "Sowing")) %>%
  mutate(Date = Sowing_Date + DAS) %>%
  left_join(weather_raw %>% select(Date, Location, Tmax), by = c("Date", "Location"))

# Generate Integrated Meteorological Plot
integrated_thermal_plot <- ggplot(weather_raw, aes(x = Date)) +
  geom_ribbon(aes(ymin = Tmin, ymax = Tmax, group = Location), fill = "grey85", alpha = 0.5) +
  geom_line(aes(y = Tmax, color = "Tmax"), linewidth = 0.8) +
  geom_line(aes(y = Tmin, color = "Tmin"), linewidth = 0.8) +
  geom_vline(data = phase_dates, aes(xintercept = as.numeric(Emergence_Date), color = Sowing), linetype = "solid", linewidth = 0.8) +
  geom_vline(data = phase_dates, aes(xintercept = as.numeric(PI_Date), color = Sowing), linetype = "dashed", linewidth = 0.8) +
  geom_point(data = micro_milestones_dates, aes(x = Date, y = Tmax, shape = Stage, fill = Sowing), size = 3, color = "black", stroke = 1) +
  facet_wrap(~ Location, ncol = 1) +
  scale_color_manual(name = "Thermal Boundaries", values = c("Tmax" = "firebrick", "Tmin" = "royalblue", "S1" = "darkorange", "S2" = "seagreen")) +
  scale_fill_manual(name = "Micro-Milestones", values = c("S1" = "darkorange", "S2" = "seagreen")) +
  scale_shape_manual(name = "Tillering Stage", values = c("T1" = 21, "T2" = 22, "T3" = 23, "Tmax" = 24)) +
  theme_bw() +
  labs(title = "Phenological Compression Under Chronic Vegetative Stress", 
       subtitle = "Lines: Emergence to PI. Points: Active tiller differentiation.", x = "Calendar Date", y = "Temperature (°C)")
ggsave("Integrated_Micro_Phenology_Plot.png", plot = integrated_thermal_plot, width = 11, height = 8, dpi = 300)


# ==========================================
# 7. TILLER CONVERSION EFFICIENCY & MOLECULAR SELECTION ####
# ==========================================
# Define Biological Environments based on GDD empirical findings
efficiency_metrics <- master_gxe %>%
  mutate(
    TCE = PT / Tmax_trait,  # Tiller Conversion Efficiency
    Stress_Regime = case_when(
      (Location == "Delhi" & Sowing == "S2") | (Location == "Dharwad" & Sowing == "S1") ~ "Optimal_Combined",
      (Location == "Delhi" & Sowing == "S1") | (Location == "Dharwad" & Sowing == "S2") ~ "Stress_Combined"
    )
  )

# Calculate Molecular Selection Matrix utilizing TCE
combined_selection <- efficiency_metrics %>%
  group_by(Genotype, Stress_Regime) %>%
  summarise(Mean_TCE = mean(TCE, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = Stress_Regime, values_from = Mean_TCE)

pop_opt <- mean(combined_selection$Optimal_Combined, na.rm = TRUE)
pop_stress <- mean(combined_selection$Stress_Combined, na.rm = TRUE)

combined_selection <- combined_selection %>%
  mutate(
    Category = case_when(
      Optimal_Combined >= pop_opt & Stress_Combined >= pop_stress ~ "Stable High Efficiency",
      Optimal_Combined >= pop_opt & Stress_Combined < pop_stress ~ "Contrasting Susceptible",
      Optimal_Combined < pop_opt & Stress_Combined >= pop_stress ~ "Contrasting Tolerant",
      Optimal_Combined < pop_opt & Stress_Combined < pop_stress ~ "Stable Low Efficiency"
    )
  )
write.csv(combined_selection, "Combined_Molecular_Panel_TCE.csv", row.names = FALSE)

print("--- Master Pipeline Execution Complete ---")

```
