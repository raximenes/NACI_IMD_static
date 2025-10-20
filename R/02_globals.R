# 02_globals.R
# Global settings and constants
# Note: Some values can be overridden by run_analysis.R

# Cohort size (can be overridden by run_analysis.R)
if (!exists("n_cohort")) {
  n_cohort <- 100L
}

# Age range
age_start <- 11L
age_end   <- 100L
n_cycles  <- as.integer(age_end - age_start)  # 89

# Discount rates
dr_cost         <- 0.015
dr_qaly         <- 0.015
v_discount_cost <- 1 / (1 + dr_cost) ^ (0:n_cycles)
v_discount_qaly <- 1 / (1 + dr_qaly) ^ (0:n_cycles)

# Analysis perspective (can be overridden by run_analysis.R)
if (!exists("perspective")) {
  perspective <- "healthcare"  # or "societal"
}

# WTP range for CEAC
v_wtp <- seq(0, 200000, by = 10000)

# Strategies
v_strats <- c("MenC", "MenACWY", "MenACWY_MenB", "MenABCWY")

# Health states
v_states <- c(
  "Healthy",
  "SeroB_Infect", "SeroC_Infect", "SeroW_Infect", "SeroY_Infect",
  "Scarring", "Single_Amput", "Multiple_Amput", "Neuro_Disability",
  "Hearing_Loss", "Renal_Failure", "Seizure", "Paralysis",
  "Dead_IMD", "Background_Mortality"
)
n_states <- length(v_states)

# Societal costs for productivity loss (only applied if perspective = "societal")
c_prod_societal <- c(
  Scarring         =  5000,
  Single_Amput     = 12000,
  Multiple_Amput   = 20000,
  Neuro_Disability = 20000,
  Hearing_Loss     =  8000,
  Renal_Failure    = 10000,
  Seizure          = 10000,
  Paralysis        = 20000
)

# Log current settings
log_info("Global parameters loaded:")
log_info(paste("  - Cohort size:", n_cohort))
log_info(paste("  - Perspective:", perspective))
log_info(paste("  - Age range:", age_start, "to", age_end, "(", n_cycles, "cycles)"))
log_info(paste("  - Discount rates: Cost =", dr_cost, ", QALY =", dr_qaly))