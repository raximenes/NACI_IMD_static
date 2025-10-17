# 02_globals.R
# Global settings and constants

n_cohort  <- 100L
age_start <- 11L
age_end   <- 100L
n_cycles  <- as.integer(age_end - age_start)  # 89

dr_cost         <- 0.015
dr_qaly         <- 0.015
v_discount_cost <- 1 / (1 + dr_cost) ^ (0:n_cycles)
v_discount_qaly <- 1 / (1 + dr_qaly) ^ (0:n_cycles)

perspective <- "healthcare"  # or "societal"
v_wtp <- seq(0, 200000, by = 10000)

v_strats <- c("MenC", "MenACWY", "MenACWY_MenB", "MenABCWY")

v_states <- c(
  "Healthy",
  "SeroB_Infect", "SeroC_Infect", "SeroW_Infect", "SeroY_Infect",
  "Scarring", "Single_Amput", "Multiple_Amput", "Neuro_Disability",
  "Hearing_Loss", "Renal_Failure", "Seizure", "Paralysis",
  "Dead_IMD", "Background_Mortality"
)
n_states <- length(v_states)

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
