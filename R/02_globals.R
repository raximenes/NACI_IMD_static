################################################################################
# 02_globals.R
# Description: Global settings, constants, and output folder organization
#
# AUTHORSHIP:
# ===========
# Raphael Ximenes (raphael.ximenes@phac-aspc.gc.ca)
# Gebremedhin Gebretekle (Gebremedhin.Gebretekle@phac-aspc.gc.ca)
#
# Modelling & Health Economics Unit
# Centre for Immunization Surveillance and Programs (CISP)
# Public Health Agency of Canada (PHAC)
#
# This script defines:
# - Simulation configuration (cohort size, PSA iterations, random seed)
# - Economic parameters (discount rates)
# - Model structure (health states, strategies, WTP thresholds)
# - Folder structure for perspective-specific outputs
#
# CRITICAL: Edit SIMULATION CONFIGURATION section to adjust analysis parameters
#
################################################################################

# ============================================================
# SIMULATION CONFIGURATION (CENTRALIZED)
# ============================================================

# Cohort size for all analyses
# CRITICAL: Larger cohorts reduce stochastic noise but increase computation time
# TYPICAL: 1,000,000 for final runs; 100,000 for testing
if (!exists("n_cohort")) {
  n_cohort <- 100000L
}

# PSA configuration
# TYPICAL: 2,500-10,000 for final analysis; 100-500 for testing
if (!exists("N_SIMULATIONS")) {
  N_SIMULATIONS <- 2500L
}

ANALYSIS_PERSPECTIVE <- "both"  # "healthcare", "societal", or "both"


# Random seed for reproducibility
if (!exists("RANDOM_SEED")) {
  RANDOM_SEED <- 2025L
}

# ============================================================
# AGE AND CYCLE CONFIGURATION
# ============================================================

# Model starting and ending ages
age_start <- 11L
age_end   <- 100L

# Number of cycles (years in model)
# Calculation: age_end - age_start = 100 - 11 = 89 years
n_cycles  <- as.integer(age_end - age_start)

# ============================================================
# DISCOUNT RATES (Economic parameters)
# ============================================================

# Annual discount rates for cost and QALY calculations
# TYPICAL: 1.5% (0.015) for Canadian analyses
# PURPOSE: Reflects time value of money/health
dr_cost         <- 0.015
dr_qaly         <- 0.015

# Discount vectors pre-computed for efficiency
# Formula: 1 / (1 + discount_rate)^cycle
v_discount_cost <- 1 / (1 + dr_cost) ^ (0:n_cycles)
v_discount_qaly <- 1 / (1 + dr_qaly) ^ (0:n_cycles)

# ============================================================
# ANALYSIS PERSPECTIVE
# ============================================================

# Default analysis perspective
# OPTIONS: "healthcare" (direct costs only) or "societal" (includes productivity loss)
if (!exists("perspective")) {
  perspective <- "healthcare"
}

# ============================================================
# WILLINGNESS-TO-PAY (WTP) RANGE FOR CEAC
# ============================================================

# Range of cost-effectiveness thresholds for CEAC plots
# TYPICAL: $0 - $200,000 CAD per QALY for Canada
# CALCULATION: seq(0, 200000, by = 10000) = 21 thresholds
v_wtp <- seq(0, 200000, by = 10000)

# ============================================================
# INTERVENTION STRATEGIES
# ============================================================

# Vaccination strategies under comparison
# B/C/W/Y = serogroups covered; MenC = existing baseline
v_strats <- c("MenC", "MenACWY", "MenACWY_MenB", "MenABCWY")

# ============================================================
# HEALTH STATES
# ============================================================

# Model health states and disease sequelae
# ORGANIZED: Healthy → Infections → Sequelae → Death
v_states <- c(
  # Healthy state
  "Healthy",
  
  # Infection states (active disease by serogroup)
  "SeroB_Infect", "SeroC_Infect", "SeroW_Infect", "SeroY_Infect",
  
  # Sequelae (long-term complications)
  "Scarring", "Single_Amput", "Multiple_Amput", "Neuro_Disability",
  "Hearing_Loss", "Renal_Failure", "Seizure", "Paralysis",
  
  # Death states
  "Dead_IMD", "Background_Mortality"
)

# Total number of states
n_states <- length(v_states)

# ============================================================
# SOCIETAL PRODUCTIVITY COSTS (if perspective = "societal")
# ============================================================

# Annual productivity losses by sequela type
# SOURCE: Cost-effectiveness analysis inputs
# UNIT: CAD per year for affected individuals
# c_prod_societal <- c(
#   Scarring         =  5000,
#   Single_Amput     = 12000,
#   Multiple_Amput   = 20000,
#   Neuro_Disability = 20000,
#   Hearing_Loss     =  8000,
#   Renal_Failure    = 10000,
#   Seizure          = 10000,
#   Paralysis        = 20000
# )

# ============================================================
# OUTPUT FOLDER ORGANIZATION BY ANALYSIS TYPE
# ============================================================

# STRUCTURE:
# outputs/
# ├── figures/
# │   ├── healthcare/
# │   │   ├── Base_Case/
# │   │   ├── PSA/
# │   │   ├── OWSA/
# │   │   └── KEY_PLOTS/
# │   └── societal/
# │       ├── Base_Case/
# │       ├── PSA/
# │       ├── OWSA/
# │       └── KEY_PLOTS/
# └── tables/
#     ├── healthcare/
#     └── societal/

# Base folders (organized by perspective FIRST)
OUT_FIG_BASE_HC <- file.path(OUT_FIG, "healthcare")
OUT_FIG_BASE_SOC <- file.path(OUT_FIG, "societal")

# Create base perspective folders
dir.create(OUT_FIG_BASE_HC, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_FIG_BASE_SOC, recursive = TRUE, showWarnings = FALSE)

# Analysis type subfolders (updated when perspective changes)
# These variables are global and shared across all scripts
OUT_FIG_PSA  <- file.path(OUT_FIG, perspective, "PSA")        # For PSA parameter distributions ONLY
OUT_FIG_OWSA <- file.path(OUT_FIG, perspective, "OWSA")       # For OWSA tornado plots
OUT_FIG_KEY  <- file.path(OUT_FIG, perspective, "KEY_PLOTS")  # For CEAC, EVPI, and key PSA outputs
OUT_FIG_DET  <- file.path(OUT_FIG, perspective, "Base_Case")  # For base case ICER frontier
#OUT_TAB_PERSP <- file.path(OUT_TAB, perspective)              # For all tables

# Create perspective-specific analysis folders
dir.create(OUT_FIG_PSA, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_FIG_OWSA, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_FIG_KEY, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_FIG_DET, recursive = TRUE, showWarnings = FALSE)
#dir.create(OUT_TAB_PERSP, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# HELPER FUNCTION: Update folders when perspective changes
# ============================================================

# Call this function after changing perspective to update all paths
update_output_folders <- function(new_perspective) {
  OUT_FIG_PSA  <<- file.path(OUT_FIG, new_perspective, "PSA")
  OUT_FIG_OWSA <<- file.path(OUT_FIG, new_perspective, "OWSA")
  OUT_FIG_KEY  <<- file.path(OUT_FIG, new_perspective, "KEY_PLOTS")
  OUT_FIG_DET  <<- file.path(OUT_FIG, new_perspective, "Base_Case")
#  OUT_TAB_PERSP <<- file.path(OUT_TAB, new_perspective)
  
  # Ensure directories exist
  dir.create(OUT_FIG_PSA, recursive = TRUE, showWarnings = FALSE)
  dir.create(OUT_FIG_OWSA, recursive = TRUE, showWarnings = FALSE)
  dir.create(OUT_FIG_KEY, recursive = TRUE, showWarnings = FALSE)
  dir.create(OUT_FIG_DET, recursive = TRUE, showWarnings = FALSE)
#  dir.create(OUT_TAB_PERSP, recursive = TRUE, showWarnings = FALSE)
  
  log_info(paste("Output folders updated for perspective:", new_perspective))
}

# ============================================================
# LOG CURRENT GLOBAL SETTINGS
# ============================================================

log_info("Global parameters loaded:")
log_info(paste("  - Cohort size:", format(n_cohort, big.mark = ",")))
log_info(paste("  - PSA simulations:", format(N_SIMULATIONS, big.mark = ",")))
log_info(paste("  - Random seed:", RANDOM_SEED))
log_info(paste("  - Perspective:", perspective))
log_info(paste("  - Age range:", age_start, "to", age_end, "(", n_cycles, "cycles)"))
log_info(paste("  - Discount rates: Cost =", dr_cost, ", QALY =", dr_qaly))
log_info(paste("  - Output folders organized by: Base_Case, PSA, OWSA, KEY_PLOTS"))
log_info(paste("  - Perspective-specific subfolders created"))