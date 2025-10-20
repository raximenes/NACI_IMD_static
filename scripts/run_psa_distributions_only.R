################################################################################
# run_psa_distributions_only.R
# Script to generate and visualize PSA parameter distributions WITHOUT running
# the full model. Useful for validating parameter distributions before 
# running expensive simulations.
################################################################################
# ============================================================
# START TIMER
# ============================================================
start_time <- Sys.time()
# ============================================================
# CONFIGURATION
# ============================================================

N_SIMULATIONS <- 10000L      # Number of PSA samples to generate
RANDOM_SEED <- 2025          # For reproducibility
COHORT_SIZE <- 100000L       # Needed for parameter loading

# ============================================================
# PACKAGE LOADING
# ============================================================

if (!require("pacman", quietly = TRUE)) {
  install.packages("pacman", repos = "https://cloud.r-project.org")
}

library(pacman)
p_load(here, readxl, ggplot2, dplyr, tidyr, purrr, stringr, tibble, readr)
p_load_gh("DARTH-git/darthtools")
p_load(dampack)

message("\n========================================")
message("PSA DISTRIBUTION GENERATION (NO MODEL RUN)")
message("========================================")
message(paste("PSA Simulations:", N_SIMULATIONS))
message(paste("Random Seed:", RANDOM_SEED))
message("========================================\n")

# Set random seed
set.seed(RANDOM_SEED)

# ============================================================
# SOURCE REQUIRED FILES ONLY
# ============================================================

message("[INFO] Loading required model files...")

# Load files in order
required_files <- c(
  "R/01_setup.R",           # Paths and logging
  "R/02_globals.R",         # Global parameters
  "R/03_utils.R",           # Utility functions
  "R/04_inputs.R",          # Parameter loading
  "R/06_psa.R",             # PSA generation
  "R/06b_psa_visualization.R"  # Visualization functions
)

for (file in required_files) {
  if (file.exists(file)) {
    message(paste("  Loading:", basename(file)))
    source(file)
  } else {
    stop(paste("Required file not found:", file))
  }
}

message("[INFO] Files loaded successfully\n")

# ============================================================
# SET PARAMETERS
# ============================================================

# Set perspective (doesn't matter for distributions, but needed for loading)
perspective <- "healthcare"
n_cohort <- as.integer(COHORT_SIZE)

# ============================================================
# LOAD BASE PARAMETERS
# ============================================================

message("\n========================================")
message("LOADING BASE PARAMETERS")
message("========================================\n")

params_bc <- get_base_params()

message("[INFO] Base parameters loaded successfully\n")

# ============================================================
# GENERATE PSA SAMPLES
# ============================================================

message("\n========================================")
message("GENERATING PSA PARAMETER SAMPLES")
message("========================================\n")

psa_samples_output <- generate_psa_samples(params_bc, n_sim = N_SIMULATIONS, seed = RANDOM_SEED)
psa_parameter_samples <- psa_samples_output$samples
psa_parameter_df <- psa_samples_output$psa_df

message("[INFO] PSA parameter samples generated successfully\n")

# ============================================================
# CREATE PSA DISTRIBUTION PLOTS
# ============================================================

message("\n========================================")
message("CREATING PSA DISTRIBUTION PLOTS")
message("========================================\n")

# 1. Key parameters panel
message("[INFO] Creating key parameters panel...")
plot_psa_key_parameters_panel(psa_parameter_df, params_bc)

# 2. Validation plot
message("[INFO] Creating PSA validation plot...")
plot_psa_validation(psa_parameter_df, params_bc)

# 3. Temporal summaries for infection probabilities
message("[INFO] Creating temporal summary plots...")
for (param in c("p_B", "p_C", "p_W", "p_Y")) {
  param_cols <- grep(paste0("^", param, "\\d+$"), names(psa_parameter_df), value = TRUE)
  if (length(param_cols) > 0) {
    plot_psa_temporal_summary(
      psa_parameter_df, params_bc, param,
      plot_title = paste("Serogroup", gsub("p_", "", param), "Infection Probability")
    )
  }
}

# 4. Temporal distributions for selected cycles (ALL INFECTION PARAMETERS)
message("[INFO] Creating temporal distribution plots...")
for (param in c("p_B", "p_C", "p_W", "p_Y")) {
  param_cols <- grep(paste0("^", param, "\\d+$"), names(psa_parameter_df), value = TRUE)
  if (length(param_cols) > 0) {
    plot_psa_temporal_distributions(
      psa_parameter_df, params_bc, param,
      cycles_to_plot = c(1, 10, 20, 40, 60, 89)
    )
  }
}

# 5. Selected scalar distributions
message("[INFO] Creating scalar distribution plots...")
important_scalars <- c(
  "c_Scarring", "c_Single_Amput", "c_Multiple_Amput",
  "u_IMD", "u_Scarring", "u_Single_Amput",
  "p_seq_overall", "c_MenABCWY", "c_MenB"
)

available_scalars <- intersect(important_scalars, names(psa_parameter_df))
if (length(available_scalars) > 0) {
  plot_psa_scalar_distributions(psa_parameter_df, params_bc, available_scalars)
}

message("[INFO] All PSA distribution plots created and saved\n")

# ============================================================
# CREATE SUMMARY STATISTICS TABLE
# ============================================================

message("\n========================================")
message("CREATING PARAMETER SUMMARY TABLE")
message("========================================\n")

# Get scalar parameters
scalar_params <- names(psa_parameter_df)[!grepl("\\d+$", names(psa_parameter_df))]

# Filter to only parameters that exist in params_bc as scalars
valid_params <- scalar_params[sapply(scalar_params, function(p) {
  exists_in_bc <- !is.null(params_bc[[p]])
  is_scalar <- exists_in_bc && length(params_bc[[p]]) == 1
  is_numeric <- is_scalar && is.numeric(params_bc[[p]])
  return(is_numeric)
})]

if (length(valid_params) > 0) {
  psa_summary <- data.frame(
    Parameter = valid_params,
    Deterministic = sapply(valid_params, function(p) params_bc[[p]]),
    PSA_Mean = sapply(valid_params, function(p) mean(psa_parameter_df[[p]], na.rm = TRUE)),
    PSA_SD = sapply(valid_params, function(p) sd(psa_parameter_df[[p]], na.rm = TRUE)),
    PSA_Q025 = sapply(valid_params, function(p) quantile(psa_parameter_df[[p]], 0.025, na.rm = TRUE)),
    PSA_Q975 = sapply(valid_params, function(p) quantile(psa_parameter_df[[p]], 0.975, na.rm = TRUE)),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::mutate(
      Pct_Diff = abs((PSA_Mean - Deterministic) / Deterministic) * 100
    ) %>%
    dplyr::arrange(desc(Pct_Diff))
  
  # Save summary
  readr::write_csv(psa_summary, file.path(OUT_TAB, "psa_parameter_summary.csv"))
  
  # Print top parameters with largest deviations
  message("\n--- TOP 10 PARAMETERS WITH LARGEST PSA MEAN DEVIATION ---")
  print(head(psa_summary[, c("Parameter", "Deterministic", "PSA_Mean", "Pct_Diff")], 10))
  
  message("\n[INFO] PSA parameter summary saved to:", file.path(OUT_TAB, "psa_parameter_summary.csv"))
  
} else {
  message("[WARN] No valid scalar parameters found for summary table")
}

# ============================================================
# TEMPORAL PARAMETERS SUMMARY
# ============================================================

message("\n========================================")
message("TEMPORAL PARAMETERS SUMMARY")
message("========================================\n")

temporal_params <- c("p_B", "p_C", "p_W", "p_Y", "p_bg_mort", 
                     "p_B_DeadIMD", "p_C_DeadIMD", "p_W_DeadIMD", "p_Y_DeadIMD",
                     "c_IMD_infection")

temporal_summary_list <- list()

for (param in temporal_params) {
  param_cols <- grep(paste0("^", param, "\\d+$"), names(psa_parameter_df), value = TRUE)
  
  if (length(param_cols) > 0) {
    # Calculate summary statistics across all cycles
    mean_across_cycles <- mean(sapply(param_cols, function(col) mean(psa_parameter_df[[col]])))
    sd_across_cycles <- mean(sapply(param_cols, function(col) sd(psa_parameter_df[[col]])))
    
    temporal_summary_list[[param]] <- data.frame(
      Parameter = param,
      N_Cycles = length(param_cols),
      Mean_Across_Cycles = mean_across_cycles,
      SD_Across_Cycles = sd_across_cycles,
      stringsAsFactors = FALSE
    )
  }
}

if (length(temporal_summary_list) > 0) {
  temporal_summary <- bind_rows(temporal_summary_list)
  
  message("\n--- TEMPORAL PARAMETERS ---")
  print(temporal_summary)
  
  readr::write_csv(temporal_summary, file.path(OUT_TAB, "psa_temporal_parameters_summary.csv"))
  message("\n[INFO] Temporal parameters summary saved")
}

# ============================================================
# SAVE PSA SAMPLES WITH METADATA
# ============================================================

message("\n========================================")
message("SAVING PSA SAMPLES FOR REUSE")
message("========================================\n")

# Add metadata for validation in run_analysis.R
psa_samples_output$metadata <- list(
  n_simulations = N_SIMULATIONS,
  cohort_size = COHORT_SIZE,
  random_seed = RANDOM_SEED,
  generation_date = Sys.time(),
  r_version = R.version.string
)

# Save the PSA samples object
psa_samples_file <- file.path(OUT_TAB, "psa_samples.rds")
saveRDS(psa_samples_output, psa_samples_file)

message("[INFO] PSA samples saved with metadata to:", psa_samples_file)
message("[INFO] Metadata included:")
message(paste("  - N simulations:", N_SIMULATIONS))
message(paste("  - Cohort size:", COHORT_SIZE))
message(paste("  - Random seed:", RANDOM_SEED))
message(paste("  - Generation date:", Sys.time()))
message("\n[INFO] To load in run_analysis.R: set REUSE_PSA_SAMPLES = TRUE")

# ============================================================
# FINAL SUMMARY
# ============================================================

message("\n========================================")
message("PSA DISTRIBUTION GENERATION COMPLETE")
message("========================================")
message(paste("Total PSA samples generated:", N_SIMULATIONS))
message(paste("Total parameters varying:", ncol(psa_parameter_df)))
message(paste("Scalar parameters:", length(valid_params)))
message(paste("Temporal parameters:", length(temporal_summary_list)))
message(paste("\nPlots saved to:", OUT_FIG))
message(paste("Tables saved to:", OUT_TAB))
message("\n========================================")
message("\nPLOTS GENERATED:")
message("1. psa_key_parameters_panel.png - Overview of key parameter distributions")
message("2. psa_validation.png - Deterministic vs PSA mean comparison")
message("3. psa_temporal_summary_[param].png - Time series with 95% CI (4 plots)")
message("4. psa_dist_[param]_cycle[N].png - Distribution histograms by cycle (24 plots)")
message("5. psa_dist_[param].png - Scalar parameter distributions (9 plots)")
message("\nTOTAL: ~38 plots generated")
message("\n========================================")
message("\nRECOMMENDATIONS:")
message("1. Check psa_validation.png - points should be on diagonal line")
message("2. Review temporal summary plots - PSA mean should track deterministic")
message("3. Examine distribution histograms - check for unexpected shapes")
message("4. If distributions look good, run full analysis with run_analysis.R")
message("========================================\n")
message("\n[INFO] Script completed successfully!\n")

# ============================================================
# CALCULATE & DISPLAY TOTAL TIME
# ============================================================
end_time <- Sys.time()
total_time <- end_time - start_time
message(paste("\n[INFO] Total execution time:", round(total_time, 2)))
# ============================================================

message("\n[INFO] Script completed successfully!\n")