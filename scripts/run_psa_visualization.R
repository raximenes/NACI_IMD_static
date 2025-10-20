# run_psa_visualization.R
# Script to visualize PSA parameter distributions
# Run this separately when you want to inspect PSA distributions

# Source all necessary files
files <- list.files("R", pattern = "^[0-9]{2}.*\\.R$", full.names = TRUE)
files <- files[order(files)]
invisible(lapply(files, source))

# Source visualization functions
source("R/06b_psa_visualization.R")

log_info("=== PSA Distribution Visualization Script ===")

# Load base parameters
params_bc <- get_base_params()

# Generate PSA samples (without running full model)
log_info("\n=== Generating PSA Samples ===")
psa_output <- generate_psa_samples(params_bc, n_sim = 1000, seed = 2025)

psa_samples <- psa_output$samples
psa_df <- psa_output$psa_df

log_info(paste("Generated", nrow(psa_df), "PSA samples"))
log_info(paste("Number of parameters:", ncol(psa_df)))

# Create all visualization plots
log_info("\n=== Creating Distribution Plots ===")
create_all_psa_distribution_plots(psa_df, params_bc)

log_info("\n=== PSA Visualization Complete ===")
log_info(paste("All plots saved to:", OUT_FIG))

# Optional: Create summary statistics table
log_info("\n=== Creating PSA Summary Statistics ===")

# Get scalar parameters
scalar_params <- names(psa_df)[!grepl("\\d+$", names(psa_df))]

# Calculate statistics
psa_summary <- data.frame(
  Parameter = scalar_params,
  Deterministic = sapply(scalar_params, function(p) params_bc[[p]]),
  PSA_Mean = sapply(scalar_params, function(p) mean(psa_df[[p]], na.rm = TRUE)),
  PSA_SD = sapply(scalar_params, function(p) sd(psa_df[[p]], na.rm = TRUE)),
  PSA_Q025 = sapply(scalar_params, function(p) quantile(psa_df[[p]], 0.025, na.rm = TRUE)),
  PSA_Q975 = sapply(scalar_params, function(p) quantile(psa_df[[p]], 0.975, na.rm = TRUE))
) %>%
  mutate(
    Pct_Diff = abs((PSA_Mean - Deterministic) / Deterministic) * 100
  ) %>%
  arrange(desc(Pct_Diff))

print(psa_summary)

# Save summary
readr::write_csv(psa_summary, file.path(OUT_TAB, "psa_parameter_summary.csv"))

log_info("PSA summary statistics saved")

# Clean up
rm(psa_samples, psa_df, psa_output)
gc()

log_info("\n=== Done! ===")