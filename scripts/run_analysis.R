################################################################################
# run_analysis.R
# Main script to run the complete IMD vaccination cost-effectiveness analysis
#
# WORKFLOW OVERVIEW:
# ==================
# This script is designed to work in conjunction with run_psa_distributions_only.R
# for efficient PSA parameter validation and model execution.
#
# STEP 1: Run complete cost-effectiveness analysis
#         Run: source("run_analysis.R")
#         - Runs base_case analysis
#         - Runs probabilistic sensitivity analysis (PSA)
#         - Runs one-way sensitivity analysis (OWSA)
#         - Generates all cost-effectiveness outputs (ICER, CEAC, EVPI, etc.)
#
# KEY CONFIGURATION PARAMETERS:
# ------------------------------
# - COHORT_SIZE: Number of individuals in the cohort
# - N_SIMULATIONS: Number of PSA iterations (use ≥1000 for final analysis)
# - ANALYSIS_PERSPECTIVE: "healthcare", "societal", or "both"
# - WTP_THRESHOLD: Willingness-to-pay threshold for cost-effectiveness
# - CREATE_OWSA_PLOTS: TRUE to generate tornado diagrams
#
# OUTPUTS:
# --------
# Tables (outputs/tables/):
#   - base_case_results_[perspective].csv
#   - base_case_icers_[perspective].csv
#   - psa_long_results_[perspective].csv
#   - owsa_output_cost_[perspective].csv
#   - owsa_output_qalys_[perspective].csv
#
# Figures (outputs/figures/):
#   - base_case_icers_frontier_[perspective].png
#   - psa_ceac_[perspective].png (cost-effectiveness acceptability curve)
#   - psa_evpi_[perspective].png (expected value of perfect information)
#   - owsa_tornado_*.png (one-way sensitivity tornado diagrams)
#
################################################################################

# ============================================================
# START TIMER
# ============================================================
start_time <- Sys.time()
# ============================================================
# PACKAGE INSTALLATION AND LOADING
# ============================================================

# Install pacman if needed
if (!require("pacman", quietly = TRUE)) {
  install.packages("pacman", repos = "https://cloud.r-project.org")
}

# Load packages
library(pacman)
p_load(here, readxl, ggplot2, dplyr, tidyr, purrr, stringr,
       tibble, doParallel, foreach, parallel, diagram, reshape2, readr)
p_load_gh("DARTH-git/darthtools")
p_load(dampack)


# ============================================================
# SOURCE ALL MODEL FILES
# ============================================================

message("[INFO] Sourcing model files...")

# Source numbered R files in order
files <- list.files("R", pattern = "^[0-9]{2}.*\\.R$", full.names = TRUE)
files <- files[order(files)]

for (file in files) {
  message(paste("  Loading:", basename(file)))
  source(file)
}

message("[INFO] All model files loaded successfully\n")
# ============================================================
# CONFIGURATION SECTION - EDIT THESE VALUES
# ============================================================

# Analysis settings
N_SIMULATIONS <- 2500L       # Number of PSA simulations
COHORT_SIZE <- n_cohort       # Number of individuals in the cohort -  defined in 02_globals.R

WTP_THRESHOLD <- 50000      # Willingness-to-pay threshold (CAD per QALY)
ANALYSIS_PERSPECTIVE <-  "both"  # "healthcare", "societal", or "both"

# OWSA Configuration 
#owsa_comparators <- c("all")  # Options: "MenC", "MenACWY", c("MenC", "MenACWY"), "all"
# Note: OWSA will automatically use the same comparator as base case
owsa_strategies <- NULL  # NULL = all except comparators, or specify: c("MenACWY", "MenABCWY")

# Visualization options
CREATE_OWSA_PLOTS <- TRUE   # Set to FALSE to skip OWSA tornado plots

# SET seed for reproducibility
RANDOM_SEED <- 2025

# ============================================================
# VALIDATE AND SET PERSPECTIVES
# ============================================================

# Validate perspective input
if (!ANALYSIS_PERSPECTIVE %in% c("healthcare", "societal", "both")) {
  stop("ANALYSIS_PERSPECTIVE must be 'healthcare', 'societal', or 'both'")
}

# Determine which perspectives to run
if (ANALYSIS_PERSPECTIVE == "both") {
  perspectives_to_run <- c("healthcare", "societal")
} else {
  perspectives_to_run <- ANALYSIS_PERSPECTIVE
}

message("\n========================================")
message("IMD VACCINATION ANALYSIS - CONFIGURATION")
message("========================================")
message(paste("Cohort Size:", COHORT_SIZE))
message(paste("PSA Simulations:", N_SIMULATIONS))
message(paste("WTP Threshold: CAD", format(WTP_THRESHOLD, big.mark = ",")))
message(paste("Perspective(s):", paste(perspectives_to_run, collapse = " & ")))
message(paste("Random Seed:", RANDOM_SEED))
message(paste("OWSA Plots:", ifelse(CREATE_OWSA_PLOTS, "YES", "NO")))
message("========================================\n")

# Set random seed
set.seed(RANDOM_SEED)


# ============================================================
# RUN ANALYSES FOR EACH PERSPECTIVE
# ============================================================

# Store all results
all_results <- list()

for (current_perspective in perspectives_to_run) {
  
  tryCatch({
    
    message("\n========================================")
    message(paste("RUNNING ANALYSIS:", toupper(current_perspective), "PERSPECTIVE"))
    message("========================================\n")
    
    
    
    # UPDATE GLOBAL PERSPECTIVE AND FOLDER STRUCTURE
    
    # Set perspective
    perspective <- current_perspective
    n_cohort <- as.integer(COHORT_SIZE)
    
    update_output_folders(current_perspective)
    
    # Log where outputs will be saved
    message("[INFO] Base Case plots will be saved to: ", OUT_FIG_DET)
    message("[INFO] Key plots (CEAC, EVPI) will be saved to: ", OUT_FIG_KEY)
    message("[INFO] OWSA plots will be saved to: ", OUT_FIG_OWSA)
    
    
    
    # Reload base parameters with correct perspective
    params_bc <- get_base_params()
    
    # ============================================================
    # BASE CASE ANALYSIS
    # ============================================================
    
    message("[INFO] Running Base Case Analysis...")
    res_det_list <- eval_all_strategies(params_bc)
    det <- summarize_det(res_det_list) # Auto select the lower cost comparator

    message("Base Case Results:")
    print(det$table)
    print(det$icers)
    print(det$incremental)
    
    # Plot ICER frontier
    p_icer <- try(plot(det$icers), silent = TRUE)
    if (inherits(p_icer, "ggplot")) {
      filename <- file.path(OUT_FIG_DET, paste0("base_case_icers_frontier_", current_perspective, ".png"))
      try(ggsave(filename, p_icer, width = 7, height = 5, dpi = 300))
      try(print(p_icer))
    }
    
    # ============================================================
    # PSA ANALYSIS
    # ============================================================
    
    message("\n[INFO] Running Probabilistic Sensitivity Analysis (PSA)...")
    message(paste("[INFO] Perspective:", current_perspective))
    
    # Run PSA using the centralized run_psa function
    # This handles parameter generation and simulation in parallel
    psa_results <- run_psa(
      params = params_bc,
      n_sim = N_SIMULATIONS,
      wtp = WTP_THRESHOLD,
      seed = RANDOM_SEED
    )
    
    # Create PSA objects for analysis (CEAC, EVPI)
    objs <- make_psa_objects(psa_results, strategies = v_strats, wtp_vec = v_wtp)
    
    # Plot CEAC
    p_ceac <- try(plot(objs$ceac), silent = TRUE)
    if (inherits(p_ceac, "ggplot")) {
      filename <- file.path(OUT_FIG_DET, paste0("psa_ceac_", current_perspective, ".png"))
      try(ggsave(filename, p_ceac, width = 7, height = 5, dpi = 300))
    }
    
    # Plot EVPI
    p_evpi <- try(plot(objs$evpi), silent = TRUE)
    if (inherits(p_evpi, "ggplot")) {
      filename <- file.path(OUT_FIG_DET, paste0("psa_evpi_", current_perspective, ".png"))
      try(ggsave(filename, p_evpi, width = 7, height = 5, dpi = 300))
    }
    
    # ============================================================
    # OWSA
    # ============================================================
    
    if (CREATE_OWSA_PLOTS) {
      
      message("\n[INFO] Running One-Way Sensitivity Analysis (OWSA)...")
      message("[INFO] Using same comparator as base case: ", det$comparator)
      message("[INFO] OWSA plots will be saved to: ", OUT_FIG_OWSA)
      
      # Use same comparator as deterministic analysis for consistency

      owsa_out <- execute_owsa(
        params_bc, 
        comparators = det$comparator,  # Auto-selected from base case
        strategies_to_test = owsa_strategies
      )
      
      # Generate tornado plots
      create_owsa_plots(
        owsa_out, 
        det, 
        strategies = owsa_strategies,
        wtp = WTP_THRESHOLD,
        plot_type = "both",
        perspective = current_perspective
      )
      
      message("[INFO] OWSA completed with comparator: ", det$comparator)
      
    } else {
      owsa_out <- NULL
      message("\n[INFO] OWSA plots skipped (CREATE_OWSA_PLOTS = FALSE)")
    }
    
    # ============================================================
    # SAVE RESULTS
    # ============================================================
    
    message("\n[INFO] Saving results...")
    
    # Add perspective suffix to filenames
    suffix <- paste0("_", current_perspective)
    
    # Save base case results table
    tryCatch({
      message("[DEBUG] Saving base case results table...")
      readr::write_csv(det$table, file.path(OUT_TAB, paste0("base_case_results", suffix, ".csv")))
      message("[DEBUG] ✓ Base Case table saved")
    }, error = function(e) {
      message("[ERROR] Failed to save Base Case table: ", e$message)
    })
    
    # Save ICERs
    tryCatch({
      message("[DEBUG] Saving ICERs...")
      icers_df <- as.data.frame(det$icers, stringsAsFactors = FALSE)
      readr::write_csv(icers_df, file.path(OUT_TAB, paste0("base_case_icers", suffix, ".csv")))
      message("[DEBUG] ✓ ICERs saved")
    }, error = function(e) {
      message("[ERROR] Failed to save ICERs: ", e$message)
      message("[INFO] Trying alternative ICER save method...")
      tryCatch({
        icers_df <- data.frame(
          Strategy = as.character(det$icers$Strategy),
          Cost = as.numeric(det$icers$Cost),
          Effect = as.numeric(det$icers$Effect),
          Inc_Cost = as.numeric(det$icers$Inc_Cost),
          Inc_Effect = as.numeric(det$icers$Inc_Effect),
          ICER = as.numeric(det$icers$ICER),
          Status = as.character(det$icers$Status),
          stringsAsFactors = FALSE
        )
        readr::write_csv(icers_df, file.path(OUT_TAB, paste0("base_case_icers", suffix, ".csv")))
        message("[DEBUG] ✓ ICERs saved (alternative method)")
      }, error = function(e2) {
        message("[ERROR] Alternative ICER save also failed: ", e2$message)
      })
    })
    
    # Save incremental epidemiological results
    tryCatch({
      message("[DEBUG] Saving incremental epi results...")
      readr::write_csv(det$incremental, file.path(OUT_TAB, paste0("base_case_incremental_epi", suffix, ".csv")))
      message("[DEBUG] ✓ Incremental epi saved")
    }, error = function(e) {
      message("[ERROR] Failed to save incremental epi: ", e$message)
      message("[DEBUG] Class of det$incremental: ", class(det$incremental))
      message("[DEBUG] Is dataframe?: ", is.data.frame(det$incremental))
    })
    
    # Save PSA results
    tryCatch({
      message("[DEBUG] Saving PSA results...")
      readr::write_csv(psa_results, file.path(OUT_TAB, paste0("psa_long_results", suffix, ".csv")))
      message("[DEBUG] ✓ PSA results saved")
    }, error = function(e) {
      message("[ERROR] Failed to save PSA results: ", e$message)
      message("[DEBUG] Class of psa_results: ", class(psa_results))
      message("[DEBUG] Is dataframe?: ", is.data.frame(psa_results))
    })
    
    # Save OWSA results
    if (!is.null(owsa_out)) {
      tryCatch({
        message("[DEBUG] Saving OWSA results...")
        # OWSA now uses single comparator (same as base case)
        comp <- det$comparator
        message("[DEBUG] Comparator used: ", comp)
        
        readr::write_csv(owsa_out$owsa_Cost, 
                         file.path(OUT_TAB, paste0("owsa_output_cost_", comp, suffix, ".csv")))
        readr::write_csv(owsa_out$owsa_QALYs, 
                         file.path(OUT_TAB, paste0("owsa_output_qalys_", comp, suffix, ".csv")))
        
        message("[DEBUG] ✓ OWSA saved (comparator: ", comp, ")")
      }, error = function(e) {
        message("[ERROR] Failed to save OWSA: ", e$message)
      })
    }
    
    # Store results
    all_results[[current_perspective]] <- list(
      base_case = det,
      psa = psa_results,
      psa_objects = objs,
      owsa = owsa_out
    )
    
    message(paste("\n[INFO] ===", toupper(current_perspective), "perspective analysis completed successfully ===\n"))
  }, error = function(e) {
    message("\n========================================")
    message(paste("ERROR:", toupper(current_perspective), "PERSPECTIVE FAILED"))
    message("========================================")
    message(paste("Error message:", e$message))
    message("\nStacktrace:")
    traceback()
    message("========================================\n")
  })
}

# ============================================================
# PRINT SUMMARY FOR ALL PERSPECTIVES
# ============================================================

message("\n========================================")
message("ANALYSIS SUMMARY - ALL PERSPECTIVES")
message("========================================")

for (persp in perspectives_to_run) {
  if (!is.null(all_results[[persp]])) {
    message("\n", strrep("=", 50))
    message(paste("PERSPECTIVE:", toupper(persp)))
    message(strrep("=", 50))
    
    # Base Case results
    message("\n--- BASE CASE RESULTS ---")
    print(all_results[[persp]]$base_case$table)
    
    message("\n--- ICERS ---")
    print(all_results[[persp]]$base_case$icers)
    
    # PSA summary
    if (!is.null(all_results[[persp]]$psa)) {
      message("\n--- PSA SUMMARY ---")
      psa_summary_by_strat <- all_results[[persp]]$psa %>%
        dplyr::group_by(Strategy) %>%
        dplyr::summarise(
          Mean_Cost = mean(Cost),
          SD_Cost = sd(Cost),
          Mean_QALYs = mean(QALYs),
          SD_QALYs = sd(QALYs),
          Mean_NMB = mean(NMB),
          SD_NMB = sd(NMB),
          .groups = "drop"
        )
      print(psa_summary_by_strat)
    }
  }
}

message("\n[INFO] Script completed successfully!\n")

# ============================================================
# CALCULATE & DISPLAY TOTAL TIME
# ============================================================
end_time <- Sys.time()
total_time <- end_time - start_time
message(paste("\n[INFO] Total execution time:", round(total_time, 2)))
# ============================================================

message("\n========================================")
message("ALL ANALYSES COMPLETE")
message("========================================")
message(paste("Results saved to:", OUT_TAB))
message(paste("Figures saved to:", OUT_FIG))
message(paste("Perspectives analyzed:", paste(perspectives_to_run, collapse = ", ")))
message("========================================\n")

# ============================================================
# CLEANUP
# ============================================================

# Optional: Clean up large objects to free memory
# Uncomment if needed
# rm(all_results)
# gc()

message("\n[INFO] Script completed. Check outputs folder for results.\n")