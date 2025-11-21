################################################################################
# 01_validate_model.R
# 
# PURPOSE: Comprehensive validation of IMD vaccination model
# - Checks model structure and parameters
# - Validates transition probabilities
# - Verifies cost and utility assignments
# - Tests PSA parameter variation
# - Generates validation report
#
# WHEN TO RUN: 
# - After any changes to model structure
# - Before running full analysis
# - To diagnose unexpected results
#
# USAGE:
#   source("scripts/01_validate_model.R")
#
# OUTPUT:
#   - Validation report (console)
#   - Validation plots in outputs/figures/[perspective]/validation/
#   - Validation tables in outputs/tables/
#
################################################################################

# ============================================================
# SETUP
# ============================================================

start_time <- Sys.time()

cat("\n")
cat(strrep("=", 80), "\n")
cat("IMD MODEL VALIDATION\n")
cat(strrep("=", 80), "\n\n")

# Load required packages
if (!require("pacman", quietly = TRUE)) {
  install.packages("pacman", repos = "https://cloud.r-project.org")
}
library(pacman)
p_load(here, readxl, ggplot2, dplyr, tidyr, purrr, tibble)

# Source all model files
cat("[INFO] Loading model files...\n")
files <- list.files("R", pattern = "^[0-9]{2}.*\\.R$", full.names = TRUE)
files <- files[order(files)]
invisible(lapply(files, source))
cat("[INFO] Model files loaded\n\n")

# ============================================================
# CONFIGURATION
# ============================================================

PERSPECTIVE <- "healthcare"  # Test perspective
N_TEST_ITERATIONS <- 10      # Number of PSA iterations to test

# Update folders for validation outputs
perspective <- PERSPECTIVE
update_output_folders(PERSPECTIVE)

# Create validation folder
validation_dir <- file.path(OUT_FIG, PERSPECTIVE, "validation")
dir.create(validation_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# VALIDATION CHECKS
# ============================================================

validation_results <- list()

cat(strrep("-", 80), "\n")
cat("VALIDATION CHECK 1: PARAMETER LOADING\n")
cat(strrep("-", 80), "\n")

# Load base parameters
params_bc <- tryCatch({
  get_base_params()
}, error = function(e) {
  cat("[ERROR] Failed to load parameters:", e$message, "\n")
  return(NULL)
})

if (!is.null(params_bc)) {
  cat("[PASS] Parameters loaded successfully\n")
  cat(sprintf("  - Number of cycles: %d\n", params_bc$n_cycles))
  cat(sprintf("  - Perspective: %s\n", params_bc$perspective))
  cat(sprintf("  - Discount rate (cost): %.3f\n", params_bc$d_c))
  cat(sprintf("  - Discount rate (effect): %.3f\n", params_bc$d_e))
  validation_results$params_loaded <- TRUE
} else {
  cat("[FAIL] Parameters not loaded\n")
  validation_results$params_loaded <- FALSE
  stop("Cannot continue without parameters")
}

# ============================================================
cat("\n")
cat(strrep("-", 80), "\n")
cat("VALIDATION CHECK 2: COST PARAMETERS\n")
cat(strrep("-", 80), "\n")

# Check all required cost parameters exist
required_costs <- c(
  "c_admin", "c_Healthy", "c_Scarring", "c_Single_Amput", 
  "c_Multiple_Amput", "c_Neuro_Disab", "c_Hearing_Loss",
  "c_Renal_Failure", "c_Seizure", "c_Paralysis", "c_Dead"
)

missing_costs <- setdiff(required_costs, names(params_bc))
if (length(missing_costs) == 0) {
  cat("[PASS] All required cost parameters present\n")
  validation_results$costs_complete <- TRUE
} else {
  cat("[FAIL] Missing cost parameters:\n")
  cat("  ", paste(missing_costs, collapse = ", "), "\n")
  validation_results$costs_complete <- FALSE
}

# Check for negative costs
negative_costs <- sapply(required_costs, function(c) {
  val <- params_bc[[c]]
  !is.null(val) && any(val < 0, na.rm = TRUE)
})

if (any(negative_costs)) {
  cat("[WARN] Negative cost values detected:\n")
  cat("  ", paste(names(negative_costs)[negative_costs], collapse = ", "), "\n")
  validation_results$no_negative_costs <- FALSE
} else {
  cat("[PASS] No negative costs\n")
  validation_results$no_negative_costs <- TRUE
}

# Check infection costs (time-varying)
if (!is.null(params_bc$c_IMD_infection)) {
  if (length(params_bc$c_IMD_infection) == params_bc$n_cycles) {
    cat("[PASS] Infection costs (time-varying) correct length\n")
    validation_results$infection_costs_ok <- TRUE
  } else {
    cat("[FAIL] Infection costs length mismatch\n")
    cat(sprintf("  Expected: %d, Got: %d\n", params_bc$n_cycles, length(params_bc$c_IMD_infection)))
    validation_results$infection_costs_ok <- FALSE
  }
} else {
  cat("[FAIL] Infection costs not found\n")
  validation_results$infection_costs_ok <- FALSE
}

# ============================================================
cat("\n")
cat(strrep("-", 80), "\n")
cat("VALIDATION CHECK 3: UTILITY PARAMETERS\n")
cat(strrep("-", 80), "\n")

required_utilities <- c(
  "u_Healthy", "u_IMD", "u_Scarring", "u_Single_Amput",
  "u_Multiple_Amput", "u_Neuro_Disability", "u_Hearing_Loss",
  "u_Renal_Failure", "u_Seizure", "u_Paralysis", "u_Dead"
)

missing_utilities <- setdiff(required_utilities, names(params_bc))
if (length(missing_utilities) == 0) {
  cat("[PASS] All required utility parameters present\n")
  validation_results$utilities_complete <- TRUE
} else {
  cat("[FAIL] Missing utility parameters:\n")
  cat("  ", paste(missing_utilities, collapse = ", "), "\n")
  validation_results$utilities_complete <- FALSE
}

# Check utility range [0, 1]
out_of_range <- sapply(required_utilities, function(u) {
  val <- params_bc[[u]]
  !is.null(val) && (any(val < 0, na.rm = TRUE) || any(val > 1, na.rm = TRUE))
})

if (any(out_of_range)) {
  cat("[WARN] Utilities outside [0,1] range:\n")
  for (u in names(out_of_range)[out_of_range]) {
    cat(sprintf("  %s: %.3f\n", u, params_bc[[u]]))
  }
  validation_results$utilities_in_range <- FALSE
} else {
  cat("[PASS] All utilities in valid range [0,1]\n")
  validation_results$utilities_in_range <- TRUE
}

# ============================================================
cat("\n")
cat(strrep("-", 80), "\n")
cat("VALIDATION CHECK 4: TRANSITION PROBABILITIES\n")
cat(strrep("-", 80), "\n")

# Build transition array for one strategy
test_strategy <- "MenC"
cat(sprintf("[INFO] Testing strategy: %s\n", test_strategy))

a_P <- tryCatch({
  build_transition_array(params_bc, test_strategy)
}, error = function(e) {
  cat("[ERROR] Failed to build transition array:", e$message, "\n")
  return(NULL)
})

if (!is.null(a_P)) {
  cat("[PASS] Transition array built successfully\n")
  
  # Check that probabilities sum to 1
  prob_sums <- apply(a_P, 3, function(m) rowSums(m))
  all_sum_to_one <- all(abs(prob_sums - 1) < 1e-6)
  
  if (all_sum_to_one) {
    cat("[PASS] All transition probabilities sum to 1\n")
    validation_results$probs_sum_to_one <- TRUE
  } else {
    cat("[FAIL] Some probabilities don't sum to 1\n")
    max_error <- max(abs(prob_sums - 1))
    cat(sprintf("  Maximum error: %.6f\n", max_error))
    validation_results$probs_sum_to_one <- FALSE
  }
  
  # Check for negative probabilities
  has_negative <- any(a_P < 0, na.rm = TRUE)
  if (has_negative) {
    cat("[FAIL] Negative probabilities detected\n")
    validation_results$no_negative_probs <- FALSE
  } else {
    cat("[PASS] No negative probabilities\n")
    validation_results$no_negative_probs <- TRUE
  }
  
  validation_results$transition_array_ok <- TRUE
} else {
  cat("[FAIL] Could not build transition array\n")
  validation_results$transition_array_ok <- FALSE
  validation_results$probs_sum_to_one <- FALSE
  validation_results$no_negative_probs <- FALSE
}

# ============================================================
cat("\n")
cat(strrep("-", 80), "\n")
cat("VALIDATION CHECK 5: MODEL EXECUTION\n")
cat(strrep("-", 80), "\n")

# Run model for all strategies
all_results <- tryCatch({
  eval_all_strategies(params_bc)
}, error = function(e) {
  cat("[ERROR] Model execution failed:", e$message, "\n")
  return(NULL)
})

if (!is.null(all_results)) {
  cat("[PASS] Model executed successfully\n")
  cat(sprintf("  - Strategies evaluated: %d\n", length(all_results)))
  
  # Check that all results have required components
  all_complete <- all(sapply(all_results, function(r) {
    all(c("cost", "qalys", "epi", "trace", "flow") %in% names(r))
  }))
  
  if (all_complete) {
    cat("[PASS] All strategies returned complete results\n")
    validation_results$model_runs <- TRUE
  } else {
    cat("[FAIL] Some results incomplete\n")
    validation_results$model_runs <- FALSE
  }
  
  # Check for reasonable values
  cat("\n[INFO] Strategy results:\n")
  for (strat in names(all_results)) {
    r <- all_results[[strat]]
    cat(sprintf("  %s:\n", strat))
    cat(sprintf("    Cost:  $%s\n", format(round(r$cost), big.mark = ",")))
    cat(sprintf("    QALYs: %.2f\n", r$qalys))
    cat(sprintf("    Infections: %.0f\n", sum(r$epi$infections)))
    cat(sprintf("    Deaths: %.0f\n", r$epi$deaths_imd))
  }
  
} else {
  cat("[FAIL] Model execution failed\n")
  validation_results$model_runs <- FALSE
}

# ============================================================
cat("\n")
cat(strrep("-", 80), "\n")
cat("VALIDATION CHECK 6: PSA PARAMETER VARIATION\n")
cat(strrep("-", 80), "\n")

cat(sprintf("[INFO] Testing with %d PSA iterations...\n", N_TEST_ITERATIONS))

psa_test <- tryCatch({
  generate_psa_samples(params_bc, n_sim = N_TEST_ITERATIONS, seed = 2025)
}, error = function(e) {
  cat("[ERROR] PSA sample generation failed:", e$message, "\n")
  return(NULL)
})

if (!is.null(psa_test)) {
  cat("[PASS] PSA samples generated\n")
  
  # NEW APPROACH: Check variation in parameter SAMPLES (not psa_df)
  samples <- psa_test$samples
  
  if (length(samples) == 0) {
    cat("[FAIL] No PSA samples returned\n")
    validation_results$psa_generates <- FALSE
    validation_results$psa_varies <- FALSE
  } else {
    validation_results$psa_generates <- TRUE
    
    # Get all parameter names from first sample
    all_param_names <- names(samples[[1]])
    
    # Remove non-numeric parameters
    numeric_params <- all_param_names[sapply(samples[[1]], function(x) {
      is.numeric(x) && length(x) == 1
    })]
    
    cat(sprintf("  - PSA iterations: %d\n", length(samples)))
    cat(sprintf("  - Total scalar parameters: %d\n", length(numeric_params)))
    
    # Check which parameters vary across iterations
    varying_params <- character(0)
    fixed_params <- character(0)
    
    for (param_name in numeric_params) {
      # Extract values across all iterations
      values <- sapply(samples, function(s) {
        val <- s[[param_name]]
        if (is.null(val) || !is.numeric(val)) return(NA)
        return(val)
      })
      
      # Skip if all NA
      if (all(is.na(values))) next
      
      # Check for variation
      param_sd <- sd(values, na.rm = TRUE)
      param_mean <- mean(values, na.rm = TRUE)
      
      # Parameter varies if SD > 0 and CV > 0.0001
      if (!is.na(param_sd) && param_sd > 1e-10) {
        if (abs(param_mean) > 1e-10) {
          cv <- param_sd / abs(param_mean)
          if (cv > 1e-4) {
            varying_params <- c(varying_params, param_name)
          } else {
            fixed_params <- c(fixed_params, param_name)
          }
        } else {
          # Mean ~0 but SD > 0, so it varies
          varying_params <- c(varying_params, param_name)
        }
      } else {
        fixed_params <- c(fixed_params, param_name)
      }
    }
    
    cat(sprintf("  - Parameters varying in PSA: %d\n", length(varying_params)))
    cat(sprintf("  - Parameters fixed in PSA: %d\n", length(fixed_params)))
    
    if (length(varying_params) > 0) {
      cat("[PASS] PSA has parameter variation\n")
      validation_results$psa_varies <- TRUE
      
      cat("\n[INFO] Examples of varying parameters:\n")
      # Show first 10 varying parameters
      show_params <- head(varying_params, 10)
      for (param in show_params) {
        values <- sapply(samples, function(s) s[[param]])
        cat(sprintf("  - %s: mean = %.4f, SD = %.4f\n", 
                    param, 
                    mean(values, na.rm = TRUE), 
                    sd(values, na.rm = TRUE)))
      }
      if (length(varying_params) > 10) {
        cat(sprintf("  ... and %d more varying parameters\n", length(varying_params) - 10))
      }
    } else {
      cat("[FAIL] No parameter variation in PSA\n")
      validation_results$psa_varies <- FALSE
    }
    
    if (length(fixed_params) > 0) {
      cat("\n[INFO] Examples of fixed parameters (expected for vaccine prices):\n")
      # Show first 10 fixed parameters
      show_fixed <- head(fixed_params, 10)
      cat("  ", paste(show_fixed, collapse = ", "), "\n")
      if (length(fixed_params) > 10) {
        cat(sprintf("  ... and %d more fixed parameters\n", length(fixed_params) - 10))
      }
    }
  }
  
} else {
  cat("[FAIL] PSA generation failed\n")
  validation_results$psa_generates <- FALSE
  validation_results$psa_varies <- FALSE
}

# ============================================================
cat("\n")
cat(strrep("-", 80), "\n")
cat("VALIDATION CHECK 7: OUTPUT FOLDER STRUCTURE\n")
cat(strrep("-", 80), "\n")

required_folders <- c(
  OUT_FIG_DET,
  OUT_FIG_OWSA,
  OUT_FIG_KEY,
  OUT_TAB
)

all_exist <- all(sapply(required_folders, dir.exists))

if (all_exist) {
  cat("[PASS] All required output folders exist\n")
  validation_results$folders_exist <- TRUE
} else {
  cat("[FAIL] Some output folders missing\n")
  missing <- required_folders[!sapply(required_folders, dir.exists)]
  cat("  Missing:\n")
  for (f in missing) {
    cat("    -", f, "\n")
  }
  validation_results$folders_exist <- FALSE
}

cat("\n[INFO] Output structure:\n")
cat("  Base Case:  ", OUT_FIG_DET, "\n")
cat("  OWSA:       ", OUT_FIG_OWSA, "\n")
cat("  Key Plots:  ", OUT_FIG_KEY, "\n")
cat("  Tables:     ", OUT_TAB, "\n")

# ============================================================
cat("\n")
cat(strrep("=", 80), "\n")
cat("VALIDATION SUMMARY\n")
cat(strrep("=", 80), "\n\n")

# Count passes and fails
all_checks <- unlist(validation_results)
n_pass <- sum(all_checks == TRUE, na.rm = TRUE)
n_fail <- sum(all_checks == FALSE, na.rm = TRUE)
n_total <- length(all_checks)

cat(sprintf("Total checks: %d\n", n_total))
cat(sprintf("  [PASS]: %d (%.0f%%)\n", n_pass, 100 * n_pass / n_total))
cat(sprintf("  [FAIL]: %d (%.0f%%)\n", n_fail, 100 * n_fail / n_total))
cat("\n")

# Print individual results
for (check_name in names(validation_results)) {
  result <- validation_results[[check_name]]
  symbol <- ifelse(result, "✓", "✗")
  status <- ifelse(result, "PASS", "FAIL")
  cat(sprintf("  [%s] %s: %s\n", symbol, check_name, status))
}

# Overall assessment
cat("\n")
if (n_fail == 0) {
  cat(strrep("=", 80), "\n")
  cat("✓✓✓ MODEL VALIDATION PASSED ✓✓✓\n")
  cat(strrep("=", 80), "\n")
  cat("\nYour model:\n")
  cat("  ✓ Has all required parameters\n")
  cat("  ✓ Uses valid transition probabilities\n")
  cat("  ✓ Executes without errors\n")
  cat("  ✓ Generates PSA samples with variation\n")
  cat("  ✓ Has correct output structure\n")
  cat("\nThe model is ready for full analysis!\n")
} else {
  cat(strrep("=", 80), "\n")
  cat("✗✗✗ VALIDATION ISSUES DETECTED ✗✗✗\n")
  cat(strrep("=", 80), "\n")
  cat("\nPlease review the failed checks above.\n")
  
  # Special note about PSA variation
  if (!validation_results$psa_varies && validation_results$psa_generates) {
    cat("\n[NOTE] PSA Parameter Variation:\n")
    cat("  The PSA generates samples successfully, but no variation was detected.\n")
    cat("  This can happen if:\n")
    cat("    1. Excel file has no 'sd' column (all parameters fixed by design)\n")
    cat("    2. All SD values are 0 or NA (intentionally fixed parameters)\n")
    cat("    3. Only vaccine prices are present (these are ALWAYS fixed)\n")
    cat("\n  If you WANT PSA variation:\n")
    cat("    - Add 'sd' column to Excel sheets\n")
    cat("    - Set SD > 0 for parameters that should vary\n")
    cat("    - Re-run this validation\n")
    cat("\n  If parameters are INTENTIONALLY fixed:\n")
    cat("    - This is OK for deterministic analysis\n")
    cat("    - PSA will run but all iterations will be identical\n")
    cat("    - This warning can be ignored\n")
  } else {
    cat("\nCommon issues:\n")
    cat("  - Missing parameters in Excel file\n")
    cat("  - Incorrect SD values (too high causing negative probs)\n")
    cat("  - Missing 'sd' column in Excel\n")
  }
}

# Save validation report
validation_file <- file.path(OUT_TAB, "model_validation_report.csv")
validation_df <- data.frame(
  Check = names(validation_results),
  Status = ifelse(unlist(validation_results), "PASS", "FAIL"),
  stringsAsFactors = FALSE
)
write.csv(validation_df, validation_file, row.names = FALSE)
cat(sprintf("\nValidation report saved to: %s\n", validation_file))

# ============================================================
# EXECUTION TIME
# ============================================================

end_time <- Sys.time()
elapsed <- difftime(end_time, start_time, units = "secs")

cat("\n")
cat(strrep("=", 80), "\n")
cat(sprintf("Validation completed in %.1f seconds\n", elapsed))
cat(strrep("=", 80), "\n\n")