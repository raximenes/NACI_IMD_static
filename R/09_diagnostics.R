################################################################################
# 09_diagnostics.R
# Model Diagnostics and Validation Functions
# 
# OVERVIEW:
# =========
# This module provides comprehensive diagnostic capabilities for validating
# the IMD vaccination model. It enables detailed inspection of model internals
# at both the deterministic and PSA levels.
#
# KEY FEATURES:
# =============
# 1. Extract and save complete eval_all_strategies() outputs
# 2. Create diagnostic summaries for trace matrices and flows
# 3. Validate state transitions and probability conservation
# 4. Generate comparison tables between strategies
# 5. Save PSA-level diagnostics for random samples
#
# USAGE:
# ======
# Source this file after the main model files, then call:
#   - save_deterministic_diagnostics(res_det_list, params, "output_folder")
#   - save_psa_diagnostics(psa_samples, n_samples, "output_folder")
#
################################################################################

## ======================
## DIAGNOSTIC UTILITIES
## ======================

# Helper function to create diagnostic output folder
create_diagnostics_folder <- function(base_path = OUT_TAB) {
  diag_path <- file.path(base_path, "diagnostics")
  dir.create(diag_path, recursive = TRUE, showWarnings = FALSE)
  return(diag_path)
}

# Helper function to save a large object with compression
save_diagnostic_object <- function(obj, name, path, compress = TRUE) {
  filepath <- file.path(path, paste0(name, ".rds"))
  
  tryCatch({
    saveRDS(obj, filepath, compress = compress)
    log_info(paste("âœ“ Saved:", name, "to", basename(filepath)))
    return(TRUE)
  }, error = function(e) {
    log_error(paste("Failed to save", name, ":", e$message))
    return(FALSE)
  })
}

## ======================
## DETERMINISTIC DIAGNOSTICS
## ======================

# Extract comprehensive information from a single strategy's results
extract_strategy_details <- function(strategy_name, strategy_result, params) {
  
  # Basic outcomes
  basic <- list(
    strategy = strategy_name,
    total_cost = strategy_result$cost,
    total_qalys = strategy_result$qalys
  )
  
  # Epidemiological outcomes
  epi <- strategy_result$epi
  epi_summary <- list(
    total_infections = sum(epi$infections),
    infections_by_serogroup = epi$infections,
    total_sequelae = sum(epi$sequelae),
    sequelae_by_type = epi$sequelae,
    deaths_imd = epi$deaths_imd
  )
  
  # Vaccination details
  vacc <- list(
    persons_vaccinated = strategy_result$vacc$persons,
    doses_administered = strategy_result$vacc$doses
  )
  
  # Trace matrix summary (state occupancy over time)
  trace <- strategy_result$trace
  trace_summary <- list(
    n_cycles = nrow(trace) - 1,
    n_states = ncol(trace),
    initial_state = trace[1, ],
    final_state = trace[nrow(trace), ],
    healthy_over_time = trace[, "Healthy"],
    dead_imd_over_time = trace[, "Dead_IMD"],
    background_mort_over_time = trace[, "Background_Mortality"]
  )
  
  # Flow matrix summary (transitions between states)
  flow <- strategy_result$flow
  
  # Calculate total flows by type
  flow_summary <- list(
    total_infections = sum(
      flow["Healthy", "SeroB_Infect", ] +
        flow["Healthy", "SeroC_Infect", ] +
        flow["Healthy", "SeroW_Infect", ] +
        flow["Healthy", "SeroY_Infect", ]
    ),
    total_deaths_imd = sum(
      flow["SeroB_Infect", "Dead_IMD", ] +
        flow["SeroC_Infect", "Dead_IMD", ] +
        flow["SeroW_Infect", "Dead_IMD", ] +
        flow["SeroY_Infect", "Dead_IMD", ]
    ),
    total_background_deaths = sum(flow[, "Background_Mortality", ]),
    total_sequelae = sum(flow[, names(params$w_seq), ])
  )
  
  # Combine all details
  details <- list(
    basic = basic,
    epidemiology = epi_summary,
    vaccination = vacc,
    trace_summary = trace_summary,
    flow_summary = flow_summary,
    full_trace = trace,
    full_flow = flow
  )
  
  return(details)
}

# Save comprehensive deterministic diagnostics
save_deterministic_diagnostics <- function(res_det_list, params, 
                                           output_folder = NULL,
                                           perspective = "healthcare") {
  
  log_info("\n=== Saving Deterministic Model Diagnostics ===")
  
  # Create output folder
  if (is.null(output_folder)) {
    output_folder <- create_diagnostics_folder()
  }
  
  # Add perspective to folder name
  output_folder <- file.path(output_folder, paste0("deterministic_", perspective))
  dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)
  
  log_info(paste("Output folder:", output_folder))
  
  # Extract details for each strategy
  all_details <- list()
  
  for (strategy_name in names(res_det_list)) {
    log_info(paste("Processing:", strategy_name))
    
    strategy_result <- res_det_list[[strategy_name]]
    details <- extract_strategy_details(strategy_name, strategy_result, params)
    
    all_details[[strategy_name]] <- details
    
    # Save individual strategy details
    save_diagnostic_object(
      details, 
      paste0("strategy_", strategy_name, "_full_details"), 
      output_folder
    )
  }
  
  # Create comparison tables across strategies
  comparison_tables <- create_comparison_tables(all_details, params)
  
  # Save comparison tables as CSV
  for (table_name in names(comparison_tables)) {
    filepath <- file.path(output_folder, paste0(table_name, ".csv"))
    tryCatch({
      readr::write_csv(comparison_tables[[table_name]], filepath)
      log_info(paste("âœ“ Saved:", table_name))
    }, error = function(e) {
      log_error(paste("Failed to save", table_name, ":", e$message))
    })
  }
  
  # Save trace matrices for all strategies
  trace_matrices <- lapply(all_details, function(x) x$full_trace)
  save_diagnostic_object(trace_matrices, "all_trace_matrices", output_folder)
  
  # Save flow arrays for all strategies
  flow_arrays <- lapply(all_details, function(x) x$full_flow)
  save_diagnostic_object(flow_arrays, "all_flow_arrays", output_folder)
  
  # Save complete diagnostics object
  diagnostics <- list(
    strategy_details = all_details,
    comparison_tables = comparison_tables,
    parameters_used = params,
    metadata = list(
      perspective = perspective,
      n_cohort = n_cohort,
      n_cycles = params$n_cycles,
      timestamp = Sys.time()
    )
  )
  
  save_diagnostic_object(diagnostics, "complete_diagnostics", output_folder)
  
  # Create validation report
  validation <- validate_deterministic_results(all_details, params)
  filepath <- file.path(output_folder, "validation_report.txt")
  writeLines(validation$report, filepath)
  log_info("âœ“ Saved: validation_report.txt")
  
  log_info("\n=== Deterministic Diagnostics Complete ===")
  log_info(paste("Files saved to:", output_folder))
  
  return(invisible(diagnostics))
}

# Create comparison tables across strategies
create_comparison_tables <- function(all_details, params) {
  
  strategies <- names(all_details)
  
  # Table 1: Basic Outcomes
  basic_outcomes <- do.call(rbind, lapply(strategies, function(s) {
    data.frame(
      Strategy = s,
      Total_Cost = all_details[[s]]$basic$total_cost,
      Total_QALYs = all_details[[s]]$basic$total_qalys,
      stringsAsFactors = FALSE
    )
  }))
  
  # Table 2: Epidemiological Outcomes
  epi_outcomes <- do.call(rbind, lapply(strategies, function(s) {
    epi <- all_details[[s]]$epidemiology
    data.frame(
      Strategy = s,
      Total_Infections = epi$total_infections,
      Infections_B = epi$infections_by_serogroup["B"],
      Infections_C = epi$infections_by_serogroup["C"],
      Infections_W = epi$infections_by_serogroup["W"],
      Infections_Y = epi$infections_by_serogroup["Y"],
      Total_Sequelae = epi$total_sequelae,
      Deaths_IMD = epi$deaths_imd,
      stringsAsFactors = FALSE
    )
  }))
  
  # Table 3: Vaccination Outcomes
  vacc_outcomes <- do.call(rbind, lapply(strategies, function(s) {
    vacc <- all_details[[s]]$vaccination
    data.frame(
      Strategy = s,
      Persons_Vaccinated = vacc$persons_vaccinated,
      Doses_Administered = vacc$doses_administered,
      stringsAsFactors = FALSE
    )
  }))
  
  # Table 4: State Occupancy at Key Timepoints
  state_occupancy <- do.call(rbind, lapply(strategies, function(s) {
    trace <- all_details[[s]]$full_trace
    n_cycles <- nrow(trace) - 1
    
    # Initial, midpoint, and final
    data.frame(
      Strategy = s,
      Healthy_Initial = trace[1, "Healthy"],
      Healthy_Midpoint = trace[n_cycles %/% 2, "Healthy"],
      Healthy_Final = trace[n_cycles + 1, "Healthy"],
      Dead_IMD_Final = trace[n_cycles + 1, "Dead_IMD"],
      Background_Mortality_Final = trace[n_cycles + 1, "Background_Mortality"],
      stringsAsFactors = FALSE
    )
  }))
  
  # Table 5: Flow Summary
  flow_summary <- do.call(rbind, lapply(strategies, function(s) {
    flow_sum <- all_details[[s]]$flow_summary
    data.frame(
      Strategy = s,
      Total_Infections_Flow = flow_sum$total_infections,
      Total_Deaths_IMD_Flow = flow_sum$total_deaths_imd,
      Total_Background_Deaths = flow_sum$total_background_deaths,
      Total_Sequelae_Flow = flow_sum$total_sequelae,
      stringsAsFactors = FALSE
    )
  }))
  
  # Table 6: Sequelae by Type
  sequelae_types <- names(params$w_seq)
  sequelae_by_type <- do.call(rbind, lapply(strategies, function(s) {
    seq_data <- all_details[[s]]$epidemiology$sequelae_by_type
    row_data <- as.data.frame(t(seq_data))
    row_data$Strategy <- s
    return(row_data[, c("Strategy", sequelae_types)])
  }))
  
  return(list(
    basic_outcomes = basic_outcomes,
    epidemiological_outcomes = epi_outcomes,
    vaccination_outcomes = vacc_outcomes,
    state_occupancy = state_occupancy,
    flow_summary = flow_summary,
    sequelae_by_type = sequelae_by_type
  ))
}

# Validate deterministic results
validate_deterministic_results <- function(all_details, params) {
  
  report_lines <- c(
    "=======================================================",
    "DETERMINISTIC MODEL VALIDATION REPORT",
    "=======================================================",
    paste("Generated:", Sys.time()),
    paste("Cohort size:", n_cohort),
    paste("Number of cycles:", params$n_cycles),
    ""
  )
  
  issues_found <- 0
  
  for (strategy_name in names(all_details)) {
    details <- all_details[[strategy_name]]
    
    report_lines <- c(
      report_lines,
      paste("\n--- STRATEGY:", strategy_name, "---")
    )
    
    # Check 1: State occupancy sums to 1 (or cohort size)
    trace <- details$full_trace
    row_sums <- rowSums(trace)
    max_deviation <- max(abs(row_sums - 1))
    
    report_lines <- c(
      report_lines,
      paste("State occupancy sum check:"),
      paste("  Max deviation from 1.0:", format(max_deviation, scientific = TRUE))
    )
    
    if (max_deviation > 1e-6) {
      report_lines <- c(report_lines, "  âš  WARNING: State probabilities do not sum to 1!")
      issues_found <- issues_found + 1
    } else {
      report_lines <- c(report_lines, "  âœ“ PASS")
    }
    
    # Check 2: No negative values in trace
    min_trace <- min(trace)
    report_lines <- c(
      report_lines,
      paste("Negative values check:"),
      paste("  Minimum value in trace:", format(min_trace, scientific = TRUE))
    )
    
    if (min_trace < -1e-10) {
      report_lines <- c(report_lines, "  âš  WARNING: Negative probabilities detected!")
      issues_found <- issues_found + 1
    } else {
      report_lines <- c(report_lines, "  âœ“ PASS")
    }
    
    # Check 3: Death states are absorbing
    dead_imd <- trace[, "Dead_IMD"]
    bg_mort <- trace[, "Background_Mortality"]
    
    dead_imd_decreasing <- any(diff(dead_imd) < -1e-10)
    bg_mort_decreasing <- any(diff(bg_mort) < -1e-10)
    
    report_lines <- c(
      report_lines,
      "Death state monotonicity check:"
    )
    
    if (dead_imd_decreasing) {
      report_lines <- c(report_lines, "  âš  WARNING: Dead_IMD state decreased over time!")
      issues_found <- issues_found + 1
    } else {
      report_lines <- c(report_lines, "  âœ“ Dead_IMD is monotonic increasing")
    }
    
    if (bg_mort_decreasing) {
      report_lines <- c(report_lines, "  âš  WARNING: Background_Mortality state decreased over time!")
      issues_found <- issues_found + 1
    } else {
      report_lines <- c(report_lines, "  âœ“ Background_Mortality is monotonic increasing")
    }
    
    # Check 4: Flow consistency with trace
    flow <- details$full_flow
    flow_sum <- details$flow_summary
    
    # Total infections from flow should match epidemiology
    epi_infections <- details$epidemiology$total_infections
    flow_infections <- flow_sum$total_infections * n_cohort
    
    infection_diff <- abs(epi_infections - flow_infections)
    
    report_lines <- c(
      report_lines,
      "Flow-Epidemiology consistency check:",
      paste("  Epi infections:", format(epi_infections, big.mark = ",")),
      paste("  Flow infections:", format(flow_infections, big.mark = ",")),
      paste("  Difference:", format(infection_diff, big.mark = ","))
    )
    
    if (infection_diff > 1) {
      report_lines <- c(report_lines, "  âš  WARNING: Flow and epi infections do not match!")
      issues_found <- issues_found + 1
    } else {
      report_lines <- c(report_lines, "  âœ“ PASS")
    }
  }
  
  # Summary
  report_lines <- c(
    report_lines,
    "\n=======================================================",
    "VALIDATION SUMMARY",
    "=======================================================",
    paste("Total issues found:", issues_found),
    ifelse(issues_found == 0, 
           "âœ“ ALL CHECKS PASSED", 
           "âš  REVIEW WARNINGS ABOVE"),
    "======================================================="
  )
  
  return(list(
    report = report_lines,
    issues_found = issues_found
  ))
}

## ======================
## PSA DIAGNOSTICS
## ======================

# Save PSA diagnostics for a sample of iterations
save_psa_diagnostics <- function(psa_samples, params_base, n_samples = 10, 
                                 output_folder = NULL, perspective = "healthcare",
                                 seed = 2025) {
  
  log_info("\n=== Saving PSA Model Diagnostics ===")
  log_info(paste("Sampling", n_samples, "iterations for detailed inspection"))
  
  # Create output folder
  if (is.null(output_folder)) {
    output_folder <- create_diagnostics_folder()
  }
  
  # Add perspective to folder name
  output_folder <- file.path(output_folder, paste0("psa_", perspective))
  dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)
  
  log_info(paste("Output folder:", output_folder))
  
  # Select random iterations to inspect
  set.seed(seed)
  total_iterations <- length(psa_samples)
  selected_iterations <- sort(sample(1:total_iterations, min(n_samples, total_iterations)))
  
  log_info(paste("Selected iterations:", paste(selected_iterations, collapse = ", ")))
  
  # Extract diagnostics for each selected iteration
  psa_diagnostics <- list()
  
  for (i in selected_iterations) {
    log_info(paste("Processing PSA iteration", i, "of", total_iterations))
    
    # Run model for this iteration
    params_i <- psa_samples[[i]]
    res_i <- eval_all_strategies(params_i)
    
    # Extract details for each strategy
    iter_details <- list()
    for (strategy_name in names(res_i)) {
      iter_details[[strategy_name]] <- extract_strategy_details(
        strategy_name, 
        res_i[[strategy_name]], 
        params_i
      )
    }
    
    # Save iteration diagnostics
    psa_diagnostics[[paste0("iteration_", i)]] <- list(
      iteration = i,
      parameters = params_i,
      results = iter_details
    )
    
    # Save individual iteration file
    save_diagnostic_object(
      psa_diagnostics[[paste0("iteration_", i)]], 
      paste0("psa_iteration_", i, "_full_details"), 
      output_folder
    )
  }
  
  # Create summary comparison across PSA iterations
  psa_summary <- create_psa_iteration_summary(psa_diagnostics)
  
  # Save summary tables
  for (table_name in names(psa_summary$tables)) {
    filepath <- file.path(output_folder, paste0("psa_", table_name, ".csv"))
    tryCatch({
      readr::write_csv(psa_summary$tables[[table_name]], filepath)
      log_info(paste("âœ“ Saved: psa_", table_name))
    }, error = function(e) {
      log_error(paste("Failed to save psa_", table_name, ":", e$message))
    })
  }
  
  # Save complete PSA diagnostics
  complete_psa_diagnostics <- list(
    selected_iterations = selected_iterations,
    iteration_details = psa_diagnostics,
    summary = psa_summary,
    base_parameters = params_base,
    metadata = list(
      perspective = perspective,
      n_cohort = n_cohort,
      total_psa_iterations = total_iterations,
      sampled_iterations = n_samples,
      timestamp = Sys.time()
    )
  )
  
  save_diagnostic_object(complete_psa_diagnostics, "complete_psa_diagnostics", output_folder)
  
  log_info("\n=== PSA Diagnostics Complete ===")
  log_info(paste("Files saved to:", output_folder))
  
  return(invisible(complete_psa_diagnostics))
}

# Create summary comparison across PSA iterations
create_psa_iteration_summary <- function(psa_diagnostics) {
  
  iterations <- names(psa_diagnostics)
  
  # Extract all strategies (assuming same across iterations)
  first_iter <- psa_diagnostics[[1]]
  strategies <- names(first_iter$results)
  
  # Table 1: Outcomes by iteration and strategy
  outcomes_table <- do.call(rbind, lapply(iterations, function(iter) {
    iter_num <- psa_diagnostics[[iter]]$iteration
    
    do.call(rbind, lapply(strategies, function(s) {
      details <- psa_diagnostics[[iter]]$results[[s]]
      
      data.frame(
        Iteration = iter_num,
        Strategy = s,
        Total_Cost = details$basic$total_cost,
        Total_QALYs = details$basic$total_qalys,
        Total_Infections = details$epidemiology$total_infections,
        Deaths_IMD = details$epidemiology$deaths_imd,
        Total_Sequelae = details$epidemiology$total_sequelae,
        stringsAsFactors = FALSE
      )
    }))
  }))
  
  # Table 2: Infections by serogroup
  infections_table <- do.call(rbind, lapply(iterations, function(iter) {
    iter_num <- psa_diagnostics[[iter]]$iteration
    
    do.call(rbind, lapply(strategies, function(s) {
      inf <- psa_diagnostics[[iter]]$results[[s]]$epidemiology$infections_by_serogroup
      
      data.frame(
        Iteration = iter_num,
        Strategy = s,
        Infections_B = inf["B"],
        Infections_C = inf["C"],
        Infections_W = inf["W"],
        Infections_Y = inf["Y"],
        stringsAsFactors = FALSE
      )
    }))
  }))
  
  # Table 3: Parameter variation summary (for debugging)
  param_variation <- do.call(rbind, lapply(iterations, function(iter) {
    iter_num <- psa_diagnostics[[iter]]$iteration
    params <- psa_diagnostics[[iter]]$parameters
    
    data.frame(
      Iteration = iter_num,
      coverage_ABCWY = params$coverage_ABCWY,
      coverage_ACWY = params$coverage_ACWY,
      p_seq_overall = params$p_seq_overall,
      u_IMD = params$u_IMD,
      c_Scarring = params$c_Scarring,
      stringsAsFactors = FALSE
    )
  }))
  
  return(list(
    tables = list(
      outcomes_by_iteration = outcomes_table,
      infections_by_serogroup = infections_table,
      parameter_variation = param_variation
    ),
    summary_stats = list(
      n_iterations = length(iterations),
      n_strategies = length(strategies),
      strategies = strategies
    )
  ))
}

## ======================
## CONVENIENCE FUNCTIONS
## ======================

# Load saved diagnostics
load_diagnostics <- function(perspective = "healthcare", 
                             type = c("deterministic", "psa"),
                             base_path = OUT_TAB) {
  
  type <- match.arg(type)
  
  # CRITICAL FIX: Match where run_diagnostics.R actually saves files
  # Both deterministic and PSA use perspective-specific subfolders
  
  if (type == "psa") {
    # PSA files are saved in psa_[perspective] subfolder
    # Actual location: outputs/tables/healthcare/psa_healthcare/complete_psa_diagnostics.rds
    filepath_primary <- file.path(base_path, perspective, paste0("psa_", perspective), "complete_psa_diagnostics.rds")
    filepath_fallback <- file.path(base_path, perspective, "diagnostics_psa.rds")
    
  } else {
    # Deterministic files are saved in deterministic_[perspective] subfolder
    # Actual location: outputs/tables/healthcare/deterministic_healthcare/complete_diagnostics.rds
    filepath_primary <- file.path(base_path, perspective, paste0("deterministic_", perspective), "complete_diagnostics.rds")
    filepath_fallback <- file.path(base_path, perspective, "diagnostics_det.rds")
  }
  
  # Try primary location first (where run_diagnostics.R actually saves)
  if (file.exists(filepath_primary)) {
    filepath <- filepath_primary
  } else if (file.exists(filepath_fallback)) {
    filepath <- filepath_fallback
  } else {
    # Neither location exists - show error
    log_error(paste("Diagnostics file not found for", type, "perspective:", perspective))
    message("\nTried locations:")
    message(paste("  1.", filepath_primary))
    message(paste("  2.", filepath_fallback))
    message("\nMake sure you've run: source('run_diagnostics.R')")
    stop(paste("Cannot find diagnostics file"))
  }
  
  log_info(paste("Loading diagnostics from:", filepath))
  diagnostics <- readRDS(filepath)
  
  return(diagnostics)
}

# Quick summary of diagnostics
summarize_diagnostics <- function(diagnostics) {
  
  if ("strategy_details" %in% names(diagnostics)) {
    # Deterministic diagnostics
    cat("\n=== DETERMINISTIC DIAGNOSTICS SUMMARY ===\n")
    cat(paste("Perspective:", diagnostics$metadata$perspective, "\n"))
    cat(paste("N Cohort:", diagnostics$metadata$n_cohort, "\n"))
    cat(paste("N Cycles:", diagnostics$metadata$n_cycles, "\n"))
    cat(paste("Strategies:", length(diagnostics$strategy_details), "\n"))
    cat(paste("  -", paste(names(diagnostics$strategy_details), collapse = ", "), "\n"))
    
  } else if ("iteration_details" %in% names(diagnostics)) {
    # PSA diagnostics
    cat("\n=== PSA DIAGNOSTICS SUMMARY ===\n")
    cat(paste("Perspective:", diagnostics$metadata$perspective, "\n"))
    cat(paste("Total PSA iterations:", diagnostics$metadata$total_psa_iterations, "\n"))
    cat(paste("Sampled iterations:", diagnostics$metadata$sampled_iterations, "\n"))
    cat(paste("Selected iterations:", paste(diagnostics$selected_iterations, collapse = ", "), "\n"))
    
  } else {
    cat("Unknown diagnostics format\n")
  }
}

## ======================
## EXPORT MESSAGE
## ======================

log_info("âœ“ Diagnostics module (09_diagnostics.R) loaded successfully")
log_info("Available functions:")
log_info("  - save_deterministic_diagnostics()")
log_info("  - save_psa_diagnostics()")
log_info("  - load_diagnostics()")
log_info("  - summarize_diagnostics()")