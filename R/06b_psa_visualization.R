## ======================
## 6b) PSA Visualization Functions (FIXED - No log_warn)
## ======================
# Functions to visualize PSA parameter distributions
# All log_warn replaced with message() for parallel compatibility

# ============================================================================
# TEMPORAL SUMMARY PLOT EXPLANATION
# ============================================================================
# 
# This plot shows how a time-varying parameter changes over cycles (age)
# with uncertainty from PSA:
#
# LINES:
#   - Blue solid line = PSA Mean (average of all PSA simulations)
#   - Red dashed line = Deterministic (fixed value from Excel)
#
# SHADED AREA:
#   - Light blue ribbon = 95% Credible Interval (CI)
#   - Shows the range where 95% of PSA simulations fall
#   - Wider area = more uncertainty in the parameter
#
# INTERPRETATION:
#   ✓ Lines should be CLOSE: If PSA Mean ≈ Deterministic, your PSA is 
#     well-calibrated and distributions are correctly specified
#   ✗ Lines far apart: May indicate problems with SD specification in Excel
#     or distribution parameterization
#
# UNCERTAINTY PATTERNS:
#   - Narrow CI = parameter is well-known with low uncertainty
#   - Wide CI = high uncertainty, need better data or more research
#   - CI width may increase at older ages due to sparse epidemiological data
#
# NOTE: With few PSA iterations (e.g., n=4), the PSA mean may be unstable.
#       Use n≥1000 for reliable PSA results.
# ============================================================================

# Plot mean and 95% CI over time for temporal parameters
plot_psa_temporal_summary <- function(psa_df, params_base,
                                      param_name = "p_B",
                                      plot_title = "Infection Probability Over Time") {
  
  log_info(paste("Creating temporal summary plot for", param_name, "..."))
  
  # Gather all cycles
  param_cols <- grep(paste0("^", param_name, "\\d+$"), names(psa_df), value = TRUE)
  
  if (length(param_cols) == 0) {
    message(paste("[WARN] No temporal parameters found for", param_name))
    return(invisible(NULL))
  }
  
  n_cycles <- length(param_cols)
  
  # Calculate statistics for each cycle
  summary_data <- data.frame(
    cycle = 1:n_cycles,
    mean = numeric(n_cycles),
    lower = numeric(n_cycles),
    upper = numeric(n_cycles),
    deterministic = params_base[[param_name]]
  )
  
  for (i in 1:n_cycles) {
    col_name <- paste0(param_name, i)
    samples <- psa_df[[col_name]]
    
    summary_data$mean[i] <- mean(samples)
    summary_data$lower[i] <- quantile(samples, 0.025)
    summary_data$upper[i] <- quantile(samples, 0.975)
  }
  
  # Create plot
  p <- ggplot(summary_data, aes(x = cycle)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), 
                fill = "#69b3a2", 
                alpha = 0.3) +
    geom_line(aes(y = mean, color = "PSA Mean"), 
              size = 1) +
    geom_line(aes(y = deterministic, color = "Deterministic"), 
              size = 1, 
              linetype = "dashed") +
    scale_color_manual(
      values = c("PSA Mean" = "#404080", "Deterministic" = "red"),
      name = ""
    ) +
    labs(
      title = plot_title,
      subtitle = "Mean with 95% Credible Interval (PSA)",
      x = "Cycle (Age)",
      y = "Parameter Value"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      legend.position = "top"
    )
  
  print(p)
  
  filename <- file.path(OUT_FIG, paste0("psa_temporal_summary_", param_name, ".png"))
  try(ggsave(filename, p, width = 8, height = 5, dpi = 300))
  
  log_info("Temporal summary plot created")
  
  return(invisible(p))
}

# ============================================================================
# TEMPORAL DISTRIBUTION PLOTS EXPLANATION (Individual Cycles)
# ============================================================================
#
# These histograms show the DISTRIBUTION of a parameter at SPECIFIC cycles
# (ages) based on PSA simulations.
#
# WHAT YOU SEE:
#   - Histogram bars = Frequency of parameter values across PSA simulations
#   - Black density curve = Smooth estimate of the distribution shape
#   - Red dashed line = Deterministic value (from Excel)
#
# INTERPRETATION:
#
# 1. DISTRIBUTION SHAPE:
#    - Bell-shaped (Normal) = Most common for costs, utilities
#    - Right-skewed (Gamma) = Common for cost parameters (can't be negative)
#    - Bounded (Beta) = For probabilities (must be between 0 and 1)
#
# 2. SPREAD (Width):
#    - Narrow distribution = Low uncertainty, precise parameter
#    - Wide distribution = High uncertainty, imprecise parameter
#    - Spread is controlled by SD in Excel
#
# 3. CENTRAL TENDENCY:
#    - Peak of distribution should be near the deterministic value (red line)
#    - If far apart, check SD specification in Excel
#
# 4. COMPARING ACROSS CYCLES:
#    Example: p_B_cycle1 vs p_B_cycle89
#    - Cycle 1 (age 11): High infection risk, narrow distribution (good data)
#    - Cycle 89 (age 100): Low infection risk, wider distribution (sparse data)
#    - Uncertainty often increases with age due to limited epidemiological data
#
# PRACTICAL USE:
#   - Identify which cycles/ages have the most uncertainty
#   - Prioritize data collection for parameters with wide distributions
#   - Validate that distributions match expected parameter ranges
#
# NOTE: With n=4 simulations, histograms will be very crude (only 4 bars).
#       Use n≥1000 for smooth, interpretable distributions.
# ============================================================================

# Plot distributions for time-varying parameters
plot_psa_temporal_distributions <- function(psa_df, params_base,
                                            param_name = "p_B",
                                            cycles_to_plot = c(1, 5, 10, 20, 40, 60, 89)) {
  
  log_info(paste("Creating temporal distribution plots for", param_name, "..."))
  
  # Gather all cycles for this parameter
  param_cols <- grep(paste0("^", param_name, "\\d+$"), names(psa_df), value = TRUE)
  
  if (length(param_cols) == 0) {
    message(paste("[WARN] No temporal parameters found for", param_name))
    return(invisible(NULL))
  }
  
  n_cycles <- length(param_cols)
  log_info(paste("Found", n_cycles, "cycles for", param_name))
  
  # Ensure cycles_to_plot are valid
  cycles_to_plot <- cycles_to_plot[cycles_to_plot <= n_cycles]
  
  plots_list <- list()
  
  # Plot distributions for selected cycles
  for (cycle in cycles_to_plot) {
    col_name <- paste0(param_name, cycle)
    
    if (!col_name %in% names(psa_df)) next
    
    samples <- psa_df[[col_name]]
    det_value <- params_base[[param_name]][cycle]
    
    p <- ggplot(data.frame(value = samples), aes(x = value)) +
      geom_histogram(aes(y = after_stat(density)), 
                     bins = 40, 
                     fill = "#e07a5f", 
                     alpha = 0.7,
                     color = "white") +
      geom_density(color = "#3d405b", size = 1.2) +
      geom_vline(xintercept = det_value, 
                 color = "red", 
                 linetype = "dashed", 
                 size = 1) +
      labs(
        title = paste0(param_name, " - Cycle ", cycle),
        subtitle = paste0("Mean = ", format(mean(samples), scientific = TRUE, digits = 3), 
                          " | SD = ", format(sd(samples), scientific = TRUE, digits = 3),
                          " | Det = ", format(det_value, scientific = TRUE, digits = 3)),
        x = "Parameter Value",
        y = "Density"
      ) +
      theme_minimal(base_size = 11) +
      theme(
        plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 9, color = "gray30")
      )
    
    plots_list[[paste0("cycle_", cycle)]] <- p
    
    filename <- file.path(OUT_FIG, paste0("psa_dist_", param_name, "_cycle", cycle, ".png"))
    try(ggsave(filename, p, width = 6, height = 4, dpi = 300))
  }
  
  log_info(paste("Created", length(plots_list), "temporal distribution plots"))
  
  return(invisible(plots_list))
}

# ============================================================================
# SCALAR DISTRIBUTION PLOTS EXPLANATION
# ============================================================================
#
# These histograms show PSA distributions for SCALAR parameters (single values
# that don't change over time, like costs, utilities, overall probabilities).
#
# SAME INTERPRETATION as temporal distributions, but for:
#   - Costs: c_Scarring, c_Single_Amput, etc.
#   - Utilities: u_IMD, u_Scarring, etc.
#   - Fixed probabilities: p_seq_overall, coverage rates
#   - Vaccine costs: c_MenABCWY, c_MenB, etc.
#
# KEY DIFFERENCES FROM TEMPORAL PLOTS:
#   - Single histogram (not one per cycle)
#   - Represents a constant value used throughout the model
#   - Often has tighter distributions (better data available)
#
# NOTE: With n=4 simulations, you'll only see 4 discrete bars.
#       Use n≥1000 for smooth, bell-shaped distributions.
# ============================================================================

# Plot distributions for scalar parameters
plot_psa_scalar_distributions <- function(psa_df, params_base, 
                                          param_names = NULL, 
                                          max_plots = 20) {
  
  log_info("Creating PSA distribution plots for scalar parameters...")
  
  # Get all scalar parameters (exclude those with numbers at end)
  all_params <- names(psa_df)[!grepl("\\d+$", names(psa_df))]
  
  # Filter to specified parameters or use all
  if (!is.null(param_names)) {
    params_to_plot <- intersect(param_names, all_params)
  } else {
    params_to_plot <- all_params
  }
  
  # Limit number of plots
  if (length(params_to_plot) > max_plots) {
    message(paste("[WARN] Too many parameters. Plotting only first", max_plots))
    params_to_plot <- params_to_plot[1:max_plots]
  }
  
  if (length(params_to_plot) == 0) {
    message("[WARN] No scalar parameters found to plot")
    return(invisible(NULL))
  }
  
  plots_list <- list()
  
  for (param in params_to_plot) {
    # Get PSA samples
    samples <- psa_df[[param]]
    
    # Get deterministic value
    det_value <- params_base[[param]]
    
    # Create histogram with density curve
    p <- ggplot(data.frame(value = samples), aes(x = value)) +
      geom_histogram(aes(y = after_stat(density)), 
                     bins = 40, 
                     fill = "#69b3a2", 
                     alpha = 0.7,
                     color = "white") +
      geom_density(color = "#404080", size = 1.2) +
      geom_vline(xintercept = det_value, 
                 color = "red", 
                 linetype = "dashed", 
                 size = 1) +
      labs(
        title = paste("PSA Distribution:", param),
        subtitle = paste0("Mean = ", round(mean(samples), 6), 
                          " | SD = ", round(sd(samples), 6),
                          " | Det = ", round(det_value, 6)),
        x = "Parameter Value",
        y = "Density"
      ) +
      theme_minimal(base_size = 11) +
      theme(
        plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 9, color = "gray30")
      )
    
    plots_list[[param]] <- p
    
    # Save individual plot
    filename <- file.path(OUT_FIG, paste0("psa_dist_", param, ".png"))
    try(ggsave(filename, p, width = 6, height = 4, dpi = 300))
  }
  
  log_info(paste("Created", length(plots_list), "scalar distribution plots"))
  
  return(invisible(plots_list))
}

# Create combined panel of distributions for key parameters
plot_psa_key_parameters_panel <- function(psa_df, params_base) {
  
  log_info("Creating key parameters panel...")
  
  # Define key parameters to show
  key_params <- c(
    "c_Scarring", "c_Single_Amput", "c_Multiple_Amput",
    "u_IMD", "u_Scarring", "u_Single_Amput",
    "p_seq_overall", "coverage_ACWY"
  )
  
  # Filter to available parameters
  available_params <- intersect(key_params, names(psa_df))
  
  if (length(available_params) == 0) {
    message("[WARN] No key parameters available for panel plot")
    return(invisible(NULL))
  }
  
  # Prepare data for plotting
  plot_data <- lapply(available_params, function(param) {
    data.frame(
      parameter = param,
      value = psa_df[[param]],
      deterministic = params_base[[param]]
    )
  }) %>% bind_rows()
  
  # Create faceted plot
  p <- ggplot(plot_data, aes(x = value)) +
    geom_histogram(aes(y = after_stat(density)), 
                   bins = 30, 
                   fill = "#69b3a2", 
                   alpha = 0.7,
                   color = "white") +
    geom_vline(aes(xintercept = deterministic), 
               color = "red", 
               linetype = "dashed", 
               size = 0.8) +
    facet_wrap(~parameter, scales = "free", ncol = 3) +
    labs(
      title = "PSA Distributions for Key Parameters",
      subtitle = "Red dashed line = deterministic value",
      x = "Parameter Value",
      y = "Density"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      strip.text = element_text(face = "bold", size = 9),
      strip.background = element_rect(fill = "gray90", color = NA)
    )
  
  print(p)
  
  filename <- file.path(OUT_FIG, "psa_key_parameters_panel.png")
  try(ggsave(filename, p, width = 12, height = 8, dpi = 300))
  
  log_info("Key parameters panel created")
  
  return(invisible(p))
}

# Compare deterministic vs PSA mean for all parameters
plot_psa_validation <- function(psa_df, params_base) {
  
  log_info("Creating PSA validation plot...")
  
  # Get scalar parameters only
  scalar_params <- names(psa_df)[!grepl("\\d+$", names(psa_df))]
  
  # Calculate PSA means (safely handle non-numeric values)
  psa_means <- sapply(scalar_params, function(param) {
    val <- mean(psa_df[[param]], na.rm = TRUE)
    if (is.finite(val)) val else NA
  })
  
  # Get deterministic values (safely handle non-numeric values)
  det_values <- sapply(scalar_params, function(param) {
    val <- params_base[[param]]
    if (is.null(val)) return(NA)
    if (is.numeric(val) && length(val) == 1) return(as.numeric(val))
    return(NA)
  })
  
  # Create comparison dataframe and remove NAs
  comparison <- data.frame(
    parameter = scalar_params,
    deterministic = det_values,
    psa_mean = psa_means,
    stringsAsFactors = FALSE
  )
  
  # Remove rows with NA or non-finite values
  comparison <- comparison[is.finite(comparison$deterministic) & 
                             is.finite(comparison$psa_mean) &
                             comparison$deterministic != 0, ]
  
  if (nrow(comparison) == 0) {
    message("[WARN] No valid parameters for validation plot")
    return(invisible(NULL))
  }
  
  # Calculate ratio
  comparison$ratio <- comparison$psa_mean / comparison$deterministic
  
  # Plot
  p <- ggplot(comparison, aes(x = deterministic, y = psa_mean)) +
    geom_abline(intercept = 0, slope = 1, 
                color = "red", linetype = "dashed", size = 1) +
    geom_point(size = 3, alpha = 0.6, color = "#404080") +
    labs(
      title = "PSA Validation: Deterministic vs PSA Mean",
      subtitle = "Points should align with diagonal line",
      x = "Deterministic Value",
      y = "PSA Mean Value"
    ) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))
  
  print(p)
  
  filename <- file.path(OUT_FIG, "psa_validation.png")
  try(ggsave(filename, p, width = 7, height = 6, dpi = 300))
  
  # Print parameters with largest deviations
  comparison <- comparison %>%
    mutate(pct_diff = abs(ratio - 1) * 100) %>%
    arrange(desc(pct_diff))
  
  log_info("Top 5 parameters with largest PSA mean deviation:")
  print(head(comparison[, c("parameter", "deterministic", "psa_mean", "pct_diff")], 5))
  
  return(invisible(comparison))
}

# Master function to create all PSA distribution plots
create_all_psa_distribution_plots <- function(psa_df, params_base) {
  
  log_info("\n=== Creating All PSA Distribution Visualizations ===")
  
  # 1. Key parameters panel
  plot_psa_key_parameters_panel(psa_df, params_base)
  
  # 2. Validation plot
  plot_psa_validation(psa_df, params_base)
  
  # 3. Temporal summaries for infection probabilities
  for (param in c("p_B", "p_C", "p_W", "p_Y")) {
    param_cols <- grep(paste0("^", param, "\\d+$"), names(psa_df), value = TRUE)
    if (length(param_cols) > 0) {
      plot_psa_temporal_summary(
        psa_df, params_base, param,
        plot_title = paste("Serogroup", gsub("p_", "", param), "Infection Probability")
      )
    }
  }
  
  # 4. Temporal distributions for selected cycles
  for (param in c("p_B", "p_C")) {
    param_cols <- grep(paste0("^", param, "\\d+$"), names(psa_df), value = TRUE)
    if (length(param_cols) > 0) {
      plot_psa_temporal_distributions(
        psa_df, params_base, param,
        cycles_to_plot = c(1, 10, 20, 40, 60, 89)
      )
    }
  }
  
  # 5. Selected scalar distributions
  important_scalars <- c(
    "c_Scarring", "c_Single_Amput", "c_Multiple_Amput",
    "u_IMD", "u_Scarring", "u_Single_Amput",
    "p_seq_overall", "c_MenABCWY", "c_MenB"
  )
  
  available_scalars <- intersect(important_scalars, names(psa_df))
  if (length(available_scalars) > 0) {
    plot_psa_scalar_distributions(psa_df, params_base, available_scalars)
  }
  
  log_info("=== All PSA distribution plots created ===\n")
  
  return(invisible(NULL))
}