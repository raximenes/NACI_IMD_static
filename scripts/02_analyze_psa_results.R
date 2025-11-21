################################################################################
# 02_analyze_psa_results.R
#
# PURPOSE: Comprehensive analysis of PSA results
# - Load and validate PSA results from run_analysis.R
# - Calculate summary statistics by strategy
# - Analyze parameter variation and impact
# - Generate diagnostic plots
# - Create correlation analysis
#
# PREREQUISITES:
#   - run_analysis.R must have been executed
#   - PSA results files must exist in outputs/tables/
#
# USAGE:
#   source("scripts/02_analyze_psa_results.R")
#
# OUTPUT:
#   - Analysis tables in outputs/tables/psa_analysis/
#   - Diagnostic plots in outputs/figures/[perspective]/PSA/
#
################################################################################

# ============================================================
# SETUP
# ============================================================

start_time <- Sys.time()

cat("\n")
cat(strrep("=", 80), "\n")
cat("PSA RESULTS ANALYSIS\n")
cat(strrep("=", 80), "\n\n")

# Load packages
if (!require("pacman", quietly = TRUE)) {
  install.packages("pacman", repos = "https://cloud.r-project.org")
}
library(pacman)
p_load(readr, dplyr, ggplot2, tidyr, corrplot, gridExtra)

# ============================================================
# CONFIGURATION
# ============================================================

PERSPECTIVES <- c("healthcare", "societal")  # Which perspectives to analyze
WTP_THRESHOLD <- 50000                        # WTP for NMB calculations

# Output directory
analysis_dir <- file.path("outputs", "tables", "psa_analysis")
dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# LOAD PSA RESULTS
# ============================================================

cat(strrep("-", 80), "\n")
cat("LOADING PSA RESULTS\n")
cat(strrep("-", 80), "\n\n")

psa_data <- list()
psa_samples_data <- NULL

for (persp in PERSPECTIVES) {
  psa_file <- file.path("outputs", "tables", paste0("psa_long_results_", persp, ".csv"))
  
  if (file.exists(psa_file)) {
    cat(sprintf("[%s] Loading PSA results...\n", toupper(persp)))
    psa_data[[persp]] <- read_csv(psa_file, show_col_types = FALSE)
    cat(sprintf("  ✓ Loaded %d rows (%d simulations)\n", 
                nrow(psa_data[[persp]]), 
                length(unique(psa_data[[persp]]$sim))))
  } else {
    cat(sprintf("[%s] PSA results not found: %s\n", toupper(persp), psa_file))
  }
}

if (length(psa_data) == 0) {
  stop("\n[ERROR] No PSA results found. Run run_analysis.R first.\n")
}

# Load PSA samples if available
psa_samples_file <- file.path("outputs", "tables", "psa_samples.rds")
if (file.exists(psa_samples_file)) {
  cat("\n[INFO] Loading PSA parameter samples...\n")
  psa_samples_data <- readRDS(psa_samples_file)
  cat("  ✓ Loaded parameter samples\n")
}

cat("\n")

# ============================================================
# ANALYSIS LOOP - ONE PERSPECTIVE AT A TIME
# ============================================================

for (persp in names(psa_data)) {
  
  cat(strrep("=", 80), "\n")
  cat(sprintf("ANALYZING: %s PERSPECTIVE\n", toupper(persp)))
  cat(strrep("=", 80), "\n\n")
  
  psa_results <- psa_data[[persp]]
  n_sims <- length(unique(psa_results$sim))
  n_strats <- length(unique(psa_results$Strategy))
  
  # Create perspective-specific plot directory
  plot_dir <- file.path("outputs", "figures", persp, "PSA")
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  
  # ==============================================================
  # 1. SUMMARY STATISTICS
  # ==============================================================
  
  cat(strrep("-", 80), "\n")
  cat("1. SUMMARY STATISTICS\n")
  cat(strrep("-", 80), "\n\n")
  
  psa_summary <- psa_results %>%
    group_by(Strategy) %>%
    summarise(
      N = n(),
      
      # Cost statistics
      Cost_Mean = mean(Cost, na.rm = TRUE),
      Cost_SD = sd(Cost, na.rm = TRUE),
      Cost_Median = median(Cost, na.rm = TRUE),
      Cost_Q025 = quantile(Cost, 0.025, na.rm = TRUE),
      Cost_Q975 = quantile(Cost, 0.975, na.rm = TRUE),
      Cost_CV = sd(Cost, na.rm = TRUE) / mean(Cost, na.rm = TRUE),
      
      # QALY statistics
      QALY_Mean = mean(QALYs, na.rm = TRUE),
      QALY_SD = sd(QALYs, na.rm = TRUE),
      QALY_Median = median(QALYs, na.rm = TRUE),
      QALY_Q025 = quantile(QALYs, 0.025, na.rm = TRUE),
      QALY_Q975 = quantile(QALYs, 0.975, na.rm = TRUE),
      QALY_CV = sd(QALYs, na.rm = TRUE) / mean(QALYs, na.rm = TRUE),
      
      # NMB statistics
      NMB_Mean = mean(NMB, na.rm = TRUE),
      NMB_SD = sd(NMB, na.rm = TRUE),
      NMB_Median = median(NMB, na.rm = TRUE),
      NMB_Q025 = quantile(NMB, 0.025, na.rm = TRUE),
      NMB_Q975 = quantile(NMB, 0.975, na.rm = TRUE),
      
      .groups = "drop"
    )
  
  print(psa_summary)
  
  # Save summary
  write_csv(psa_summary, 
            file.path(analysis_dir, paste0("psa_summary_", persp, ".csv")))
  cat(sprintf("\n✓ Saved: %s\n\n", 
              file.path(analysis_dir, paste0("psa_summary_", persp, ".csv"))))
  
  # ==============================================================
  # 2. INCREMENTAL ANALYSIS
  # ==============================================================
  
  cat(strrep("-", 80), "\n")
  cat("2. INCREMENTAL ANALYSIS (vs MenC)\n")
  cat(strrep("-", 80), "\n\n")
  
  # Calculate incremental outcomes for each iteration
  psa_incremental <- psa_results %>%
    select(sim, Strategy, Cost, QALYs, NMB) %>%
    pivot_wider(names_from = Strategy, values_from = c(Cost, QALYs, NMB)) %>%
    mutate(
      # Incremental vs MenC
      Inc_Cost_MenACWY = Cost_MenACWY - Cost_MenC,
      Inc_Cost_MenACWY_MenB = Cost_MenACWY_MenB - Cost_MenC,
      Inc_Cost_MenABCWY = Cost_MenABCWY - Cost_MenC,
      
      Inc_QALY_MenACWY = QALYs_MenACWY - QALYs_MenC,
      Inc_QALY_MenACWY_MenB = QALYs_MenACWY_MenB - QALYs_MenC,
      Inc_QALY_MenABCWY = QALYs_MenABCWY - QALYs_MenC,
      
      # ICERs
      ICER_MenACWY = Inc_Cost_MenACWY / Inc_QALY_MenACWY,
      ICER_MenACWY_MenB = Inc_Cost_MenACWY_MenB / Inc_QALY_MenACWY_MenB,
      ICER_MenABCWY = Inc_Cost_MenABCWY / Inc_QALY_MenABCWY
    )
  
  # Summary of incremental results
  incremental_summary <- psa_incremental %>%
    select(sim, starts_with("Inc_"), starts_with("ICER_")) %>%
    pivot_longer(cols = -sim, names_to = "Metric", values_to = "Value") %>%
    group_by(Metric) %>%
    summarise(
      Mean = mean(Value, na.rm = TRUE),
      SD = sd(Value, na.rm = TRUE),
      Median = median(Value, na.rm = TRUE),
      Q025 = quantile(Value, 0.025, na.rm = TRUE),
      Q975 = quantile(Value, 0.975, na.rm = TRUE),
      .groups = "drop"
    )
  
  print(incremental_summary)
  
  write_csv(incremental_summary, 
            file.path(analysis_dir, paste0("psa_incremental_", persp, ".csv")))
  cat(sprintf("\n✓ Saved: %s\n\n", 
              file.path(analysis_dir, paste0("psa_incremental_", persp, ".csv"))))
  
  # ==============================================================
  # 3. PROBABILITY OF COST-EFFECTIVENESS
  # ==============================================================
  
  cat(strrep("-", 80), "\n")
  cat("3. PROBABILITY OF COST-EFFECTIVENESS\n")
  cat(strrep("-", 80), "\n\n")
  
  # Calculate probability each strategy has highest NMB
  prob_ce <- psa_results %>%
    group_by(sim) %>%
    mutate(Best = Strategy[which.max(NMB)]) %>%
    ungroup() %>%
    group_by(Best) %>%
    summarise(
      N_Best = n(),
      Probability = N_Best / n_sims,
      .groups = "drop"
    ) %>%
    arrange(desc(Probability))
  
  cat(sprintf("At WTP = $%s per QALY:\n", format(WTP_THRESHOLD, big.mark = ",")))
  print(prob_ce)
  
  write_csv(prob_ce, 
            file.path(analysis_dir, paste0("psa_prob_ce_", persp, ".csv")))
  cat(sprintf("\n✓ Saved: %s\n\n", 
              file.path(analysis_dir, paste0("psa_prob_ce_", persp, ".csv"))))
  
  # ==============================================================
  # 4. DIAGNOSTIC PLOTS
  # ==============================================================
  
  cat(strrep("-", 80), "\n")
  cat("4. GENERATING DIAGNOSTIC PLOTS\n")
  cat(strrep("-", 80), "\n\n")
  
  # Plot 1: Cost distributions
  p1 <- ggplot(psa_results, aes(x = Cost, fill = Strategy)) +
    geom_density(alpha = 0.5) +
    scale_x_continuous(labels = scales::comma) +
    labs(
      title = sprintf("PSA Cost Distributions (%s perspective)", persp),
      subtitle = sprintf("n = %d iterations", n_sims),
      x = "Total Cost (CAD)",
      y = "Density"
    ) +
    theme_bw() +
    theme(legend.position = "bottom")
  
  ggsave(file.path(plot_dir, paste0("psa_cost_distributions_", persp, ".png")),
         p1, width = 10, height = 6, dpi = 300)
  cat("  ✓ Cost distributions\n")
  
  # Plot 2: QALY distributions
  p2 <- ggplot(psa_results, aes(x = QALYs, fill = Strategy)) +
    geom_density(alpha = 0.5) +
    scale_x_continuous(labels = scales::comma) +
    labs(
      title = sprintf("PSA QALY Distributions (%s perspective)", persp),
      subtitle = sprintf("n = %d iterations", n_sims),
      x = "Total QALYs",
      y = "Density"
    ) +
    theme_bw() +
    theme(legend.position = "bottom")
  
  ggsave(file.path(plot_dir, paste0("psa_qaly_distributions_", persp, ".png")),
         p2, width = 10, height = 6, dpi = 300)
  cat("  ✓ QALY distributions\n")
  
  # Plot 3: Cost-effectiveness plane
  p3 <- ggplot(psa_results, aes(x = QALYs, y = Cost, color = Strategy)) +
    geom_point(alpha = 0.3, size = 1) +
    stat_ellipse(level = 0.95, size = 1) +
    scale_y_continuous(labels = scales::comma) +
    scale_x_continuous(labels = scales::comma) +
    labs(
      title = sprintf("Cost-Effectiveness Plane (%s perspective)", persp),
      subtitle = sprintf("n = %d iterations | Ellipses = 95%% CI", n_sims),
      x = "QALYs",
      y = "Cost (CAD)"
    ) +
    theme_bw() +
    theme(legend.position = "bottom")
  
  ggsave(file.path(plot_dir, paste0("psa_ce_plane_", persp, ".png")),
         p3, width = 10, height = 8, dpi = 300)
  cat("  ✓ Cost-effectiveness plane\n")
  
  # Plot 4: NMB distributions
  p4 <- ggplot(psa_results, aes(x = NMB, fill = Strategy)) +
    geom_density(alpha = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
    scale_x_continuous(labels = scales::comma) +
    labs(
      title = sprintf("Net Monetary Benefit Distributions (%s perspective)", persp),
      subtitle = sprintf("WTP = $%s | n = %d iterations", 
                         format(WTP_THRESHOLD, big.mark = ","), n_sims),
      x = "Net Monetary Benefit (CAD)",
      y = "Density"
    ) +
    theme_bw() +
    theme(legend.position = "bottom")
  
  ggsave(file.path(plot_dir, paste0("psa_nmb_distributions_", persp, ".png")),
         p4, width = 10, height = 6, dpi = 300)
  cat("  ✓ NMB distributions\n")
  
  # Plot 5: ICER distributions
  icer_data <- psa_incremental %>%
    select(sim, starts_with("ICER_")) %>%
    pivot_longer(cols = starts_with("ICER_"), 
                 names_to = "Comparison", 
                 values_to = "ICER") %>%
    mutate(Comparison = gsub("ICER_", "", Comparison))
  
  # Remove infinite ICERs for plotting
  icer_data_finite <- icer_data %>%
    filter(is.finite(ICER) & abs(ICER) < 500000)
  
  p5 <- ggplot(icer_data_finite, aes(x = ICER, fill = Comparison)) +
    geom_histogram(alpha = 0.6, bins = 50) +
    geom_vline(xintercept = WTP_THRESHOLD, 
               linetype = "dashed", color = "red", size = 1) +
    facet_wrap(~Comparison, scales = "free_y", ncol = 1) +
    scale_x_continuous(labels = scales::comma) +
    labs(
      title = sprintf("ICER Distributions vs MenC (%s perspective)", persp),
      subtitle = sprintf("WTP threshold (red line) = $%s", 
                         format(WTP_THRESHOLD, big.mark = ",")),
      x = "ICER (CAD per QALY)",
      y = "Frequency"
    ) +
    theme_bw() +
    theme(legend.position = "none")
  
  ggsave(file.path(plot_dir, paste0("psa_icer_distributions_", persp, ".png")),
         p5, width = 10, height = 10, dpi = 300)
  cat("  ✓ ICER distributions\n")
  
  cat(sprintf("\n✓ All plots saved to: %s\n\n", plot_dir))
  
  # ==============================================================
  # 5. VARIATION ANALYSIS
  # ==============================================================
  
  cat(strrep("-", 80), "\n")
  cat("5. COEFFICIENT OF VARIATION (CV) ANALYSIS\n")
  cat(strrep("-", 80), "\n\n")
  
  cv_analysis <- psa_summary %>%
    select(Strategy, Cost_CV, QALY_CV) %>%
    mutate(
      Cost_CV_Status = case_when(
        Cost_CV < 0.05 ~ "Very Low",
        Cost_CV < 0.15 ~ "Low-Moderate",
        Cost_CV < 0.30 ~ "Moderate (typical)",
        Cost_CV < 0.50 ~ "High",
        TRUE ~ "Very High"
      ),
      QALY_CV_Status = case_when(
        QALY_CV < 0.01 ~ "Very Low (typical for QALYs)",
        QALY_CV < 0.05 ~ "Low",
        QALY_CV < 0.15 ~ "Moderate",
        TRUE ~ "High"
      )
    )
  
  cat("Coefficient of Variation indicates the relative uncertainty:\n")
  cat("  CV < 0.15 = Low uncertainty (typical for well-defined parameters)\n")
  cat("  CV 0.15-0.30 = Moderate uncertainty (expected for costs)\n")
  cat("  CV > 0.30 = High uncertainty (may need investigation)\n\n")
  
  print(cv_analysis)
  
  write_csv(cv_analysis, 
            file.path(analysis_dir, paste0("psa_cv_analysis_", persp, ".csv")))
  cat(sprintf("\n✓ Saved: %s\n\n", 
              file.path(analysis_dir, paste0("psa_cv_analysis_", persp, ".csv"))))
  
} # End perspective loop

# ============================================================
# PARAMETER CORRELATION ANALYSIS (if samples available)
# ============================================================

if (!is.null(psa_samples_data)) {
  
  cat(strrep("=", 80), "\n")
  cat("PARAMETER CORRELATION ANALYSIS\n")
  cat(strrep("=", 80), "\n\n")
  
  # This analysis uses the first perspective's results
  first_persp <- names(psa_data)[1]
  psa_results <- psa_data[[first_persp]]
  
  cat(sprintf("[INFO] Using %s perspective data\n", first_persp))
  cat("[INFO] Extracting parameter values from PSA samples...\n\n")
  
  # Extract scalar parameter values
  psa_samples <- psa_samples_data$samples
  n_sims <- length(psa_samples)
  
  # Get parameter names
  param_names <- names(psa_samples[[1]])
  
  # Extract only scalar, numeric parameters that vary
  param_values <- list()
  for (param in param_names) {
    values <- sapply(1:n_sims, function(i) {
      val <- psa_samples[[i]][[param]]
      if (is.null(val) || !is.numeric(val)) return(NA)
      if (length(val) > 1) return(val[1])  # Take first value for vectors
      return(val)
    })
    
    # Only include if parameter varies
    # Add robust SD check
    param_sd <- sd(values, na.rm = TRUE)
    
    # Check if SD is valid and > 0
    if (!is.na(param_sd) && !is.nan(param_sd) && is.finite(param_sd) && param_sd > 1e-10) {
      param_values[[param]] <- values
    }
  }
  
  if (length(param_values) > 0) {
    cat(sprintf("✓ Extracted %d varying parameters\n\n", length(param_values)))
    
    # Merge with results for one strategy
    example_strategy <- "MenACWY"
    
    param_df <- as.data.frame(param_values)
    param_df$sim <- 1:n_sims
    
    merged_data <- psa_results %>%
      filter(Strategy == example_strategy) %>%
      left_join(param_df, by = "sim")
    
    # Select key parameters for correlation
    key_params <- intersect(
      c("coverage_ABCWY", "coverage_ACWY", "coverage_C", 
        "p_seq_overall", "u_Healthy", "u_IMD"),
      names(param_values)
    )
    
    if (length(key_params) >= 3) {
      cor_vars <- c(key_params[1:min(8, length(key_params))], "Cost", "QALYs", "NMB")
      cor_data <- merged_data %>%
        select(all_of(cor_vars)) %>%
        na.omit()
      
      # Calculate correlation matrix
      cor_matrix <- cor(cor_data)
      
      # Save correlation plot
      png(file.path(analysis_dir, "parameter_correlations.png"),
          width = 10, height = 10, units = "in", res = 300)
      corrplot(cor_matrix, method = "color", type = "upper",
               tl.col = "black", tl.srt = 45,
               title = sprintf("Parameter Correlations - %s", example_strategy),
               mar = c(0,0,2,0))
      dev.off()
      
      cat("✓ Correlation analysis complete\n")
      cat(sprintf("✓ Saved: %s\n\n", 
                  file.path(analysis_dir, "parameter_correlations.png")))
    }
  }
}

# ============================================================
# SUMMARY
# ============================================================

end_time <- Sys.time()
elapsed <- difftime(end_time, start_time, units = "secs")

cat("\n")
cat(strrep("=", 80), "\n")
cat("ANALYSIS COMPLETE\n")
cat(strrep("=", 80), "\n\n")

cat("Files created:\n")
cat(sprintf("  - Summary tables: %s\n", analysis_dir))
for (persp in names(psa_data)) {
  cat(sprintf("  - %s plots: outputs/figures/%s/PSA/\n", persp, persp))
}

cat(sprintf("\nTotal time: %.1f seconds\n", elapsed))
cat(strrep("=", 80), "\n\n")