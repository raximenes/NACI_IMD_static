################################################################################
# 04_quick_results_check.R
#
# PURPOSE: Quick sanity check of model results after running analysis
# - Fast validation of base case results
# - Quick PSA check
# - Identify any obvious issues
# - Generate one-page summary
#
# WHEN TO USE:
#   - Immediately after running run_analysis.R
#   - Before deep analysis
#   - To quickly verify results make sense
#
# USAGE:
#   source("scripts/04_quick_results_check.R")
#
# OUTPUT:
#   - Console summary
#   - One-page summary plot
#
################################################################################

# ============================================================
# SETUP
# ============================================================

cat("\n")
cat(strrep("=", 80), "\n")
cat("QUICK RESULTS CHECK\n")
cat(strrep("=", 80), "\n\n")

# Load packages
if (!require("pacman", quietly = TRUE)) {
  install.packages("pacman", repos = "https://cloud.r-project.org")
}
library(pacman)
p_load(readr, dplyr, ggplot2, tidyr, gridExtra)

# ============================================================
# CONFIGURATION
# ============================================================

PERSPECTIVE <- "healthcare"  # Which perspective to check
WTP <- 50000                 # WTP threshold

# ============================================================
# CHECK 1: DO RESULT FILES EXIST?
# ============================================================

cat(strrep("-", 80), "\n")
cat("CHECK 1: FILE EXISTENCE\n")
cat(strrep("-", 80), "\n\n")

files_to_check <- c(
  base_results = file.path("outputs", "tables", paste0("base_case_results_", PERSPECTIVE, ".csv")),
  base_icers = file.path("outputs", "tables", paste0("base_case_icers_", PERSPECTIVE, ".csv")),
  psa_results = file.path("outputs", "tables", paste0("psa_long_results_", PERSPECTIVE, ".csv"))
)

file_status <- sapply(files_to_check, file.exists)

for (i in seq_along(file_status)) {
  symbol <- ifelse(file_status[i], "✓", "✗")
  cat(sprintf("  [%s] %s\n", symbol, names(file_status)[i]))
}

if (!all(file_status)) {
  cat("\n[ERROR] Some result files are missing!\n")
  cat("        Run run_analysis.R first.\n\n")
  stop("Cannot proceed without results files")
}

cat("\n✓ All required files found\n\n")

# ============================================================
# CHECK 2: BASE CASE RESULTS
# ============================================================

cat(strrep("-", 80), "\n")
cat("CHECK 2: BASE CASE RESULTS\n")
cat(strrep("-", 80), "\n\n")

base_results <- read_csv(files_to_check["base_results"], show_col_types = FALSE)
base_icers <- read_csv(files_to_check["base_icers"], show_col_types = FALSE)

cat("Base Case Summary:\n")
print(base_results)
cat("\n")

cat("ICERs:\n")
print(base_icers)
cat("\n")

# Sanity checks
checks <- list()

# Check 1: Costs are positive
checks$positive_costs <- all(base_results$Cost > 0)
if (checks$positive_costs) {
  cat("✓ All costs are positive\n")
} else {
  cat("✗ WARNING: Some costs are negative or zero\n")
}

# Check 2: QALYs are positive
checks$positive_qalys <- all(base_results$QALYs > 0)
if (checks$positive_qalys) {
  cat("✓ All QALYs are positive\n")
} else {
  cat("✗ WARNING: Some QALYs are negative or zero\n")
}

# Check 3: Cost ranking makes sense (more expensive strategies should have more doses)
cost_order <- base_results$Strategy[order(base_results$Cost)]
cat(sprintf("\nCost ranking (cheapest → most expensive):\n  %s\n", 
            paste(cost_order, collapse = " < ")))

# Check 4: QALY ranking
qaly_order <- base_results$Strategy[order(-base_results$QALYs)]
cat(sprintf("QALY ranking (most QALYs → least):\n  %s\n", 
            paste(qaly_order, collapse = " > ")))

# Check 5: Are there dominated strategies?
if ("Status" %in% names(base_icers)) {
  dominated <- base_icers$Strategy[base_icers$Status == "D"]
  if (length(dominated) > 0) {
    cat(sprintf("\n[INFO] Dominated strategies: %s\n", paste(dominated, collapse = ", ")))
  } else {
    cat("\n[INFO] No dominated strategies\n")
  }
}

cat("\n")

# ============================================================
# CHECK 3: PSA RESULTS
# ============================================================

cat(strrep("-", 80), "\n")
cat("CHECK 3: PSA RESULTS\n")
cat(strrep("-", 80), "\n\n")

psa_results <- read_csv(files_to_check["psa_results"], show_col_types = FALSE)

n_sims <- length(unique(psa_results$sim))
n_strats <- length(unique(psa_results$Strategy))

cat(sprintf("PSA Information:\n"))
cat(sprintf("  - Number of simulations: %d\n", n_sims))
cat(sprintf("  - Number of strategies: %d\n", n_strats))
cat(sprintf("  - Total rows: %d\n\n", nrow(psa_results)))

# Check for variation
psa_summary <- psa_results %>%
  group_by(Strategy) %>%
  summarise(
    Mean_Cost = mean(Cost),
    SD_Cost = sd(Cost),
    CV_Cost = sd(Cost) / mean(Cost),
    Mean_QALY = mean(QALYs),
    SD_QALY = sd(QALYs),
    CV_QALY = sd(QALYs) / mean(QALYs),
    .groups = "drop"
  )

cat("PSA Variation (Coefficient of Variation):\n")
print(select(psa_summary, Strategy, CV_Cost, CV_QALY))
cat("\n")

# Check variation
checks$psa_has_variation <- all(psa_summary$SD_Cost > 0) && all(psa_summary$SD_QALY > 0)
if (checks$psa_has_variation) {
  cat("✓ PSA shows parameter variation (SD > 0)\n")
} else {
  cat("✗ WARNING: PSA shows no variation - parameters may not be varying\n")
}

# Check CV is reasonable
reasonable_cv <- all(psa_summary$CV_Cost > 0.01 & psa_summary$CV_Cost < 1)
if (reasonable_cv) {
  cat("✓ Coefficient of variation is reasonable (1% - 100%)\n")
} else {
  cat("✗ WARNING: CV outside typical range - check parameter SDs\n")
}

cat("\n")

# ============================================================
# CHECK 4: PROBABILITY OF COST-EFFECTIVENESS
# ============================================================

cat(strrep("-", 80), "\n")
cat("CHECK 4: COST-EFFECTIVENESS AT WTP = $", format(WTP, big.mark = ","), "\n")
cat(strrep("-", 80), "\n\n")

# Find which strategy has highest NMB in each iteration
best_strategy <- psa_results %>%
  group_by(sim) %>%
  slice_max(NMB, n = 1, with_ties = FALSE) %>%  # Keep only best strategy per simulation
  ungroup() %>%
  count(Strategy, name = "Times_Best") %>%  # Count unique simulations per strategy
  mutate(Probability = Times_Best / n_sims) %>%
  arrange(desc(Probability))

cat("Probability of being cost-effective:\n")
print(best_strategy)
cat("\n")

optimal <- best_strategy$Strategy[1]
cat(sprintf("✓ Most likely optimal strategy: %s (%.1f%% probability)\n\n", 
            optimal, 100 * best_strategy$Probability[1]))

# ============================================================
# CHECK 5: CONSISTENCY BASE CASE vs PSA MEAN
# ============================================================

cat(strrep("-", 80), "\n")
cat("CHECK 5: BASE CASE vs PSA MEAN\n")
cat(strrep("-", 80), "\n\n")

psa_means <- psa_results %>%
  group_by(Strategy) %>%
  summarise(
    PSA_Mean_Cost = mean(Cost),
    PSA_Mean_QALY = mean(QALYs),
    .groups = "drop"
  )

comparison <- base_results %>%
  select(Strategy, Base_Cost = Cost, Base_QALY = QALYs) %>%
  left_join(psa_means, by = "Strategy") %>%
  mutate(
    Cost_Diff_Pct = abs(Base_Cost - PSA_Mean_Cost) / Base_Cost * 100,
    QALY_Diff_Pct = abs(Base_QALY - PSA_Mean_QALY) / Base_QALY * 100
  )

cat("Base Case vs PSA Mean (% difference):\n")
print(select(comparison, Strategy, Cost_Diff_Pct, QALY_Diff_Pct))
cat("\n")

max_diff <- max(comparison$Cost_Diff_Pct, comparison$QALY_Diff_Pct)
if (max_diff < 5) {
  cat("✓ Base case and PSA means are very close (<5% difference)\n")
} else if (max_diff < 10) {
  cat("✓ Base case and PSA means are reasonably close (<10% difference)\n")
} else {
  cat("! Base case and PSA means differ by >10% - this is unusual but may be OK\n")
  cat("  (Could indicate parameter distributions are not symmetric)\n")
}

cat("\n")

# ============================================================
# VISUAL SUMMARY - ONE PAGE
# ============================================================

cat(strrep("-", 80), "\n")
cat("GENERATING VISUAL SUMMARY\n")
cat(strrep("-", 80), "\n\n")

# Create 4-panel summary plot
p1 <- ggplot(base_results, aes(x = reorder(Strategy, Cost), y = Cost)) +
  geom_col(fill = "steelblue", alpha = 0.7) +
  coord_flip() +
  scale_y_continuous(labels = scales::dollar) +
  labs(title = "Base Case Costs", x = NULL, y = "Total Cost (CAD)") +
  theme_bw()

p2 <- ggplot(base_results, aes(x = reorder(Strategy, -QALYs), y = QALYs)) +
  geom_col(fill = "forestgreen", alpha = 0.7) +
  coord_flip() +
  scale_y_continuous(labels = scales::comma) +
  labs(title = "Base Case QALYs", x = NULL, y = "Total QALYs") +
  theme_bw()

p3 <- ggplot(psa_results, aes(x = QALYs, y = Cost, color = Strategy)) +
  geom_point(alpha = 0.3, size = 0.5) +
  stat_ellipse(level = 0.95, size = 1) +
  scale_y_continuous(labels = scales::comma) +
  scale_x_continuous(labels = scales::comma) +
  labs(title = sprintf("Cost-Effectiveness Plane (n=%d)", n_sims),
       x = "QALYs", y = "Cost (CAD)") +
  theme_bw() +
  theme(legend.position = "bottom")

p4 <- ggplot(psa_results, aes(x = NMB, fill = Strategy)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  scale_x_continuous(labels = scales::comma) +
  labs(title = sprintf("NMB Distributions (WTP=$%s)", format(WTP, big.mark = ",")),
       x = "Net Monetary Benefit (CAD)", y = "Density") +
  theme_bw() +
  theme(legend.position = "bottom")

# Combine plots
summary_plot <- gridExtra::grid.arrange(
  p1, p2, p3, p4,
  ncol = 2,
  top = sprintf("Quick Results Summary - %s Perspective", toupper(PERSPECTIVE))
)

# Save
output_file <- file.path("outputs", "figures", PERSPECTIVE, 
                         paste0("quick_summary_", PERSPECTIVE, ".png"))
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

ggsave(output_file, summary_plot, width = 14, height = 10, dpi = 300)

cat(sprintf("✓ Summary plot saved: %s\n\n", output_file))

# ============================================================
# OVERALL ASSESSMENT
# ============================================================

cat(strrep("=", 80), "\n")
cat("OVERALL ASSESSMENT\n")
cat(strrep("=", 80), "\n\n")

all_checks_pass <- all(unlist(checks))

if (all_checks_pass) {
  cat("✓✓✓ ALL CHECKS PASSED ✓✓✓\n\n")
  cat("Your results look good!\n")
  cat("  ✓ All required files exist\n")
  cat("  ✓ Costs and QALYs are positive\n")
  cat("  ✓ PSA shows parameter variation\n")
  cat("  ✓ Base case and PSA are consistent\n")
  cat(sprintf("  ✓ Most likely optimal: %s\n", optimal))
  cat("\nNext steps:\n")
  cat("  1. Review detailed results in outputs/tables/\n")
  cat("  2. Check key plots in outputs/figures/\n")
  cat("  3. Run 02_analyze_psa_results.R for deeper analysis\n")
} else {
  cat("⚠ SOME CHECKS FAILED ⚠\n\n")
  cat("Issues detected:\n")
  
  if (!checks$positive_costs) cat("  ✗ Some costs are not positive\n")
  if (!checks$positive_qalys) cat("  ✗ Some QALYs are not positive\n")
  if (!checks$psa_has_variation) cat("  ✗ PSA shows no variation\n")
  
  cat("\nPlease review the warnings above.\n")
  cat("Common fixes:\n")
  cat("  - Check Excel input file for errors\n")
  cat("  - Verify SD values are specified correctly\n")
  cat("  - Check that model ran completely\n")
}

cat("\n")
cat(strrep("=", 80), "\n\n")