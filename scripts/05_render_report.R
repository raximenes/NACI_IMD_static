################################################################################
# 05_render_report.R
#
# PURPOSE: Render the RMarkdown report for specified perspectives.
#
# USAGE:
#   source("scripts/05_render_report.R")
#
# OUTPUT:
#   - Word/HTML reports in reports/
#
################################################################################

# ============================================================
# SETUP
# ============================================================

if (!require("pacman", quietly = TRUE)) {
  install.packages("pacman", repos = "https://cloud.r-project.org")
}
library(pacman)
p_load(rmarkdown, tools)

# ============================================================
# CONFIGURATION
# ============================================================

# Perspectives to generate reports for
PERSPECTIVES <- c("healthcare", "societal")

# Report template
RMD_FILE <- file.path("reports", "imd_analysis_report.Rmd")

# ============================================================
# RENDER LOOP
# ============================================================

if (!file.exists(RMD_FILE)) {
  stop("Report template not found: ", RMD_FILE)
}

cat("\n========================================\n")
cat("RENDERING REPORTS\n")
cat("========================================\n")

for (persp in PERSPECTIVES) {

  cat(sprintf("\nProcessing perspective: %s...\n", toupper(persp)))

  # Check if results exist
  results_file <- file.path("outputs", "tables", paste0("base_case_results_", persp, ".csv"))
  if (!file.exists(results_file)) {
    cat(sprintf("[WARN] Results not found for %s. Skipping report generation.\n", persp))
    next
  }

  # Define output filename
  output_file <- paste0("IMD_Report_", tools::toTitleCase(persp), "_", format(Sys.Date(), "%Y-%m-%d"), ".docx")

  tryCatch({
    render(
      input = RMD_FILE,
      output_file = output_file,
      output_dir = "reports",
      params = list(perspective = persp),
      quiet = TRUE
    )
    cat(sprintf("[PASS] Report generated: reports/%s\n", output_file))
  }, error = function(e) {
    cat(sprintf("[FAIL] Failed to render report for %s: %s\n", persp, e$message))
  })
}

cat("\n========================================\n")
cat("REPORT GENERATION COMPLETE\n")
cat("========================================\n")
