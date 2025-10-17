# 01_setup.R
# Packages, theming, logging, and project paths

if (!require("pacman")) install.packages("pacman", repos = "https://cloud.r-project.org")
library(pacman)

# Core libs
p_load(here, readxl, ggplot2, dplyr, tidyr, purrr, stringr, tibble,
       doParallel, foreach, parallel, diagram, reshape2, readr)
p_load_gh("DARTH-git/darthtools")
p_load(dampack)

theme_set(ggplot2::theme_minimal())

# ---- Logging helpers ----
log_info <- function(msg)  message(paste0("[INFO] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " - ", msg))
log_warn <- function(msg)  warning(paste0("[WARN] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " - ", msg), call. = FALSE, immediate. = TRUE)
log_error <- function(msg) stop(paste0("[ERROR] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " - ", msg), call. = FALSE)

# ---- Project paths (auto-detected from .Rproj) ----
PROJECT_DIR <- here::here()
IN_RAW  <- file.path(PROJECT_DIR, "data", "raw")
IN_PROC <- file.path(PROJECT_DIR, "data", "processed")
OUT_FIG <- file.path(PROJECT_DIR, "outputs", "figures")
OUT_TAB <- file.path(PROJECT_DIR, "outputs", "tables")

dir.create(IN_PROC, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_FIG, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_TAB, recursive = TRUE, showWarnings = FALSE)

log_info(paste("Project root:", PROJECT_DIR))
log_info(paste("R version:", R.version.string))
