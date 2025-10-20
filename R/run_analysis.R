# scripts/run_analysis.R
# Source all numbered R files and run main()

install.packages(c("pacman","here","readxl","ggplot2","dplyr","tidyr","purrr","stringr",
                   "tibble","doParallel","foreach","parallel","diagram","reshape2","readr"))
pacman::p_load_gh("DARTH-git/darthtools"); pacman::p_load(dampack)


files <- list.files("R", pattern = "^[0-9]{2}_.*\\.R$", full.names = TRUE)
files <- files[order(files)]  # ensure 01..08 order
invisible(lapply(files, source))

tryCatch({
  results <- main()
  message("\n[INFO] === SUCCESS: All analyses completed without errors ===")
}, error = function(e) {
  message("\n[ERROR] Analysis failed: ", e$message)
  traceback()
})




