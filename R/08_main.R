## ======================
## 8) Main Analysis Execution (UPDATED)
## ======================

main <- function() {
  log_info("=== Starting IMD Vaccination Cost-Utility Model ===")
  log_info(paste("Working directory:", getwd()))
  
  # Load base-case parameters
  params_bc <- get_base_params()
  
  ## --- BASE CASE ANALYSIS ---
  log_info("\n=== Running Base Case Analysis ===")
  res_det_list <- eval_all_strategies(params_bc)
  det <- summarize_det(res_det_list, comparator = "MenC")
  
  log_info("\nBase Case Results:")
  print(det$table)
  print(det$icers)
  print(det$incremental)
  
  # Plot ICER frontier
  p_icer <- try(plot(det$icers), silent = TRUE)
  if (inherits(p_icer, "ggplot")) {
    try(ggsave(file.path(OUT_FIG_DET, "base_case_icers_frontier.png"), p_icer, width = 7, height = 5, dpi = 300))
    try(print(p_icer))
  }
  
  ## --- PROBABILISTIC SENSITIVITY ANALYSIS (PSA) ---
  log_info("\n=== Running Probabilistic Sensitivity Analysis (PSA) ===")
  psa_output <- run_psa(params_bc, n_sim = 1000, wtp = 50000)
  psa_long <- psa_output$results
  psa_df <- psa_output$psa_df
  
  objs <- make_psa_objects(psa_long, strategies = v_strats, wtp_vec = v_wtp)
  
  # Plot CEAC
  p_ceac <- try(plot(objs$ceac), silent = TRUE)
  if (inherits(p_ceac, "ggplot")) try(ggsave(file.path(OUT_FIG_KEY, "psa_ceac.png"), p_ceac, width = 7, height = 5, dpi = 300))
  
  
  # Plot EVPI
  p_evpi <- try(plot(objs$evpi), silent = TRUE)
  if (inherits(p_evpi, "ggplot")) try(ggsave(file.path(OUT_FIG_KEY, "psa_evpi.png"), p_evpi, width = 7, height = 5, dpi = 300))
  
  ## --- ONE-WAY SENSITIVITY ANALYSIS (OWSA) ---
  log_info("\n=== Running One-Way Sensitivity Analysis (OWSA) ===")
  owsa_out <- execute_owsa(params_bc, comparator = "MenC")
  
  # Generate tornado plots (choose plot_type: "incremental", "nmb", "icer", or "both")
  create_owsa_plots(
    owsa_out, 
    det, 
    strategies = v_strats,
    wtp = 50000,
    plot_type = "both"  # Creates both incremental and NMB plots
  )
 
  ## --- SAVE RESULTS ---
  log_info("\n=== Saving Results ===")
  readr::write_csv(det$table,        file.path(OUT_TAB, "base_case_results.csv"))
  readr::write_csv(as.data.frame(det$icers), file.path(OUT_TAB, "base_case_icers.csv"))
  readr::write_csv(det$incremental,  file.path(OUT_TAB, "base_case_incremental_epi.csv"))
  readr::write_csv(psa_long,         file.path(OUT_TAB, "psa_long_results.csv"))
  if (is.list(owsa_out)) {
    readr::write_csv(owsa_out$owsa_Cost,  file.path(OUT_TAB, "owsa_output_cost.csv"))
    readr::write_csv(owsa_out$owsa_QALYs, file.path(OUT_TAB, "owsa_output_qalys.csv"))
  }
  
  log_info("\n=== Analysis Complete ===")
  log_info(paste("Figures:", OUT_FIG))
  log_info(paste("Tables:",  OUT_TAB))
  
  invisible(list(base_case = det, psa = psa_long, psa_parameters = psa_df, psa_objects = objs, owsa = owsa_out))
}