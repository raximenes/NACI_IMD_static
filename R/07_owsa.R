
## ======================
## 7) OWSA Functions (CORRECTED)
## ======================
# FUNCTION 1: Define parameter ranges for the OWSA.
create_owsa_ranges <- function(params) {
  rng <- list()
  # Helper function to add a new parameter range
  add_range <- function(param_name, min_val, max_val) {
    rng[[length(rng) + 1]] <<- tibble::tibble(pars = param_name, min = min_val, max = max_val)
  }
  
  # Add all parameter ranges to test here
  add_range("coverage_ABCWY", max(0, params$coverage_ABCWY*0.8), min(1, params$coverage_ABCWY*1.2))
  add_range("coverage_ACWY",  max(0, params$coverage_ACWY*0.8),  min(1, params$coverage_ACWY*1.2))
  add_range("coverage_C",     max(0, params$coverage_C*0.8),     min(1, params$coverage_C*1.2))
  add_range("coverage_B",     max(0, params$coverage_B*0.8),     min(1, params$coverage_B*1.2))
  add_range("mult_p_B", 0.5, 1.5); add_range("mult_p_C", 0.5, 1.5); add_range("mult_p_W", 0.5, 1.5); add_range("mult_p_Y", 0.5, 1.5)
  add_range("mult_cfr_B", 0.5, 1.5); add_range("mult_cfr_C", 0.5, 1.5); add_range("mult_cfr_W", 0.5, 1.5); add_range("mult_cfr_Y", 0.5, 1.5)
  add_range("ve_scale_B", 0.80, 1.20); add_range("ve_scale_C", 0.80, 1.20); add_range("ve_scale_W", 0.80, 1.20); add_range("ve_scale_Y", 0.80, 1.20)
  add_range("mult_bg_mort", 0.8, 1.2)
  add_range("mult_c_IMD",   0.5, 1.5)
  add_range("p_seq_overall", 0.20, 0.30)
  for (nm in c("u_IMD","u_Scarring","u_Single_Amput","u_Multiple_Amput","u_Neuro_Disability", "u_Hearing_Loss","u_Renal_Failure","u_Seizure","u_Paralysis")) {
    base <- params[[nm]]; add_range(nm, max(0.01, base-0.05), min(0.99, base+0.05))
  }
  for (nm in c("c_Scarring","c_Single_Amput","c_Multiple_Amput","c_Neuro_Disab", "c_Hearing_Loss","c_Renal_Failure","c_Seizure","c_Paralysis")) {
    base <- params[[nm]]; add_range(nm, max(0, base*0.7), base*1.3)
  }
  add_range("c_admin", max(0, params$c_admin*0.8), params$c_admin*1.2)
  
  # Return a standard data.frame
  return(as.data.frame(dplyr::bind_rows(rng)))
}


# FUNCTION 2: The model function that we will now call.
owsa_model_function <- function(param_name, param_value, all_params) {
  p <- all_params
  p[[param_name]] <- param_value
  results_list <- eval_all_strategies(p)
  costs <- sapply(results_list, function(x) x$cost)
  qalys <- sapply(results_list, function(x) x$qalys)
  # Return a simple data frame with the results
  return(data.frame(Strategy = names(costs), Cost = costs, QALYs = qalys))
}

# Execute OWSA with incremental outcomes
execute_owsa <- function(params_basecase, comparator = "MenC", nsamp = 25) {
  log_info("Running One-Way Sensitivity Analysis (OWSA)...")
  
  ranges_df <- create_owsa_ranges(params_basecase)
  
  # Store results for each parameter
  full_owsa_results <- list()
  
  # Loop through each parameter
  for (i in 1:nrow(ranges_df)) {
    param_name <- ranges_df$pars[i]
    min_val <- ranges_df$min[i]
    max_val <- ranges_df$max[i]
    
    log_info(paste("Testing parameter:", param_name))
    
    # Run model at minimum value
    res_min <- owsa_model_function(param_name, min_val, params_basecase)
    
    # Run model at maximum value
    res_max <- owsa_model_function(param_name, max_val, params_basecase)
    
    # Store results
    full_owsa_results[[param_name]] <- list(min = res_min, max = res_max)
  }
  
  # Format output with INCREMENTAL values vs comparator
  owsa_cost_df <- do.call(rbind, lapply(names(full_owsa_results), function(name) {
    min_data <- full_owsa_results[[name]]$min
    max_data <- full_owsa_results[[name]]$max
    
    # Get comparator costs
    comp_cost_min <- min_data$Cost[min_data$Strategy == comparator]
    comp_cost_max <- max_data$Cost[max_data$Strategy == comparator]
    
    data.frame(
      pars = name,
      par_value_min = ranges_df$min[ranges_df$pars == name],
      par_value_max = ranges_df$max[ranges_df$pars == name],
      strategy = min_data$Strategy,
      Cost_min = min_data$Cost,
      Cost_max = max_data$Cost,
      # INCREMENTAL COSTS
      Inc_Cost_min = min_data$Cost - comp_cost_min,
      Inc_Cost_max = max_data$Cost - comp_cost_max,
      stringsAsFactors = FALSE
    )
  }))
  
  owsa_qaly_df <- do.call(rbind, lapply(names(full_owsa_results), function(name) {
    min_data <- full_owsa_results[[name]]$min
    max_data <- full_owsa_results[[name]]$max
    
    # Get comparator QALYs
    comp_qaly_min <- min_data$QALYs[min_data$Strategy == comparator]
    comp_qaly_max <- max_data$QALYs[max_data$Strategy == comparator]
    
    data.frame(
      pars = name,
      par_value_min = ranges_df$min[ranges_df$pars == name],
      par_value_max = ranges_df$max[ranges_df$pars == name],
      strategy = min_data$Strategy,
      QALYs_min = min_data$QALYs,
      QALYs_max = max_data$QALYs,
      # INCREMENTAL QALYs
      Inc_QALYs_min = min_data$QALYs - comp_qaly_min,
      Inc_QALYs_max = max_data$QALYs - comp_qaly_max,
      stringsAsFactors = FALSE
    )
  }))
  
  log_info("OWSA complete")
  
  return(list(
    owsa_Cost = owsa_cost_df,
    owsa_QALYs = owsa_qaly_df,
    comparator = comparator
  ))
}

# Create OWSA tornado plots with INCREMENTAL values and NMB option
create_owsa_plots <- function(owsa_out, det_results, strategies = v_strats, 
                              threshold = 0.05, top_n = 10, 
                              wtp = 50000, plot_type = "incremental") {
  log_info("Generating OWSA tornado plots...")
  
  comparator <- owsa_out$comparator
  
  for (strategy_to_plot in strategies) {
    
    # Skip if plotting the comparator itself (incremental values would be zero)
    if (strategy_to_plot == comparator && plot_type == "incremental") {
      log_info(paste("Skipping tornado plot for comparator:", strategy_to_plot))
      next
    }
    
    # Base values for the selected strategy
    base_cost <- det_results$table$Cost[det_results$table$Strategy == strategy_to_plot]
    base_qaly <- det_results$table$QALYs[det_results$table$Strategy == strategy_to_plot]
    
    # Base values for comparator
    comp_cost <- det_results$table$Cost[det_results$table$Strategy == comparator]
    comp_qaly <- det_results$table$QALYs[det_results$table$Strategy == comparator]
    
    # Incremental base values
    base_inc_cost <- base_cost - comp_cost
    base_inc_qaly <- base_qaly - comp_qaly
    base_icer <- ifelse(base_inc_qaly != 0, base_inc_cost / base_inc_qaly, NA)
    base_nmb <- base_inc_qaly * wtp - base_inc_cost
    
    ## --- OPTION 1: INCREMENTAL COST TORNADO ---
    if (plot_type %in% c("incremental", "both")) {
      owsa_cost_data <- owsa_out$owsa_Cost %>%
        dplyr::filter(strategy == strategy_to_plot) %>%
        dplyr::mutate(range = abs(Inc_Cost_max - Inc_Cost_min))
      
      max_range_cost <- max(owsa_cost_data$range, na.rm = TRUE)
      
      owsa_cost_data <- owsa_cost_data %>%
        dplyr::filter(range >= threshold * max_range_cost) %>%
        dplyr::arrange(desc(range)) %>%
        dplyr::slice_head(n = top_n) %>%
        dplyr::mutate(pars = factor(pars, levels = pars))
      
      p_owsa_inc_cost <- ggplot(owsa_cost_data, aes(y = pars)) +
        geom_segment(
          aes(x = Inc_Cost_min, xend = Inc_Cost_max, yend = pars),
          color = "#0072B2",
          size = 4,
          alpha = 0.8
        ) +
        geom_vline(
          xintercept = base_inc_cost,
          linetype = "dashed",
          color = "red",
          size = 1
        ) +
        geom_vline(xintercept = 0, linetype = "solid", color = "black", size = 0.5) +
        geom_point(aes(x = base_inc_cost), color = "red", size = 3) +
        scale_y_discrete(limits = rev) +
        labs(
          title = paste0("OWSA: Incremental Cost\n", strategy_to_plot, " vs ", comparator),
          x = "Incremental Cost (CAD)",
          y = "Parameter"
        ) +
        theme_minimal(base_size = 12)
      
      print(p_owsa_inc_cost)
      filename <- file.path(OUT_FIG, paste0("owsa_tornado_inc_cost_", strategy_to_plot, "_vs_", comparator, ".png"))
      try(ggsave(filename, p_owsa_inc_cost, width = 8, height = 7, dpi = 300))
      
      ## --- INCREMENTAL QALYs TORNADO ---
      owsa_qaly_data <- owsa_out$owsa_QALYs %>%
        dplyr::filter(strategy == strategy_to_plot) %>%
        dplyr::mutate(range = abs(Inc_QALYs_max - Inc_QALYs_min))
      
      max_range_qaly <- max(owsa_qaly_data$range, na.rm = TRUE)
      
      owsa_qaly_data <- owsa_qaly_data %>%
        dplyr::filter(range >= threshold * max_range_qaly) %>%
        dplyr::arrange(desc(range)) %>%
        dplyr::slice_head(n = top_n) %>%
        dplyr::mutate(pars = factor(pars, levels = pars))
      
      p_owsa_inc_qaly <- ggplot(owsa_qaly_data, aes(y = pars)) +
        geom_segment(
          aes(x = Inc_QALYs_min, xend = Inc_QALYs_max, yend = pars),
          color = "#009E73",
          size = 4,
          alpha = 0.8
        ) +
        geom_vline(
          xintercept = base_inc_qaly,
          linetype = "dashed",
          color = "red",
          size = 1
        ) +
        geom_vline(xintercept = 0, linetype = "solid", color = "black", size = 0.5) +
        geom_point(aes(x = base_inc_qaly), color = "red", size = 3) +
        scale_y_discrete(limits = rev) +
        labs(
          title = paste0("OWSA: Incremental QALYs\n", strategy_to_plot, " vs ", comparator),
          x = "Incremental QALYs",
          y = "Parameter"
        ) +
        theme_minimal(base_size = 12)
      
      print(p_owsa_inc_qaly)
      filename <- file.path(OUT_FIG, paste0("owsa_tornado_inc_qalys_", strategy_to_plot, "_vs_", comparator, ".png"))
      try(ggsave(filename, p_owsa_inc_qaly, width = 8, height = 7, dpi = 300))
    }
    
    ## --- OPTION 2: NET MONETARY BENEFIT (NMB) TORNADO ---
    if (plot_type %in% c("nmb", "both")) {
      owsa_nmb_data <- owsa_out$owsa_Cost %>%
        dplyr::filter(strategy == strategy_to_plot) %>%
        dplyr::left_join(
          owsa_out$owsa_QALYs %>% 
            dplyr::filter(strategy == strategy_to_plot) %>%
            dplyr::select(pars, Inc_QALYs_min, Inc_QALYs_max),
          by = "pars"
        ) %>%
        dplyr::mutate(
          NMB_min = Inc_QALYs_min * wtp - Inc_Cost_min,
          NMB_max = Inc_QALYs_max * wtp - Inc_Cost_max,
          range = abs(NMB_max - NMB_min)
        )
      
      max_range_nmb <- max(owsa_nmb_data$range, na.rm = TRUE)
      
      owsa_nmb_data <- owsa_nmb_data %>%
        dplyr::filter(range >= threshold * max_range_nmb) %>%
        dplyr::arrange(desc(range)) %>%
        dplyr::slice_head(n = top_n) %>%
        dplyr::mutate(pars = factor(pars, levels = pars))
      
      p_owsa_nmb <- ggplot(owsa_nmb_data, aes(y = pars)) +
        geom_segment(
          aes(x = NMB_min, xend = NMB_max, yend = pars),
          color = "#D55E00",
          size = 4,
          alpha = 0.8
        ) +
        geom_vline(
          xintercept = base_nmb,
          linetype = "dashed",
          color = "red",
          size = 1
        ) +
        geom_vline(xintercept = 0, linetype = "solid", color = "black", size = 0.5) +
        geom_point(aes(x = base_nmb), color = "red", size = 3) +
        scale_y_discrete(limits = rev) +
        labs(
          title = paste0("OWSA: Net Monetary Benefit (WTP = $", 
                         format(wtp, big.mark = ","), ")\n",
                         strategy_to_plot, " vs ", comparator),
          x = "Net Monetary Benefit (CAD)",
          y = "Parameter",
          caption = "NMB > 0 favors intervention"
        ) +
        theme_minimal(base_size = 12)
      
      print(p_owsa_nmb)
      filename <- file.path(OUT_FIG, paste0("owsa_tornado_nmb_", strategy_to_plot, "_vs_", comparator, ".png"))
      try(ggsave(filename, p_owsa_nmb, width = 8, height = 7, dpi = 300))
    }
    
    ## --- OPTION 3: ICER TORNADO (Advanced) ---
    if (plot_type == "icer") {
      owsa_icer_data <- owsa_out$owsa_Cost %>%
        dplyr::filter(strategy == strategy_to_plot) %>%
        dplyr::left_join(
          owsa_out$owsa_QALYs %>% 
            dplyr::filter(strategy == strategy_to_plot) %>%
            dplyr::select(pars, Inc_QALYs_min, Inc_QALYs_max),
          by = "pars"
        ) %>%
        dplyr::mutate(
          ICER_min = ifelse(Inc_QALYs_min != 0, Inc_Cost_min / Inc_QALYs_min, NA),
          ICER_max = ifelse(Inc_QALYs_max != 0, Inc_Cost_max / Inc_QALYs_max, NA),
          range = abs(ICER_max - ICER_min)
        ) %>%
        dplyr::filter(!is.na(ICER_min) & !is.na(ICER_max) & is.finite(ICER_min) & is.finite(ICER_max))
      
      if (nrow(owsa_icer_data) > 0) {
        max_range_icer <- max(owsa_icer_data$range, na.rm = TRUE)
        
        owsa_icer_data <- owsa_icer_data %>%
          dplyr::filter(range >= threshold * max_range_icer) %>%
          dplyr::arrange(desc(range)) %>%
          dplyr::slice_head(n = top_n) %>%
          dplyr::mutate(pars = factor(pars, levels = pars))
        
        p_owsa_icer <- ggplot(owsa_icer_data, aes(y = pars)) +
          geom_segment(
            aes(x = ICER_min, xend = ICER_max, yend = pars),
            color = "#CC79A7",
            size = 4,
            alpha = 0.8
          ) +
          geom_vline(
            xintercept = base_icer,
            linetype = "dashed",
            color = "red",
            size = 1
          ) +
          geom_vline(xintercept = wtp, linetype = "dotted", color = "blue", size = 1) +
          geom_point(aes(x = base_icer), color = "red", size = 3) +
          scale_y_discrete(limits = rev) +
          scale_x_continuous(labels = scales::comma) +
          labs(
            title = paste0("OWSA: ICER\n", strategy_to_plot, " vs ", comparator),
            x = "ICER (CAD per QALY)",
            y = "Parameter",
            caption = paste0("Blue line = WTP threshold ($", format(wtp, big.mark = ","), "/QALY)")
          ) +
          theme_minimal(base_size = 12)
        
        print(p_owsa_icer)
        filename <- file.path(OUT_FIG, paste0("owsa_tornado_icer_", strategy_to_plot, "_vs_", comparator, ".png"))
        try(ggsave(filename, p_owsa_icer, width = 8, height = 7, dpi = 300))
      }
    }
    
    log_info(paste("OWSA plots created for strategy:", strategy_to_plot))
  }
}
