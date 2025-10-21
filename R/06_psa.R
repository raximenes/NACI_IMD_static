## ======================
## 6) PSA Functions (UPDATED - FIXED)
## ======================
# - Using SD values directly from Excel (no fixed values)
# - Supporting time-varying SD vectors
# - Vaccine costs now included with their SDs from Excel
# - Fixed: Replaced all log_warn with message() or log_info()

# Generate PSA samples using Excel SD values
generate_psa_samples <- function(params, n_sim = 1000, seed = 2025) {
  set.seed(seed)
  log_info(paste("Generating", n_sim, "PSA samples using Excel SD values..."))
  
  flat_means <- list()
  flat_sds <- list()
  dists <- list()
  
  # Helper to add scalar parameter
  # Only adds to PSA if SD is valid (not NULL, NA, or <= 0)
  # Otherwise parameter stays FIXED at its deterministic value
  .add_scalar <- function(nm, mean_val, sd_val, dist) {
    # Check if parameter should be fixed (no SD)
    if (is.null(sd_val) || is.na(sd_val) || sd_val <= 0) {
      return()
    }
    
    if (is.null(mean_val) || is.na(mean_val)) {
      return()
    }
    if (dist == "beta") {
      if (mean_val < 1e-5) {
        # For very small probabilities, use gamma instead
        adj <- .gamma_safe(mean_val, sd_val)
        if (is.null(adj)) return()
        
        flat_means[[nm]] <<- adj$mean
        flat_sds[[nm]] <<- adj$sd
        dists[[nm]] <<- "gamma"
      } else {
        adj <- .beta_safe(mean_val, sd_val)
        if (is.null(adj)) return()
        
        flat_means[[nm]] <<- adj$mean
        flat_sds[[nm]] <<- adj$sd
        dists[[nm]] <<- dist
      }
    } else if (dist == "gamma") {
      adj <- .gamma_safe(mean_val, sd_val)
      if (is.null(adj)) return()
      
      flat_means[[nm]] <<- adj$mean
      flat_sds[[nm]] <<- adj$sd
      dists[[nm]] <<- dist
    }
  }
  
  # Helper to add vector parameter
  # Only adds elements that have valid SDs
  # Elements without SD remain FIXED at deterministic values
  .add_vector <- function(nm, mean_vec, sd_vec, dist) {
    if (is.null(mean_vec)) return()
    if (is.null(sd_vec)) {
      return()
    }
    
    L <- length(mean_vec)
    if (length(sd_vec) != L) {
      sd_vec <- rep(sd_vec, length.out = L)
    }
    
    for (i in seq_len(L)) {
      if (!is.na(sd_vec[i]) && sd_vec[i] > 0) {
        .add_scalar(paste0(nm, i), mean_vec[i], sd_vec[i], dist)
      }
    }
  }
  
  ## ---- COSTS (Gamma distribution) ----
  
  # Sequelae costs (scalar)
  cost_seq_params <- c(
    "c_Healthy", "c_Scarring", "c_Single_Amput", "c_Multiple_Amput",
    "c_Neuro_Disab", "c_Hearing_Loss", "c_Renal_Failure",
    "c_Seizure", "c_Paralysis", "c_Dead"
  )
  
  for (nm in cost_seq_params) {
    sd_nm <- paste0("sd_", nm)
    .add_scalar(nm, params[[nm]], params[[sd_nm]], "gamma")
  }
  
  # Infection costs (time-varying vector)
  .add_vector("c_IMD_infection", params$c_IMD_infection, params$sd_c_IMD_infection, "gamma")
  
  # Vaccine costs (now included in PSA with Excel SDs)
  if (exists("file_price") && file.exists(file_price)) {
    pr <- readxl::read_excel(file_price)
    
    get_price_sd <- function(nm) {
      if (!"sd" %in% names(pr)) return(NA)
      v <- pr$sd[pr$Name == nm]
      if (length(v) == 0) return(NA)
      as.numeric(v[1])
    }
    
    .add_scalar("c_MenABCWY", params$c_MenABCWY, get_price_sd("c_MenABCWY"), "gamma")
    .add_scalar("c_MenACWY",  params$c_MenACWY,  get_price_sd("c_MenACWY"),  "gamma")
    .add_scalar("c_MenC",     params$c_MenC,     get_price_sd("c_MenC"),     "gamma")
    .add_scalar("c_MenB",     params$c_MenB,     get_price_sd("c_MenB"),     "gamma")
  }
  
  ## ---- PROBABILITIES (Beta distribution) ----
  
  # Infection probabilities (time-varying vectors)
  prob_inf_params <- c("p_B", "p_C", "p_W", "p_Y")
  for (nm in prob_inf_params) {
    sd_nm <- paste0("sd_", nm)
    .add_vector(nm, params[[nm]], params[[sd_nm]], "beta")
  }
  
  # Case fatality rates (time-varying vectors)
  cfr_params <- c("p_B_DeadIMD", "p_C_DeadIMD", "p_W_DeadIMD", "p_Y_DeadIMD")
  for (nm in cfr_params) {
    sd_nm <- paste0("sd_", nm)
    .add_vector(nm, params[[nm]], params[[sd_nm]], "beta")
  }
  
  # Background mortality (time-varying vector)
  .add_vector("p_bg_mort", params$p_bg_mort, params$sd_p_bg_mort, "beta")
  
  # Overall sequelae probability (scalar)
  .add_scalar("p_seq_overall", params$p_seq_overall, params$sd_p_seq_overall, "beta")
  
  # Sequelae shares (scalar, independent Beta for each)
  for (nm in names(params$w_seq)) {
    .add_scalar(paste0("share_", nm), params$w_seq[[nm]], params$sd_w_seq[[nm]], "beta")
  }
  
  # Coverage (scalar)
  for (nm in c("coverage_ABCWY", "coverage_ACWY", "coverage_C", "coverage_B")) {
    sd_nm <- paste0("sd_", nm)
    if (!is.null(params[[sd_nm]])) {
      .add_scalar(nm, params[[nm]], params[[sd_nm]], "beta")
    }
  }
  
  ## ---- UTILITIES (Beta distribution) ----
  
  util_params <- c(
    "u_Healthy", "u_IMD", "u_Scarring", "u_Single_Amput",
    "u_Multiple_Amput", "u_Neuro_Disability", "u_Hearing_Loss",
    "u_Renal_Failure", "u_Seizure", "u_Paralysis", "u_Dead"
  )
  
  for (nm in util_params) {
    sd_nm <- paste0("sd_", nm)
    .add_scalar(nm, params[[nm]], params[[sd_nm]], "beta")
  }
  
  ## ---- VACCINE EFFECTIVENESS (Beta distribution) ----
  # VE base values with SD from Excel
  if (!is.null(params$ve_data)) {
    for (i in seq_len(nrow(params$ve_data))) {
      vac_name <- params$ve_data$vaccine[i]
      ve_mean <- params$ve_data$effectiveness[i]
      ve_sd <- params$ve_data$sd_effectiveness[i]
      
      .add_scalar(paste0("ve_base_", vac_name), ve_mean, ve_sd, "beta")
    }
  }
  
  ## ---- BUILD DAMPACK INPUT ----
  
  param_names <- names(flat_means)
  dist_types  <- unlist(dists)
  means_vec   <- unlist(flat_means)
  sds_vec     <- unlist(flat_sds)
  
  dists_params <- lapply(seq_along(param_names), function(i) {
    c(mean = means_vec[i], sd = sds_vec[i])
  })
  names(dists_params) <- param_names
  
  param_types <- rep("mean, sd", length(param_names))
  names(param_types) <- param_names
  
  # Generate PSA samples using dampack
  psa_df <- dampack::gen_psa_samp(
    params = param_names,
    dists  = dist_types,
    parameterization_types = param_types,
    dists_params = dists_params,
    nsamp  = n_sim
  )
  
  # Remove nsamp column if present
  if ("nsamp" %in% names(psa_df)) {
    psa_df$nsamp <- NULL
  }
  
  ## ---- CONVERT TO PARAMETER SETS ----
  
  # Helper to gather vectors from column names
  gather_vec <- function(prefix, df_row) {
    cols <- grep(paste0("^", prefix, "\\d+$"), names(df_row), value = TRUE)
    if (length(cols) == 0) return(NULL)
    as.numeric(df_row[, cols, drop = TRUE])
  }
  
  # Convert to list of parameter sets
  samples <- lapply(seq_len(n_sim), function(i) {
    row <- psa_df[i, , drop = FALSE]
    s <- params  # Start with base parameters (FIXED values preserved)
    
    ## Update scalar cost parameters (only if they were in PSA)
    for (nm in cost_seq_params) {
      if (nm %in% names(row)) {
        s[[nm]] <- as.numeric(row[[nm]])
      }
    }
    
    # Update vaccine costs (only if they were in PSA)
    for (nm in c("c_MenABCWY", "c_MenACWY", "c_MenC", "c_MenB")) {
      if (nm %in% names(row)) {
        s[[nm]] <- as.numeric(row[[nm]])
      }
    }
    
    ## Update vector cost parameters (only elements that were in PSA)
    vcost_psa <- gather_vec("c_IMD_infection", row)
    if (!is.null(vcost_psa)) {
      vcost_final <- params$c_IMD_infection
      psa_indices <- as.numeric(gsub("c_IMD_infection", "", names(row)[grep("^c_IMD_infection\\d+$", names(row))]))
      for (idx in psa_indices) {
        if (idx <= length(vcost_final)) {
          vcost_final[idx] <- pmax(vcost_psa[which(psa_indices == idx)], 0)
        }
      }
      s$c_IMD_infection <- vcost_final
    }
    
    ## Update infection probabilities (only elements that were in PSA)
    for (nm in prob_inf_params) {
      vec_psa <- gather_vec(nm, row)
      if (!is.null(vec_psa)) {
        vec_final <- params[[nm]]
        psa_indices <- as.numeric(gsub(nm, "", names(row)[grep(paste0("^", nm, "\\d+$"), names(row))]))
        for (idx in psa_indices) {
          if (idx <= length(vec_final)) {
            vec_final[idx] <- .clamp(vec_psa[which(psa_indices == idx)], 0, 0.999)
          }
        }
        s[[nm]] <- vec_final
      }
    }
    
    ## Update CFR vectors (only elements that were in PSA)
    for (nm in cfr_params) {
      vec_psa <- gather_vec(nm, row)
      if (!is.null(vec_psa)) {
        vec_final <- params[[nm]]
        psa_indices <- as.numeric(gsub(nm, "", names(row)[grep(paste0("^", nm, "\\d+$"), names(row))]))
        for (idx in psa_indices) {
          if (idx <= length(vec_final)) {
            vec_final[idx] <- .clamp(vec_psa[which(psa_indices == idx)], 0, 0.999)
          }
        }
        s[[nm]] <- vec_final
      }
    }
    
    ## Update background mortality (only elements that were in PSA)
    vec_bg_psa <- gather_vec("p_bg_mort", row)
    if (!is.null(vec_bg_psa)) {
      vec_bg_final <- params$p_bg_mort
      psa_indices <- as.numeric(gsub("p_bg_mort", "", names(row)[grep("^p_bg_mort\\d+$", names(row))]))
      for (idx in psa_indices) {
        if (idx <= length(vec_bg_final)) {
          vec_bg_final[idx] <- .clamp(vec_bg_psa[which(psa_indices == idx)], 0, 0.999)
        }
      }
      s$p_bg_mort <- vec_bg_final
    }
    
    ## Update sequelae probability (only if in PSA)
    if ("p_seq_overall" %in% names(row)) {
      s$p_seq_overall <- .clamp(as.numeric(row[["p_seq_overall"]]), 0, 0.999)
    }
    
    ## Update and normalize sequelae shares (only those that were in PSA)
    shares_now <- sapply(names(params$w_seq), function(nm) {
      share_nm <- paste0("share_", nm)
      if (share_nm %in% names(row)) {
        as.numeric(row[[share_nm]])
      } else {
        params$w_seq[[nm]]
      }
    })
    s$w_seq <- .normalize_shares(shares_now, default_share = params$w_seq)
    
    ## Update coverage (only if in PSA)
    for (nm in c("coverage_ABCWY", "coverage_ACWY", "coverage_C", "coverage_B")) {
      if (nm %in% names(row)) {
        s[[nm]] <- .clamp(as.numeric(row[[nm]]), 0, 1)
      }
    }
    
    ## Update utilities (only if in PSA)
    for (nm in util_params) {
      if (nm %in% names(row)) {
        s[[nm]] <- .clamp(as.numeric(row[[nm]]), 0, 1)
      }
    }
    
    ## Update VE base values and recalculate VE vectors (only if in PSA)
    if (!is.null(params$ve_data)) {
      ve_data_updated <- params$ve_data
      ve_changed <- FALSE
      
      for (j in seq_len(nrow(ve_data_updated))) {
        vac_name <- ve_data_updated$vaccine[j]
        ve_nm <- paste0("ve_base_", vac_name)
        
        if (ve_nm %in% names(row)) {
          ve_data_updated$effectiveness[j] <- .clamp(as.numeric(row[[ve_nm]]), 0, 1)
          ve_changed <- TRUE
        }
      }
      
      # Only recalculate if VE values changed
      if (ve_changed) {
        ve_all <- calculate_all_ve(ve_data_updated, n_cycles = params$n_cycles)
        
        s$ve_MenABCWY_forACWY <- ve_all$Effectiveness[ve_all$Vaccine == "MenABCWY_for_SeroACWY"] * s$coverage_ABCWY
        s$ve_MenABCWY_forB    <- ve_all$Effectiveness[ve_all$Vaccine == "MenABCWY_for_SeroB"]    * s$coverage_ABCWY
        s$ve_MenACWY          <- ve_all$Effectiveness[ve_all$Vaccine == "MenACWY"]               * s$coverage_ACWY
        s$ve_MenC             <- ve_all$Effectiveness[ve_all$Vaccine == "MenC"]                  * s$coverage_C
        s$ve_MenB             <- ve_all$Effectiveness[ve_all$Vaccine == "MenB"]                  * s$coverage_B
      }
    }
    
    return(s)
  })
  
  log_info(paste("Total parameters varying in PSA:", length(flat_means)))
  log_info("PSA samples generated successfully with Excel SD values")
  
  # Return both samples and the raw PSA dataframe for visualization
  return(list(
    samples = samples,
    psa_df = psa_df
  ))
}

# Run PSA with parallel processing
run_psa <- function(params, n_sim = 1000, wtp = 50000) {
  # Setup parallel cluster
  n_cores <- parallel::detectCores() - 1
  if (is.na(n_cores) || n_cores < 1) n_cores <- 1
  
  cl <- parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cl)
  
  log_info(paste("Running", n_sim, "PSA simulations on", n_cores, "cores..."))
  
  # Generate all parameter samples
  psa_output <- generate_psa_samples(params, n_sim = n_sim)
  samp <- psa_output$samples
  psa_df <- psa_output$psa_df
  
  # Run simulations in parallel
  res <- foreach::foreach(
    i = 1:n_sim,
    .combine = 'rbind',
    .packages = c('dplyr', 'tibble', 'darthtools'),
    .export = c(
      # Functions
      "eval_all_strategies", "eval_strategy", "build_transition_array",
      "run_markov", "get_serogroup_ve", "shift_ve", ".clamp",
      "get_state_costs", "get_state_utils", "vaccination_costs",
      # Logging functions (CRITICAL!)
      "log_info", "log_warn", "log_error",
      # Global variables
      "v_strats", "v_states", "n_states", "c_prod_societal",
      "v_discount_cost", "v_discount_qaly", "n_cohort"
    )
  ) %dopar% {
    # Evaluate all strategies for this parameter set
    r_i <- eval_all_strategies(samp[[i]])
    
    # Return results as tibble
    tibble::tibble(
      sim = i,
      Strategy = names(r_i),
      Cost  = vapply(r_i, function(x) x$cost, numeric(1)),
      QALYs = vapply(r_i, function(x) x$qalys, numeric(1))
    )
  }
  
  # Stop parallel cluster
  parallel::stopCluster(cl)
  
  # Calculate NMB
  res <- res %>%
    dplyr::mutate(NMB = QALYs * wtp - Cost)
  
  log_info("PSA completed successfully")
  
  # Return both results and parameter distributions
  return(list(
    results = res,
    psa_df = psa_df
  ))
}

# Create PSA objects for dampack
make_psa_objects <- function(psa_df_long, strategies, wtp_vec = seq(0, 200000, by = 10000)) {
  log_info("Creating PSA objects for CEA analysis...")
  
  sims <- sort(unique(psa_df_long$sim))
  
  # Build cost dataframe
  cost_df <- lapply(strategies, function(s) {
    psa_df_long %>%
      dplyr::filter(Strategy == s) %>%
      dplyr::arrange(sim) %>%
      dplyr::pull(Cost)
  }) %>%
    do.call(cbind, .) %>%
    as.data.frame()
  
  colnames(cost_df) <- strategies
  rownames(cost_df) <- sims
  
  # Build effectiveness dataframe
  eff_df <- lapply(strategies, function(s) {
    psa_df_long %>%
      dplyr::filter(Strategy == s) %>%
      dplyr::arrange(sim) %>%
      dplyr::pull(QALYs)
  }) %>%
    do.call(cbind, .) %>%
    as.data.frame()
  
  colnames(eff_df) <- strategies
  rownames(eff_df) <- sims
  
  # Create dampack objects
  psa_obj <- dampack::make_psa_obj(
    cost = cost_df,
    effectiveness = eff_df,
    strategies = strategies
  )
  
  ceac_obj <- dampack::ceac(wtp = wtp_vec, psa = psa_obj)
  evpi_obj <- dampack::calc_evpi(psa = psa_obj, wtp = wtp_vec, pop = n_cohort)
  
  log_info("PSA objects created successfully")
  
  return(list(
    psa = psa_obj,
    ceac = ceac_obj,
    evpi = evpi_obj
  ))
}