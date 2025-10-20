## ======================
## 6) PSA Functions
## ======================

# Build standard deviation values for PSA
build_sd_values <- function(params) {
  sd <- list()
  
  # Costs (Gamma) — exclude vaccine prices
  rel_sd_cost <- 0.20
  cost_params <- c(
    "c_Healthy", "c_Scarring", "c_Single_Amput", "c_Multiple_Amput",
    "c_Neuro_Disab", "c_Hearing_Loss", "c_Renal_Failure",
    "c_Seizure", "c_Paralysis", "c_Dead"
  )
  
  for (nm in cost_params) {
    sd[[nm]] <- abs(params[[nm]] * rel_sd_cost)
  }
  
  # Infection costs (vector)
  mu_inf <- as.numeric(params$c_IMD_infection)
  sd[["c_IMD_infection"]] <- ifelse(mu_inf == 0, NA_real_, pmax(mu_inf * rel_sd_cost, 0))
  
  # Probabilities (Beta) — per-element capped SD
  sd_prob <- 0.10
  prob_params <- c(
    "p_B", "p_C", "p_W", "p_Y",
    "p_B_DeadIMD", "p_C_DeadIMD", "p_W_DeadIMD", "p_Y_DeadIMD",
    "p_bg_mort"
  )
  
  for (nm in prob_params) {
    mu_vec <- as.numeric(params[[nm]])
    sd_vec <- sapply(mu_vec, function(mu) {
      mu2 <- min(1 - 1e-6, max(1e-6, mu))
      min(sd_prob, sqrt(mu2 * (1 - mu2)) / 2)
    })
    sd[[nm]] <- sd_vec
  }
  
  sd[["p_seq_overall"]] <- sd_prob / 4
  
  # Sequela shares (independent Beta)
  for (nm in names(params$w_seq)) {
    sd[[paste0("share_", nm)]] <- 0.05
  }
  
  # Coverage (Beta)
  for (nm in c("coverage_ABCWY", "coverage_ACWY", "coverage_C", "coverage_B")) {
    sd[[nm]] <- 0.05
  }
  
  # Utilities (Beta) — fix u_Healthy and u_Dead
  sd_util <- 0.02
  util_params <- c(
    "u_Healthy", "u_IMD", "u_Scarring", "u_Single_Amput",
    "u_Multiple_Amput", "u_Neuro_Disability", "u_Hearing_Loss",
    "u_Renal_Failure", "u_Seizure", "u_Paralysis", "u_Dead"
  )
  
  for (nm in util_params) {
    sd[[nm]] <- sd_util
  }
  
  # Fix perfect health and death utilities
  sd[["u_Healthy"]] <- NA
  sd[["u_Dead"]] <- NA
  
  # Multipliers/scalers: keep fixed in PSA (handled in OWSA)
  fixed_params <- c(
    "mult_p_B", "mult_p_C", "mult_p_W", "mult_p_Y",
    "mult_cfr_B", "mult_cfr_C", "mult_cfr_W", "mult_cfr_Y",
    "mult_bg_mort", "mult_c_IMD",
    "ve_scale_B", "ve_scale_C", "ve_scale_W", "ve_scale_Y"
  )
  
  for (nm in fixed_params) {
    sd[[nm]] <- NA
  }
  
  # Vaccine prices excluded from PSA
  for (nm in c("c_MenABCWY", "c_MenACWY", "c_MenC", "c_MenB")) {
    sd[[nm]] <- NA
  }
  
  return(sd)
}

# Generate PSA samples
generate_psa_samples <- function(params, n_sim = 1000, seed = 2025) {
  set.seed(seed)
  log_info(paste("Generating", n_sim, "PSA samples..."))
  
  sdv <- build_sd_values(params)
  
  flat_means <- list()
  flat_sds <- list()
  dists <- list()
  
  # Helper to add scalar parameter
  .add_scalar <- function(nm, mean_val, sd_val, dist) {
    if (is.null(mean_val) || is.null(sd_val) || is.na(sd_val)) return()
    
    if (dist == "beta") {
      adj <- .beta_safe(mean_val, sd_val)
      if (is.null(adj)) return()
      flat_means[[nm]] <<- adj$mean
      flat_sds[[nm]] <<- adj$sd
      dists[[nm]] <<- dist
      
    } else if (dist == "gamma") {
      adj <- .gamma_safe(mean_val, sd_val)
      if (is.null(adj)) return()
      flat_means[[nm]] <<- adj$mean
      flat_sds[[nm]] <<- adj$sd
      dists[[nm]] <<- dist
    }
  }
  
  # Helper to add vector parameter
  .add_vector <- function(nm, mean_vec, sd_vec, dist) {
    if (is.null(mean_vec) || is.null(sd_vec)) return()
    
    L <- length(mean_vec)
    if (length(sd_vec) != L) {
      sd_vec <- rep(sd_vec, length.out = L)
    }
    
    for (i in seq_len(L)) {
      .add_scalar(paste0(nm, i), mean_vec[i], sd_vec[i], dist)
    }
  }
  
  # Add costs (Gamma)
  cost_params <- c(
    "c_Healthy", "c_Scarring", "c_Single_Amput", "c_Multiple_Amput",
    "c_Neuro_Disab", "c_Hearing_Loss", "c_Renal_Failure",
    "c_Seizure", "c_Paralysis", "c_Dead"
  )
  
  for (nm in cost_params) {
    .add_scalar(nm, params[[nm]], sdv[[nm]], "gamma")
  }
  
  .add_vector("c_IMD_infection", params$c_IMD_infection, sdv$c_IMD_infection, "gamma")
  
  # Add probabilities (Beta)
  prob_params <- c(
    "p_B", "p_C", "p_W", "p_Y",
    "p_B_DeadIMD", "p_C_DeadIMD", "p_W_DeadIMD", "p_Y_DeadIMD",
    "p_bg_mort"
  )
  
  for (nm in prob_params) {
    .add_vector(nm, params[[nm]], sdv[[nm]], "beta")
  }
  
  .add_scalar("p_seq_overall", params$p_seq_overall, sdv$p_seq_overall, "beta")
  
  # Add shares (Beta)
  for (nm in names(params$w_seq)) {
    .add_scalar(paste0("share_", nm), params$w_seq[[nm]], sdv[[paste0("share_", nm)]], "beta")
  }
  
  # Add coverage (Beta)
  for (nm in c("coverage_ABCWY", "coverage_ACWY", "coverage_C", "coverage_B")) {
    .add_scalar(nm, params[[nm]], sdv[[nm]], "beta")
  }
  
  # Add utilities (Beta)
  util_params <- c(
    "u_Healthy", "u_IMD", "u_Scarring", "u_Single_Amput",
    "u_Multiple_Amput", "u_Neuro_Disability", "u_Hearing_Loss",
    "u_Renal_Failure", "u_Seizure", "u_Paralysis", "u_Dead"
  )
  
  for (nm in util_params) {
    .add_scalar(nm, params[[nm]], sdv[[nm]], "beta")
  }
  
  # Build gen_psa_samp inputs
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
  
  # Helper to gather vectors from column names
  gather_vec <- function(prefix, df_row) {
    cols <- grep(paste0("^", prefix, "\\d+$"), names(df_row), value = TRUE)
    if (length(cols) == 0) return(NULL)
    as.numeric(df_row[, cols, drop = TRUE])
  }
  
  # Convert to list of parameter sets (more efficient than loop)
  samples <- lapply(seq_len(n_sim), function(i) {
    row <- psa_df[i, , drop = FALSE]
    s <- params  # Start with base parameters
    
    # Update scalar parameters
    scalar_params <- c(
      cost_params, "p_seq_overall",
      paste0("coverage_", c("ABCWY", "ACWY", "C", "B")),
      util_params,
      paste0("share_", names(params$w_seq))
    )
    
    for (nm in scalar_params) {
      if (nm %in% names(row)) {
        s[[nm]] <- as.numeric(row[[nm]])
      }
    }
    
    # Update vector parameters
    for (nm in prob_params) {
      vec <- gather_vec(nm, row)
      if (!is.null(vec)) {
        s[[nm]] <- .clamp(vec, 0, 0.999)
      }
    }
    
    vcost <- gather_vec("c_IMD_infection", row)
    if (!is.null(vcost)) {
      s$c_IMD_infection <- pmax(vcost, 0)
    }
    
    # Normalize shares
    shares_now <- sapply(names(params$w_seq), function(nm) {
      s[[paste0("share_", nm)]]
    })
    s$w_seq <- .normalize_shares(shares_now, default_share = params$w_seq)
    
    # Clamp coverage and rescale VE accordingly
    s$coverage_ABCWY <- .clamp(s$coverage_ABCWY, 0, 1)
    s$coverage_ACWY  <- .clamp(s$coverage_ACWY,  0, 1)
    s$coverage_C     <- .clamp(s$coverage_C,     0, 1)
    s$coverage_B     <- .clamp(s$coverage_B,     0, 1)
    
    # Rescale VE vectors by new coverage
    s$ve_MenABCWY_forACWY <- (params$ve_MenABCWY_forACWY / params$coverage_ABCWY) * s$coverage_ABCWY
    s$ve_MenABCWY_forB    <- (params$ve_MenABCWY_forB    / params$coverage_ABCWY) * s$coverage_ABCWY
    s$ve_MenACWY          <- (params$ve_MenACWY          / params$coverage_ACWY)  * s$coverage_ACWY
    s$ve_MenC             <- (params$ve_MenC             / params$coverage_C)     * s$coverage_C
    s$ve_MenB             <- (params$ve_MenB             / params$coverage_B)     * s$coverage_B
    
    return(s)
  })
  
  log_info("PSA samples generated successfully")
  return(samples)
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
  samp <- generate_psa_samples(params, n_sim = n_sim)
  
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
  return(res)
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

