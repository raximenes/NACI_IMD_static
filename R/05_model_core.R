## ======================
## 5) Model Core Functions
## ======================

# Strategy-specific VE per serogroup (with validation)
get_serogroup_ve <- function(params, strategy, t) {
  # Get VE scalers
  sc <- c(
    B = params$ve_scale_B,
    C = params$ve_scale_C,
    W = params$ve_scale_W,
    Y = params$ve_scale_Y
  )
  
  # Validate scalers
  if (any(sc < 0) || any(sc > 2)) {
    log_error("VE scalers must be between 0 and 2")
  }
  
  # Initialize VE vector
  ve <- c(B = 0, C = 0, W = 0, Y = 0)
  
  # Apply strategy-specific VE
  if (strategy == "MenC") {
    ve["C"] <- params$ve_MenC[t] * sc["C"]
    
  } else if (strategy == "MenACWY") {
    ve[c("C", "W", "Y")] <- params$ve_MenACWY[t] * sc[c("C", "W", "Y")]
    
  } else if (strategy == "MenACWY_MenB") {
    ve[c("C", "W", "Y")] <- params$ve_MenACWY[t] * sc[c("C", "W", "Y")]
    # Booster at age 16 (5 years after initial dose at age 11)
    ve_b_shift <- shift_ve(params$ve_MenB, 5)
    ve["B"] <- ve_b_shift[t] * sc["B"]
    
  } else if (strategy == "MenABCWY") {
    ve["B"] <- params$ve_MenABCWY_forB[t] * sc["B"]
    ve[c("C", "W", "Y")] <- params$ve_MenABCWY_forACWY[t] * sc[c("C", "W", "Y")]
  }
  
  # Clamp VE to valid range
  ve_clamped <- .clamp(ve, 0, 1)
  names(ve_clamped) <- names(ve)  
  
  # Warn if clamping occurred
  if (any(ve != ve_clamped)) {
    log_warn(sprintf("VE clamped for strategy %s at cycle %d", strategy, t))
  }
  
  return(ve_clamped)
}

# Build transition probability array
build_transition_array <- function(params, strategy) {
  S <- length(v_states)
  Tt <- params$n_cycles
  
  a_P <- array(0, 
               dim = c(S, S, Tt),
               dimnames = list(v_states, v_states, paste0("t", 1:Tt)))
  
  for (t in 1:Tt) {
    # Background mortality
    pbg <- .clamp(params$p_bg_mort[t] * params$mult_bg_mort, 0, 0.999)
    
    # Vaccine effectiveness
    ve <- get_serogroup_ve(params, strategy, t)
    
    # Infection risks after VE and multipliers
    p_inf_B <- .clamp(params$p_B[t] * params$mult_p_B * (1 - ve[["B"]]), 0, 0.999)
    p_inf_C <- .clamp(params$p_C[t] * params$mult_p_C * (1 - ve[["C"]]), 0, 0.999)
    p_inf_W <- .clamp(params$p_W[t] * params$mult_p_W * (1 - ve[["W"]]), 0, 0.999)
    p_inf_Y <- .clamp(params$p_Y[t] * params$mult_p_Y * (1 - ve[["Y"]]), 0, 0.999)
    
    p_inf_total <- p_inf_B + p_inf_C + p_inf_W + p_inf_Y
    p_stay_healthy <- 1 - pbg - p_inf_total
    
    # Check for invalid probabilities
    if (is.na(p_stay_healthy)) {
      log_error(sprintf("NA probability for staying healthy at cycle %d", t))
    }
    
    # Handle probability overflow
    if (p_stay_healthy < 0) {
      log_warn(sprintf(
        "Negative stay-healthy probability at cycle %d: %.4f. Scaling infection probabilities.",
        t, p_stay_healthy
      ))
      
      # Scale down infection probabilities proportionally
      scale <- ifelse(p_inf_total > 0, (1 - pbg) / p_inf_total, 0)
      p_inf_B <- p_inf_B * scale
      p_inf_C <- p_inf_C * scale
      p_inf_W <- p_inf_W * scale
      p_inf_Y <- p_inf_Y * scale
      p_stay_healthy <- 1 - pbg - (p_inf_B + p_inf_C + p_inf_W + p_inf_Y)
    }
    
    # Healthy state transitions
    a_P["Healthy", "Healthy", t]              <- p_stay_healthy
    a_P["Healthy", "SeroB_Infect", t]         <- p_inf_B
    a_P["Healthy", "SeroC_Infect", t]         <- p_inf_C
    a_P["Healthy", "SeroW_Infect", t]         <- p_inf_W
    a_P["Healthy", "SeroY_Infect", t]         <- p_inf_Y
    a_P["Healthy", "Background_Mortality", t] <- pbg
    
    # Case fatality rates with multipliers
    cfrB <- .clamp(params$p_B_DeadIMD[t] * params$mult_cfr_B, 0, 0.999)
    cfrC <- .clamp(params$p_C_DeadIMD[t] * params$mult_cfr_C, 0, 0.999)
    cfrW <- .clamp(params$p_W_DeadIMD[t] * params$mult_cfr_W, 0, 0.999)
    cfrY <- .clamp(params$p_Y_DeadIMD[t] * params$mult_cfr_Y, 0, 0.999)
    pseq <- .clamp(params$p_seq_overall, 0, 0.999)
    
    # Fill infection state transitions - SeroB
    a_P["SeroB_Infect", "Dead_IMD", t] <- cfrB
    a_P["SeroB_Infect", "Background_Mortality", t] <- (1 - cfrB) * pbg
    palive_no_bg_B <- (1 - cfrB) * (1 - pbg)
    for (sq in names(params$w_seq)) {
      a_P["SeroB_Infect", sq, t] <- palive_no_bg_B * pseq * params$w_seq[[sq]]
    }
    a_P["SeroB_Infect", "Healthy", t] <- palive_no_bg_B * (1 - pseq)
    
    # Fill infection state transitions - SeroC
    a_P["SeroC_Infect", "Dead_IMD", t] <- cfrC
    a_P["SeroC_Infect", "Background_Mortality", t] <- (1 - cfrC) * pbg
    palive_no_bg_C <- (1 - cfrC) * (1 - pbg)
    for (sq in names(params$w_seq)) {
      a_P["SeroC_Infect", sq, t] <- palive_no_bg_C * pseq * params$w_seq[[sq]]
    }
    a_P["SeroC_Infect", "Healthy", t] <- palive_no_bg_C * (1 - pseq)
    
    # Fill infection state transitions - SeroW
    a_P["SeroW_Infect", "Dead_IMD", t] <- cfrW
    a_P["SeroW_Infect", "Background_Mortality", t] <- (1 - cfrW) * pbg
    palive_no_bg_W <- (1 - cfrW) * (1 - pbg)
    for (sq in names(params$w_seq)) {
      a_P["SeroW_Infect", sq, t] <- palive_no_bg_W * pseq * params$w_seq[[sq]]
    }
    a_P["SeroW_Infect", "Healthy", t] <- palive_no_bg_W * (1 - pseq)
    
    # Fill infection state transitions - SeroY
    a_P["SeroY_Infect", "Dead_IMD", t] <- cfrY
    a_P["SeroY_Infect", "Background_Mortality", t] <- (1 - cfrY) * pbg
    palive_no_bg_Y <- (1 - cfrY) * (1 - pbg)
    for (sq in names(params$w_seq)) {
      a_P["SeroY_Infect", sq, t] <- palive_no_bg_Y * pseq * params$w_seq[[sq]]
    }
    a_P["SeroY_Infect", "Healthy", t] <- palive_no_bg_Y * (1 - pseq)
    
    # Sequelae states: absorbing with background mortality
    for (sq in names(params$w_seq)) {
      a_P[sq, sq, t] <- 1 - pbg
      a_P[sq, "Background_Mortality", t] <- pbg
    }
    
    # Death states: absorbing
    a_P["Dead_IMD", "Dead_IMD", t] <- 1
    a_P["Background_Mortality", "Background_Mortality", t] <- 1
  }
  
  # Validate transition array
  try({
    darthtools::check_transition_probability(a_P, err_stop = FALSE, verbose = FALSE)
  }, silent = TRUE)
  
  try({
    darthtools::check_sum_of_transition_array(
      a_P,
      n_cycles = params$n_cycles,
      err_stop = FALSE,
      verbose = FALSE
    )
  }, silent = TRUE)
  
  return(a_P)
}

# Run Markov cohort simulation with flow tracking
run_markov <- function(a_P, params) {
  Tt <- params$n_cycles
  
  # Initialize trace matrix
  m_M <- matrix(0,
                nrow = Tt + 1,
                ncol = n_states,
                dimnames = list(0:Tt, v_states))
  
  # Initial state vector (everyone starts healthy)
  v_m0 <- rep(0, n_states)
  names(v_m0) <- v_states
  v_m0["Healthy"] <- 1
  m_M[1, ] <- v_m0
  
  # Initialize flow array
  a_flow <- array(0,
                  dim = c(n_states, n_states, Tt),
                  dimnames = list(v_states, v_states, paste0("t", 1:Tt)))
  
  # Markov simulation
  for (t in 1:Tt) {
    # Calculate flows for this cycle
    a_flow[, , t] <- diag(m_M[t, ]) %*% a_P[, , t]
    
    # Update state occupancy
    m_M[t + 1, ] <- m_M[t, ] %*% a_P[, , t]
  }
  
  return(list(trace = m_M, flow = a_flow))
}

# Get state costs (with societal perspective adjustment)
# get_state_costs <- function(params, t) {
#   # Sequelae costs
#   c_seq <- c(
#     Scarring         = params$c_Scarring,
#     Single_Amput     = params$c_Single_Amput,
#     Multiple_Amput   = params$c_Multiple_Amput,
#     Neuro_Disability = params$c_Neuro_Disab,
#     Hearing_Loss     = params$c_Hearing_Loss,
#     Renal_Failure    = params$c_Renal_Failure,
#     Seizure          = params$c_Seizure,
#     Paralysis        = params$c_Paralysis
#   )
#   
#   # Add productivity costs if societal perspective
#   if (identical(params$perspective, "societal")) {
#     # Build societal costs from individual parameters
#     c_prod <- c(
#       Scarring         = params$c_prod_Scarring,
#       Single_Amput     = params$c_prod_Single_Amput,
#       Multiple_Amput   = params$c_prod_Multiple_Amput,
#       Neuro_Disability = params$c_prod_Neuro_Disability,
#       Hearing_Loss     = params$c_prod_Hearing_Loss,
#       Renal_Failure    = params$c_prod_Renal_Failure,
#       Seizure          = params$c_prod_Seizure,
#       Paralysis        = params$c_prod_Paralysis
#     )
#     c_seq <- c_seq + c_prod[names(c_seq)]
#   }
#   
#   # Infection costs (with multiplier)
#   c_inf <- c(
#     SeroB_Infect = params$c_IMD_infection[t],
#     SeroC_Infect = params$c_IMD_infection[t],
#     SeroW_Infect = params$c_IMD_infection[t],
#     SeroY_Infect = params$c_IMD_infection[t]
#   ) * params$mult_c_IMD
#   
#   # Other state costs
#   c_other <- c(
#     Healthy = params$c_Healthy,
#     Dead_IMD = params$c_Dead,
#     Background_Mortality = params$c_Dead
#   )
#   
#   # Combine in correct state order
#   costs <- setNames(
#     c(c_other["Healthy"], c_inf, c_seq, c_other[c("Dead_IMD", "Background_Mortality")]),
#     v_states
#   )
#   
#   return(costs)
# }
get_state_costs <- function(params, t) {
  
  # Base healthcare costs (must exist)
  costs <- c(
    Healthy = params$c_Healthy,
    Scarring = params$c_Scarring,
    Single_Amput = params$c_Single_Amput,
    Multiple_Amput = params$c_Multiple_Amput,
    Neuro_Disability = params$c_Neuro_Disab,
    Hearing_Loss = params$c_Hearing_Loss,
    Renal_Failure = params$c_Renal_Failure,
    Seizure = params$c_Seizure,
    Paralysis = params$c_Paralysis,
    Death = 0
  )
  
  # Verify all costs are numeric
  if (any(is.na(costs)) || any(!is.numeric(costs))) {
    na_states <- names(costs)[is.na(costs)]
    stop(sprintf("Base healthcare costs have NA values for: %s", 
                 paste(na_states, collapse = ", ")))
  }
  
  # Add societal productivity costs if perspective is societal
  if (!is.null(params$perspective) && params$perspective == "societal") {
    
    # Define sequelae states (exclude Healthy and Death)
    sequelae_states <- c("Scarring", "Single_Amput", "Multiple_Amput",
                         "Neuro_Disability", "Hearing_Loss", "Renal_Failure",
                         "Seizure", "Paralysis")
    
    # Add productivity costs to each sequelae state
    for (state in sequelae_states) {
      prod_cost_param <- paste0("c_prod_", state)
      prod_cost <- params[[prod_cost_param]]
      
      # Safety check: only add if productivity cost exists and is valid
      if (!is.null(prod_cost) && !is.na(prod_cost) && is.numeric(prod_cost)) {
        costs[state] <- costs[state] + prod_cost
      } else {
        warning(sprintf("Productivity cost %s is NULL/NA or invalid. Using 0.", prod_cost_param))
        # costs[state] remains unchanged (no productivity cost added)
      }
    }
  }
  # Final verification: ensure no NA values in the result
  if (any(is.na(costs))) {
    na_states <- names(costs)[is.na(costs)]
    stop(sprintf("FINAL costs vector has NA values for states: %s. Check Excel data.", 
                 paste(na_states, collapse = ", ")))
  }
  return(costs)
}
# Get state utilities
get_state_utils <- function(params) {
  utils <- c(
    Healthy          = params$u_Healthy,
    SeroB_Infect     = params$u_IMD,
    SeroC_Infect     = params$u_IMD,
    SeroW_Infect     = params$u_IMD,
    SeroY_Infect     = params$u_IMD,
    Scarring         = params$u_Scarring,
    Single_Amput     = params$u_Single_Amput,
    Multiple_Amput   = params$u_Multiple_Amput,
    Neuro_Disability = params$u_Neuro_Disability,
    Hearing_Loss     = params$u_Hearing_Loss,
    Renal_Failure    = params$u_Renal_Failure,
    Seizure          = params$u_Seizure,
    Paralysis        = params$u_Paralysis,
    Dead_IMD         = params$u_Dead,
    Background_Mortality = params$u_Dead
  )
  
  return(utils)
}

# Calculate vaccination costs by cycle
vaccination_costs <- function(params, strategy) {
  v <- rep(0, params$n_cycles + 1)
  
  if (strategy == "MenC") {
    v[1] <- params$coverage_C * (params$c_MenC + params$c_admin)
    
  } else if (strategy == "MenACWY") {
    v[1] <- params$coverage_ACWY * (params$c_MenACWY + params$c_admin)
    
  } else if (strategy == "MenACWY_MenB") {
    # Initial dose at age 11
    v[1] <- params$coverage_ACWY * (params$c_MenACWY + params$c_admin)
    # Booster at age 16 (cycle 6)
    if (params$n_cycles >= 6) {
      v[6 + 1] <- params$coverage_B * (params$c_MenB + params$c_admin)
    }
    
  } else if (strategy == "MenABCWY") {
    v[1] <- params$coverage_ABCWY * (params$c_MenABCWY + params$c_admin)
  }
  
  return(v)
}

# Evaluate a single strategy
eval_strategy <- function(params, strategy) {
  # Build transition array
  a_P <- build_transition_array(params, strategy)
  
  # Run Markov simulation
  sim <- run_markov(a_P, params)
  m_M <- sim$trace
  a_flow <- sim$flow
  
  # Get utilities
  v_util <- get_state_utils(params)
  
  # Get cycle correction weights
  v_wcc <- darthtools::gen_wcc(params$n_cycles, method = "Simpson1/3")
  
  # Initialize cost and QALY vectors
  v_cost_cycle <- rep(0, params$n_cycles + 1)
  v_qaly_cycle <- rep(0, params$n_cycles + 1)
  v_vacc_cycle <- vaccination_costs(params, strategy)
  
  # Calculate costs and QALYs per cycle
  for (t in 1:params$n_cycles) {
    c_state <- get_state_costs(params, t)
    v_cost_cycle[t] <- sum(m_M[t, ] * c_state) + v_vacc_cycle[t]
    v_qaly_cycle[t] <- sum(m_M[t, ] * v_util)
  }
  
  # Final cycle
  c_state_final <- get_state_costs(params, params$n_cycles)
  v_cost_cycle[params$n_cycles + 1] <- sum(m_M[params$n_cycles + 1, ] * c_state_final) + 
    v_vacc_cycle[params$n_cycles + 1]
  v_qaly_cycle[params$n_cycles + 1] <- sum(m_M[params$n_cycles + 1, ] * v_util)
  
  # Discount and sum
  tot_cost <- sum(v_cost_cycle * v_discount_cost * v_wcc) * n_cohort
  tot_qaly <- sum(v_qaly_cycle * v_discount_qaly * v_wcc) * n_cohort
  
  # Epidemiological outcomes: infections by serogroup
  n_inf <- c(
    B = sum(a_flow["Healthy", "SeroB_Infect", ]),
    C = sum(a_flow["Healthy", "SeroC_Infect", ]),
    W = sum(a_flow["Healthy", "SeroW_Infect", ]),
    Y = sum(a_flow["Healthy", "SeroY_Infect", ])
  ) * n_cohort
  
  # Sequelae by type
  seq_types <- names(params$w_seq)
  n_seq <- setNames(rep(0, length(seq_types)), seq_types)
  
  for (sq in seq_types) {
    n_seq[sq] <- sum(
      a_flow["SeroB_Infect", sq, ] +
        a_flow["SeroC_Infect", sq, ] +
        a_flow["SeroW_Infect", sq, ] +
        a_flow["SeroY_Infect", sq, ]
    ) * n_cohort
  }
  
  # IMD deaths
  n_deaths_imd <- sum(
    a_flow["SeroB_Infect", "Dead_IMD", ] +
      a_flow["SeroC_Infect", "Dead_IMD", ] +
      a_flow["SeroW_Infect", "Dead_IMD", ] +
      a_flow["SeroY_Infect", "Dead_IMD", ]
  ) * n_cohort
  
  # Number vaccinated (persons and doses)
  n_persons <- switch(
    strategy,
    "MenC"          = params$coverage_C,
    "MenACWY"       = params$coverage_ACWY,
    "MenACWY_MenB"  = 1 - (1 - params$coverage_ACWY) * (1 - params$coverage_B),
    "MenABCWY"      = params$coverage_ABCWY,
    0
  ) * n_cohort
  
  n_doses <- switch(
    strategy,
    "MenC"          = params$coverage_C,
    "MenACWY"       = params$coverage_ACWY,
    "MenACWY_MenB"  = params$coverage_ACWY + params$coverage_B,
    "MenABCWY"      = params$coverage_ABCWY,
    0
  ) * n_cohort
  
  return(list(
    cost = tot_cost,
    qalys = tot_qaly,
    epi = list(
      infections = n_inf,
      sequelae = n_seq,
      deaths_imd = n_deaths_imd
    ),
    vacc = list(
      persons = n_persons,
      doses = n_doses
    ),
    trace = m_M,
    flow = a_flow
  ))
}

# Evaluate all strategies
eval_all_strategies <- function(params) {
  res_list <- lapply(v_strats, function(s) {
    eval_strategy(params, s)
  })
  names(res_list) <- v_strats
  
  return(res_list)
}

# Summarize deterministic results
summarize_det <- function(res_list, comparator = "MenC") {
  # Extract costs and QALYs
  df <- tibble::tibble(
    Strategy = names(res_list),
    Cost     = vapply(res_list, function(x) x$cost, numeric(1)),
    QALYs    = vapply(res_list, function(x) x$qalys, numeric(1))
  )
  
  # Calculate ICERs
  icers <- dampack::calculate_icers(
    cost = df$Cost,
    effect = df$QALYs,
    strategies = df$Strategy
  )
  
  # Incremental epidemiological outcomes vs comparator
  comp <- res_list[[comparator]]
  
  sum_inf <- function(x) sum(x$epi$infections)
  sum_seq <- function(x) sum(x$epi$sequelae)
  
  df_inc <- tibble::tibble(
    Strategy = names(res_list),
    Infections = vapply(res_list, sum_inf, numeric(1)),
    Sequelae   = vapply(res_list, sum_seq, numeric(1)),
    Deaths_IMD = vapply(res_list, function(x) x$epi$deaths_imd, numeric(1)),
    Vacc_Persons = vapply(res_list, function(x) x$vacc$persons, numeric(1)),
    Vacc_Doses   = vapply(res_list, function(x) x$vacc$doses, numeric(1))
  ) %>%
    dplyr::mutate(
      Infections_prevented = Infections[Strategy == comparator][1] - Infections,
      Sequelae_prevented   = Sequelae[Strategy == comparator][1] - Sequelae,
      Deaths_prevented     = Deaths_IMD[Strategy == comparator][1] - Deaths_IMD,
      NNV_infection_person = ifelse(Infections_prevented > 0, 
                                    Vacc_Persons / Infections_prevented, 
                                    NA_real_),
      NNV_infection_dose   = ifelse(Infections_prevented > 0, 
                                    Vacc_Doses / Infections_prevented, 
                                    NA_real_),
      NNV_death_person     = ifelse(Deaths_prevented > 0, 
                                    Vacc_Persons / Deaths_prevented, 
                                    NA_real_),
      NNV_death_dose       = ifelse(Deaths_prevented > 0, 
                                    Vacc_Doses / Deaths_prevented, 
                                    NA_real_)
    )
  
  return(list(
    table = df,
    icers = icers,
    incremental = df_inc
  ))
}