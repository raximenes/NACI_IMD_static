###############################################################################
# 06_psa_functions.R - FULLY COMMENTED VERSION
#
# Purpose: Probabilistic Sensitivity Analysis (PSA) implementation
#
# Key Components:
#   1. Local helper functions for Beta/Gamma parameterization
#   2. PSA sample generation using Excel SD values
#   3. Parallel PSA execution
#   4. CEAC/EVPI object creation for dampack
#
# Critical Features:
#   - Uses SD values directly from Excel (no arbitrary assumptions)
#   - Supports time-varying parameters (vectors)
#   - Vaccine costs ALWAYS fixed (excluded from PSA)
#   - VE-coverage rescaling for consistency
#   - Parallel processing for speed
###############################################################################

## ======================
## 6) PSA Functions
## ======================

# ============================================================
# LOCAL HELPER FUNCTIONS (Beta & Gamma Parameterization)
# ============================================================
# These functions convert (mean, SD) to distribution parameters
# Required because R's rbeta() and rgamma() need shape parameters,
# not mean and SD
#
# WHY LOCAL? To avoid darthtools dependency and have full control

# Convert mean and SD to Beta distribution parameters (alpha, beta)
beta_parm_local <- function(mean_val, sd_val) {
  # Beta distribution: mean = alpha / (alpha + beta)
  # Variance = (alpha * beta) / ((alpha + beta)^2 * (alpha + beta + 1))
  # SD = sqrt(variance)
  
  if (mean_val <= 0 || mean_val >= 1 || sd_val <= 0) {
    return(NULL)
  }
  
  # Method of moments (solve for alpha and beta from mean and variance):
  # alpha + beta = mean * (1 - mean) / sd^2 - 1
  # alpha = mean * (alpha + beta)
  # beta = (1 - mean) * (alpha + beta)
  
  sum_ab <- mean_val * (1 - mean_val) / (sd_val^2) - 1
  
  if (sum_ab <= 0) {
    return(NULL)
  }
  
  alpha <- mean_val * sum_ab
  beta <- (1 - mean_val) * sum_ab
  
  if (alpha <= 0 || beta <= 0) {
    return(NULL)
  }
  
  list(alpha = alpha, beta = beta)
}

# Convert mean and SD to Gamma distribution parameters (shape, rate)
gamma_parm_local <- function(mean_val, sd_val) {
  # Gamma distribution: mean = shape / rate = alpha / beta
  # Variance = shape / rate^2 = alpha / beta^2
  # SD = sqrt(variance)
  
  if (mean_val <= 0 || sd_val <= 0) {
    return(NULL)
  }
  
  # Method of moments (solve for shape and rate from mean and variance):
  # shape = mean^2 / variance = mean^2 / sd^2
  # rate = mean / variance = mean / sd^2
  
  shape <- mean_val^2 / (sd_val^2)
  rate <- mean_val / (sd_val^2)
  
  if (shape <= 0 || rate <= 0) {
    return(NULL)
  }
  
  list(alpha = shape, beta = rate)  # R uses alpha/beta naming for rgamma
}

# Safe wrapper for Beta parameterization with validation
beta_safe_local <- function(mean_val, sd_val) {
  result <- beta_parm_local(mean_val, sd_val)
  
  if (is.null(result)) {
    return(NULL)
  }
  
  # Check if parameters are reasonable
  if (is.na(result$alpha) || is.na(result$beta) || 
      result$alpha <= 0 || result$beta <= 0) {
    return(NULL)
  }
  
  # Return in format expected by rbeta
  list(mean = mean_val, sd = sd_val, alpha = result$alpha, beta = result$beta)
}

# Safe wrapper for Gamma parameterization
gamma_safe_local <- function(mean_val, sd_val) {
  result <- gamma_parm_local(mean_val, sd_val)
  
  if (is.null(result)) {
    return(NULL)
  }
  
  # Check if parameters are reasonable
  if (is.na(result$alpha) || is.na(result$beta) || 
      result$alpha <= 0 || result$beta <= 0) {
    return(NULL)
  }
  
  # Return in format expected by rgamma
  list(mean = mean_val, sd = sd_val, alpha = result$alpha, beta = result$beta)
}

# Clamp values to range [min, max]
.clamp <- function(x, min_val = 0, max_val = 1) {
  pmax(min_val, pmin(x, max_val))
}

# ============================================================
# MAIN PSA FUNCTION: Generate Parameter Samples
# ============================================================
#
# PURPOSE: Create n_sim sets of parameter values by sampling from
#          distributions defined by (mean, SD) pairs from Excel
#
# THREE-LIST ARCHITECTURE:
# ───────────────────────────────────────────────────────────
# This function builds 3 parallel lists that work together:
#
# 1. flat_means: Shape/location parameters (alpha for Beta/Gamma)
#    - For Beta: this is shape1 (alpha)
#    - For Gamma: this is shape (alpha)
#
# 2. flat_sds: Scale parameters (beta for both)
#    - For Beta: this is shape2 (beta)
#    - For Gamma: this is rate (beta)
#
# 3. dists: Distribution type ("beta" or "gamma")
#    - Tells us which R function to call during sampling
#
# WHY 3 LISTS?
# ────────────
# - R's rbeta() and rgamma() require SHAPE PARAMETERS, not (mean, SD)
# - We convert (mean, SD) → (shape params) ONCE before sampling
# - During sampling (loop), we look up the shape params + dist type
# - This is MUCH faster than re-converting for each iteration
#
# EXAMPLE:
# ────────
# Parameter: coverage_ABCWY with mean=0.90, SD=0.05
# 
# Step 1: Convert to Beta parameters
#   beta_safe_local(0.90, 0.05) → alpha=153, beta=17
#
# Step 2: Store in lists
#   flat_means[["coverage_ABCWY"]] ← 153
#   flat_sds[["coverage_ABCWY"]]   ← 17
#   dists[["coverage_ABCWY"]]      ← "beta"
#
# Step 3: Sample (in loop, 1000 times)
#   for (i in 1:1000) {
#     value <- rbeta(1, flat_means[["coverage_ABCWY"]], 
#                       flat_sds[["coverage_ABCWY"]])
#   }
#
# VECTOR PARAMETERS:
# ──────────────────
# For time-varying params (e.g., p_B with 89 elements):
#   - We store as p_B1, p_B2, ..., p_B89
#   - Each gets its own entry in the 3 lists
#   - After sampling, we reconstruct the vector

generate_psa_samples <- function(params, n_sim = 1000, seed = 2025) {
  set.seed(seed)
  log_info(paste("Generating", n_sim, "PSA samples using Excel SD values..."))
  
  flat_means <- list() # Will store shape1/shape parameters
  flat_sds <- list() # Will store shape2/rate parameters
  dists <- list() # Will store "beta" or "gamma"
  
  # This function:
  # 1. Checks if SD is valid (exists, not NA, > 0)
  # 2. If no valid SD → parameter stays FIXED at mean value
  # 3. If valid SD → converts (mean, SD) to shape params
  # 4. Adds to the 3 lists
  #
  # <<super>> syntax: modifies parent environment lists
  
  .add_scalar <- function(nm, mean_val, sd_val, dist) {
    # Check if parameter should be fixed (no SD)
    if (is.null(sd_val) || is.na(sd_val) || sd_val <= 0) {
      return()  # Parameter will remain FIXED at mean_val throughout all PSA iterations
    }
    
    if (is.null(mean_val) || is.na(mean_val)) {
      return()
    }
    
    if (dist == "beta") {
      if (mean_val < 1e-5) {
        # For very small probabilities, use gamma instead
        adj <- gamma_safe_local(mean_val, sd_val)
        if (is.null(adj)) return()
        
        flat_means[[nm]] <<- adj$alpha
        flat_sds[[nm]] <<- adj$beta
        dists[[nm]] <<- "gamma"
      } else {
        adj <- beta_safe_local(mean_val, sd_val)
        if (is.null(adj)) return()
        
        flat_means[[nm]] <<- adj$alpha
        flat_sds[[nm]] <<- adj$beta
        dists[[nm]] <<- "beta"
      }
    } else if (dist == "gamma") {
      adj <- gamma_safe_local(mean_val, sd_val)
      if (is.null(adj)) return()
      
      flat_means[[nm]] <<- adj$alpha
      flat_sds[[nm]] <<- adj$beta
      dists[[nm]] <<- "gamma"
    }
  }
  ## ---- SOCIETAL PRODUCTIVITY COSTS (Gamma distribution) ----
  
  # Only include if base case has societal costs defined
  if (!is.null(params$c_prod_Scarring)) {
    societal_cost_params <- c(
      "c_prod_Scarring", "c_prod_Single_Amput", "c_prod_Multiple_Amput",
      "c_prod_Neuro_Disability", "c_prod_Hearing_Loss", "c_prod_Renal_Failure",
      "c_prod_Seizure", "c_prod_Paralysis"
    )
    
    for (nm in societal_cost_params) {
      sd_nm <- paste0("sd_", nm)
      # Only add to PSA if SD exists and is > 0
      # Otherwise, societal costs remain FIXED at deterministic values
      if (!is.null(params[[sd_nm]])) {
        .add_scalar(nm, params[[nm]], params[[sd_nm]], "gamma")
      }
    }
  }
  
  ## ---- SOCIETAL CAREGIVER COSTS (Gamma distribution) ----  
  
  # Only include if base case has caregiver costs defined
  if (!is.null(params$c_caregiver_infection)) {
    # Caregiver cost during infection (one-time)
    .add_scalar("c_caregiver_infection", 
                params$c_caregiver_infection, 
                params$sd_c_caregiver_infection, 
                "gamma")
    
    # Caregiver cost for sequelae (annual recurring)
    .add_scalar("c_caregiver_sequelae", 
                params$c_caregiver_sequelae, 
                params$sd_c_caregiver_sequelae, 
                "gamma")
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
  
  # Administration cost
  .add_scalar("c_admin", params$c_admin, params$sd_c_admin, "gamma")
  
  # Infection costs (time-varying vector)
  .add_vector("c_IMD_infection", params$c_IMD_infection, params$sd_c_IMD_infection, "gamma")
  
  # IMPORTANT: Vaccine prices are NEVER included in PSA
  # They are always fixed at their deterministic values
  # Even if SD values exist in Excel, they are ignored
  log_info("Vaccine prices (c_MenABCWY, c_MenACWY, c_MenC, c_MenB) are FIXED in PSA")
  
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
  
  # VE parameters are usually provided as vectors (one per cycle)
  # If scalar, expand to vector
  ve_params <- c("ve_MenC", "ve_MenB", "ve_MenACWY", 
                 "ve_MenABCWY_forB", "ve_MenABCWY_forACWY")
  
  for (nm in ve_params) {
    if (!is.null(params[[nm]])) {
      sd_nm <- paste0("sd_", nm)
      .add_vector(nm, params[[nm]], params[[sd_nm]], "beta")
    }
  }
  
  ## ---- MULTIPLIER PARAMETERS (Gamma for ratios) ----
  
  # Infection rate multipliers
  for (sg in c("B", "C", "W", "Y")) {
    nm <- paste0("mult_p_", sg)
    sd_nm <- paste0("sd_", nm)
    if (!is.null(params[[nm]]) && !is.null(params[[sd_nm]])) {
      .add_scalar(nm, params[[nm]], params[[sd_nm]], "gamma")
    }
  }
  
  # CFR multipliers
  for (sg in c("B", "C", "W", "Y")) {
    nm <- paste0("mult_cfr_", sg)
    sd_nm <- paste0("sd_", nm)
    if (!is.null(params[[nm]]) && !is.null(params[[sd_nm]])) {
      .add_scalar(nm, params[[nm]], params[[sd_nm]], "gamma")
    }
  }
  
  # VE scalers
  for (sg in c("B", "C", "W", "Y")) {
    nm <- paste0("ve_scale_", sg)
    sd_nm <- paste0("sd_", nm)
    if (!is.null(params[[nm]]) && !is.null(params[[sd_nm]])) {
      .add_scalar(nm, params[[nm]], params[[sd_nm]], "gamma")
    }
  }
  
  # Other multipliers
  if (!is.null(params$mult_bg_mort)) {
    .add_scalar("mult_bg_mort", params$mult_bg_mort, params$sd_mult_bg_mort, "gamma")
  }
  if (!is.null(params$mult_c_IMD)) {
    .add_scalar("mult_c_IMD", params$mult_c_IMD, params$sd_mult_c_IMD, "gamma")
  }
  
  # Generate all samples for each iteration
  log_info(paste("Sampling", length(flat_means), "parameters across", n_sim, "iterations..."))

  # ============================================================
  # SAMPLING LOOP: Generate n_sim parameter sets 
  # ============================================================
  
  # For each iteration:
  # 1. Sample from distributions using the 3 lists
  # 2. Copy fixed parameters (those not in flat_means)
  # 3. Reconstruct vectors from individual samples
  # 4. Store complete parameter set
  samples <- vector("list", n_sim)
  
  for (i in seq_len(n_sim)) {
    if (i %% 100 == 0) log_info(paste("  PSA iteration", i, "of", n_sim))
    
    s <- list() # This iteration's parameter set
    
    # ── Step 1: Sample from distributions ──
    # This is where the 3 lists are used together!
    for (nm in names(flat_means)) {
      d <- dists[[nm]]            # Get distribution type
      param1 <- flat_means[[nm]]  # Get shape1/shape parameter
      param2 <- flat_sds[[nm]]    # Get shape2/rate parameter
      
      if (d == "beta") {
        # Beta distribution: rbeta(n, shape1, shape2)
        s[[nm]] <- rbeta(1, param1, param2) # used for probabilities/utilities
      } else if (d == "gamma") {
        # Gamma distribution: rgamma(n, shape, rate)
        s[[nm]] <- rgamma(1, param1, param2) # used for costs/multipliers
      }
    }
    
    # ── Step 2: Copy fixed parameters ──
    # Parameters NOT in flat_means stay at deterministic values
    fixed_scalar_params <- setdiff(names(params), names(flat_means))
    for (nm in fixed_scalar_params) {
      if (!is.list(params[[nm]])) {
        s[[nm]] <- params[[nm]]
      }
    }
    
    # ── Step 3: Ensure societal costs are present ──
    # Even if not sampled (no SD), they need to exist in the parameter set
    societal_cost_params <- c(
      "c_prod_Scarring", "c_prod_Single_Amput", "c_prod_Multiple_Amput",
      "c_prod_Neuro_Disability", "c_prod_Hearing_Loss", "c_prod_Renal_Failure",
      "c_prod_Seizure", "c_prod_Paralysis"
    )
    
    for (nm in societal_cost_params) {
      if (is.null(s[[nm]]) && !is.null(params[[nm]])) {
        # If not sampled and exists in base params, copy it
        s[[nm]] <- params[[nm]]   # Use deterministic value
      } else if (is.null(s[[nm]])) {
        # If doesn't exist at all, set to 0
        s[[nm]] <- 0
      }
    }
    
    samples[[i]] <- s
  }
  
  log_info(paste(n_sim, "PSA samples generated successfully"))
  
  # ============================================================
  # VE-COVERAGE RESCALING 
  # ============================================================
  # 
  # PURPOSE: Ensure population-level VE reflects coverage changes
  # 
  # RATIONALE: When coverage increases/decreases in PSA, the population-level
  # vaccine effectiveness should scale proportionally. Without this adjustment,
  # the model assumes same VE regardless of coverage level.
  #
  # EXAMPLE:
  # - Base case: coverage_ABCWY = 0.95, ve_MenABCWY = 0.85
  # - PSA sample: coverage_ABCWY = 0.80 (5 percentage points lower)
  # - Rescaled VE: 0.85 * (0.80/0.95) = 0.716 (proportionally lower)
  #
  # This aligns the modular model with the reference model's behavior
  # and ensures consistency with policy scenarios.
  #
  log_info("Applying VE-coverage rescaling to PSA samples...")
  
  for (i in seq_len(n_sim)) {
    s <- samples[[i]]
    
    # Get base parameters for rescaling calculations
    base_cov_ABCWY <- params$coverage_ABCWY
    base_cov_ACWY <- params$coverage_ACWY
    base_cov_C <- params$coverage_C
    base_cov_B <- params$coverage_B
    
    # Guard against division by zero
    eps <- 1e-8
    
    # ---- Rescale ACWY-based strategies ----
    
    # MenACWY strategy: rescale based on ACWY coverage
    if (!is.null(s$coverage_ACWY) && base_cov_ACWY > eps && 
        !is.null(params$ve_MenACWY)) {
      # For vector VE, rescale element-wise
      if (is.vector(params$ve_MenACWY) && length(params$ve_MenACWY) > 1) {
        # Time-varying VE: rescale each cycle
        s$ve_MenACWY <- params$ve_MenACWY * (s$coverage_ACWY / base_cov_ACWY)
      } else if (is.numeric(params$ve_MenACWY) && length(params$ve_MenACWY) == 1) {
        # Scalar VE: rescale the single value
        s$ve_MenACWY <- params$ve_MenACWY * (s$coverage_ACWY / base_cov_ACWY)
      }
      
      # Clamp to [0, 1]
      if (!is.null(s$ve_MenACWY)) {
        s$ve_MenACWY <- .clamp(s$ve_MenACWY, 0, 1)
      }
    }
    
    # MenABCWY strategy for ACWY component: rescale based on ABCWY coverage
    if (!is.null(s$coverage_ABCWY) && base_cov_ABCWY > eps && 
        !is.null(params$ve_MenABCWY_forACWY)) {
      if (is.vector(params$ve_MenABCWY_forACWY) && length(params$ve_MenABCWY_forACWY) > 1) {
        s$ve_MenABCWY_forACWY <- params$ve_MenABCWY_forACWY * (s$coverage_ABCWY / base_cov_ABCWY)
      } else if (is.numeric(params$ve_MenABCWY_forACWY) && length(params$ve_MenABCWY_forACWY) == 1) {
        s$ve_MenABCWY_forACWY <- params$ve_MenABCWY_forACWY * (s$coverage_ABCWY / base_cov_ABCWY)
      }
      
      if (!is.null(s$ve_MenABCWY_forACWY)) {
        s$ve_MenABCWY_forACWY <- .clamp(s$ve_MenABCWY_forACWY, 0, 1)
      }
    }
    
    # ---- Rescale C-only strategy ----
    
    # MenC strategy: rescale based on C coverage
    if (!is.null(s$coverage_C) && base_cov_C > eps && 
        !is.null(params$ve_MenC)) {
      if (is.vector(params$ve_MenC) && length(params$ve_MenC) > 1) {
        s$ve_MenC <- params$ve_MenC * (s$coverage_C / base_cov_C)
      } else if (is.numeric(params$ve_MenC) && length(params$ve_MenC) == 1) {
        s$ve_MenC <- params$ve_MenC * (s$coverage_C / base_cov_C)
      }
      
      if (!is.null(s$ve_MenC)) {
        s$ve_MenC <- .clamp(s$ve_MenC, 0, 1)
      }
    }
    
    # ---- Rescale B-containing strategies ----
    
    # MenB strategy: rescale based on B coverage
    if (!is.null(s$coverage_B) && base_cov_B > eps && 
        !is.null(params$ve_MenB)) {
      if (is.vector(params$ve_MenB) && length(params$ve_MenB) > 1) {
        s$ve_MenB <- params$ve_MenB * (s$coverage_B / base_cov_B)
      } else if (is.numeric(params$ve_MenB) && length(params$ve_MenB) == 1) {
        s$ve_MenB <- params$ve_MenB * (s$coverage_B / base_cov_B)
      }
      
      if (!is.null(s$ve_MenB)) {
        s$ve_MenB <- .clamp(s$ve_MenB, 0, 1)
      }
    }
    
    # MenABCWY strategy for B component: rescale based on ABCWY coverage
    if (!is.null(s$coverage_ABCWY) && base_cov_ABCWY > eps && 
        !is.null(params$ve_MenABCWY_forB)) {
      if (is.vector(params$ve_MenABCWY_forB) && length(params$ve_MenABCWY_forB) > 1) {
        s$ve_MenABCWY_forB <- params$ve_MenABCWY_forB * (s$coverage_ABCWY / base_cov_ABCWY)
      } else if (is.numeric(params$ve_MenABCWY_forB) && length(params$ve_MenABCWY_forB) == 1) {
        s$ve_MenABCWY_forB <- params$ve_MenABCWY_forB * (s$coverage_ABCWY / base_cov_ABCWY)
      }
      
      if (!is.null(s$ve_MenABCWY_forB)) {
        s$ve_MenABCWY_forB <- .clamp(s$ve_MenABCWY_forB, 0, 1)
      }
    }
    
    samples[[i]] <- s
  }
  
  log_info("VE-coverage rescaling complete")
  
  # ============================================================
  # END VE-COVERAGE RESCALING
  # ============================================================
  
  list(samples = samples, psa_df = NULL)
}

# ============================================================
# PARALLEL PSA EXECUTION
# ============================================================

# Run PSA with parallel processing (cleaned up & reproducible)
run_psa <- function(params, n_sim = 1000, wtp = 50000, seed = 2025) {
  # Determine cores (at least 1; cap at n_sim to avoid idle workers)
  n_cores <- parallel::detectCores()
  if (is.na(n_cores) || n_cores < 1) n_cores <- 1L
  n_cores <- max(1L, min(n_cores - 1L, n_sim))
  
  # Create & register cluster
  cl <- parallel::makeCluster(n_cores)
  on.exit({ try(parallel::stopCluster(cl), silent = TRUE) }, add = TRUE)
  doParallel::registerDoParallel(cl)
  
  # Reproducible RNG across workers
  parallel::clusterSetRNGStream(cl, iseed = seed)
  
  log_info(paste("Running", n_sim, "PSA simulations on", n_cores, "cores..."))
  
  # Generate parameter samples (also reproducible)
  samp <- generate_psa_samples(params, n_sim = n_sim, seed = seed)
  
  # One-time export of everything workers need
  export_names <- c(
    # parameter samples used by foreach body
    "samp",
    # model functions
    "eval_all_strategies", "eval_strategy", "build_transition_array",
    "run_markov", "get_serogroup_ve", "shift_ve", ".clamp",
    "get_state_costs", "get_state_utils", "vaccination_costs",
    # logging helpers
    "log_info", "log_warn", "log_error",
    # globals
    "v_strats", "v_states", "n_states",
    "v_discount_cost", "v_discount_qaly", "n_cohort"
  )
  parallel::clusterExport(cl, unique(export_names), envir = environment())
  
  # Ensure required packages are loaded on workers (for nested function calls)
  parallel::clusterEvalQ(cl, {
    suppressPackageStartupMessages({
      library(dplyr); library(tibble)
    })
    NULL
  })
  
  # Parallel simulation (no .export here → no duplicate-export warnings)
  res <- foreach::foreach(
    i = seq_len(n_sim),
    .combine = "rbind"
  ) %dopar% {
    r_i <- eval_all_strategies(samp$samples[[i]]) # here the sample is used
    tibble::tibble(
      sim = i,
      Strategy = names(r_i),
      Cost  = vapply(r_i, function(x) x$cost,  numeric(1)),
      QALYs = vapply(r_i, function(x) x$qalys, numeric(1))
    )
  }
  
  # Calculate NMB (base R to avoid extra dependencies here)
  res$NMB <- res$QALYs * wtp - res$Cost
  
  log_info("PSA completed successfully")
  
  # Return PSA results
  res
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