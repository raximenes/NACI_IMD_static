## ======================
## 3) Utility Functions
## ======================

# Clamp values between bounds
.clamp <- function(x, lo = 0, hi = 1) {
  pmin(pmax(x, lo), hi)
}

# Safe Beta parameterization for PSA
.beta_safe <- function(mu, sd, eps = 1e-6) {
  mu2 <- min(1 - eps, max(eps, as.numeric(mu)))
  sd_max <- sqrt(mu2 * (1 - mu2)) * 0.99
  sd2 <- min(as.numeric(sd), sd_max)
  if (!is.finite(sd2) || sd2 <= 0) return(NULL)
  list(mean = mu2, sd = sd2)
}

# Safe Gamma parameterization for PSA
.gamma_safe <- function(mu, sd) {
  mu <- as.numeric(mu)
  sd <- as.numeric(sd)
  if (!is.finite(mu) || !is.finite(sd) || mu <= 0 || sd <= 0) return(NULL)
  list(mean = mu, sd = sd)
}

# Linear VE decay calculation
calculate_all_ve <- function(data, n_cycles = 89, waning_type = "linear") {
  years <- 1:n_cycles
  
  out_list <- lapply(seq_len(nrow(data)), function(i) {
    vac      <- data$vaccine[i]
    base_ve  <- as.numeric(data$effectiveness[i])
    max_year <- as.integer(data$duration[i])
    
    # Compute effectiveness decay
    ve_vec <- vapply(years, function(yr) {
      if (yr >= 1 && yr <= max_year) {
        if (waning_type == "linear") {
          # Linear interpolation between two points (x0, y0) and (x1, y1)
          x0 <- 1; x1 <- max_year; y0 <- base_ve; y1 <- 0
          y0 + (y1 - y0) / (x1 - x0) * (yr - x0)
        } else if (waning_type == "exponential") {
          # Exponential decay (alternative model)
          decay_rate <- -log(0.01) / max_year  # 99% reduction at max_year
          base_ve * exp(-decay_rate * (yr - 1))
        } else {
          stop("Unknown waning_type: ", waning_type)
        }
      } else {
        0
      }
    }, numeric(1))
    
    tibble::tibble(Cycle = years, Vaccine = vac, Effectiveness = ve_vec)
  })
  
  dplyr::bind_rows(out_list)
}

# Shift VE vector by k years (for booster scenarios)
shift_ve <- function(ve_vec, k_shift) {
  if (k_shift <= 0) {
    return(ve_vec)
  } else {
    # Add zeros at the start and remove last k_shift elements
    c(rep(0, k_shift), head(ve_vec, -k_shift))
  }
}

# Normalize shares to sum to 1 (handles NAs properly)
.normalize_shares <- function(x, default_share) {
  # Remove NAs first
  x_clean <- x[!is.na(x)]
  
  if (length(x_clean) == 0) {
    log_warn("All shares are NA, using default")
    return(default_share)
  }
  
  s <- sum(x_clean)
  
  if (s <= 0) {
    log_warn("Sum of shares is zero or negative, using default")
    return(default_share)
  }
  
  # Detect if values are percentages (sum > 1.5 suggests percentages)
  if (s > 1.5) {
    x_clean <- x_clean / 100
    log_info("Shares converted from percentages to proportions")
  }
  
  # Normalize
  normalized <- x_clean / sum(x_clean)
  
  # Return full vector with NAs replaced
  result <- x
  result[!is.na(x)] <- normalized
  
  return(result)
}