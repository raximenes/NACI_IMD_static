# 04_inputs.R (UPDATED - COMPLETE VERSION)
# Read Excel inputs and build base-case parameters
# - Removed get_col, len_ok, conv functions
# - Using direct column names from Excel
# - Using darthtools::rate_to_prob for conversions
# - Reading SD values directly from Excel
# 
# IMPORTANT: Parameters WITHOUT SD values in Excel will remain FIXED in PSA
# EXCEPTION: Vaccine prices are ALWAYS FIXED regardless of SD presence
# - To make a parameter vary in PSA: provide a valid SD > 0
# - To keep a parameter fixed in PSA: set SD to 0, NA, or leave blank
# - This applies to both scalar and time-varying (vector) parameters


# - Discount rates loaded from Excel (d_c, d_e in model_settings)
# - Societal costs loaded from Excel (Societal_Costs sheet)
# - Vaccine prices ALWAYS FIXED in PSA (no SD used even if present)
# 





file_main  <- file.path(IN_RAW, "IMD Data.xls")
file_price <- file.path(IN_RAW, "vaccine_costs.xlsx")

get_base_params <- function() {
  
  # ============================================================
  # STEP 1: Read Excel file and load all sheets
  # ============================================================
  
  if (!file.exists(file_main)) {
    log_error(paste("Excel file not found:", file_main))
  }
  
  log_info(paste("Reading data from:", file_main))
  sheet_names <- readxl::excel_sheets(file_main)
  data_list <- lapply(sheet_names, function(s) readxl::read_excel(file_main, sheet = s))
  names(data_list) <- sheet_names
  
  # ============================================================
  # STEP 2: Helper functions (NO DEFAULTS - strict validation)
  # ============================================================
  
  # Get single value from name-value table
  getv <- function(df, nm, sheet_name = "unknown") {
    if (is.null(df) || !all(c("Name", "Value") %in% names(df))) {
      log_error(paste0("Sheet '", sheet_name, "' is missing or incorrectly formatted"))
    }
    v <- df$Value[df$Name == nm]
    if (length(v) == 0) {
      log_error(paste0("Parameter '", nm, "' not found in sheet '", sheet_name, "'"))
    }
    as.numeric(v[1])
  }
  
  # Get SD value from name-value table (returns NA if not found - for optional SDs)
  get_sd <- function(df, nm, sheet_name = "unknown") {
    if (is.null(df) || !"sd" %in% names(df)) {
      log_warn(paste0("Sheet '", sheet_name, "' does not have 'sd' column"))
      return(NA)
    }
    v <- df$sd[df$Name == nm]
    if (length(v) == 0) {
      log_warn(paste0("SD for '", nm, "' not found in sheet '", sheet_name, "'"))
      return(NA)
    }
    as.numeric(v[1])
  }
  
  # ============================================================
  # STEP 3: Model Settings
  # ============================================================
  
  ms <- data_list[["model_settings"]]
  if (is.null(ms)) {
    log_error("Sheet 'model_settings' not found in Excel file")
  }
  
  nC  <- getv(ms, "n_cycles", "model_settings")
  cyc <- getv(ms, "cycle_len", "model_settings")
  
  # Discount rates from Excel (d_c and d_e)
  dc  <- getv(ms, "d_c", "model_settings")
  de  <- getv(ms, "d_e", "model_settings")
  
  # Coverage rates
  cov_ABCWY <- getv(ms, "coverage_ABCWY", "model_settings")
  cov_ACWY  <- getv(ms, "coverage_ACWY", "model_settings")
  cov_C     <- getv(ms, "coverage_C", "model_settings")
  cov_B     <- getv(ms, "coverage_B", "model_settings")
  
  # Coverage SDs (for PSA)
  sd_cov_ABCWY <- get_sd(ms, "coverage_ABCWY", "model_settings")
  sd_cov_ACWY  <- get_sd(ms, "coverage_ACWY", "model_settings")
  sd_cov_C     <- get_sd(ms, "coverage_C", "model_settings")
  sd_cov_B     <- get_sd(ms, "coverage_B", "model_settings")
  
  log_info(paste("Model settings: n_cycles =", nC, ", cycle_length =", cyc))
  
  # ============================================================
  # STEP 4: Societal Productivity Costs (from Excel)
  # ============================================================
  
  soc_costs <- data_list[["Societal_Costs"]]
  if (is.null(soc_costs)) {
    log_warn("Sheet 'Societal_Costs' not found - using zero productivity costs")
    c_prod_societal <- c(
      Scarring         = 0,
      Single_Amput     = 0,
      Multiple_Amput   = 0,
      Neuro_Disability = 0,
      Hearing_Loss     = 0,
      Renal_Failure    = 0,
      Seizure          = 0,
      Paralysis        = 0
    )
  } else {
    if (!all(c("Name", "Value") %in% names(soc_costs))) {
      log_error("Sheet 'Societal_Costs' must have 'Name' and 'Value' columns")
    }
    
    getsc <- function(nm) {
      v <- soc_costs$Value[soc_costs$Name == nm]
      if (length(v) == 0) {
        log_warn(paste0("Societal cost for '", nm, "' not found - using 0"))
        return(0)
      }
      as.numeric(v[1])
    }
    
    c_prod_societal <- c(
      Scarring         = getsc("Scarring"),
      Single_Amput     = getsc("Single_Amput"),
      Multiple_Amput   = getsc("Multiple_Amput"),
      Neuro_Disability = getsc("Neuro_Disability"),
      Hearing_Loss     = getsc("Hearing_Loss"),
      Renal_Failure    = getsc("Renal_Failure"),
      Seizure          = getsc("Seizure"),
      Paralysis        = getsc("Paralysis")
    )
    
    log_info("Societal productivity costs loaded from Excel")
  }
  # Read SD for societal costs (if they should vary in PSA)
  # If no SD column exists, these costs remain FIXED
  get_sd_societal <- function(nm) {
    if (!"sd" %in% names(soc_costs)) {
      return(NA)
    }
    v <- soc_costs$sd[soc_costs$Name == nm]
    if (length(v) == 0) {
      return(NA)
    }
    as.numeric(v[1])
  }
  
  sd_c_prod_societal <- c(
    Scarring         = get_sd_societal("Scarring"),
    Single_Amput     = get_sd_societal("Single_Amput"),
    Multiple_Amput   = get_sd_societal("Multiple_Amput"),
    Neuro_Disability = get_sd_societal("Neuro_Disability"),
    Hearing_Loss     = get_sd_societal("Hearing_Loss"),
    Renal_Failure    = get_sd_societal("Renal_Failure"),
    Seizure          = get_sd_societal("Seizure"),
    Paralysis        = get_sd_societal("Paralysis")
  )
  # ============================================================
  # STEP 5: Vaccine Prices
  # ============================================================
  
  if (!file.exists(file_price)) {
    log_error(paste("Vaccine price file not found:", file_price))
  }
  
  log_info(paste("Reading vaccine prices from:", file_price))
  pr <- readxl::read_excel(file_price)
  
  if (!all(c("Name", "Value") %in% names(pr))) {
    log_error("vaccine_costs.xlsx must have 'Name' and 'Value' columns")
  }
  
  getp <- function(nm) {
    v <- pr$Value[pr$Name == nm]
    if (length(v) == 0) {
      log_error(paste0("Vaccine price '", nm, "' not found in vaccine_costs.xlsx"))
    }
    as.numeric(v[1])
  }
  
  get_price_sd <- function(nm) {
    if (!"sd" %in% names(pr)) {
      log_warn("No 'sd' column in vaccine_costs.xlsx - vaccine prices will not vary in PSA")
      return(NA)
    }
    v <- pr$sd[pr$Name == nm]
    if (length(v) == 0) {
      log_warn(paste0("SD for '", nm, "' not found - this vaccine price will not vary in PSA"))
      return(NA)
    }
    as.numeric(v[1])
  }
  
  c_MenABCWY <- getp("c_MenABCWY")
  c_MenACWY  <- getp("c_MenACWY")
  c_MenC     <- getp("c_MenC")
  c_MenB     <- getp("c_MenB")
  
  # IMPORTANT: Vaccine prices are ALWAYS FIXED in PSA
  # Set SD to NA regardless of what's in Excel
  sd_c_MenABCWY <- NA
  sd_c_MenACWY  <- NA
  sd_c_MenC     <- NA
  sd_c_MenB     <- NA
  
  log_info("Vaccine prices loaded (ALWAYS FIXED in PSA)")
  
  # ============================================================
  # STEP 6: Other Costs (admin, sequelae)
  # ============================================================
  
  costs_tbl <- data_list[["costs"]]
  if (is.null(costs_tbl)) {
    log_error("Sheet 'costs' not found in Excel file")
  }
  
  getc <- function(nm) {
    if (!all(c("Name", "Value") %in% names(costs_tbl))) {
      log_error("Sheet 'costs' is missing 'Name' or 'Value' columns")
    }
    v <- costs_tbl$Value[costs_tbl$Name == nm]
    if (length(v) == 0) {
      log_error(paste0("Cost parameter '", nm, "' not found in 'costs' sheet"))
    }
    as.numeric(v[1])
  }
  
  get_sd_c <- function(nm) {
    if (!"sd" %in% names(costs_tbl)) {
      log_warn("Sheet 'costs' does not have 'sd' column - costs will not vary in PSA")
      return(NA)
    }
    v <- costs_tbl$sd[costs_tbl$Name == nm]
    if (length(v) == 0) {
      log_warn(paste0("SD for cost '", nm, "' not found"))
      return(NA)
    }
    as.numeric(v[1])
  }
  
  c_admin          <- getc("c_admin")
  c_Healthy        <- getc("c_Healthy")
  c_Scarring       <- getc("c_Scarring")
  c_Single_Amput   <- getc("c_Single_Amput")
  c_Multiple_Amput <- getc("c_Multiple_Amput")
  c_Neuro_Disab    <- getc("c_Neuro_Disab")
  c_Hearing_Loss   <- getc("c_Hearing_Loss")
  c_Renal_Failure  <- getc("c_Renal_Failure")
  c_Seizure        <- getc("c_Seizure")
  c_Paralysis      <- getc("c_Paralysis")
  c_Dead           <- getc("c_Dead")
  
  sd_c_admin          <- get_sd_c("c_admin")
  sd_c_Healthy        <- get_sd_c("c_Healthy")
  sd_c_Scarring       <- get_sd_c("c_Scarring")
  sd_c_Single_Amput   <- get_sd_c("c_Single_Amput")
  sd_c_Multiple_Amput <- get_sd_c("c_Multiple_Amput")
  sd_c_Neuro_Disab    <- get_sd_c("c_Neuro_Disab")
  sd_c_Hearing_Loss   <- get_sd_c("c_Hearing_Loss")
  sd_c_Renal_Failure  <- get_sd_c("c_Renal_Failure")
  sd_c_Seizure        <- get_sd_c("c_Seizure")
  sd_c_Paralysis      <- get_sd_c("c_Paralysis")
  sd_c_Dead           <- get_sd_c("c_Dead")
  
  # ============================================================
  # STEP 7: Infection Costs (time-varying vector)
  # ============================================================
  
  cost_imd <- data_list[["cost_IMD"]]
  if (is.null(cost_imd)) {
    log_error("Sheet 'cost_IMD' not found in Excel file")
  }
  
  if (!"cost" %in% names(cost_imd)) {
    log_error("Sheet 'cost_IMD' must have 'cost' column")
  }
  
  c_IMD_infection <- as.numeric(cost_imd$cost)
  
  if (length(c_IMD_infection) != nC) {
    log_error(paste0("cost_IMD has ", length(c_IMD_infection), 
                     " rows but n_cycles is ", nC))
  }
  
  if ("sd" %in% names(cost_imd)) {
    sd_c_IMD_infection <- as.numeric(cost_imd$sd)
  } else {
    log_warn("Sheet 'cost_IMD' does not have 'sd' column - infection costs will not vary in PSA")
    sd_c_IMD_infection <- rep(NA, nC)
  }
  
  # ============================================================
  # STEP 8: Infection Probabilities (time-varying vectors)
  # ============================================================
  
  inf <- data_list[["infection"]]
  if (is.null(inf)) {
    log_error("Sheet 'infection' not found in Excel file")
  }
  
  required_inf_cols <- c(
    "SerogroupB_Infection", "sd_SerogroupB_Infection",
    "SerogroupC_Infection", "sd_SerogroupC_Infection",
    "SerogroupW_Infection", "sd_SerogroupW_Infection",
    "SerogroupY_Infection", "sd_SerogroupY_Infection"
  )
  
  missing_inf_cols <- setdiff(required_inf_cols, names(inf))
  if (length(missing_inf_cols) > 0) {
    log_error(paste0("Missing required columns in 'infection' sheet: ", 
                     paste(missing_inf_cols, collapse = ", ")))
  }
  
  # Read rates
  p_B_raw <- as.numeric(inf$SerogroupB_Infection)
  p_C_raw <- as.numeric(inf$SerogroupC_Infection)
  p_W_raw <- as.numeric(inf$SerogroupW_Infection)
  p_Y_raw <- as.numeric(inf$SerogroupY_Infection)
  
  # Convert rates to probabilities using darthtools
  p_B <- darthtools::rate_to_prob(p_B_raw, t = cyc)
  p_C <- darthtools::rate_to_prob(p_C_raw, t = cyc)
  p_W <- darthtools::rate_to_prob(p_W_raw, t = cyc)
  p_Y <- darthtools::rate_to_prob(p_Y_raw, t = cyc)
  
  # Validate probabilities
  if (any(p_B < 0 | p_B > 1, na.rm = TRUE)) {
    log_error("Invalid probabilities in p_B after conversion from rates")
  }
  if (any(p_C < 0 | p_C > 1, na.rm = TRUE)) {
    log_error("Invalid probabilities in p_C after conversion from rates")
  }
  if (any(p_W < 0 | p_W > 1, na.rm = TRUE)) {
    log_error("Invalid probabilities in p_W after conversion from rates")
  }
  if (any(p_Y < 0 | p_Y > 1, na.rm = TRUE)) {
    log_error("Invalid probabilities in p_Y after conversion from rates")
  }
  
  # Clamp to safe range for model stability
  p_B <- .clamp(p_B, 0, 0.999)
  p_C <- .clamp(p_C, 0, 0.999)
  p_W <- .clamp(p_W, 0, 0.999)
  p_Y <- .clamp(p_Y, 0, 0.999)
  
  # Read SD values
  sd_p_B <- as.numeric(inf$sd_SerogroupB_Infection)
  sd_p_C <- as.numeric(inf$sd_SerogroupC_Infection)
  sd_p_W <- as.numeric(inf$sd_SerogroupW_Infection)
  sd_p_Y <- as.numeric(inf$sd_SerogroupY_Infection)
  
  log_info("Infection probabilities loaded and converted from rates")
  
  # ============================================================
  # STEP 9: Mortality (background + CFRs)
  # ============================================================
  
  mort <- data_list[["mortality"]]
  if (is.null(mort)) {
    log_error("Sheet 'mortality' not found in Excel file")
  }
  
  required_mort_cols <- c(
    "Background_Mortality", "sd_Background_Mortality",
    "SerogroupB_Dead", "sd_SerogroupB_Dead",
    "SerogroupC_Dead", "sd_SerogroupC_Dead",
    "SerogroupW_Dead", "sd_SerogroupW_Dead",
    "SerogroupY_Dead", "sd_SerogroupY_Dead"
  )
  
  missing_mort_cols <- setdiff(required_mort_cols, names(mort))
  if (length(missing_mort_cols) > 0) {
    log_error(paste0("Missing required columns in 'mortality' sheet: ", 
                     paste(missing_mort_cols, collapse = ", ")))
  }
  
  # Background mortality (convert from rate)
  bg_raw <- as.numeric(mort$Background_Mortality)
  p_bg <- darthtools::rate_to_prob(bg_raw, t = cyc)
  
  if (any(p_bg < 0 | p_bg > 1, na.rm = TRUE)) {
    log_error("Invalid probabilities in background mortality after conversion")
  }
  
  p_bg <- .clamp(p_bg, 0, 0.999)
  sd_p_bg <- as.numeric(mort$sd_Background_Mortality)
  
  # Case fatality rates (already probabilities)
  p_B_CFR <- as.numeric(mort$SerogroupB_Dead)
  p_C_CFR <- as.numeric(mort$SerogroupC_Dead)
  p_W_CFR <- as.numeric(mort$SerogroupW_Dead)
  p_Y_CFR <- as.numeric(mort$SerogroupY_Dead)
  
  # Validate CFRs
  if (any(p_B_CFR < 0 | p_B_CFR > 1, na.rm = TRUE)) {
    log_error("Invalid CFR values for Serogroup B")
  }
  if (any(p_C_CFR < 0 | p_C_CFR > 1, na.rm = TRUE)) {
    log_error("Invalid CFR values for Serogroup C")
  }
  if (any(p_W_CFR < 0 | p_W_CFR > 1, na.rm = TRUE)) {
    log_error("Invalid CFR values for Serogroup W")
  }
  if (any(p_Y_CFR < 0 | p_Y_CFR > 1, na.rm = TRUE)) {
    log_error("Invalid CFR values for Serogroup Y")
  }
  
  p_B_CFR <- .clamp(p_B_CFR, 0, 0.999)
  p_C_CFR <- .clamp(p_C_CFR, 0, 0.999)
  p_W_CFR <- .clamp(p_W_CFR, 0, 0.999)
  p_Y_CFR <- .clamp(p_Y_CFR, 0, 0.999)
  
  # Read SD values
  sd_p_B_CFR <- as.numeric(mort$sd_SerogroupB_Dead)
  sd_p_C_CFR <- as.numeric(mort$sd_SerogroupC_Dead)
  sd_p_W_CFR <- as.numeric(mort$sd_SerogroupW_Dead)
  sd_p_Y_CFR <- as.numeric(mort$sd_SerogroupY_Dead)
  
  log_info("Mortality parameters loaded")
  
  # ============================================================
  # STEP 10: Sequelae Probabilities and Shares
  # ============================================================
  
  seq_tbl <- data_list[["sequelae_probp_IMD"]]
  if (is.null(seq_tbl)) {
    log_error("Sheet 'sequelae_probp_IMD' not found in Excel file")
  }
  
  if (!all(c("Name", "Value") %in% names(seq_tbl))) {
    log_error("Sheet 'sequelae_probp_IMD' must have 'Name' and 'Value' columns")
  }
  
  gv_seq <- function(nm) {
    v <- seq_tbl$Value[seq_tbl$Name == nm]
    if (length(v) == 0) {
      log_error(paste0("Parameter '", nm, "' not found in 'sequelae_probp_IMD' sheet"))
    }
    as.numeric(v[1])
  }
  
  get_sd_seq <- function(nm) {
    if (!"sd" %in% names(seq_tbl)) {
      log_warn("Sheet 'sequelae_probp_IMD' does not have 'sd' column")
      return(NA)
    }
    v <- seq_tbl$sd[seq_tbl$Name == nm]
    if (length(v) == 0) {
      log_warn(paste0("SD for '", nm, "' not found in 'sequelae_probp_IMD' sheet"))
      return(NA)
    }
    as.numeric(v[1])
  }
  
  # Overall sequelae probability - WITH ENHANCED VALIDATION
  p_seq_overall <- gv_seq("p_IMD_overall")
  
  # CRITICAL VALIDATION: Check for NA first
  if (is.na(p_seq_overall)) {
    log_error(
      "CRITICAL ERROR: p_IMD_overall is NA! ",
      "Check that Excel sheet 'sequelae_probp_IMD' contains parameter 'p_IMD_overall' with valid value"
    )
  }
  
  # Range validation
  if (p_seq_overall < 0 || p_seq_overall > 1) {
    log_error(
      sprintf(
        "CRITICAL ERROR: p_IMD_overall=%.6f is outside valid range [0,1]",
        p_seq_overall
      )
    )
  }
  
  # WARNING: If sequelae probability is zero
  if (p_seq_overall == 0) {
    log_warn(
      "WARNING: p_IMD_overall is EXACTLY ZERO. ",
      "NO SEQUELAE CASES will be generated. ",
      "Check if this is intentional."
    )
  }
  
  sd_p_seq_overall <- get_sd_seq("p_IMD_overall")
  
  # Sequelae shares - WITH ENHANCED VALIDATION
  shares_raw <- c(
    Scarring         = gv_seq("p_IMD_Scarring"),
    Single_Amput     = gv_seq("p_IMD_Single_Amput"),
    Multiple_Amput   = gv_seq("p_IMD_Multiple_Amput"),
    Neuro_Disability = gv_seq("p_IMD_Neuro_Disability"),
    Hearing_Loss     = gv_seq("p_IMD_Hearing_Loss"),
    Renal_Failure    = gv_seq("p_IMD_Renal_Failure"),
    Seizure          = gv_seq("p_IMD_Seizure"),
    Paralysis        = gv_seq("p_IMD_Paralysis")
  )
  
  # Check for NA values in shares
  if (any(is.na(shares_raw))) {
    na_types <- names(shares_raw)[is.na(shares_raw)]
    log_error(
      sprintf(
        "CRITICAL ERROR: Sequelae shares contain NA values: %s. ",
        paste(na_types, collapse = ", ")
      ),
      "Check Excel sheet 'sequelae_probp_IMD' for missing entries."
    )
  }
  
  # Check that sum of shares is positive
  shares_sum <- sum(shares_raw)
  if (shares_sum <= 0) {
    log_error(
      sprintf(
        "CRITICAL ERROR: Sum of sequelae shares is %.6f (must be > 0)",
        shares_sum
      )
    )
  }
  
  # Normalize shares to sum to 1 - SAFE VERSION
  w_seq <- shares_raw / shares_sum
  
  # Final validation: check no NA values remain
  if (any(is.na(w_seq))) {
    log_error("CRITICAL ERROR: w_seq contains NA values after normalization!")
  }
  
  # Log sequelae distribution for audit trail
  log_info("Sequelae distribution (normalized):")
  for (sq_name in names(w_seq)) {
    log_info(sprintf("  %s: %.4f", sq_name, w_seq[sq_name]))
  }
  
  # Read SD for sequelae shares
  sd_shares <- c(
    Scarring         = get_sd_seq("p_IMD_Scarring"),
    Single_Amput     = get_sd_seq("p_IMD_Single_Amput"),
    Multiple_Amput   = get_sd_seq("p_IMD_Multiple_Amput"),
    Neuro_Disability = get_sd_seq("p_IMD_Neuro_Disability"),
    Hearing_Loss     = get_sd_seq("p_IMD_Hearing_Loss"),
    Renal_Failure    = get_sd_seq("p_IMD_Renal_Failure"),
    Seizure          = get_sd_seq("p_IMD_Seizure"),
    Paralysis        = get_sd_seq("p_IMD_Paralysis")
  )
  
  log_info("Sequelae parameters loaded and validated")
  
  # ============================================================
  # STEP 11: Utilities
  # ============================================================
  
  utab <- data_list[["utilities"]]
  if (is.null(utab)) {
    log_error("Sheet 'utilities' not found in Excel file")
  }
  
  if (!all(c("Name", "Value") %in% names(utab))) {
    log_error("Sheet 'utilities' must have 'Name' and 'Value' columns")
  }
  
  getu <- function(nm) {
    v <- utab$Value[utab$Name == nm]
    if (length(v) == 0) {
      log_error(paste0("Utility parameter '", nm, "' not found in 'utilities' sheet"))
    }
    as.numeric(v[1])
  }
  
  get_sd_u <- function(nm) {
    if (!"sd" %in% names(utab)) {
      log_warn("Sheet 'utilities' does not have 'sd' column - utilities will not vary in PSA")
      return(NA)
    }
    v <- utab$sd[utab$Name == nm]
    if (length(v) == 0) {
      log_warn(paste0("SD for utility '", nm, "' not found"))
      return(NA)
    }
    as.numeric(v[1])
  }
  
  u_Healthy          <- getu("u_Healthy")
  u_IMD              <- getu("u_IMD")
  u_Scarring         <- getu("u_Scarring")
  u_Single_Amput     <- getu("u_Single_Amput")
  u_Multiple_Amput   <- getu("u_Multiple_Amput")
  u_Neuro_Disability <- getu("u_Neuro_Disability")
  u_Hearing_Loss     <- getu("u_Hearing_Loss")
  u_Renal_Failure    <- getu("u_Renal_Failure")
  u_Seizure          <- getu("u_Seizure")
  u_Paralysis        <- getu("u_Paralysis")
  u_Dead             <- getu("u_Dead")
  
  # Validate utilities
  util_values <- c(u_Healthy, u_IMD, u_Scarring, u_Single_Amput, u_Multiple_Amput,
                   u_Neuro_Disability, u_Hearing_Loss, u_Renal_Failure, u_Seizure, 
                   u_Paralysis, u_Dead)
  if (any(util_values < 0 | util_values > 1, na.rm = TRUE)) {
    log_error("All utility values must be between 0 and 1")
  }
  
  # Read SD (NA for fixed values like u_Healthy and u_Dead)
  sd_u_Healthy          <- get_sd_u("u_Healthy")
  sd_u_IMD              <- get_sd_u("u_IMD")
  sd_u_Scarring         <- get_sd_u("u_Scarring")
  sd_u_Single_Amput     <- get_sd_u("u_Single_Amput")
  sd_u_Multiple_Amput   <- get_sd_u("u_Multiple_Amput")
  sd_u_Neuro_Disability <- get_sd_u("u_Neuro_Disability")
  sd_u_Hearing_Loss     <- get_sd_u("u_Hearing_Loss")
  sd_u_Renal_Failure    <- get_sd_u("u_Renal_Failure")
  sd_u_Seizure          <- get_sd_u("u_Seizure")
  sd_u_Paralysis        <- get_sd_u("u_Paralysis")
  sd_u_Dead             <- get_sd_u("u_Dead")
  
  log_info("Utility parameters loaded")
  
  # ============================================================
  # STEP 12: Vaccine Effectiveness
  # ============================================================
  
  ve_df <- data_list[["vac_effect"]]
  if (is.null(ve_df)) {
    log_error("Sheet 'vac_effect' not found in Excel file")
  }
  
  required_ve_cols <- c("Name", "Value", "sd", "duration")
  missing_ve_cols <- setdiff(required_ve_cols, names(ve_df))
  if (length(missing_ve_cols) > 0) {
    log_error(paste0("Missing required columns in 'vac_effect' sheet: ", 
                     paste(missing_ve_cols, collapse = ", ")))
  }
  
  vaccine_names <- c("MenABCWY_for_SeroACWY", "MenABCWY_for_SeroB", "MenACWY", "MenC", "MenB")
  
  # Check all vaccines present
  missing_vaccines <- setdiff(vaccine_names, ve_df$Name)
  if (length(missing_vaccines) > 0) {
    log_error(paste0("Missing vaccines in 'vac_effect' sheet: ", 
                     paste(missing_vaccines, collapse = ", ")))
  }
  
  ve_data <- tibble::tibble(
    vaccine          = vaccine_names,
    effectiveness    = as.numeric(ve_df$Value[match(vaccine_names, ve_df$Name)]),
    sd_effectiveness = as.numeric(ve_df$sd[match(vaccine_names, ve_df$Name)]),
    duration         = as.integer(ve_df$duration[match(vaccine_names, ve_df$Name)])
  )
  
  # Check for NAs
  if (any(is.na(ve_data$effectiveness))) {
    log_error("NA values found in vaccine effectiveness")
  }
  if (any(is.na(ve_data$duration))) {
    log_error("NA values found in vaccine duration")
  }
  
  # Calculate VE over time with linear waning
  ve_all <- calculate_all_ve(ve_data, n_cycles = nC)
  
  # Extract VE vectors and adjust by coverage
  ve_MenABCWY_forACWY <- ve_all$Effectiveness[ve_all$Vaccine == "MenABCWY_for_SeroACWY"] * cov_ABCWY
  ve_MenABCWY_forB    <- ve_all$Effectiveness[ve_all$Vaccine == "MenABCWY_for_SeroB"]    * cov_ABCWY
  ve_MenACWY          <- ve_all$Effectiveness[ve_all$Vaccine == "MenACWY"]               * cov_ACWY
  ve_MenC             <- ve_all$Effectiveness[ve_all$Vaccine == "MenC"]                  * cov_C
  ve_MenB             <- ve_all$Effectiveness[ve_all$Vaccine == "MenB"]                  * cov_B
  
  log_info("Vaccine effectiveness parameters loaded")
  
  # ============================================================
  # STEP 13: Build Complete Parameter List
  # ============================================================
  
  params <- list(
    # Model settings
    n_cycles = nC, 
    cycle_length = cyc, 
    d_c = dc, 
    d_e = de,
    perspective = perspective,
    
    # Coverage
    coverage_ABCWY = cov_ABCWY, 
    coverage_ACWY = cov_ACWY, 
    coverage_C = cov_C, 
    coverage_B = cov_B,
    sd_coverage_ABCWY = sd_cov_ABCWY,
    sd_coverage_ACWY = sd_cov_ACWY,
    sd_coverage_C = sd_cov_C,
    sd_coverage_B = sd_cov_B,
    
    # Vaccine costs
    c_MenABCWY = c_MenABCWY, 
    c_MenACWY = c_MenACWY, 
    c_MenC = c_MenC, 
    c_MenB = c_MenB,
    sd_c_MenABCWY = sd_c_MenABCWY,
    sd_c_MenACWY = sd_c_MenACWY,
    sd_c_MenC = sd_c_MenC,
    sd_c_MenB = sd_c_MenB,
    
    # Other costs
    c_admin = c_admin,
    sd_c_admin = sd_c_admin,
    c_Healthy = c_Healthy,
    sd_c_Healthy = sd_c_Healthy,
    c_Scarring = c_Scarring,
    sd_c_Scarring = sd_c_Scarring,
    c_Single_Amput = c_Single_Amput,
    sd_c_Single_Amput = sd_c_Single_Amput,
    c_Multiple_Amput = c_Multiple_Amput,
    sd_c_Multiple_Amput = sd_c_Multiple_Amput,
    c_Neuro_Disab = c_Neuro_Disab,
    sd_c_Neuro_Disab = sd_c_Neuro_Disab,
    c_Hearing_Loss = c_Hearing_Loss,
    sd_c_Hearing_Loss = sd_c_Hearing_Loss,
    c_Renal_Failure = c_Renal_Failure,
    sd_c_Renal_Failure = sd_c_Renal_Failure,
    c_Seizure = c_Seizure,
    sd_c_Seizure = sd_c_Seizure,
    c_Paralysis = c_Paralysis,
    sd_c_Paralysis = sd_c_Paralysis,
    c_Dead = c_Dead,
    sd_c_Dead = sd_c_Dead,
    
    # Infection costs (time-varying)
    c_IMD_infection = c_IMD_infection,
    sd_c_IMD_infection = sd_c_IMD_infection,
    
    # Infection probabilities (time-varying)
    p_B = p_B,
    sd_p_B = sd_p_B,
    p_C = p_C,
    sd_p_C = sd_p_C,
    p_W = p_W,
    sd_p_W = sd_p_W,
    p_Y = p_Y,
    sd_p_Y = sd_p_Y,
    
    # Case fatality rates (time-varying)
    p_B_DeadIMD = p_B_CFR,
    sd_p_B_DeadIMD = sd_p_B_CFR,
    p_C_DeadIMD = p_C_CFR,
    sd_p_C_DeadIMD = sd_p_C_CFR,
    p_W_DeadIMD = p_W_CFR,
    sd_p_W_DeadIMD = sd_p_W_CFR,
    p_Y_DeadIMD = p_Y_CFR,
    sd_p_Y_DeadIMD = sd_p_Y_CFR,
    
    # Background mortality (time-varying)
    p_bg_mort = p_bg,
    sd_p_bg_mort = sd_p_bg,
    
    # Sequelae
    p_seq_overall = p_seq_overall,
    sd_p_seq_overall = sd_p_seq_overall,
    w_seq = w_seq,
    sd_w_seq = sd_shares,
    
    # Utilities
    u_Healthy = u_Healthy,
    sd_u_Healthy = sd_u_Healthy,
    u_IMD = u_IMD,
    sd_u_IMD = sd_u_IMD,
    u_Scarring = u_Scarring,
    sd_u_Scarring = sd_u_Scarring,
    u_Single_Amput = u_Single_Amput,
    sd_u_Single_Amput = sd_u_Single_Amput,
    u_Multiple_Amput = u_Multiple_Amput,
    sd_u_Multiple_Amput = sd_u_Multiple_Amput,
    u_Neuro_Disability = u_Neuro_Disability,
    sd_u_Neuro_Disability = sd_u_Neuro_Disability,
    u_Hearing_Loss = u_Hearing_Loss,
    sd_u_Hearing_Loss = sd_u_Hearing_Loss,
    u_Renal_Failure = u_Renal_Failure,
    sd_u_Renal_Failure = sd_u_Renal_Failure,
    u_Seizure = u_Seizure,
    sd_u_Seizure = sd_u_Seizure,
    u_Paralysis = u_Paralysis,
    sd_u_Paralysis = sd_u_Paralysis,
    u_Dead = u_Dead,
    sd_u_Dead = sd_u_Dead,
    
    # Vaccine effectiveness (vectors) - INTERNAL NAMING
    ve_MenABCWY_forACWY = ve_MenABCWY_forACWY,
    ve_MenABCWY_forB = ve_MenABCWY_forB,
    ve_MenACWY = ve_MenACWY,
    ve_MenC = ve_MenC,
    ve_MenB = ve_MenB,
    
    # Vaccine effectiveness (vectors) - EXTERNAL NAMING (for PSA/downstream code)
    ve_base_MenABCWY_for_SeroACWY = ve_MenABCWY_forACWY,
    ve_base_MenABCWY_for_SeroB = ve_MenABCWY_forB,
    ve_base_MenACWY = ve_MenACWY,
    ve_base_MenC = ve_MenC,
    ve_base_MenB = ve_MenB,
    
    # SDs for VE parameters (for PSA uncertainty)
    sd_ve_base_MenABCWY_for_SeroACWY = ve_data$sd_effectiveness[ve_data$vaccine == "MenABCWY_for_SeroACWY"],
    sd_ve_base_MenABCWY_for_SeroB = ve_data$sd_effectiveness[ve_data$vaccine == "MenABCWY_for_SeroB"],
    sd_ve_base_MenACWY = ve_data$sd_effectiveness[ve_data$vaccine == "MenACWY"],
    sd_ve_base_MenC = ve_data$sd_effectiveness[ve_data$vaccine == "MenC"],
    sd_ve_base_MenB = ve_data$sd_effectiveness[ve_data$vaccine == "MenB"],
    
    # Store ve_data for reference
    ve_data = ve_data,  # Store for PSA recalculation
    
    # Societal productivity costs (from Excel)
    c_prod_Scarring = c_prod_societal["Scarring"],
    c_prod_Single_Amput = c_prod_societal["Single_Amput"],
    c_prod_Multiple_Amput = c_prod_societal["Multiple_Amput"],
    c_prod_Neuro_Disability = c_prod_societal["Neuro_Disability"],
    c_prod_Hearing_Loss = c_prod_societal["Hearing_Loss"],
    c_prod_Renal_Failure = c_prod_societal["Renal_Failure"],
    c_prod_Seizure = c_prod_societal["Seizure"],
    c_prod_Paralysis = c_prod_societal["Paralysis"],
    
    # SDs for societal productivity costs (for PSA)
    sd_c_prod_Scarring = sd_c_prod_societal["Scarring"],
    sd_c_prod_Single_Amput = sd_c_prod_societal["Single_Amput"],
    sd_c_prod_Multiple_Amput = sd_c_prod_societal["Multiple_Amput"],
    sd_c_prod_Neuro_Disability = sd_c_prod_societal["Neuro_Disability"],
    sd_c_prod_Hearing_Loss = sd_c_prod_societal["Hearing_Loss"],
    sd_c_prod_Renal_Failure = sd_c_prod_societal["Renal_Failure"],
    sd_c_prod_Seizure = sd_c_prod_societal["Seizure"],
    sd_c_prod_Paralysis = sd_c_prod_societal["Paralysis"],
    
    # Multipliers for sensitivity analysis (OWSA)
    mult_p_B = 1.0,
    mult_p_C = 1.0,
    mult_p_W = 1.0,
    mult_p_Y = 1.0,
    mult_cfr_B = 1.0,
    mult_cfr_C = 1.0,
    mult_cfr_W = 1.0,
    mult_cfr_Y = 1.0,
    mult_bg_mort = 1.0,
    mult_c_IMD = 1.0,
    ve_scale_B = 1.0,
    ve_scale_C = 1.0,
    ve_scale_W = 1.0,
    ve_scale_Y = 1.0
  )
  
  log_info("Base parameters loaded successfully (with SDs from Excel)")
  log_info(paste("  - Parameters with valid SD will vary in PSA"))
  log_info(paste("  - Parameters without SD or SD=0/NA will remain FIXED in PSA"))
  
  return(params)
}