# 04_inputs.R
# Read Excel inputs and build base-case parameters

file_main  <- file.path(IN_RAW, "IMD Data.xls")
file_price <- file.path(IN_RAW, "vaccine_costs.xlsx")

get_base_params <- function() {
  if (!file.exists(file_main))  log_error(paste("Excel file not found:", file_main))
  
  log_info(paste("Reading data from:", file_main))
  sheet_names <- readxl::excel_sheets(file_main)
  data_list <- lapply(sheet_names, function(s) readxl::read_excel(file_main, sheet = s))
  names(data_list) <- sheet_names
  
  getv <- function(df, nm, default) {
    if (is.null(df) || !all(c("Name", "Value") %in% names(df))) return(default)
    v <- df$Value[df$Name == nm]; if (length(v) == 0) default else as.numeric(v[1])
  }
  
  # ---- Model settings ----
  ms  <- data_list[["model_settings"]]
  nC  <- getv(ms, "n_cycles", n_cycles)
  cyc <- getv(ms, "cycle_len", 1)
  dc  <- getv(ms, "d_c", 0.015)
  de  <- getv(ms, "d_e", 0.015)
  
  cov_ABCWY <- getv(ms, "coverage_ABCWY", 0.90)
  cov_ACWY  <- getv(ms, "coverage_ACWY",  0.94)
  cov_C     <- getv(ms, "coverage_C",     0.93)
  cov_B     <- getv(ms, "coverage_B",     0.60)
  
  # ---- Costs + prices ----
  costs_tbl <- data_list[["costs"]]
  c_admin <- getv(costs_tbl, "c_admin", 14.70)
  
  if (file.exists(file_price)) {
    log_info(paste("Reading vaccine prices from:", file_price))
    pr <- readxl::read_excel(file_price)
    getp <- function(nm, default) { v <- pr$Value[pr$Name == nm]; if (length(v) == 0) default else as.numeric(v[1]) }
    c_MenABCWY <- getp("c_MenABCWY", 92)
    c_MenACWY  <- getp("c_MenACWY",  32)
    c_MenC     <- getp("c_MenC",     13)
    c_MenB     <- getp("c_MenB",     82)
  } else {
    c_MenABCWY <- 92; c_MenACWY <- 32; c_MenC <- 13; c_MenB <- 82
    log_warn("vaccine_costs.xlsx not found. Using default prices (CAD).")
  }
  
  
  ## ---- Infection rates/probabilities ----
  inf <- data_list[["infection"]]
  
  # Helper to get column by multiple possible names
  get_col <- function(df, candidates) {
    idx <- which(names(df) %in% candidates)
    if (length(idx) == 0) return(NULL)
    as.numeric(df[[idx[1]]])
  }
  
  # Improved conversion function for rates/probabilities
  conv <- function(x) {
    if (is.null(x)) return(rep(0, nC))
    
    # Check if values look like percentages (all > 1 but < 100)
    if (all(x > 1 & x <= 100, na.rm = TRUE)) {
      log_info("Converting percentages to proportions")
      x <- x / 100
    }
    
    # If still > 1, treat as rates and convert
    if (any(x > 1, na.rm = TRUE)) {
      log_info("Converting rates to probabilities")
      x <- darthtools::rate_to_prob(x, t = cyc)
    }
    
    .clamp(x, 0, 1)
  }
  
  # Ensure vectors are correct length
  len_ok <- function(v) {
    if (length(v) < nC) {
      rep(v, length.out = nC)
    } else {
      head(v, nC)
    }
  }
  
  p_B <- len_ok(conv(get_col(inf, c("SerogroupB_Infection", "r_SerogroupB_Infection", "p_SerogroupB_Infection"))))
  p_C <- len_ok(conv(get_col(inf, c("SerogroupC_Infection", "r_SerogroupC_Infection", "p_SerogroupC_Infection"))))
  p_W <- len_ok(conv(get_col(inf, c("SerogroupW_Infection", "r_SerogroupW_Infection", "p_SerogroupW_Infection"))))
  p_Y <- len_ok(conv(get_col(inf, c("SerogroupY_Infection", "r_SerogroupY_Infection", "p_SerogroupY_Infection"))))
  
  ## ---- Mortality ----
  mort <- data_list[["mortality"]]
  bg <- get_col(mort, c("Background_Mortality", "p_Background_Mortality", "r_Background_Mortality"))
  p_bg <- len_ok(conv(bg))
  
  to_prob <- function(x, default) {
    if (is.null(x)) return(rep(default, nC))
    conv(x)
  }
  
  p_B_CFR <- len_ok(to_prob(get_col(mort, c("SerogroupB_Dead", "p_SerogroupB_Dead", "r_SerogroupB_Dead")), 0.10))
  p_C_CFR <- len_ok(to_prob(get_col(mort, c("SerogroupC_Dead", "p_SerogroupC_Dead", "r_SerogroupC_Dead")), 0.15))
  p_W_CFR <- len_ok(to_prob(get_col(mort, c("SerogroupW_Dead", "p_SerogroupW_Dead", "r_SerogroupW_Dead")), 0.12))
  p_Y_CFR <- len_ok(to_prob(get_col(mort, c("SerogroupY_Dead", "p_SerogroupY_Dead", "r_SerogroupY_Dead")), 0.11))
  
  ## ---- Sequelae shares and overall probability ----
  seq_tbl <- data_list[["sequelae_probp_IMD"]]
  gv <- function(nm, default) {
    v <- seq_tbl$Value[seq_tbl$Name == nm]
    if (length(v) == 0) default else as.numeric(v[1])
  }
  
  shares_raw <- c(
    Scarring         = gv("p_IMD_Scarring",         0.10),
    Single_Amput     = gv("p_IMD_Single_Amput",     0.08),
    Multiple_Amput   = gv("p_IMD_Multiple_Amput",   0.04),
    Neuro_Disability = gv("p_IMD_Neuro_Disability", 0.20),
    Hearing_Loss     = gv("p_IMD_Hearing_Loss",     0.25),
    Renal_Failure    = gv("p_IMD_Renal_Failure",    0.05),
    Seizure          = gv("p_IMD_Seizure",          0.18),
    Paralysis        = gv("p_IMD_Paralysis",        0.10)
  )
  
  w_seq <- .normalize_shares(shares_raw, default_share = shares_raw / sum(shares_raw))
  p_seq_overall <- .clamp(gv("p_IMD_overall", 0.25), 0, 1)
  
  ## ---- Vaccine effectiveness (VE) data with linear waning ----
  ve_df <- data_list[["vac_effect"]]
  
  if (is.null(ve_df)) {
    log_warn("vac_effect sheet not found. Using default VE values.")
    ve_df <- tibble::tibble(
      Name = c("MenABCWY_for_SeroACWY", "MenABCWY_for_SeroB", "MenACWY", "MenC", "MenB"),
      Value = c(0.85, 0.85, 0.85, 0.85, 0.75),
      duration = c(10, 5, 10, 10, 5)
    )
  }
  
  # FIXED: Define vaccine_names before using it
  vaccine_names <- c("MenABCWY_for_SeroACWY", "MenABCWY_for_SeroB", "MenACWY", "MenC", "MenB")
  
  ve_data <- tibble::tibble(
    vaccine       = vaccine_names,
    effectiveness = as.numeric(ve_df$Value[match(vaccine_names, ve_df$Name)]),
    duration      = as.integer(ve_df$duration[match(vaccine_names, ve_df$Name)])
  )
  
  # Replace any NAs with defaults
  ve_data$effectiveness[is.na(ve_data$effectiveness)] <- c(0.85, 0.85, 0.85, 0.85, 0.75)[is.na(ve_data$effectiveness)]
  ve_data$duration[is.na(ve_data$duration)] <- c(10, 5, 10, 10, 5)[is.na(ve_data$duration)]
  
  ve_all <- calculate_all_ve(ve_data, n_cycles = nC)
  
  # Extract VE vectors and adjust by coverage
  ve_MenABCWY_forACWY <- ve_all$Effectiveness[ve_all$Vaccine == "MenABCWY_for_SeroACWY"] * cov_ABCWY
  ve_MenABCWY_forB    <- ve_all$Effectiveness[ve_all$Vaccine == "MenABCWY_for_SeroB"]    * cov_ABCWY
  ve_MenACWY          <- ve_all$Effectiveness[ve_all$Vaccine == "MenACWY"]               * cov_ACWY
  ve_MenC             <- ve_all$Effectiveness[ve_all$Vaccine == "MenC"]                  * cov_C
  ve_MenB             <- ve_all$Effectiveness[ve_all$Vaccine == "MenB"]                  * cov_B
  
  ## ---- Costs (annual sequelae; infection cost vector) ----
  getc <- function(nm, default, df = costs_tbl) {
    if (is.null(df) || !all(c("Name", "Value") %in% names(df))) {
      return(default)
    }
    v <- df$Value[df$Name == nm]
    if (length(v) == 0) default else as.numeric(v[1])
  }
  
  c_Healthy        <- getc("c_Healthy", 0)
  c_Scarring       <- getc("c_Scarring", 13617)
  c_Single_Amput   <- getc("c_Single_Amput", 38187)
  c_Multiple_Amput <- getc("c_Multiple_Amput", 66000)
  c_Neuro_Disab    <- getc("c_Neuro_Disab", 15644)
  c_Hearing_Loss   <- getc("c_Hearing_Loss", 25093)
  c_Renal_Failure  <- getc("c_Renal_Failure", 18149)
  c_Seizure        <- getc("c_Seizure", 13454)
  c_Paralysis      <- getc("c_Paralysis", 10000)
  c_Dead           <- getc("c_Dead", 0)
  
  cost_imd <- data_list[["cost_IMD"]]
  if (!is.null(cost_imd) && "cost" %in% names(cost_imd)) {
    c_IMD_infection <- as.numeric(cost_imd$cost)
  } else {
    c_IMD_infection <- rep(0, nC)
  }
  
  if (length(c_IMD_infection) < nC) {
    c_IMD_infection <- rep(c_IMD_infection, length.out = nC)
  }
  
  ## ---- Utilities ----
  utab <- data_list[["utilities"]]
  getu <- function(nm, default) {
    if (is.null(utab) || !all(c("Name", "Value") %in% names(utab))) {
      return(default)
    }
    v <- utab$Value[utab$Name == nm]
    if (length(v) == 0) default else as.numeric(v[1])
  }
  
  u_Healthy          <- getu("u_Healthy", 1.00)
  u_IMD              <- getu("u_IMD",     0.91)
  u_Scarring         <- getu("u_Scarring", 0.95)
  u_Single_Amput     <- getu("u_Single_Amput", 0.70)
  u_Multiple_Amput   <- getu("u_Multiple_Amput", 0.61)
  u_Neuro_Disability <- getu("u_Neuro_Disability", 0.06)
  u_Hearing_Loss     <- getu("u_Hearing_Loss", 0.72)
  u_Renal_Failure    <- getu("u_Renal_Failure", 0.82)
  u_Seizure          <- getu("u_Seizure", 0.83)
  u_Paralysis        <- getu("u_Paralysis", 0.26)
  u_Dead             <- getu("u_Dead", 0.00)
  
  ## ---- Build parameter list ----
  params <- list(
    n_cycles = nC, cycle_length = cyc, d_c = dc, d_e = de,
    coverage_ABCWY = cov_ABCWY, coverage_ACWY = cov_ACWY, coverage_C = cov_C, coverage_B = cov_B,
    c_MenABCWY = c_MenABCWY, c_MenACWY = c_MenACWY, c_MenC = c_MenC, c_MenB = c_MenB, c_admin = c_admin,
    p_B = p_B, p_C = p_C, p_W = p_W, p_Y = p_Y,
    p_B_DeadIMD = p_B_CFR, p_C_DeadIMD = p_C_CFR, p_W_DeadIMD = p_W_CFR, p_Y_DeadIMD = p_Y_CFR,
    p_bg_mort = p_bg, p_seq_overall = p_seq_overall, w_seq = w_seq,
    ve_MenABCWY_forACWY = ve_MenABCWY_forACWY, ve_MenABCWY_forB = ve_MenABCWY_forB,
    ve_MenACWY = ve_MenACWY, ve_MenC = ve_MenC, ve_MenB = ve_MenB,
    c_IMD_infection = c_IMD_infection, c_Healthy = c_Healthy, c_Scarring = c_Scarring,
    c_Single_Amput = c_Single_Amput, c_Multiple_Amput = c_Multiple_Amput, c_Neuro_Disab = c_Neuro_Disab,
    c_Hearing_Loss = c_Hearing_Loss, c_Renal_Failure = c_Renal_Failure, c_Seizure = c_Seizure,
    c_Paralysis = c_Paralysis, c_Dead = c_Dead,
    u_Healthy = u_Healthy, u_IMD = u_IMD, u_Scarring = u_Scarring, u_Single_Amput = u_Single_Amput,
    u_Multiple_Amput = u_Multiple_Amput, u_Neuro_Disability = u_Neuro_Disability,
    u_Hearing_Loss = u_Hearing_Loss, u_Renal_Failure = u_Renal_Failure, u_Seizure = u_Seizure,
    u_Paralysis = u_Paralysis, u_Dead = u_Dead,
    perspective = perspective,
    mult_p_B = 1.0, mult_p_C = 1.0, mult_p_W = 1.0, mult_p_Y = 1.0,
    mult_cfr_B = 1.0, mult_cfr_C = 1.0, mult_cfr_W = 1.0, mult_cfr_Y = 1.0,
    mult_bg_mort = 1.0, mult_c_IMD = 1.0,
    ve_scale_B = 1.0, ve_scale_C = 1.0, ve_scale_W = 1.0, ve_scale_Y = 1.0
  )
  
  log_info("Base parameters loaded successfully")
  return(params)
}
