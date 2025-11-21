################################################################################
# 03_model_diagnostics.R
#
# PURPOSE: Deep dive into model mechanics and intermediate outputs
# - Trace analysis (state occupancy over time)
# - Flow analysis (transitions between states)
# - Epidemiological outcomes by serogroup
# - Vaccine effectiveness application
# - Cost and QALY accumulation patterns
#
# WHEN TO USE:
#   - To understand HOW the model works
#   - To debug unexpected results
#   - To verify model logic
#   - For presentations/publications showing model dynamics
#
# USAGE:
#   source("scripts/03_model_diagnostics.R")
#
# OUTPUT:
#   - Diagnostic plots in outputs/figures/[perspective]/diagnostics/
#   - Diagnostic tables in outputs/tables/diagnostics/
#
################################################################################

# ============================================================
# SETUP
# ============================================================

start_time <- Sys.time()

cat("\n")
cat(strrep("=", 80), "\n")
cat("MODEL DIAGNOSTICS - DETAILED MECHANICS\n")
cat(strrep("=", 80), "\n\n")

# Load packages
if (!require("pacman", quietly = TRUE)) {
  install.packages("pacman", repos = "https://cloud.r-project.org")
}
library(pacman)
p_load(here, readxl, ggplot2, dplyr, tidyr, purrr, tibble, viridis, reshape2)

# Source all model files
cat("[INFO] Loading model files...\n")
files <- list.files("R", pattern = "^[0-9]{2}.*\\.R$", full.names = TRUE)
files <- files[order(files)]
invisible(lapply(files, source))
cat("[INFO] Model files loaded\n\n")

# ============================================================
# CONFIGURATION
# ============================================================

PERSPECTIVE <- "healthcare"
STRATEGIES_TO_ANALYZE <- c("MenC", "MenACWY", "MenABCWY")  # Select strategies

# Update folders
perspective <- PERSPECTIVE
update_output_folders(PERSPECTIVE)

# Create diagnostics folders
diag_dir <- file.path(OUT_FIG, PERSPECTIVE, "diagnostics")
dir.create(diag_dir, recursive = TRUE, showWarnings = FALSE)

diag_table_dir <- file.path(OUT_TAB, "diagnostics")
dir.create(diag_table_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# LOAD PARAMETERS AND RUN MODEL
# ============================================================

cat(strrep("-", 80), "\n")
cat("RUNNING MODEL FOR DIAGNOSTICS\n")
cat(strrep("-", 80), "\n\n")

params_bc <- get_base_params()
all_results <- eval_all_strategies(params_bc)

cat(sprintf("✓ Model executed for %d strategies\n\n", length(all_results)))

# ============================================================
# 1. STATE OCCUPANCY ANALYSIS
# ============================================================

cat(strrep("-", 80), "\n")
cat("1. STATE OCCUPANCY OVER TIME\n")
cat(strrep("-", 80), "\n\n")

# Extract trace matrices for selected strategies
for (strat in STRATEGIES_TO_ANALYZE) {
  
  if (!strat %in% names(all_results)) next
  
  cat(sprintf("[%s] Analyzing state occupancy...\n", strat))
  
  trace <- all_results[[strat]]$trace
  n_cycles <- nrow(trace) - 1
  
  # Convert to long format
  trace_long <- as.data.frame(trace) %>%
    mutate(Cycle = 0:n_cycles) %>%
    pivot_longer(cols = -Cycle, names_to = "State", values_to = "Proportion")
  
  # Plot: All states
  p1 <- ggplot(trace_long, aes(x = Cycle, y = Proportion, color = State)) +
    geom_line(size = 0.8) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    labs(
      title = sprintf("State Occupancy Over Time - %s", strat),
      subtitle = sprintf("Cohort of %s individuals", format(n_cohort, big.mark = ",")),
      x = "Cycle (Years)",
      y = "Proportion of Cohort"
    ) +
    theme_bw() +
    theme(legend.position = "right")
  
  ggsave(file.path(diag_dir, paste0("trace_all_states_", strat, ".png")),
         p1, width = 12, height = 7, dpi = 300)
  
  # Plot: Focus on disease states (exclude Healthy and deaths)
  disease_states <- trace_long %>%
    filter(!State %in% c("Healthy", "Dead_IMD", "Background_Mortality"))
  
  if (nrow(disease_states) > 0) {
    p2 <- ggplot(disease_states, aes(x = Cycle, y = Proportion, color = State)) +
      geom_line(size = 1) +
      scale_y_continuous(labels = scales::percent) +
      labs(
        title = sprintf("Disease State Occupancy - %s", strat),
        subtitle = "Excluding Healthy and Death states",
        x = "Cycle (Years)",
        y = "Proportion of Cohort"
      ) +
      theme_bw() +
      theme(legend.position = "right")
    
    ggsave(file.path(diag_dir, paste0("trace_disease_states_", strat, ".png")),
           p2, width = 12, height = 7, dpi = 300)
  }
  
  # Save state occupancy table (first 20 cycles)
  trace_summary <- as.data.frame(trace[1:min(21, nrow(trace)), ]) %>%
    mutate(Cycle = 0:min(20, n_cycles)) %>%
    select(Cycle, everything())
  
  write.csv(trace_summary, 
            file.path(diag_table_dir, paste0("state_occupancy_", strat, ".csv")),
            row.names = FALSE)
  
  cat(sprintf("  ✓ Plots and table saved for %s\n", strat))
}

cat("\n")

# ============================================================
# 2. FLOW ANALYSIS (TRANSITIONS)
# ============================================================

cat(strrep("-", 80), "\n")
cat("2. FLOW ANALYSIS - STATE TRANSITIONS\n")
cat(strrep("-", 80), "\n\n")

for (strat in STRATEGIES_TO_ANALYZE) {
  
  if (!strat %in% names(all_results)) next
  
  cat(sprintf("[%s] Analyzing flows...\n", strat))
  
  flow <- all_results[[strat]]$flow
  
  # Calculate cumulative flows from Healthy to infection states
  infection_flows <- data.frame(
    Cycle = 1:params_bc$n_cycles,
    To_SeroB = flow["Healthy", "SeroB_Infect", ],
    To_SeroC = flow["Healthy", "SeroC_Infect", ],
    To_SeroW = flow["Healthy", "SeroW_Infect", ],
    To_SeroY = flow["Healthy", "SeroY_Infect", ]
  ) %>%
    mutate(
      Total_Infections = To_SeroB + To_SeroC + To_SeroW + To_SeroY,
      Cum_SeroB = cumsum(To_SeroB) * n_cohort,
      Cum_SeroC = cumsum(To_SeroC) * n_cohort,
      Cum_SeroW = cumsum(To_SeroW) * n_cohort,
      Cum_SeroY = cumsum(To_SeroY) * n_cohort
    )
  
  # Plot: Cumulative infections by serogroup
  infection_long <- infection_flows %>%
    select(Cycle, starts_with("Cum_")) %>%
    pivot_longer(cols = starts_with("Cum_"), 
                 names_to = "Serogroup", 
                 values_to = "Cumulative_Cases") %>%
    mutate(Serogroup = gsub("Cum_", "", Serogroup))
  
  p3 <- ggplot(infection_long, aes(x = Cycle, y = Cumulative_Cases, color = Serogroup)) +
    geom_line(size = 1.2) +
    scale_y_continuous(labels = scales::comma) +
    labs(
      title = sprintf("Cumulative Infections by Serogroup - %s", strat),
      subtitle = sprintf("Cohort of %s individuals", format(n_cohort, big.mark = ",")),
      x = "Cycle (Years)",
      y = "Cumulative Infections"
    ) +
    theme_bw() +
    theme(legend.position = "right")
  
  ggsave(file.path(diag_dir, paste0("cumulative_infections_", strat, ".png")),
         p3, width = 10, height = 6, dpi = 300)
  
  # Plot: Flow FROM infections TO sequelae
  sequela_states <- c("Scarring", "Single_Amput", "Multiple_Amput", 
                      "Neuro_Disability", "Hearing_Loss", "Renal_Failure", 
                      "Seizure", "Paralysis")
  
  sequela_flows <- data.frame(
    Cycle = 1:params_bc$n_cycles
  )
  
  for (seq in sequela_states) {
    if (seq %in% dimnames(flow)[[2]]) {
      # Sum flows from all infection states to this sequela
      sequela_flows[[seq]] <- (
        flow["SeroB_Infect", seq, ] +
          flow["SeroC_Infect", seq, ] +
          flow["SeroW_Infect", seq, ] +
          flow["SeroY_Infect", seq, ]
      ) * n_cohort
    }
  }
  
  if (ncol(sequela_flows) > 1) {
    sequela_long <- sequela_flows %>%
      pivot_longer(cols = -Cycle, names_to = "Sequela", values_to = "New_Cases")
    
    p4 <- ggplot(sequela_long, aes(x = Cycle, y = New_Cases, color = Sequela)) +
      geom_line(size = 0.8) +
      scale_y_continuous(labels = scales::comma) +
      labs(
        title = sprintf("Incident Sequelae Cases per Cycle - %s", strat),
        x = "Cycle (Years)",
        y = "New Sequelae Cases"
      ) +
      theme_bw() +
      theme(legend.position = "right")
    
    ggsave(file.path(diag_dir, paste0("incident_sequelae_", strat, ".png")),
           p4, width = 12, height = 7, dpi = 300)
  }
  
  # Save flow summary table
  write.csv(infection_flows,
            file.path(diag_table_dir, paste0("infection_flows_", strat, ".csv")),
            row.names = FALSE)
  
  cat(sprintf("  ✓ Flow plots and table saved for %s\n", strat))
}

cat("\n")

# ============================================================
# 3. EPIDEMIOLOGICAL OUTCOMES COMPARISON
# ============================================================

cat(strrep("-", 80), "\n")
cat("3. EPIDEMIOLOGICAL OUTCOMES COMPARISON\n")
cat(strrep("-", 80), "\n\n")

# Extract epidemiological outcomes
epi_summary <- data.frame(
  Strategy = names(all_results),
  Total_Infections = sapply(all_results, function(r) sum(r$epi$infections)),
  Infections_B = sapply(all_results, function(r) r$epi$infections["B"]),
  Infections_C = sapply(all_results, function(r) r$epi$infections["C"]),
  Infections_W = sapply(all_results, function(r) r$epi$infections["W"]),
  Infections_Y = sapply(all_results, function(r) r$epi$infections["Y"]),
  Deaths_IMD = sapply(all_results, function(r) r$epi$deaths_imd),
  Total_Sequelae = sapply(all_results, function(r) sum(r$epi$sequelae))
)

# Calculate infections prevented vs MenC
menc_infections <- epi_summary$Total_Infections[epi_summary$Strategy == "MenC"]
epi_summary$Infections_Prevented <- menc_infections - epi_summary$Total_Infections

print(epi_summary)
cat("\n")

# Plot: Infections by serogroup
epi_long <- epi_summary %>%
  select(Strategy, Infections_B, Infections_C, Infections_W, Infections_Y) %>%
  pivot_longer(cols = starts_with("Infections_"), 
               names_to = "Serogroup", 
               values_to = "Cases") %>%
  mutate(Serogroup = gsub("Infections_", "", Serogroup))

p5 <- ggplot(epi_long, aes(x = Strategy, y = Cases, fill = Serogroup)) +
  geom_col(position = "stack") +
  scale_y_continuous(labels = scales::comma) +
  labs(
    title = "Total Infections by Serogroup and Strategy",
    subtitle = sprintf("Cohort of %s individuals over lifetime", 
                       format(n_cohort, big.mark = ",")),
    x = "Strategy",
    y = "Total Infections"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(diag_dir, "infections_by_serogroup.png"),
       p5, width = 10, height = 6, dpi = 300)

# Save epidemiological summary
write.csv(epi_summary,
          file.path(diag_table_dir, "epidemiological_summary.csv"),
          row.names = FALSE)

cat("✓ Epidemiological comparison complete\n\n")

# ============================================================
# 4. VACCINE EFFECTIVENESS APPLICATION
# ============================================================

cat(strrep("-", 80), "\n")
cat("4. VACCINE EFFECTIVENESS OVER TIME\n")
cat(strrep("-", 80), "\n\n")

# Calculate VE for each strategy over time
ve_data <- list()

for (strat in STRATEGIES_TO_ANALYZE) {
  ve_by_cycle <- data.frame(
    Cycle = 1:params_bc$n_cycles,
    Strategy = strat
  )
  
  # Get VE for each serogroup at each cycle
  for (t in 1:params_bc$n_cycles) {
    ve <- get_serogroup_ve(params_bc, strat, t)
    ve_by_cycle$VE_B[t] <- ve["B"]
    ve_by_cycle$VE_C[t] <- ve["C"]
    ve_by_cycle$VE_W[t] <- ve["W"]
    ve_by_cycle$VE_Y[t] <- ve["Y"]
  }
  
  ve_data[[strat]] <- ve_by_cycle
}

ve_combined <- bind_rows(ve_data) %>%
  pivot_longer(cols = starts_with("VE_"), 
               names_to = "Serogroup", 
               values_to = "Effectiveness") %>%
  mutate(Serogroup = gsub("VE_", "", Serogroup))

# Plot VE over time
p6 <- ggplot(ve_combined, aes(x = Cycle, y = Effectiveness, 
                              color = Serogroup, linetype = Strategy)) +
  geom_line(size = 1) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(
    title = "Vaccine Effectiveness by Serogroup Over Time",
    subtitle = "Shows waning immunity and booster effects",
    x = "Cycle (Years from vaccination)",
    y = "Vaccine Effectiveness"
  ) +
  theme_bw() +
  theme(legend.position = "right")

ggsave(file.path(diag_dir, "vaccine_effectiveness_time.png"),
       p6, width = 12, height = 7, dpi = 300)

cat("✓ Vaccine effectiveness analysis complete\n\n")

# ============================================================
# 5. COST AND QALY ACCUMULATION
# ============================================================

cat(strrep("-", 80), "\n")
cat("5. COST AND QALY ACCUMULATION PATTERNS\n")
cat(strrep("-", 80), "\n\n")

# This requires re-running model with cycle-by-cycle tracking
# We'll approximate using state occupancy and per-cycle costs/QALYs

for (strat in STRATEGIES_TO_ANALYZE[1:2]) {  # Just show 2 for brevity
  
  if (!strat %in% names(all_results)) next
  
  trace <- all_results[[strat]]$trace
  
  # Get per-cycle costs and QALYs
  v_util <- get_state_utils(params_bc)
  v_wcc <- darthtools::gen_wcc(params_bc$n_cycles, method = "Simpson1/3")
  
  cycle_data <- data.frame(
    Cycle = 1:params_bc$n_cycles
  )
  
  for (t in 1:params_bc$n_cycles) {
    c_state <- get_state_costs(params_bc, t)
    
    # Undiscounted per-cycle values
    cycle_data$Cost_Undiscounted[t] <- sum(trace[t, ] * c_state) * n_cohort
    cycle_data$QALY_Undiscounted[t] <- sum(trace[t, ] * v_util) * n_cohort
    
    # Discounted
    cycle_data$Cost_Discounted[t] <- cycle_data$Cost_Undiscounted[t] * v_discount_cost[t] * v_wcc[t]
    cycle_data$QALY_Discounted[t] <- cycle_data$QALY_Undiscounted[t] * v_discount_qaly[t] * v_wcc[t]
  }
  
  # Add cumulative
  cycle_data$Cum_Cost <- cumsum(cycle_data$Cost_Discounted)
  cycle_data$Cum_QALY <- cumsum(cycle_data$QALY_Discounted)
  
  # Plot cumulative costs and QALYs
  p7 <- ggplot(cycle_data, aes(x = Cycle)) +
    geom_line(aes(y = Cum_Cost / 1e6), color = "red", size = 1) +
    scale_y_continuous(
      name = "Cumulative Discounted Cost (Millions CAD)",
      labels = scales::dollar,
      sec.axis = sec_axis(~.*1e6/max(cycle_data$Cum_Cost)*max(cycle_data$Cum_QALY), 
                          name = "Cumulative Discounted QALYs",
                          labels = scales::comma)
    ) +
    geom_line(aes(y = Cum_QALY * max(Cum_Cost) / max(Cum_QALY) / 1e6), 
              color = "blue", size = 1) +
    labs(
      title = sprintf("Cost and QALY Accumulation - %s", strat),
      subtitle = "Red = Cost | Blue = QALYs (discounted)",
      x = "Cycle (Years)"
    ) +
    theme_bw()
  
  ggsave(file.path(diag_dir, paste0("cost_qaly_accumulation_", strat, ".png")),
         p7, width = 10, height = 6, dpi = 300)
}

cat("✓ Accumulation analysis complete\n\n")

# ============================================================
# SUMMARY
# ============================================================

end_time <- Sys.time()
elapsed <- difftime(end_time, start_time, units = "secs")

cat("\n")
cat(strrep("=", 80), "\n")
cat("DIAGNOSTICS COMPLETE\n")
cat(strrep("=", 80), "\n\n")

cat("Outputs created:\n")
cat(sprintf("  - Diagnostic plots: %s\n", diag_dir))
cat(sprintf("  - Diagnostic tables: %s\n", diag_table_dir))
cat(sprintf("\nTotal time: %.1f seconds\n", elapsed))
cat(strrep("=", 80), "\n\n")
