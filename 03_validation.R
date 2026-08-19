# 03_validation.R
# Requires: source("00_parameters.R"), source("00b_le_transitions.R"), source("02_calibration.R")
# F2 progression validation + Hagstrom time-to-severe-disease

#############################################################
############### VALIDATION: F2 PROGRESSION ##################
############################################################# 

N_F2_placebo    <- 81
obs_progress_F2 <- 29    

predict_progression_F2 <- function(p_prog_local) {
  final_F2 <- run_single_start(p_prog_local, "F2",
                                cycle_offset = cycle_offset_essence)
  p_prog   <- sum(final_F2[c("F3", "F4_CC", "DCC", "HCC",
                              "LT_Y1", "Post_LT")])

  expected_n <- p_prog * N_F2_placebo
  ci_lo <- qbinom(0.025, N_F2_placebo, p_prog)
  ci_hi <- qbinom(0.975, N_F2_placebo, p_prog)

  list(p_prog = p_prog, expected_n = expected_n,
       ci_lo = ci_lo, ci_hi = ci_hi)
}

pp_chosen   <- build_p_prog_from_annual(candidate_sets[[best]], p_prog_month)
val_F2_prog <- predict_progression_F2(pp_chosen)

cat("\n=== Validation: F2 Progression in ESSENCE Placebo ===\n")
cat("Source:    Table S6, Sanyal et al. NEJM 2025 supplement\n")
cat("Observed:  ", obs_progress_F2, "/", N_F2_placebo, "\n", sep = "")
cat("Predicted: ", round(val_F2_prog$expected_n, 1),
    " (95% CI: ", val_F2_prog$ci_lo, "-", val_F2_prog$ci_hi, ")\n", sep = "")


###############################################################
####### TIME-TO-SEVERE-DISEASE VALIDATION OF ALL STATES #######
###############################################################

target_hagstrom <- data.frame(
  start_stage  = c("F0", "F1", "F2", "F3", "F4_CC"),
  point_estimates = c(30.5, 35.6, 19.4, 6.0, 5.6),
  time_low     = c(21.5,  25.6,  9.3,  2.3,  0.9),
  time_high    = c(39.6,  45.4, 29.5,  9.6, 10.3),
  endpoint     = c("Severe Disease", "Severe Disease", "Severe Disease", "Severe Disease", "Decompensation")
)

age_hagstrom     <- 48.2
cycle_offset_hagstrom <- round((age_hagstrom - age_start) / cycle_length)

time_to_severe <- function(start_state, threshold = 0.10) {
  
  v_init_val <- setNames(rep(0, n_states), v_states)
  v_init_val[start_state] <- 1
  
  aP_val <- build_a_P(
    rr_reg             = c(LSM = 1, Semaglutide = 1),
    rr_prog            = c(LSM = 1, Semaglutide = 1),
    p_prog_month_local = p_prog_month,
    treat_dur_cycles   = c(LSM = 0L, Semaglutide = 0L),
    treat_start_cycles = c(LSM = 1L, Semaglutide = 1L)
  )
  
  max_cyc <- n_cycles - cycle_offset_hagstrom
  
  m_tr <- matrix(0, nrow = max_cyc + 1, ncol = n_states,
               dimnames = list(NULL, v_states))
  m_tr[1, ] <- v_init_val
  
  # Make severe states absorbing
  for (s in c("DCC", "HCC", "LT_Y1", "Post_LT")) {
    for (cyc in 1:dim(aP_val)[3]) {
      aP_val[s, , cyc, "LSM"] <- 0
      aP_val[s, s, cyc, "LSM"] <- 1
    }
  }

  # Run the Markov trace
  for (t in 1:max_cyc) {
    m_tr[t + 1, ] <- m_tr[t, ] %*% aP_val[,,cycle_offset_hagstrom + t,"LSM"]
  }

  severe_states <- c("DCC", "HCC", "LT_Y1", "Post_LT")
  pct_severe <- rowSums(m_tr[1:(max_cyc + 1), severe_states])

  cycle_at_threshold <- which(pct_severe >= threshold)[1]
  
  if (is.na(cycle_at_threshold)) {
    return(NA_real_)
  }
  
  (cycle_at_threshold - 1) * cycle_length
}

# ---- Run Hagstrom validation for each starting stage ----
val2_results <- data.frame(
  Stage      = c("F0", "F1", "F2", "F3", "F4_CC"),
  Years_to_10pct_severe = sapply(c("F0", "F1", "F2", "F3", "F4_CC"),
                                  time_to_severe, threshold = 0.10)
)

val2_results <- merge(val2_results, target_hagstrom,
                       by.x = "Stage", by.y = "start_stage", all.x = TRUE)

val2_results$Pass <- with(val2_results,
  ifelse(is.na(Years_to_10pct_severe), "Not reached",
    ifelse(Years_to_10pct_severe >= time_low &
           Years_to_10pct_severe <= time_high, "PASS", "FAIL"))
)

cat("\n==========================================================\n")
cat(" VALIDATION 2: Hagstrom 2017 Time-to-Severe-Disease\n")
cat("==========================================================\n")
cat(" Severe disease = DCC, HCC, LT, or post-LT (liver outcomes)\n")
cat(" Threshold      = 10% of cohort reaches severe disease\n\n")
print(val2_results, row.names = FALSE)
cat("==========================================================\n\n")


# Plot: Time-to-severe trajectories for each starting stage
plot_severe_trajectory <- function() {
  
  starting_stages <- c("F0", "F1", "F2", "F3", "F4_CC")
  
  trajectories <- lapply(starting_stages, function(stg) {
    v_init_val <- setNames(rep(0, n_states), v_states)
    v_init_val[stg] <- 1
    
    aP_val <- build_a_P(
      rr_reg             = c(LSM = 1, Semaglutide = 1),
      rr_prog            = c(LSM = 1, Semaglutide = 1),
      p_prog_month_local = p_prog_month,
      treat_dur_cycles   = c(LSM = 0L, Semaglutide = 0L),
      treat_start_cycles = c(LSM = 1L, Semaglutide = 1L)
    )
    
    max_cyc <- n_cycles - cycle_offset_hagstrom
    m_tr <- matrix(0, nrow = max_cyc + 1, ncol = n_states,
                   dimnames = list(NULL, v_states))
    m_tr[1, ] <- v_init_val
    
    # Make severe states absorbing
    for (s in c("DCC", "HCC", "LT_Y1", "Post_LT")) {
      for (cyc in 1:dim(aP_val)[3]) {
        aP_val[s, , cyc, "LSM"] <- 0
        aP_val[s, s, cyc, "LSM"] <- 1
      }
    }

    # Run the Markov trace
    for (t in 1:max_cyc) {
      m_tr[t + 1, ] <- m_tr[t, ] %*% aP_val[,,cycle_offset_hagstrom + t,"LSM"]
    }

    severe_states <- c("DCC", "HCC", "LT_Y1", "Post_LT")
    pct_severe <- rowSums(m_tr[1:(max_cyc + 1), severe_states])
    
    data.frame(
      Stage      = stg,
      Year       = (0:max_cyc) * cycle_length,
      PctSevere  = pct_severe
    )
  })
  
  df_traj <- do.call(rbind, trajectories) %>%
    filter(Year <= 40)
  
  ggplot(df_traj, aes(x = Year, y = PctSevere, color = Stage)) +
    geom_line(linewidth = 1) +
    geom_hline(yintercept = 0.10, linetype = "dashed", color = "gray40") +
    annotate("text", x = 35, y = 0.105, label = "10% threshold",
             hjust = 0, size = 3.5, color = "gray40") +
    scale_x_continuous("Years from baseline", breaks = seq(0, 40, by = 5)) +
    scale_y_continuous("Cumulative proportion reaching severe disease",
                       labels = percent_format(), limits = c(0, 0.5)) +
    scale_color_manual(
      values = c(F0 = "#56B4E9", F1 = "#0072B2", F2 = "#009E73",
                 F3 = "#E69F00", F4_CC = "#D55E00")
    ) +
    labs(title = "Time-to-Severe-Disease by Baseline Fibrosis Stage",
         subtitle = paste0("Cohort starting at age ", age_hagstrom, ", no treatment"),
         color = "Baseline Stage") +
    theme_bw(base_size = 13) +
    theme(legend.position = "right")
}

plot_severe_trajectory()

cat("03_validation.R complete.\n")
