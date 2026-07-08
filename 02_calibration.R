# 02_calibration.R
# Requires: source("00_parameters.R"), source("00b_le_transitions.R"), source("01_model_functions.R")
# Calibrates fibrosis transition probabilities to ESSENCE placebo data.
# Outputs: p_prog_month (updated), candidate_sets, best, build_p_prog_from_annual,
#          p_prog_month_le_low, p_prog_month_le_high 

#############################################################
####################### CALIBRATION #########################
#############################################################

# ---- 1. ESSENCE trial parameters ----
N_trial_biopsy <- 222
obs_regress    <- 50
weeks_trial    <- 72
n_cyc_trial    <- round(weeks_trial / 52 / cycle_length)

age_essence     <- 56.0
cycle_offset_essence <- round((age_essence - age_start) / cycle_length)

v_init_essence <- c(F0 = 0, F1 = 0, F2 = 0.313, F3 = 0.687, F4_CC = 0,
                    DCC = 0, HCC = 0, LT_Y1 = 0, LT_Y1_P = 0,
                    Post_LT = 0, Dead = 0)

# ---- 2. HCC/DCC transitions (single source: nonfib_annual, 00_parameters.R) ----
hcc_dcc_shared <- c(
  F3_HCC  = nonfib_annual$F3_HCC,
  F4_DCC  = nonfib_annual$F4_DCC,
  F4_HCC  = nonfib_annual$F4_HCC,
  DCC_HCC = nonfib_annual$DCC_HCC
)

# ---- 3. Fibrosis transition sets (calculated in 001b) ----

candidate_sets <- lapply(
  list(obs_low    = le_obs_low,
       obs        = le_obs_mean,
       obs_high   = le_obs_high,
       trial_low  = le_trial_low,
       trial      = le_trial_mean,
       trial_high = le_trial_high),
  function(x) c(x, hcc_dcc_shared)
)

# ---- 4. Convert annual probabilities -> monthly probabilities ----
build_p_prog_from_annual <- function(vals_annual, p_prog_template) {
  pp <- p_prog_template  
  for (nm in names(vals_annual)) {
    pp[[nm]] <- rate_to_prob(prob_to_rate(vals_annual[[nm]]), cycle_length)
  }
  pp
}

# ---- 5. Run ESSENCE cohort from a single starting state ----
run_single_start <- function(p_prog_local, start_state, 
                             cycle_offset = cycle_offset_essence) {
  v0 <- setNames(rep(0, n_states), v_states)
  v0[start_state] <- 1

  aP <- build_a_P(
    rr_reg             = c(LSM = 1, Semaglutide = 1),
    rr_prog            = c(LSM = 1, Semaglutide = 1),
    p_prog_month_local = p_prog_local,
    treat_dur_cycles   = c(LSM = 0L, Semaglutide = 0L),
    treat_start_cycles = c(LSM = 1L, Semaglutide = 1L)
  )

  m_tr <- matrix(0, nrow = n_cyc_trial + 1, ncol = n_states,
                 dimnames = list(NULL, v_states))
  m_tr[1, ] <- v0

  for (t in 1:n_cyc_trial) {
    m_tr[t + 1, ] <- m_tr[t, ] %*% aP[,,cycle_offset + t,"LSM"]
  }
  m_tr[n_cyc_trial + 1, ]
}

# ---- 6. Predicted proportion with fibrosis regression at week 72 ----
predict_regression <- function(p_prog_local) {
  final_F2 <- run_single_start(p_prog_local, "F2")
  final_F3 <- run_single_start(p_prog_local, "F3")

  p_reg_F2 <- sum(final_F2[c("F0", "F1")])           
  p_reg_F3 <- sum(final_F3[c("F0", "F1", "F2")])    

  p_reg_total <- 0.313 * p_reg_F2 + 0.687 * p_reg_F3
  expected_n  <- p_reg_total * N_trial_biopsy

  ci_lo <- qbinom(0.025, N_trial_biopsy, p_reg_total)
  ci_hi <- qbinom(0.975, N_trial_biopsy, p_reg_total)

  list(p_reg_total = p_reg_total,
       expected_n  = expected_n,
       ci_lo       = ci_lo,
       ci_hi       = ci_hi,
       p_reg_F2    = p_reg_F2,
       p_reg_F3    = p_reg_F3)
}

# ---- 7. Calibration loop: fibrosis regression ----
results <- data.frame(
  set        = character(),
  p_reg_F2   = numeric(),
  p_reg_F3   = numeric(),
  p_reg_mix  = numeric(),
  expected_n = numeric(),
  ci_lo      = integer(),
  ci_hi      = integer(),
  hits_50    = logical(),
  diff_50    = numeric(),
  stringsAsFactors = FALSE
)

for (nm in names(candidate_sets)) {
  pp  <- build_p_prog_from_annual(candidate_sets[[nm]], p_prog_month)
  out <- predict_regression(pp)
  results <- rbind(results, data.frame(
    set        = nm,
    p_reg_F2   = round(out$p_reg_F2, 3),
    p_reg_F3   = round(out$p_reg_F3, 3),
    p_reg_mix  = round(out$p_reg_total, 3),
    expected_n = round(out$expected_n, 1),
    ci_lo      = out$ci_lo,
    ci_hi      = out$ci_hi,
    hits_50    = (50 >= out$ci_lo) & (50 <= out$ci_hi),
    diff_50    = round(out$expected_n - 50, 1)
  ))
}

# ---- 8. Pick the best-fitting parameter set ----
eligible <- results[results$set %in% c("obs", "trial"), ]
best <- eligible$set[which.min(abs(eligible$diff_50))]
cat("\nBest-fitting parameter set:", best, "\n")

# ---- 9. Plot: regression calibration ----
results$set <- factor(results$set,
                      levels = c("obs_low", "obs", "obs_high",
                                 "trial_low", "trial", "trial_high"))

ggplot(results, aes(x = set, y = expected_n)) +
  geom_point(size = 3, color = "#0072B2") +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.2,
                color = "#0072B2", linewidth = 0.8) +
  geom_hline(yintercept = 50, linetype = "dashed", linewidth = 0.8) +
  annotate("text", x = 6.3, y = 50, label = "Observed = 50",
           hjust = 0, size = 3.5) +
  labs(x     = "Transition Probabilities Source",
       y     = "Predicted Regressions (N = 222)",
       title = "ESSENCE Calibration: Biopsy-Confirmed Fibrosis Regression") +
  theme_bw(base_size = 13) +
  theme(plot.margin = margin(5, 60, 5, 5))

# ---- 10. Commit the best candidate set to p_prog_month ----
best_annual <- candidate_sets[[best]]
p_prog_month <- build_p_prog_from_annual(best_annual, p_prog_month)

cat("\nFibrosis transitions in p_prog_month now set to '", best,
    "' from the calibration.\n", sep = "")

#############################################################
########### OWSA bounds: Le et al. trial range ##############
#############################################################

p_prog_month_le_low  <- build_p_prog_from_annual(candidate_sets$trial_low,  p_prog_month)
p_prog_month_le_high <- build_p_prog_from_annual(candidate_sets$trial_high, p_prog_month)

cat("02_calibration.R complete.\n")
