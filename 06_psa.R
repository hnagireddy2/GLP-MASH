# 06_psa.R
# Requires: source("00_parameters.R"), source("00b_le_transitions.R"), source("02_calibration.R")
# PSA parameter generation, PSA iterator, and simulation loop
   
#############################################################
########## Define a generate_psa_params() function ##########
#############################################################

generate_psa_params <- function(n_sim) {

  ## Helpers
  gamma_params <- function(mean, se) {
    list(shape = (mean/se)^2, scale = se^2/mean)
  }
  beta_params <- function(mean, se) {
    # clamp se so alpha/beta stay positive
    se   <- min(se, sqrt(mean * (1 - mean)) * 0.99)
    alpha <- mean * (mean*(1-mean)/se^2 - 1)
    beta  <- (1-mean) * (mean*(1-mean)/se^2 - 1)
    list(shape1 = max(alpha, 0.01), shape2 = max(beta, 0.01))
  }
  ## SE implied by a reported 95% CI, assuming it was built as mean +/-
  ## 1.96*SE (the standard normal-approximation relationship). Verified
  ## exactly against this model's literature-sourced hazard CIs -- e.g.
  ## h_DCC_Death's (0.1216, 0.2784) bounds reconstruct only when SE = 0.04,
  ## which is exactly what this formula returns.
  se_from_ci95 <- function(lo, hi) (hi - lo) / (2 * 1.96)
  rgamma_ci <- function(n, mean, lo, hi) {
    p <- gamma_params(mean, se_from_ci95(lo, hi))
    rgamma(n, p$shape, scale = p$scale)
  }
  rbeta_se <- function(n, mean, se) {
    p <- beta_params(mean, se)
    rbeta(n, p$shape1, p$shape2)
  }

  ## Fibrosis-transition mean/low/high (annual), tied to whichever Le et al.
  ## candidate set calibration selected (02_calibration.R's `best`), so PSA's
  ## central estimate always matches the calibrated base case instead of
  ## drifting from it.
  fib_mean <- candidate_sets[[best]]
  fib_low  <- candidate_sets[[paste0(best, "_low")]]
  fib_high <- candidate_sets[[paste0(best, "_high")]]

  df <- data.frame(

    ## TREATMENT EFFECTS (log-normal, SEs from ESSENCE trial cell counts)
    rr_sema_regress  = rlnorm(n_sim,
                              meanlog = log(1.6637),
                              sdlog   = 0.1280),
    rr_sema_progress = rlnorm(n_sim,
                              meanlog = log(0.5785),
                              sdlog   = 0.2116),

    ## FIBROSIS TRANSITION HAZARDS (annual, gamma)
    h_F0_F1 = rgamma_ci(n_sim, fib_mean["F0_F1"], fib_low["F0_F1"], fib_high["F0_F1"]),
    h_F1_F0 = rgamma_ci(n_sim, fib_mean["F1_F0"], fib_low["F1_F0"], fib_high["F1_F0"]),
    h_F1_F2 = rgamma_ci(n_sim, fib_mean["F1_F2"], fib_low["F1_F2"], fib_high["F1_F2"]),
    h_F2_F1 = rgamma_ci(n_sim, fib_mean["F2_F1"], fib_low["F2_F1"], fib_high["F2_F1"]),
    h_F2_F3 = rgamma_ci(n_sim, fib_mean["F2_F3"], fib_low["F2_F3"], fib_high["F2_F3"]),
    h_F3_F2 = rgamma_ci(n_sim, fib_mean["F3_F2"], fib_low["F3_F2"], fib_high["F3_F2"]),
    h_F3_F4 = rgamma_ci(n_sim, fib_mean["F3_F4"], fib_low["F3_F4"], fib_high["F3_F4"]),
    h_F4_F3 = rgamma_ci(n_sim, fib_mean["F4_F3"], fib_low["F4_F3"], fib_high["F4_F3"]),

    ## Advanced-disease hazards (annual). Means = base case (nonfib_annual);
    ## ranges = 95% CI from Kim S1 SEs, or Kim's 20%-of-mean rule for the
    ## Rustgi-sourced LT transitions. se_from_ci95() recovers the exact SE.
    h_F3_HCC       = rgamma_ci(n_sim, 0.0034, 0.0021, 0.0047),
    h_F4_HCC       = rgamma_ci(n_sim, 0.0378, 0.0213, 0.0543),
    h_F4_DCC       = rgamma_ci(n_sim, 0.0659, 0.0400, 0.0918),
    h_DCC_HCC      = rgamma_ci(n_sim, 0.0378, 0.0213, 0.0543),
    h_DCC_LT       = rgamma_ci(n_sim, 0.0230, 0.0140, 0.0320),
    h_DCC_Death    = rgamma_ci(n_sim, 0.20,   0.1216, 0.2784),
    h_HCC_LT       = rgamma_ci(n_sim, 0.0300, 0.0182, 0.0418),
    h_HCC_Death    = rgamma_ci(n_sim, 0.1305, 0.1049, 0.1561),
    h_LT_Death     = rgamma_ci(n_sim, 0.0400, 0.0243, 0.0557),
    h_PostLT_Death = rgamma_ci(n_sim, 0.0820, 0.0499, 0.1141),

    ## STATE COSTS (gamma) -- means/ranges from costs_base/costs_low/costs_high
    ## (00_parameters.R) so PSA always tracks the base-case cost inputs.
    cost_F0_F2 = rgamma_ci(n_sim, costs_base["F0"], costs_low["F0"], costs_high["F0"]),

    ## STATE COSTS (gamma) — drawn then rank-ordered to enforce F3 < F4_CC < HCC < DCC
    cost_F3_raw    = rgamma_ci(n_sim, costs_base["F3"],    costs_low["F3"],    costs_high["F3"]),
    cost_F4_CC_raw = rgamma_ci(n_sim, costs_base["F4_CC"], costs_low["F4_CC"], costs_high["F4_CC"]),
    cost_HCC_raw   = rgamma_ci(n_sim, costs_base["HCC"],   costs_low["HCC"],   costs_high["HCC"]),
    cost_DCC_raw   = rgamma_ci(n_sim, costs_base["DCC"],   costs_low["DCC"],   costs_high["DCC"]),

    ## LT procedure cost (gamma)
    cost_LT = rgamma_ci(n_sim, costs_base["LT"], costs_low["LT"], costs_high["LT"]),

    ## HEALTH STATE UTILITIES (decrement, beta)
    # SE is capped inside beta_params() to keep alpha/beta positive
    qdec_F0_F2  = rbeta_se(n_sim, 0.016, 0.016 * 0.10),
    qdec_F3     = rbeta_se(n_sim, 0.145, 0.145 * 0.10),
    qdec_F4_CC  = rbeta_se(n_sim, 0.145, 0.145 * 0.10),
    qdec_DCC    = rbeta_se(n_sim, 0.155, 0.155 * 0.10),
    qdec_HCC    = rbeta_se(n_sim, 0.165, 0.165 * 0.10),
    qdec_PostLT = rbeta_se(n_sim, 0.036, 0.036 * 0.10)

  )
# Enforce rank order: F3 < F4_CC < HCC < DCC
  cost_ordered <- t(apply(cbind(df$cost_F3_raw, df$cost_F4_CC_raw,
                                df$cost_HCC_raw, df$cost_DCC_raw), 1, sort))
  df$cost_F3    <- cost_ordered[, 1]
  df$cost_F4_CC <- cost_ordered[, 2]
  df$cost_HCC   <- cost_ordered[, 3]
  df$cost_DCC   <- cost_ordered[, 4]
  df$cost_F3_raw <- df$cost_F4_CC_raw <- df$cost_HCC_raw <- df$cost_DCC_raw <- NULL

  return(df)
}

#############################################################
##  run_model_psa_iter_all: PSA iterator for 3 strategies  ##
##  LSM, Sema 72w (Age 12), Sema 72w (Age 18)              ##
##  Single cohort: F2/F3 adolescents at age 12             ##
#############################################################

run_model_psa_iter_all <- function(psa_row) {
  ## ---- Build all PSA-sampled parameters ----
  rr_reg_psa  <- rr_regress;  rr_reg_psa["Semaglutide"]  <- psa_row$rr_sema_regress
  rr_prog_psa <- rr_progress; rr_prog_psa["Semaglutide"] <- psa_row$rr_sema_progress
  
  ## All sampled hazards are ANNUAL probabilities; annual_to_month() does the
  ## prob -> rate -> cycle conversion, same as how the base case is built,
  ## so the PSA is centered on base. pmin() guards prob_to_rate at p<1.
  amh <- function(p) annual_to_month(pmin(p, 0.999))
  
  p_cycle_psa <- p_prog_month          # seed from base case, then overwrite sampled cells
  p_cycle_psa$F0_F1     <- amh(psa_row$h_F0_F1)
  p_cycle_psa$F1_F0     <- amh(psa_row$h_F1_F0)
  p_cycle_psa$F1_F2     <- amh(psa_row$h_F1_F2)
  p_cycle_psa$F2_F1     <- amh(psa_row$h_F2_F1)
  p_cycle_psa$F2_F3     <- amh(psa_row$h_F2_F3)
  p_cycle_psa$F3_F2     <- amh(psa_row$h_F3_F2)
  p_cycle_psa$F3_F4     <- amh(psa_row$h_F3_F4)
  p_cycle_psa$F4_F3     <- amh(psa_row$h_F4_F3)
  p_cycle_psa$F3_HCC    <- amh(psa_row$h_F3_HCC)
  p_cycle_psa$F4_HCC    <- amh(psa_row$h_F4_HCC)
  p_cycle_psa$F4_DCC    <- amh(psa_row$h_F4_DCC)
  p_cycle_psa$DCC_HCC   <- amh(psa_row$h_DCC_HCC)
  p_cycle_psa$DCC_LT    <- amh(psa_row$h_DCC_LT)
  p_cycle_psa$DCC_Death <- amh(psa_row$h_DCC_Death)
  p_cycle_psa$HCC_LT    <- amh(psa_row$h_HCC_LT)
  p_cycle_psa$HCC_Death <- amh(psa_row$h_HCC_Death)
  p_cycle_psa$LT_Death  <- amh(psa_row$h_LT_Death)
  p_cycle_psa$PostLT_Death <- amh(psa_row$h_PostLT_Death)
  
  cost_vec_psa          <- costs_base
  cost_vec_psa["F0"]    <- psa_row$cost_F0_F2
  cost_vec_psa["F1"]    <- psa_row$cost_F0_F2
  cost_vec_psa["F2"]    <- psa_row$cost_F0_F2
  cost_vec_psa["F3"]    <- psa_row$cost_F3
  cost_vec_psa["F4_CC"] <- psa_row$cost_F4_CC
  cost_vec_psa["DCC"]   <- psa_row$cost_DCC
  cost_vec_psa["HCC"]   <- psa_row$cost_HCC
  cost_vec_psa["LT"]    <- psa_row$cost_LT

  # Drug cost fixed at base case (no PSA uncertainty)
  drug_psa <- drug_cost

  qdec_psa             <- qaly_dec_base
  qdec_psa["F0"]       <- psa_row$qdec_F0_F2
  qdec_psa["F1"]       <- psa_row$qdec_F0_F2
  qdec_psa["F2"]       <- psa_row$qdec_F0_F2
  qdec_psa["F3"]       <- psa_row$qdec_F3
  qdec_psa["F4_CC"]    <- psa_row$qdec_F4_CC
  qdec_psa["DCC"]      <- psa_row$qdec_DCC
  qdec_psa["HCC"]      <- psa_row$qdec_HCC
  qdec_psa["LT"]       <- psa_row$qdec_F4_CC
  qdec_psa["Post_LT"]  <- psa_row$qdec_PostLT

  util_mat_psa <- build_util_matrix(v_util_age_base, qdec_psa)

  ## ---- LSM / Age12 / Age18, all three strategies at once ----
  ## (LSM always has treat_dur_cycles == 0, so its trace already comes
  ## out of the "Age12" run below — no need for a separate LSM-only build.)
  run_three_strategies(rr_reg_psa, rr_prog_psa,
                       p_prog_month_local   = p_cycle_psa,
                       util_matrix           = util_mat_psa,
                       cost_vector           = cost_vec_psa,
                       drug_cost_vec         = drug_psa,
                       treat_dur_cycles_vec  = treat_dur_72w_cycles)
}

#############################################################
##  Multi-Scenario PSA Loop: All Ages × All Durations     ##
##  NOTE: Each sim runs 8 model calls (2 ages × 4 durations)
##  n_sim_all is set to 1000                               ##
##  n_sim = 1000 from PSA above for VOI                   ##
#############################################################

n_sim_all <- 1000

# Uncomment to prevent overwriting of cache:
# if (!file.exists("psa_results_all.rds")) {

  df_psa_input_all <- generate_psa_params(n_sim_all)
  df_c_all <- as.data.frame(matrix(0, nrow = n_sim_all,
                                   ncol = length(all_strat_labels)))
  df_e_all <- as.data.frame(matrix(0, nrow = n_sim_all,
                                   ncol = length(all_strat_labels)))
  colnames(df_c_all) <- colnames(df_e_all) <- all_strat_labels

  t_start_all <- Sys.time()
  for (i in 1:n_sim_all) {
    res_i <- run_model_psa_iter_all(df_psa_input_all[i, ])
    for (lbl in all_strat_labels) {
      df_c_all[i, lbl] <- res_i[[lbl]]["Cost"]
      df_e_all[i, lbl] <- res_i[[lbl]]["QALY"]
    }
    if (i %% 50 == 0) cat(i, "/", n_sim_all, "\n")
  }
  cat("Multi-scenario PSA runtime:",
      round(difftime(Sys.time(), t_start_all, units = "mins"), 1), "min\n")

  saveRDS(list(df_c_all = df_c_all, df_e_all = df_e_all,
               df_psa_input_all = df_psa_input_all),
          "psa_results_all.rds")

  # Comment out the closing bracket and the else block:
  # } else {
  #   all_psa          <- readRDS("psa_results_all.rds")
  #   df_c_all         <- all_psa$df_c_all
  #   df_e_all         <- all_psa$df_e_all
  #   df_psa_input_all <- all_psa$df_psa_input_all
  #   cat("Loaded all-strategy PSA results from cache.\n")
  # }

test_res <- run_model_psa_iter_all(df_psa_input_all[1, ])
cat("Labels returned by function:\n")
print(names(test_res))
cat("\nLabels expected (all_strat_labels):\n")
print(all_strat_labels)
cat("\nMismatches:\n")
print(setdiff(all_strat_labels, names(test_res)))
