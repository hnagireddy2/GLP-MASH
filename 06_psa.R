# 06_psa.R
# Requires: 00_parameters.R, 01_model_functions.R, 02_calibration.R
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
  se_from_range <- function(lo, hi) (hi - lo) / 4

  df <- data.frame(

    ## TREATMENT EFFECTS (log-normal, SEs from ESSENCE trial cell counts)
    rr_sema_regress  = rlnorm(n_sim,
                              meanlog = log(1.6637),
                              sdlog   = 0.1280),
    rr_sema_progress = rlnorm(n_sim,
                              meanlog = log(0.5785),
                              sdlog   = 0.2116),

    ## FIBROSIS TRANSITION HAZARDS (annual, gamma) 
    h_F0_F1 = rgamma(n_sim,
                     gamma_params(0.100, se_from_range(0.077, 0.132))$shape,
                     scale = gamma_params(0.100, se_from_range(0.077, 0.132))$scale),
    h_F1_F0 = rgamma(n_sim,
                     gamma_params(0.025, se_from_range(0.014, 0.046))$shape,
                     scale = gamma_params(0.025, se_from_range(0.014, 0.046))$scale),
    h_F1_F2 = rgamma(n_sim,
                     gamma_params(0.097, se_from_range(0.068, 0.137))$shape,
                     scale = gamma_params(0.097, se_from_range(0.068, 0.137))$scale),
    h_F2_F1 = rgamma(n_sim,
                     gamma_params(0.049, se_from_range(0.023, 0.105))$shape,
                     scale = gamma_params(0.049, se_from_range(0.023, 0.105))$scale),
    h_F2_F3 = rgamma(n_sim,
                     gamma_params(0.075, se_from_range(0.061, 0.092))$shape,
                     scale = gamma_params(0.075, se_from_range(0.061, 0.092))$scale),
    h_F3_F2 = rgamma(n_sim,
                     gamma_params(0.080, se_from_range(0.046, 0.138))$shape,
                     scale = gamma_params(0.080, se_from_range(0.046, 0.138))$scale),
    h_F3_F4 = rgamma(n_sim,
                     gamma_params(0.045, se_from_range(0.033, 0.060))$shape,
                     scale = gamma_params(0.045, se_from_range(0.033, 0.060))$scale),
    h_F4_F3 = rgamma(n_sim,
                     gamma_params(0.047, se_from_range(0.021, 0.108))$shape,
                     scale = gamma_params(0.047, se_from_range(0.021, 0.108))$scale),

    ## Advanced state hazards (annual) 
    h_F4_DCC    = rgamma(n_sim,
                         gamma_params(0.0659, se_from_range(0.03, 0.12))$shape,
                         scale = gamma_params(0.0659, se_from_range(0.03, 0.12))$scale),
    h_DCC_Death = rgamma(n_sim,
                         gamma_params(0.32, se_from_range(0.15, 0.40))$shape,
                         scale = gamma_params(0.32, se_from_range(0.15, 0.40))$scale),
    h_HCC_Death = rgamma(n_sim,
                         gamma_params(0.19, se_from_range(0.10, 0.38))$shape,
                         scale = gamma_params(0.19, se_from_range(0.10, 0.38))$scale),
    h_LT_Death  = rgamma(n_sim,
                         gamma_params(0.067, se_from_range(0.065, 0.10))$shape,
                         scale = gamma_params(0.067, se_from_range(0.065, 0.10))$scale),
    h_PostLT_Death = rgamma(n_sim,
                            gamma_params(0.036, se_from_range(0.03, 0.05))$shape,
                            scale = gamma_params(0.036, se_from_range(0.03, 0.05))$scale),

    ## STATE COSTS (gamma) 
    cost_F0_F2 = rgamma(n_sim,
                    gamma_params(8698, se_from_range(6958, 10436))$shape,
                    scale = gamma_params(8698, se_from_range(6958, 10436))$scale),
    
    ## STATE COSTS (gamma) — drawn then rank-ordered to enforce F3 < F4_CC < HCC < DCC
    cost_F3_raw    = rgamma(n_sim,
                       gamma_params(10372,  se_from_range(8297,   12447))$shape,
                       scale = gamma_params(10372,  se_from_range(8297,   12447))$scale),
    cost_F4_CC_raw = rgamma(n_sim,
                       gamma_params(42207,  se_from_range(33766,  50650))$shape,
                       scale = gamma_params(42207,  se_from_range(33766,  50650))$scale),
    cost_HCC_raw   = rgamma(n_sim,
                       gamma_params(141615, se_from_range(113292, 169939))$shape,
                       scale = gamma_params(141615, se_from_range(113292, 169939))$scale),
    cost_DCC_raw   = rgamma(n_sim,
                       gamma_params(195156, se_from_range(156125, 234187))$shape,
                       scale = gamma_params(195156, se_from_range(156125, 234187))$scale),
   
    ## LT complication additional cost (gamma)
    cost_lt_add = rgamma(n_sim,
                         gamma_params(lt_add_cost_base,
                                      se_from_range(lt_add_cost_low,
                                                    lt_add_cost_high))$shape,
                         scale = gamma_params(lt_add_cost_base,
                                              se_from_range(lt_add_cost_low,
                                                          lt_add_cost_high))$scale),

    ## HEALTH STATE UTILITIES (decrement, beta) 
    # SE is capped inside beta_params() to keep alpha/beta positive
    qdec_F0_F2  = rbeta(n_sim,
                        beta_params(0.016, 0.016 * 0.10)$shape1,
                        beta_params(0.016, 0.016 * 0.10)$shape2),
    qdec_F3     = rbeta(n_sim,
                        beta_params(0.145, 0.145 * 0.10)$shape1,
                        beta_params(0.145, 0.145 * 0.10)$shape2),
    qdec_F4_CC  = rbeta(n_sim,
                        beta_params(0.145, 0.145 * 0.10)$shape1,
                        beta_params(0.145, 0.145 * 0.10)$shape2),
    qdec_DCC    = rbeta(n_sim,
                        beta_params(0.155, 0.155 * 0.10)$shape1,
                        beta_params(0.155, 0.155 * 0.10)$shape2),
    qdec_HCC    = rbeta(n_sim,
                        beta_params(0.165, 0.165 * 0.10)$shape1,
                        beta_params(0.165, 0.165 * 0.10)$shape2),
    qdec_PostLT = rbeta(n_sim,
                        beta_params(0.036, 0.036 * 0.10)$shape1,
                        beta_params(0.036, 0.036 * 0.10)$shape2)

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

  p_cycle_psa <- p_prog_month
  p_cycle_psa$F0_F1        <- rate_to_prob(psa_row$h_F0_F1,        cycle_length)
  p_cycle_psa$F1_F0        <- rate_to_prob(psa_row$h_F1_F0,        cycle_length)
  p_cycle_psa$F1_F2        <- rate_to_prob(psa_row$h_F1_F2,        cycle_length)
  p_cycle_psa$F2_F1        <- rate_to_prob(psa_row$h_F2_F1,        cycle_length)
  p_cycle_psa$F2_F3        <- rate_to_prob(psa_row$h_F2_F3,        cycle_length)
  p_cycle_psa$F3_F2        <- rate_to_prob(psa_row$h_F3_F2,        cycle_length)
  p_cycle_psa$F3_F4        <- rate_to_prob(psa_row$h_F3_F4,        cycle_length)
  p_cycle_psa$F4_F3        <- rate_to_prob(psa_row$h_F4_F3,        cycle_length)
  p_cycle_psa$F4_DCC       <- rate_to_prob(psa_row$h_F4_DCC,       cycle_length)
  p_cycle_psa$DCC_Death    <- rate_to_prob(psa_row$h_DCC_Death,    cycle_length)
  p_cycle_psa$HCC_Death    <- rate_to_prob(psa_row$h_HCC_Death,    cycle_length)
  p_cycle_psa$LT_Death     <- rate_to_prob(psa_row$h_LT_Death,     cycle_length)
  p_cycle_psa$PostLT_Death <- rate_to_prob(psa_row$h_PostLT_Death, cycle_length)

  cost_vec_psa          <- costs_base
  cost_vec_psa["F0"]    <- psa_row$cost_F0_F2
  cost_vec_psa["F1"]    <- psa_row$cost_F0_F2
  cost_vec_psa["F2"]    <- psa_row$cost_F0_F2
  cost_vec_psa["F3"]    <- psa_row$cost_F3
  cost_vec_psa["F4_CC"] <- psa_row$cost_F4_CC
  cost_vec_psa["DCC"]   <- psa_row$cost_DCC
  cost_vec_psa["HCC"]   <- psa_row$cost_HCC
  for (st in c("LT_Y1", "LT_Y1_P"))
    cost_vec_psa[st] <- 452682 + psa_row$cost_lt_add

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
  qdec_psa["LT_Y1"]    <- psa_row$qdec_F4_CC
  qdec_psa["LT_Y1_P"]  <- psa_row$qdec_F4_CC
  qdec_psa["Post_LT"]  <- psa_row$qdec_PostLT

  util_mat_psa <- build_util_matrix(v_util_age_base, qdec_psa)

  ## ---- Strategy 1: LSM (no treatment, natural history) ----
  aP_lsm <- build_a_P(rr_reg_psa, rr_prog_psa,
                      p_prog_month_local = p_cycle_psa,
                      treat_dur_cycles   = c(LSM = 0L, Semaglutide = 0L),
                      treat_start_cycles = treat_start_immediate)

  traces_lsm <- lapply(treatments, function(stg)
    run_markov(aP_lsm[,,,stg], v_init))
  names(traces_lsm) <- treatments

  res_lsm <- summarize_strategies(
    traces_lsm, "LSM",
    util_mat_psa, v_background_cost_cycle,
    cost_vec_psa, drug_psa,
    treat_dur_cycles_vec = c(LSM = 0L, Semaglutide = 0L)
  )

  ## ---- Strategy 2: Treat at 12 (72 weeks starting immediately) ----
  aP_12 <- build_a_P(rr_reg_psa, rr_prog_psa,
                     p_prog_month_local = p_cycle_psa,
                     treat_dur_cycles   = treat_dur_72w_cycles,
                     treat_start_cycles = treat_start_immediate)

  traces_12 <- lapply(treatments, function(stg)
    run_markov(aP_12[,,,stg], v_init))
  names(traces_12) <- treatments

  res_12 <- summarize_strategies(
    traces_12, "Age12",
    util_mat_psa, v_background_cost_cycle,
    cost_vec_psa, drug_psa,
    treat_dur_cycles_vec = treat_dur_72w_cycles
  )

  ## ---- Strategy 3: Wait until 18 (natural history age 12-18, then 72 weeks) ----
  aP_18 <- build_a_P(rr_reg_psa, rr_prog_psa,
                     p_prog_month_local = p_cycle_psa,
                     treat_dur_cycles   = treat_dur_72w_cycles,
                     treat_start_cycles = treat_start_age18)

  traces_18 <- lapply(treatments, function(stg)
    run_markov(aP_18[,,,stg], v_init))
  names(traces_18) <- treatments

  res_18 <- summarize_strategies(
    traces_18, "Age18",
    util_mat_psa, v_background_cost_cycle,
    cost_vec_psa, drug_psa,
    treat_dur_cycles_vec = treat_dur_72w_cycles
  )

  ## ---- Return named list with 3 strategies ----
  list(
    "LSM"               = c(Cost = res_lsm$Cost[res_lsm$Strategy == "LSM"],
                            QALY = res_lsm$QALY[res_lsm$Strategy == "LSM"]),
    "Sema 72w (Age 12)" = c(Cost = res_12$Cost[res_12$Strategy == "Semaglutide"],
                            QALY = res_12$QALY[res_12$Strategy == "Semaglutide"]),
    "Sema 72w (Age 18)" = c(Cost = res_18$Cost[res_18$Strategy == "Semaglutide"],
                            QALY = res_18$QALY[res_18$Strategy == "Semaglutide"])
  )
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
