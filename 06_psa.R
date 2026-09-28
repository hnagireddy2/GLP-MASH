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

  se_from_ci95 <- function(lo, hi) (hi - lo) / (2 * 1.96)
  rgamma_ci <- function(n, mean, lo, hi) {
    p <- gamma_params(mean, se_from_ci95(lo, hi))
    rgamma(n, p$shape, scale = p$scale)
  }
  rbeta_se <- function(n, mean, se) {
    p <- beta_params(mean, se)
    rbeta(n, p$shape1, p$shape2)
  }

  ## Rank-order-preserving correlated sampling
  induce_rank_order <- function(U, tol = 0.001, rho_step = 0.01,
                                max_inflate = 0.4, inflate_step = 0.05) {
    n <- nrow(U); k <- ncol(U)

    ## Step 1: minimum sufficient correlation for each pair, independently
    pair_rho <- function(i, j) {
      best <- 0
      for (rho in seq(1 - rho_step, 0, by = -rho_step)) {
        L <- chol(matrix(c(1, rho, rho, 1), 2, 2))
        Y <- matrix(rnorm(n * 2), n, 2) %*% L
        ui <- U[, i]; uj <- U[, j]
        ui[order(Y[, 1])] <- sort(U[, i])
        uj[order(Y[, 2])] <- sort(U[, j])
        if (mean(ui >= uj) <= tol) best <- rho else break
      }
      best
    }

    Sigma_raw <- diag(k)
    for (i in 1:(k - 1)) for (j in (i + 1):k) {
      Sigma_raw[i, j] <- Sigma_raw[j, i] <- pair_rho(i, j)
    }

    ## Steps 2-3: nearest valid correlation matrix, then joint sampling
    adjust_and_sample <- function(Sigma) {
      eig        <- eigen(Sigma, symmetric = TRUE)
      vals_fixed <- pmax(eig$values, 1e-8)             # clip negative eigenvalues
      Sigma_psd  <- eig$vectors %*% diag(vals_fixed) %*% t(eig$vectors)
      d          <- sqrt(diag(Sigma_psd))
      Sigma_adj  <- Sigma_psd / outer(d, d); diag(Sigma_adj) <- 1

      Y <- matrix(rnorm(n * k), n, k) %*% chol(Sigma_adj)
      U_out <- U
      for (j in 1:k) U_out[order(Y[, j]), j] <- sort(U[, j])
      U_out
    }

    U_out <- adjust_and_sample(Sigma_raw)
    viol  <- mean(apply(U_out, 1, is.unsorted))

    ## Step 4: realized-violation safety net 
    inflate <- 0
    while (viol > tol && inflate < max_inflate) {
      inflate <- inflate + inflate_step
      Sigma_try <- Sigma_raw
      Sigma_try[upper.tri(Sigma_try)] <- pmin(Sigma_raw[upper.tri(Sigma_raw)] + inflate, 0.999)
      Sigma_try[lower.tri(Sigma_try)] <- t(Sigma_try)[lower.tri(Sigma_try)]
      U_out <- adjust_and_sample(Sigma_try)
      viol  <- mean(apply(U_out, 1, is.unsorted))
    }

    if (viol > tol) {
      warning("induce_rank_order(): violation rate ", round(viol, 4),
              " still above tolerance after inflation retries.")
    }
    U_out
  }

  ## Fibrosis-transition mean/low/high (annual), tied to the Le et al.
  ## candidate set calibration selected 
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

    ## Advanced-disease hazards (annual)
    h_F3_HCC       = rgamma_ci(n_sim, 0.0034, 0.0021, 0.0047),
    h_F4_HCC       = rgamma_ci(n_sim, 0.0378, 0.0213, 0.0543),
    h_F4_DCC       = rgamma_ci(n_sim, 0.0659, 0.0400, 0.0918),
    h_DCC_HCC      = rgamma_ci(n_sim, 0.0378, 0.0213, 0.0543),
    h_DCC_LT       = rgamma_ci(n_sim, 0.0230, 0.0140, 0.0320),
    h_DCC_Death    = rgamma_ci(n_sim, 0.20,   0.1216, 0.2784),
    h_HCC_LT       = rgamma_ci(n_sim, 0.0300, 0.0182, 0.0418),
    h_HCC_Death    = rgamma_ci(n_sim, 0.1305, 0.1049, 0.1561),
    h_LT_Death     = rgamma_ci(n_sim, 0.0157, 0.0095, 0.0218),
    h_PostLT_Death = rgamma_ci(n_sim, 0.040365, 0.031592, 0.049058),

    ## STATE COSTS (gamma) 
    cost_F0_F2 = rgamma_ci(n_sim, costs_base["F0"], costs_low["F0"], costs_high["F0"]),
    cost_F3_raw    = rgamma_ci(n_sim, costs_base["F3"],    costs_low["F3"],    costs_high["F3"]),
    cost_F4_CC_raw = rgamma_ci(n_sim, costs_base["F4_CC"], costs_low["F4_CC"], costs_high["F4_CC"]),
    cost_HCC_raw   = rgamma_ci(n_sim, costs_base["HCC"],   costs_low["HCC"],   costs_high["HCC"]),
    cost_DCC_raw   = rgamma_ci(n_sim, costs_base["DCC"],   costs_low["DCC"],   costs_high["DCC"]),

    ## LT procedure cost (gamma)
    cost_LT = rgamma_ci(n_sim, costs_base["LT"], costs_low["LT"], costs_high["LT"]),

    ## HEALTH STATE UTILITIES (decrement, beta)
    qdec_F0_F2  = rbeta_se(n_sim, qaly_dec_base["F0"],      qaly_dec_base["F0"]      * 0.10),
    qdec_F3     = rbeta_se(n_sim, qaly_dec_base["F3"],      qaly_dec_base["F3"]      * 0.10),
    qdec_F4_CC  = rbeta_se(n_sim, qaly_dec_base["F4_CC"],   qaly_dec_base["F4_CC"]   * 0.10),
    qdec_DCC    = rbeta_se(n_sim, qaly_dec_base["DCC"],     qaly_dec_base["DCC"]     * 0.10),
    qdec_HCC    = rbeta_se(n_sim, qaly_dec_base["HCC"],     qaly_dec_base["HCC"]     * 0.10),
    qdec_LT     = rbeta_se(n_sim, qaly_dec_base["LT"],      qaly_dec_base["LT"]      * 0.10),
    qdec_PostLT = rbeta_se(n_sim, qaly_dec_base["Post_LT"], qaly_dec_base["Post_LT"] * 0.10)

  )
# Induce rank order (F3 < F4_CC < HCC < DCC) via correlated resampling
  cost_ordered <- induce_rank_order(cbind(df$cost_F3_raw, df$cost_F4_CC_raw,
                                          df$cost_HCC_raw, df$cost_DCC_raw))
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
  # LT_Death is a one-time cumulative probability 
  p_cycle_psa$LT_Death  <- pmin(psa_row$h_LT_Death, 0.999)
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
  qdec_psa["LT"]       <- psa_row$qdec_LT
  qdec_psa["Post_LT"]  <- psa_row$qdec_PostLT

  util_mat_psa <- build_util_matrix(v_util_age_base, qdec_psa)

  ## ---- LSM / Age12 / Age18, all three strategies at once ----
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

test_res <- run_model_psa_iter_all(df_psa_input_all[1, ])
cat("Labels returned by function:\n")
print(names(test_res))
cat("\nLabels expected (all_strat_labels):\n")
print(all_strat_labels)
cat("\nMismatches:\n")
print(setdiff(all_strat_labels, names(test_res)))
