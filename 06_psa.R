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

  ## Rank-order-preserving correlated sampling, per Goldhaber-Fiebert JD,
  ## Jalal H. "Some Health States Are Better Than Others: Using Health
  ## State Rank Order to Improve Probabilistic Analyses." Med Decis Making.
  ## 2016;36(8):927-940.
  ##
  ## PROBLEM: F3 < F4_CC < HCC < DCC must hold in every PSA draw, but each
  ## state's cost also has its own literature-sourced marginal distribution
  ## (gamma, from costs_base/costs_low/costs_high) that should be preserved.
  ## Simply sorting 4 independent draws (the old approach) guarantees the
  ## ordering but destroys each state's marginal: the value assigned to F3
  ## becomes the MIN of 4 independent draws (not a draw from F3's own
  ## gamma distribution), DCC's becomes the MAX, etc. -- exactly the
  ## distortion the paper's abstract flags as a tradeoff of naive
  ## rank-preserving methods.
  ##
  ## METHOD (per-pair correlation matrix + eigenvalue correction + rank
  ## matching, following the paper's general algorithm):
  ## 1. For EACH of the 6 pairs among the 4 states independently, find the
  ##    smallest sufficient correlation (search from near 1 downward) that
  ##    keeps THAT PAIR's own violation rate below tolerance, using only
  ##    those two columns. States whose distributions are already far
  ##    apart (e.g. F3 ~$9k vs DCC ~$172k) typically need rho=0; only
  ##    genuinely close pairs (e.g. HCC ~$125k vs DCC ~$172k) need real
  ##    correlation. This avoids over-correlating pairs that don't need it
  ##    -- the failure mode of using one shared correlation for all pairs.
  ## 2. Assemble these 6 pairwise values into a single 4x4 matrix. Because
  ##    each was optimized independently, the assembled matrix isn't
  ##    guaranteed to be a valid (positive-semi-definite) correlation
  ##    matrix. Fix this via eigenvalue decomposition: clip any negative
  ##    eigenvalues up to a small positive floor, reconstruct, and
  ##    rescale to unit diagonal -- the nearest valid correlation matrix
  ##    to the pairwise-optimal one.
  ## 3. Generate ONE joint multivariate-normal reference from this final
  ##    matrix (not pairwise), then re-sort each column's own values into
  ##    the row positions given by that column's reference rank (Iman-
  ##    Conover rank matching) -- exactly as before, so marginals are
  ##    still preserved exactly.
  ## 4. Check the REALIZED violation rate across the full row. Note this
  ##    is NOT a separate check needed for non-adjacent pairs (e.g. F3 vs
  ##    HCC): if F3<F4_CC, F4_CC<HCC, and HCC<DCC all hold for a given
  ##    row's actual values, F3<HCC etc. hold automatically by
  ##    transitivity of real-number ordering -- no realized draw can
  ##    satisfy every adjacent comparison while failing a non-adjacent
  ##    one. The re-check IS still necessary for a different reason: step
  ##    2's eigenvalue correction can shrink the pairwise-optimal
  ##    correlations, and Iman-Conover only approximately induces the
  ##    target matrix in finite samples -- so the realized rate can drift
  ##    above tolerance even when every pairwise target was individually
  ##    sufficient. If it does, uniformly inflate the raw matrix's
  ##    off-diagonal entries and retry.
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

    ## Step 4: realized-violation safety net (see comment above -- this
    ## guards against eigenvalue-correction/finite-sample drift, not
    ## against non-adjacent pairs, which can't fail once adjacency holds)
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
    ## LT_Death: Sharma et al. 2018 (UNOS, n=64,977) 90-day all-cause
    ## mortality, background-netted. Source reports no CI, so the +/-39.2%
    ## relative width from the prior Rustgi-sourced estimate is carried
    ## forward, centered on the new mean.
    ## PostLT_Death: Bezinover et al. 2023, NASH/CC-specific AYA (15-39yo)
    ## curve, exponential fit to all 3 digitized points (years 1/3/5), with
    ## a proper binomial-SE-based 95% CI on the fitted annual hazard --
    ## see nonfib_annual$PostLT_Death in 00_parameters.R for the full
    ## derivation (lambda, delta-method CI, background-netting).
    ## NOTE: PostLT_Death uses the AYA-specific curve as a single flat rate
    ## for all Post_LT cycles, regardless of the patient's actual age at
    ## that point -- Bezinover also reports a 40-65yo curve that would be
    ## more age-appropriate for patients well past the AYA band. Flagged
    ## for review; not implemented (would need an age-indexed lookup, not
    ## just a different constant -- see build_a_P()'s per-cycle background
    ## mortality lookup for the analogous mechanism).
    h_LT_Death     = rgamma_ci(n_sim, 0.0157, 0.0095, 0.0218),
    h_PostLT_Death = rgamma_ci(n_sim, 0.040365, 0.031592, 0.049058),

    ## STATE COSTS (gamma) -- means/ranges from costs_base/costs_low/costs_high
    ## (00_parameters.R) so PSA always tracks the base-case cost inputs.
    cost_F0_F2 = rgamma_ci(n_sim, costs_base["F0"], costs_low["F0"], costs_high["F0"]),

    ## STATE COSTS (gamma) -- independent draws from each state's own
    ## literature-sourced marginal. Rank order (F3 < F4_CC < HCC < DCC) is
    ## induced afterward via induce_rank_order() -- see below -- rather
    ## than by sorting these draws directly.
    cost_F3_raw    = rgamma_ci(n_sim, costs_base["F3"],    costs_low["F3"],    costs_high["F3"]),
    cost_F4_CC_raw = rgamma_ci(n_sim, costs_base["F4_CC"], costs_low["F4_CC"], costs_high["F4_CC"]),
    cost_HCC_raw   = rgamma_ci(n_sim, costs_base["HCC"],   costs_low["HCC"],   costs_high["HCC"]),
    cost_DCC_raw   = rgamma_ci(n_sim, costs_base["DCC"],   costs_low["DCC"],   costs_high["DCC"]),

    ## LT procedure cost (gamma)
    cost_LT = rgamma_ci(n_sim, costs_base["LT"], costs_low["LT"], costs_high["LT"]),

    ## HEALTH STATE UTILITIES (decrement, beta)
    # SE is capped inside beta_params() to keep alpha/beta positive
    qdec_F0_F2  = rbeta_se(n_sim, qaly_dec_base["F0"],      qaly_dec_base["F0"]      * 0.10),
    qdec_F3     = rbeta_se(n_sim, qaly_dec_base["F3"],      qaly_dec_base["F3"]      * 0.10),
    qdec_F4_CC  = rbeta_se(n_sim, qaly_dec_base["F4_CC"],   qaly_dec_base["F4_CC"]   * 0.10),
    qdec_DCC    = rbeta_se(n_sim, qaly_dec_base["DCC"],     qaly_dec_base["DCC"]     * 0.10),
    qdec_HCC    = rbeta_se(n_sim, qaly_dec_base["HCC"],     qaly_dec_base["HCC"]     * 0.10),
    qdec_LT     = rbeta_se(n_sim, qaly_dec_base["LT"],      qaly_dec_base["LT"]      * 0.10),
    qdec_PostLT = rbeta_se(n_sim, qaly_dec_base["Post_LT"], qaly_dec_base["Post_LT"] * 0.10)

  )
# Induce rank order (F3 < F4_CC < HCC < DCC) via correlated resampling
  # (see induce_rank_order() above) instead of sorting -- preserves each
  # state's own marginal cost distribution exactly.
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
  # LT_Death is a one-time cumulative probability applied to LT's single
  # cycle, not a recurring rate -- bypass amh()'s annual_to_month(), same as
  # the base case (00_parameters.R).
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
