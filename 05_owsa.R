# 05_owsa.R
# Requires: source("00_parameters.R"), source("00b_le_transitions.R"), source("02_calibration.R")
# One-way sensitivity analysis + threshold sensitivity analysis 

#############################################################
############# One-Way Sensitivity Analysis ##################
#############################################################

owsa_params <- list()

###### 1. Fibrosis progression (Le et al. 95% CI)
icer_fibprog_lo <- run_owsa_icer(p_prog_cycle_local = p_cycle_fibprog_low)
icer_fibprog_hi <- run_owsa_icer(p_prog_cycle_local = p_cycle_fibprog_high)
owsa_params[["fib_progression"]] <- c(low = icer_fibprog_lo, high = icer_fibprog_hi)

###### 2. Fibrosis regression (Le et al. 95% CI)
icer_fibreg_lo <- run_owsa_icer(p_prog_cycle_local = p_cycle_fibreg_low)
icer_fibreg_hi <- run_owsa_icer(p_prog_cycle_local = p_cycle_fibreg_high)
owsa_params[["fib_regression"]] <- c(low = icer_fibreg_lo, high = icer_fibreg_hi)

###### 3. Semaglutide regression RR (ESSENCE 95% CI)
rr_reg_vec_lo <- rr_regress; rr_reg_vec_lo["Semaglutide"] <- RR_reg_lo
rr_reg_vec_hi <- rr_regress; rr_reg_vec_hi["Semaglutide"] <- RR_reg_hi

icer_rrreg_lo <- run_owsa_icer(rr_regress_vec = rr_reg_vec_lo)
icer_rrreg_hi <- run_owsa_icer(rr_regress_vec = rr_reg_vec_hi)
owsa_params[["sema_rr_regression"]] <- c(low = icer_rrreg_lo, high = icer_rrreg_hi)

###### 4. Semaglutide progression RR (ESSENCE 95% CI)
rr_prog_vec_lo <- rr_progress; rr_prog_vec_lo["Semaglutide"] <- RR_prog_lo
rr_prog_vec_hi <- rr_progress; rr_prog_vec_hi["Semaglutide"] <- RR_prog_hi

icer_rrprog_lo <- run_owsa_icer(rr_progress_vec = rr_prog_vec_lo)
icer_rrprog_hi <- run_owsa_icer(rr_progress_vec = rr_prog_vec_hi)
owsa_params[["sema_rr_progression"]] <- c(low = icer_rrprog_lo, high = icer_rrprog_hi)

###### 5. State-specific medical costs (0.75x vs 1.25x)
cost_lo <- costs_base * 0.75
cost_hi <- costs_base * 1.25

icer_cost_lo <- run_owsa_icer(cost_vector = cost_lo)
icer_cost_hi <- run_owsa_icer(cost_vector = cost_hi)
owsa_params[["state_costs"]] <- c(low = icer_cost_lo,
                                     high = icer_cost_hi)

###### 6. LT procedure cost (±20%)
costs_lt_lo <- costs_base
costs_lt_hi <- costs_base
costs_lt_lo["LT"] <- costs_low["LT"]
costs_lt_hi["LT"] <- costs_high["LT"]

icer_lt_lo <- run_owsa_icer(cost_vector = costs_lt_lo)
icer_lt_hi <- run_owsa_icer(cost_vector = costs_lt_hi)
owsa_params[["lt_cost"]] <- c(low = icer_lt_lo, high = icer_lt_hi)

###### 7. Annual cost of semaglutide
lo_drug <- drug_cost; lo_drug["Semaglutide"] <- cost_sema_low
hi_drug <- drug_cost; hi_drug["Semaglutide"] <- cost_sema_high

icer_treat_lo <- run_owsa_icer(drug_cost_vec = lo_drug)
icer_treat_hi <- run_owsa_icer(drug_cost_vec = hi_drug)
owsa_params[["sema_cost"]] <- c(low = icer_treat_lo,
                                   high = icer_treat_hi)

###### 8. Utility scale (QALYs)
icer_u_lo <- run_owsa_icer(util_matrix = m_util_base * 0.9)
icer_u_hi <- run_owsa_icer(util_matrix = pmin(m_util_base * 1.1, 1))
owsa_params[["qalys"]] <- c(low = icer_u_lo,
                               high = icer_u_hi)

base_icer <- run_owsa_icer()

#############################################################
## DAMPACK OWSA — NMB for THREE strategies
#############################################################

wtp_threshold    <- 150000

calculate_ce_out_mash <- function(l_params, n_wtp = wtp_threshold) {

  v_dwc_saved <- v_dwc
  v_dwu_saved <- v_dwu
  on.exit({
    v_dwc <<- v_dwc_saved
    v_dwu <<- v_dwu_saved
  })
  disc <- l_params[["Discount rate"]]
  v_dwc <<- 1/(1+disc)^((0:n_cycles)*cycle_length)
  v_dwu <<- 1/(1+disc)^((0:n_cycles)*cycle_length)

  rr_reg_loc  <- rr_regress
  rr_reg_loc["Semaglutide"]  <- l_params[["Semaglutide regression RR"]]

  rr_prog_loc <- rr_progress
  rr_prog_loc["Semaglutide"] <- l_params[["Semaglutide progression RR"]]

  cost_loc  <- costs_base * l_params[["State medical costs (0.75x-1.25x)"]]

  drug_loc  <- drug_cost
  drug_loc["Semaglutide"] <- l_params[["Semaglutide annual cost"]]

  qdec_loc     <- pmin(qaly_dec_base * l_params[["Utility decrements (0.9x-1.1x)"]], 1)
  util_mat_loc <- build_util_matrix(v_util_age_base, qdec_loc)

  ## Fibrosis progression / regression rate multipliers
  prog_sc <- l_params[["Fibrosis progression (0.75x-1.25x)"]]
  reg_sc  <- l_params[["Fibrosis regression (0.75x-1.25x)"]]
  p_prog_loc <- p_prog_cycle
  for (nm in fib_prog_names) {
    p_prog_loc[[nm]] <- rate_to_prob(prob_to_rate(p_prog_cycle[[nm]], cycle_length) * prog_sc,
                                     cycle_length)
  }
  for (nm in fib_reg_names) {
    p_prog_loc[[nm]] <- rate_to_prob(prob_to_rate(p_prog_cycle[[nm]], cycle_length) * reg_sc,
                                     cycle_length)
  }
  p_prog_loc$DCC_Death <- annual_to_cycle(l_params[["DCC->Death (annual)"]])
  p_prog_loc$HCC_Death <- annual_to_cycle(l_params[["HCC->Death (annual)"]])
  p_prog_loc$F4_DCC    <- annual_to_cycle(l_params[["F4->DCC (annual)"]])

  results <- run_three_strategies(rr_reg_loc, rr_prog_loc,
                                  p_prog_cycle_local  = p_prog_loc,
                                  util_matrix          = util_mat_loc,
                                  cost_vector           = cost_loc,
                                  drug_cost_vec         = drug_loc,
                                  treat_dur_cycles_vec  = treat_dur_72w_cycles)

  df_out <- data.frame(
    Strategy = names(results),
    Cost     = sapply(results, `[[`, "Cost"),
    Effect   = sapply(results, `[[`, "QALY"),
    stringsAsFactors = FALSE
  )
  df_out$NMB <- df_out$Effect * n_wtp - df_out$Cost

  df_out
}

df_params_owsa <- data.frame(
  pars = c("Semaglutide regression RR",
           "Semaglutide progression RR",
           "Fibrosis progression (0.75x-1.25x)",
           "Fibrosis regression (0.75x-1.25x)",
           "State medical costs (0.75x-1.25x)",
           "Semaglutide annual cost",
           "Utility decrements (0.9x-1.1x)",
           "Discount rate",
           "DCC->Death (annual)", 
           "HCC->Death (annual)", 
           "F4->DCC (annual)"),
  min = c(RR_reg_lo, RR_prog_lo, 0.75, 0.75, 0.75, cost_sema_low,  0.90, 0.00, 0.1216, 0.1049, 0.0400),
  max = c(RR_reg_hi, RR_prog_hi, 1.25, 1.25, 1.25, cost_sema_high, 1.10, 0.05, 0.2784, 0.1561, 0.0918)
)

l_params_basecase <- list(
  "Semaglutide regression RR"              = RR_regress,
  "Semaglutide progression RR"             = RR_progress,
  "Fibrosis progression (0.75x-1.25x)"     = 1.0,
  "Fibrosis regression (0.75x-1.25x)"      = 1.0,
  "State medical costs (0.75x-1.25x)"      = 1.0,
  "Semaglutide annual cost"                = cost_sema_base,
  "Utility decrements (0.9x-1.1x)"         = 1.0,
  "Discount rate"                          = 0.03,
  "DCC->Death (annual)"                    = 0.20,
  "HCC->Death (annual)"                    = 0.1305,
  "F4->DCC (annual)"                       = 0.0659
)

owsa_nmb <- run_owsa_det(
  params_range    = df_params_owsa,
  params_basecase = l_params_basecase,
  nsamp           = 100,
  FUN             = calculate_ce_out_mash,
  outcomes        = "NMB",
  strategies      = all_strat_labels,
  n_wtp           = wtp_threshold
)

test_out <- calculate_ce_out_mash(l_params_basecase)
cat("Strategies returned by function:\n")
print(test_out$Strategy)
cat("\nall_strat_labels passed to run_owsa_det:\n")
print(all_strat_labels)
cat("\nMatch:", identical(sort(test_out$Strategy), sort(all_strat_labels)), "\n")

# NMB OWSA plot
plot(owsa_nmb, txtsize = 11, n_x_ticks = 4, facet_scales = "free") +
  labs(y = "NMB ($ per patient)",
       title    = "One-Way Sensitivity Analysis: NMB by Parameter",
       subtitle = paste0("WTP = $150,000/QALY -- Strategies vs. LSM")) +
  theme(legend.position = "bottom")

# Optimal strategy by parameter 
owsa_opt_strat(owsa = owsa_nmb, txtsize = 11)

#############################################################
######### THRESHOLD SENSITIVITY ANALYSIS ####################
#############################################################

find_threshold <- function(param_name, base_val,
                           lo_val, hi_val,
                           run_fn,
                           wtp   = wtp_threshold,
                           tol   = 1) {

  icer_lo <- run_fn(lo_val)
  icer_hi <- run_fn(hi_val)

  crossing_exists <- !is.na(icer_lo) & !is.na(icer_hi) &
                     ((icer_lo - wtp) * (icer_hi - wtp) < 0)

  if (!crossing_exists) {
    return(data.frame(
      Parameter  = param_name,
      Base_Value = base_val,
      Threshold  = NA_real_,
      Direction  = ifelse(icer_lo > wtp & icer_hi > wtp,
                          "Always dominated",
                          "Always cost-effective"),
      Base_ICER  = round(run_fn(base_val))
    ))
  }

  a <- lo_val; b <- hi_val
  for (iter in 1:60) {
    mid      <- (a + b) / 2
    icer_mid <- run_fn(mid)
    if (is.na(icer_mid)) break
    if (abs(icer_mid - wtp) < tol) break
    if ((icer_mid - wtp) * (icer_lo - wtp) < 0) b <- mid else a <- mid
  }

  data.frame(
    Parameter  = param_name,
    Base_Value = base_val,
    Threshold  = round(mid, 4),
    Direction  = ifelse(run_fn(lo_val) < wtp,
                        "CE at low; flips at threshold",
                        "Not CE at low; flips at threshold"),
    Base_ICER  = round(run_fn(base_val))
  )
}

frontier_strats <- list(
  "Semaglutide 72w Age 12" = list(dur = treat_dur_72w_cycles,
                           start = treat_start_immediate),
  "Semaglutide 72w Age 18" = list(dur = treat_dur_72w_cycles,
                           start = treat_start_age18)
)

thresh_all_strats <- lapply(names(frontier_strats), function(strat_name) {
  strat_info <- frontier_strats[[strat_name]]
  dur_cyc    <- strat_info$dur
  start_cyc  <- strat_info$start

  results <- list()

  results[["sema_cost"]] <- find_threshold(
    "Semaglutide annual cost ($)", cost_sema_base, 1000, 80000,
    run_fn = function(x) {
      dv <- drug_cost; dv["Semaglutide"] <- x
      run_owsa_icer(drug_cost_vec = dv, treat_dur_cycles_vec = dur_cyc, treat_start_cycles = start_cyc)
    })

  results[["rr_regress"]] <- find_threshold(
    "Semaglutide regression RR", RR_regress, 1.0, RR_regress,
    run_fn = function(x) {
      rv <- rr_regress; rv["Semaglutide"] <- x
      run_owsa_icer(rr_regress_vec = rv, treat_dur_cycles_vec = dur_cyc, treat_start_cycles = start_cyc)
    })

  results[["rr_progress"]] <- find_threshold(
    "Semaglutide progression RR", RR_progress, RR_progress, 1.0,
    run_fn = function(x) {
      pv <- rr_progress; pv["Semaglutide"] <- x
      run_owsa_icer(rr_progress_vec = pv, treat_dur_cycles_vec = dur_cyc, treat_start_cycles = start_cyc)
    })

  results[["p_F3_F4"]] <- find_threshold(
    "F3->F4 annual probability", candidate_sets[[best]][["F3_F4"]], 0.01, 0.30,
    run_fn = function(x) {
      pm <- p_prog_cycle; pm$F3_F4 <- annual_to_cycle(x)
      run_owsa_icer(p_prog_cycle_local = pm, treat_dur_cycles_vec = dur_cyc, treat_start_cycles = start_cyc)
    })

  results[["state_cost"]] <- find_threshold(
    "State cost multiplier", 1.0, 0.5, 3.0,
    run_fn = function(x) {
      run_owsa_icer(cost_vector = costs_base * x, treat_dur_cycles_vec = dur_cyc, treat_start_cycles = start_cyc)
    })

  df <- dplyr::bind_rows(results)
  df$Strategy <- strat_name
  df
})

df_thresh_all <- dplyr::bind_rows(thresh_all_strats)

df_thresh_all_plot <- df_thresh_all %>%
  filter(!is.na(Threshold)) %>%
  mutate(
    pct_change = round((Threshold - Base_Value) / abs(Base_Value) * 100, 1),
    label_text = paste0("Base: ", signif(Base_Value, 3),
                        "\nThreshold: ", signif(Threshold, 3),
                        " (", ifelse(pct_change > 0, "+", ""), pct_change, "%)")
  )

df_thresh_all_plot <- df_thresh_all_plot %>%
  group_by(Strategy) %>%
  mutate(param_ordered = reorder(paste0(Strategy, "__", Parameter), abs(pct_change))) %>%
  ungroup()

ggplot(df_thresh_all_plot,
       aes(x = pct_change, y = param_ordered, fill = pct_change > 0)) +
  geom_col(alpha = 0.85, width = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.7, color = "grey40") +
  geom_text(aes(label = label_text, hjust = ifelse(pct_change >= 0, -0.05, 1.05)), size = 2.8) +
  facet_wrap(~ Strategy, ncol = 1, scales = "free_y") +
  scale_y_discrete(labels = function(x) sub("^.*?__", "", x)) +
  scale_fill_manual(
    values = c(`TRUE` = "#D55E00", `FALSE` = "#0072B2"),
    labels = c(`TRUE` = "Must increase to flip", `FALSE` = "Must decrease to flip"),
    name = "Direction"
  ) +
  scale_x_continuous(labels = function(x) paste0(x, "%"), expand = expansion(mult = c(0.45, 0.45))) +
  labs(title    = "Threshold Sensitivity Analysis \u2014 Semaglutide Treatment Strategies",
       subtitle = paste0("% change required to flip CE decision at WTP = $150K/QALY"),
       x = "% change from base value", y = NULL) +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

cat("05_owsa.R complete.\n")
