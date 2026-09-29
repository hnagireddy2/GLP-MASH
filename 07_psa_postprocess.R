# 07_psa_postprocess.R
# Requires: source("06_psa.R") (df_c_all, df_e_all, df_psa_input_all)
# PSA parameter distributions, CE scatter, CEAC, expected loss
  
#############################################################
##########  PSA Parameter Distributions Plot  ###############
#############################################################

txtsize <- 13
v_wtp <- seq(0, 600000, by = 10000)

# Display names (transitions = annual probabilities; LT death = per cycle)
param_display_names <- c(
  rr_sema_regress  = "Sema RR: Regression",
  rr_sema_progress = "Sema RR: Progression",
  p_F0_F1          = "F0\u2192F1",
  p_F1_F0          = "F1\u2192F0",
  p_F1_F2          = "F1\u2192F2",
  p_F2_F1          = "F2\u2192F1",
  p_F2_F3          = "F2\u2192F3",
  p_F3_F2          = "F3\u2192F2",
  p_F3_F4          = "F3\u2192F4",
  p_F4_F3          = "F4\u2192F3",
  p_F4_DCC         = "F4\u2192DCC",
  p_F3_HCC         = "F3\u2192HCC",
  p_F4_HCC         = "F4\u2192HCC",
  p_DCC_HCC        = "DCC\u2192HCC",
  p_DCC_LT         = "DCC\u2192LT",
  p_HCC_LT         = "HCC\u2192LT",
  p_DCC_Death      = "DCC\u2192Death",
  p_HCC_Death      = "HCC\u2192Death",
  p_LT_Death       = "LT\u2192Death (perioperative)",
  p_PostLT_Death   = "Post-LT\u2192Death",
  qdec_F0_F2       = "Util dec: F0\u2013F2",
  qdec_F3_F4       = "Util dec: F3/F4",
  qdec_DCC         = "Util dec: DCC",
  qdec_HCC         = "Util dec: HCC",
  qdec_LT          = "Util dec: LT",
  qdec_PostLT      = "Util dec: Post-LT"
)

df_param_dist_long <- df_psa_input_all[, names(param_display_names)] %>%
  reshape2::melt(variable.name = "Parameter", value.name = "Value") %>%
  mutate(Parameter = factor(
    param_display_names[as.character(Parameter)],
    levels = unname(param_display_names)
  ))

ggplot(df_param_dist_long, aes(x = Value)) +
  geom_histogram(aes(y = after_stat(density)), 
                 col = "black", fill = "gray80", bins = 60) +
  geom_density(color = "#E69F00", linewidth = 1) +
  facet_wrap(~ Parameter, scales = "free", ncol = 3) +
  labs(
    x     = "Parameter Value",
    y     = "",
    title = "PSA: Parameter Distributions — All Sampled Parameters"
  ) +
  scale_y_continuous("", breaks = NULL) +
  theme_bw(base_size = txtsize)
#############################################################
#########  PSA Cloud and Table: Strategies ##############
#############################################################

## ---- PSA object with strategies --------------------------------------
l_psa_all <- make_psa_obj(
  cost          = df_c_all,
  effectiveness = df_e_all,
  parameters    = df_psa_input_all,
  strategies    = all_strat_labels
)

## Restore readable strategy names 
l_psa_all$strategies <- all_strat_labels
colnames(l_psa_all$cost) <- colnames(l_psa_all$effectiveness) <- all_strat_labels

## ---- Compute PSA means and run calculate_icers -----------------------------
psa_means_all <- summary(l_psa_all)

icer_psa_all <- calculate_icers(
  cost       = psa_means_all$meanCost,
  effect     = psa_means_all$meanEffect,
  strategies = psa_means_all$Strategy
)

## ---- Updated CEA table ---------------------------------------------
tbl_psa_frontier <- icer_psa_all %>%
  transmute(
    Strategy         = Strategy,
    `Cost ($)`       = scales::dollar(round(Cost)),
    `QALYs`          = round(Effect, 3),
    `Delta Cost ($)` = ifelse(is.na(Inc_Cost),
                              "\u2014", scales::dollar(round(Inc_Cost))),
    `Delta QALYs`    = ifelse(is.na(Inc_Effect),
                              "\u2014", as.character(round(Inc_Effect, 3))),
    `ICER ($/QALY)`  = dplyr::case_when(
      is.na(ICER)    ~ "\u2014",
      Status == "D"  ~ "Dominated",
      Status == "ED" ~ "Ext. Dominated",
      TRUE           ~ scales::dollar(round(ICER))
    ),
    Dominance        = dplyr::case_when(
      Status == "ND" ~ "Non-dominated",
      Status == "D"  ~ "Dominated",
      Status == "ED" ~ "Weakly dominated"
    )
  )

## ---- CE Scatter (PSA cloud for strategies) --------------------------
plot(l_psa_all) +
  labs(
    title    = "PSA: Cost-Effectiveness Scatter — All Strategies",
    subtitle = paste0("Each point = 1 PSA simulation  |  n = ", n_sim_all)
  ) +
  scale_y_continuous("Discounted Cost (USD)", labels = label_dollar(scale = 1)) +
  xlab("Discounted QALYs")
#############################################################
####### All-Strategy CEAC & Expected Loss — Post-PSA ########
#############################################################

okabe_ito_3 <- c("#999999", "#0072B2", "#D55E00")
names(okabe_ito_3) <- all_strat_labels

v_wtp_all <- seq(0, 600000, by = 10000)

# ── 2. Post-PSA Cost-Effectiveness Plane (mean estimates) ───
df_psa_all_means <- data.frame(
  Strategy = all_strat_labels,
  Cost     = sapply(all_strat_labels, function(s) mean(df_c_all[[s]])),
  QALY     = sapply(all_strat_labels, function(s) mean(df_e_all[[s]]))
) %>%
  arrange(QALY)

df_frontier_all <- calculate_icers(
  cost       = df_psa_all_means$Cost,
  effect     = df_psa_all_means$QALY,
  strategies = df_psa_all_means$Strategy
)

df_psa_all_means <- df_psa_all_means %>%
  mutate(
    age_grp = dplyr::case_when(
      Strategy == "LSM"                   ~ "LSM",
      grepl("Age 12", Strategy)           ~ "Treat Age 12",
      grepl("Age 18", Strategy)           ~ "Treat Age 18"
    ),
    dom_status = dplyr::case_when(
      Strategy %in% filter(df_frontier_all, Status == "D")$Strategy  ~ "Dominated",
      Strategy %in% filter(df_frontier_all, Status == "ED")$Strategy ~ "Weakly dominated",
      TRUE ~ "Non-dominated"
    )
  )

frontier_nd <- df_frontier_all %>% filter(Status == "ND") %>% arrange(Effect)

ggplot(df_psa_all_means, aes(x = QALY, y = Cost,
                              shape = dom_status, color = age_grp)) +
  geom_line(data        = frontier_nd,
            mapping     = aes(x = Effect, y = Cost),
            inherit.aes = FALSE,
            color = "grey35", linetype = "dashed", linewidth = 0.9) +
  geom_point(size = 4, stroke = 1.2) +
  geom_text_repel(aes(label = Strategy),
                  size = 3.1, max.overlaps = 25, box.padding = 0.45,
                  segment.color = NA) +
  scale_shape_manual(
    name   = "Dominance status",
    values = c("Non-dominated"    = 16,
               "Dominated"        = 17,
               "Weakly dominated" = 15)
  ) +
  scale_color_manual(
    name   = "Strategy",
    values = c("LSM"          = "grey40",
               "Treat Age 12" = "#0072B2",
               "Treat Age 18" = "#D55E00")
  ) +
  scale_x_continuous("QALYs",       breaks = pretty_breaks(n = 6)) +
  scale_y_continuous("Cost (USD)",
                     labels = label_dollar(scale = 1),
                     breaks = pretty_breaks(n = 6)) +
  labs(
    title    = "Cost-Effectiveness Frontier",
    subtitle = paste0(
      "Points = costs & QALYs (n = ", n_sim_all, " simulations)\n",
      "\u25cf Non-dominated  \u25b2 Dominated"
    )
  ) +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom", legend.box = "vertical",
        panel.grid.major = element_line(color = "grey90"),
        panel.grid.minor = element_blank())

print(df_frontier_all)

# ── 3. CEAC — All Strategies ──────────────────────────────────
ceac_all <- ceac(wtp = v_wtp_all, psa = l_psa_all)

plot(ceac_all) +
  scale_color_manual(values = unname(okabe_ito_3)) +
  labs(
    title    = "Cost-Effectiveness Acceptability Curve — All Strategies",
    subtitle = "Probability each strategy is cost-effective at each WTP threshold",
    x        = "Willingness-to-Pay Threshold ($/QALY)",
    y        = "Probability Cost-Effective"
  ) +
  theme_bw(base_size = 14) +
  theme(legend.position  = "right",
        legend.text      = element_text(size = 8),
        legend.key.width = unit(1.5, "cm"))

# ── 4. Expected Loss Curve — All Strategies ───────────────────
elc_all <- calc_exp_loss(wtp = v_wtp_all, psa = l_psa_all)

plot(elc_all) +
  scale_color_manual(values = unname(okabe_ito_3)) +
  labs(
    title    = "Expected Loss Curve — All Strategies",
    subtitle = "Expected opportunity loss (foregone net benefit) at each WTP threshold",
    x        = "Willingness-to-Pay Threshold ($/QALY)",
    y        = "Expected Loss ($)"
  ) +
  theme_bw(base_size = 14) +
  theme(legend.position  = "right",
        legend.text      = element_text(size = 8),
        legend.key.width = unit(1.5, "cm"))


#############################################################
######  95% Uncertainty Interval for ICER  ##################
######  Sema Age 12 vs LSM (non-dominated pair)  ############
#############################################################
# Requires df_c_all, df_e_all from 06_psa.R (one row per PSA draw)

# ---- 0. Exact strategy labels ----
s_trt <- "Sema 72w (Age 12)"   
s_ref <- "LSM"

# ---- 1. Per-draw incrementals = parameter-uncertainty cloud ----
d_cost <- df_c_all[[s_trt]] - df_c_all[[s_ref]]
d_qaly <- df_e_all[[s_trt]] - df_e_all[[s_ref]]

ok     <- is.finite(d_cost) & is.finite(d_qaly)
if (any(!ok)) cat(sprintf("Dropped %d of %d draws with NA/Inf\n", sum(!ok), length(ok)))
d_cost <- d_cost[ok]
d_qaly <- d_qaly[ok]

# ---- 2. Point ICER  ----
icer_point <- mean(d_cost) / mean(d_qaly)

# ---- 3. Diagnostic: fraction of draws with ΔQALY <= 0 ----
#    If this is non-trivial (>~1-2%), the ratio interval is unstable;
#    report INMB/CEAC instead.
frac_dE_le0 <- mean(d_qaly <= 0)

# ---- 4. HEADLINE: 95% uncertainty interval, percentile method ----
icer_draws <- d_cost / d_qaly
icer_ui    <- quantile(icer_draws, c(.025, .975), na.rm = TRUE)

# ---- 5. INMB-based interval + P(cost-effective) at WTP thresholds ----
# wtp = willingness-to-pay threshold ($/QALY); base-case WTP = $150,000
wtp_vec <- c(50000, 100000, 150000)
inmb_tbl <- do.call(rbind, lapply(wtp_vec, function(wtp) {
  inmb <- wtp * d_qaly - d_cost                # INMB = WTP x dQALY - dCost
  data.frame(
    WTP        = wtp,
    INMB_mean  = mean(inmb),
    INMB_lower = quantile(inmb, .025),
    INMB_upper = quantile(inmb, .975),
    P_CE       = mean(inmb > 0)                 # P(cost-effective vs LSM) at this WTP
  )
}))

# ---- 6. Fieller-type interval (robustness check) ----
fieller_ci <- function(dc, de, level = .95) {
  ok <- is.finite(dc) & is.finite(de)        # drop NA/NaN/Inf pairs
  dc <- dc[ok]; de <- de[ok]
  if (length(dc) < 2) return(c(lower = NA, upper = NA))
  
  mc <- mean(dc); me <- mean(de)
  v11 <- var(dc); v22 <- var(de); v12 <- cov(dc, de)
  z   <- qnorm(1 - (1 - level) / 2)
  a   <- me^2 - z^2 * v22
  b   <- -2 * (mc * me - z^2 * v12)
  cc  <- mc^2 - z^2 * v11
  disc <- b^2 - 4 * a * cc
  
  # guard against any residual NA in the condition
  if (!is.finite(a) || !is.finite(disc) || isTRUE(a <= 0) || isTRUE(disc < 0))
    return(c(lower = NA, upper = NA))
  sort((-b + c(-1, 1) * sqrt(disc)) / (2 * a))
}

icer_ci_fieller <- fieller_ci(d_cost, d_qaly)

# ---- 7. Print summary ----
cat(sprintf("\nSema Age 12 vs LSM  (n = %d draws)\n", length(d_cost)))
cat(sprintf("  Point ICER          : $%s/QALY\n", format(round(icer_point), big.mark = ",")))
cat(sprintf("  95%% UI (percentile) : $%s to $%s/QALY\n",
            format(round(icer_ui[1]), big.mark = ","),
            format(round(icer_ui[2]), big.mark = ",")))
cat(sprintf("  95%% Fieller-type interval (robustness): $%s to $%s/QALY\n",
            format(round(icer_ci_fieller[1]), big.mark = ","),
            format(round(icer_ci_fieller[2]), big.mark = ",")))
cat(sprintf("  Frac draws ΔQALY<=0  : %.3f%%\n\n", 100 * frac_dE_le0))
print(inmb_tbl, row.names = FALSE)

