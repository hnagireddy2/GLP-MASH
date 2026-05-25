# 07_psa_postprocess.R
# Requires: 06_psa.R (df_c_all, df_e_all, df_psa_input_all)
# PSA parameter distributions, CE scatter, CEAC, expected loss

#############################################################
##########  PSA Parameter Distributions Plot  ###############
#############################################################

txtsize <- 13
v_wtp <- seq(0, 600000, by = 10000)

# Rename columns for readable facet labels
param_display_names <- c(
  rr_sema_regress  = "Sema RR: Regression",
  rr_sema_progress = "Sema RR: Progression",
  h_F0_F1          = "Hazard F0\u2192F1",
  h_F1_F0          = "Hazard F1\u2192F0",
  h_F1_F2          = "Hazard F1\u2192F2",
  h_F2_F1          = "Hazard F2\u2192F1",
  h_F2_F3          = "Hazard F2\u2192F3",
  h_F3_F2          = "Hazard F3\u2192F2",
  h_F3_F4          = "Hazard F3\u2192F4",
  h_F4_F3          = "Hazard F4\u2192F3",
  h_F4_DCC         = "Hazard F4\u2192DCC",
  h_DCC_Death      = "DCC\u2192Death",
  h_HCC_Death      = "HCC\u2192Death",
  h_LT_Death       = "LT\u2192Death (yr 1)",
  h_PostLT_Death   = "Post-LT\u2192Death",
  qdec_F0_F2       = "Util dec: F0\u2013F2",
  qdec_F3          = "Util dec: F3",
  qdec_F4_CC       = "Util dec: F4/CC",
  qdec_DCC         = "Util dec: DCC",
  qdec_HCC         = "Util dec: HCC",
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

## ---- Compute PSA means and run calculate_icers -----------------------------
psa_means_all <- summary(l_psa_all)

# CEAC and ELC using strategies
ceac_obj <- ceac(wtp = v_wtp, psa = l_psa_all)
elc_obj  <- calc_exp_loss(wtp = v_wtp, psa = l_psa_all)

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

okabe_ito_9 <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442",
  "#0072B2", "#D55E00", "#CC79A7", "#999999", "#000000"
)
names(okabe_ito_9) <- all_strat_labels

# ── 1. Build PSA object ─────────────────────────────────────
l_psa_all <- make_psa_obj(
  cost          = df_c_all[, all_strat_labels],
  effectiveness = df_e_all[, all_strat_labels],
  parameters    = df_psa_input_all,
  strategies    = all_strat_labels
)
 
v_wtp_all <- seq(0, 600000, by = 10000)

# ── 2. Post-PSA Cost-Effectiveness Plane (mean estimates) ───
df_psa_all_means <- data.frame(
  Strategy = all_strat_labels,
  Cost     = sapply(all_strat_labels, function(s) mean(df_c_all[[s]])),
  QALY     = sapply(all_strat_labels, function(s) mean(df_e_all[[s]]))
) %>%
  arrange(QALY)

# Define df_frontier_all FIRST so mutate() can reference it
df_frontier_all <- calculate_icers(
  cost       = df_psa_all_means$Cost,
  effect     = df_psa_all_means$QALY,
  strategies = df_psa_all_means$Strategy
)

df_psa_all_means <- df_psa_all_means %>%
  mutate(
    age_grp = dplyr::case_when(
      Strategy == "LSM"                   ~ "LSM",
      grepl("Start Age 12", Strategy)     ~ "Treat Age 12",
      grepl("Start Age 18", Strategy)     ~ "Treat Age 18"
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
               "Weakly dominated" = 1,
               "Dominated"        = 17)
  ) +
  scale_color_manual(
    name   = "Strategy",
    values = c("LSM"          = "grey40",
               "Treat Age 12" = "#0072B2",
               "Treat Age 18" = "#D55E00")
  ) +
  scale_x_continuous("Mean PSA QALYs",       breaks = pretty_breaks(n = 6)) +
  scale_y_continuous("Mean PSA Cost (USD)",
                     labels = label_dollar(scale = 1),
                     breaks = pretty_breaks(n = 6)) +
  labs(
    title    = "Post-PSA Combined CEA Frontier: All Strategies",
    subtitle = paste0(
      "Points = PSA mean costs & QALYs (n = ", n_sim_all, " simulations)\n",
      "\u25cf Non-dominated  \u25cb Weakly dominated  \u25b2 Dominated"
    )
  ) +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom", legend.box = "vertical",
        panel.grid.major = element_line(color = "grey90"),
        panel.grid.minor = element_blank())

# Print full ICER table
print(df_frontier_all)

# ── 3. CEAC — All Strategies ──────────────────────────────────
ceac_all <- ceac(wtp = v_wtp_all, psa = l_psa_all)

plot(ceac_all) +
  scale_color_manual(values = unname(okabe_ito_9)) +
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
  scale_color_manual(values = unname(okabe_ito_9)) +
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

