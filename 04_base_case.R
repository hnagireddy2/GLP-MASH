# 04_base_case.R
# Requires: 00_parameters.R, 01_model_functions.R, 02_calibration.R, 03_validation.R
# CEA frontier + external validation (uses traces from base case runs)

#############################################################
########## COMBINED CEA FRONTIER: TWO STRATEGIES ############
#############################################################

res_imm_full <- run_strategy(treat_start_immediate, "Treat at 12")
res_del_full <- run_strategy(treat_start_age18,     "Wait until 18")
res_imm <- res_imm_full$summary
res_del <- res_del_full$summary
traces_imm <- res_imm_full$traces
traces_del <- res_del_full$traces

## LSM is identical in both runs — deduplicate
res_combined <- bind_rows(res_imm, res_del) %>%
  mutate(StrategyLabel = case_when(
    Strategy == "LSM"         ~ "LSM",
    Scenario == "Treat at 12" ~ "Sema 72w (Age 12)",
    Scenario == "Wait until 18" ~ "Sema 72w (Age 18)"
  )) %>%
  distinct(StrategyLabel, .keep_all = TRUE) %>%
  arrange(Cost)

icer_res <- calculate_icers(
  cost       = res_combined$Cost,
  effect     = res_combined$QALY,
  strategies = res_combined$StrategyLabel
)

tbl_frontier <- icer_res %>%
  transmute(
    Strategy         = Strategy,
    `Cost ($)`       = scales::dollar(round(Cost)),
    `QALYs`          = round(Effect, 3),
    `Delta Cost ($)` = ifelse(is.na(Inc_Cost),
                              "\u2014", scales::dollar(round(Inc_Cost))),
    `Delta QALYs`    = ifelse(is.na(Inc_Effect),
                              "\u2014", as.character(round(Inc_Effect, 3))),
    `ICER ($/QALY)`  = case_when(
      is.na(ICER)    ~ "\u2014",
      Status == "D"  ~ "Dominated",
      Status == "ED" ~ "Ext. Dominated",
      TRUE           ~ scales::dollar(round(ICER))
    ),
    Dominance        = case_when(
      Status == "ND" ~ "Non-dominated",
      Status == "D"  ~ "Dominated",
      Status == "ED" ~ "Weakly dominated"
    )
  )

cat("\n====================================================\n")
cat(" CEA FRONTIER — Two Strategies\n")
cat(" Cohort: F2/F3 adolescents at age 12\n")
cat(" Treatment duration: 72 weeks (both arms)\n")
cat("====================================================\n")
print(tbl_frontier, row.names = FALSE)
cat("====================================================\n")

ggplot(icer_res, aes(x = Effect, y = Cost)) +
  geom_line(data = icer_res %>% filter(Status == "ND"),
            color = "grey35", linetype = "dashed", linewidth = 0.9) +
  geom_point(aes(color = Strategy), size = 5) +
  geom_text_repel(aes(label = Strategy), size = 3.5, box.padding = 0.5) +
  scale_color_manual(values = c("LSM" = "grey40",
                                "Sema 72w (Age 12)" = "#0072B2",
                                "Sema 72w (Age 18)" = "#D55E00")) +
  scale_x_continuous("Discounted QALYs", breaks = pretty_breaks(n = 6)) +
  scale_y_continuous("Discounted Cost (USD)",
                     labels = label_dollar(), breaks = pretty_breaks(n = 6)) +
  labs(title    = "CEA Frontier: Treat at 12 vs. Wait until 18",
       subtitle = "Cohort: F2/F3 adolescents at age 12; treatment duration = 72 weeks") +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom")


# ============================================================
# EXTERNAL VALIDATION
# ============================================================

txtsize <- 13

target_hcc_incidence  <- list(
  annual_low  = 0.005,
  annual_high = 0.030,
  source      = "Kanwal et al. Gastroenterology 2022; Sanyal et al. NEJM 2021"
)

target_dcc_rate <- list(
  annual_low  = 0.027,
  annual_high = 0.066,
  source      = "Sanyal et al. NEJM 2021; Le et al. CGH 2023"
)

target_dcc_mortality <- list(
  annual_low  = 0.10,
  annual_high = 0.25,
  source      = "Ng et al. CGH 2023; D'Amico et al. J Hepatol 2014"
)

trace_lsm  <- traces_imm[["LSM"]]
ages_trace <- age_start + (0:n_cycles) * cycle_length

# ---- Target 1: Annual HCC incidence from F4 ----
f4_window <- 100:200
hcc_monthly <- mean(
  diff(trace_lsm[f4_window, "HCC"]) /
    pmax(trace_lsm[f4_window[-length(f4_window)], "F4_CC"], 1e-8)
)
hcc_incidence_annual <- rate_to_prob(prob_to_rate(hcc_monthly) * 12)

cat("=== Target 1: Annual HCC incidence from F4 ===\n")
cat("  Model predicted:  ", round(hcc_incidence_annual * 100, 2), "% /year\n")
cat("  Literature range: ", target_hcc_incidence$annual_low * 100, "-",
    target_hcc_incidence$annual_high * 100, "% /year\n")
cat("  Source:", target_hcc_incidence$source, "\n\n")

# ---- Target 2: Annual DCC rate from F4 ----
dcc_monthly <- mean(
  diff(trace_lsm[f4_window, "DCC"]) /
    pmax(trace_lsm[f4_window[-length(f4_window)], "F4_CC"], 1e-8)
)
dcc_rate_annual <- rate_to_prob(prob_to_rate(dcc_monthly) * 12)

cat("=== Target 2: Annual DCC rate from F4 ===\n")
cat("  Model predicted:  ", round(dcc_rate_annual * 100, 2), "% /year\n")
cat("  Literature range: ", target_dcc_rate$annual_low * 100, "-",
    target_dcc_rate$annual_high * 100, "% /year\n")
cat("  Source:", target_dcc_rate$source, "\n\n")

# ---- Target 3: Annual mortality from DCC ----
dcc_mort_annual <- rate_to_prob(prob_to_rate(p_prog_month$DCC_Death) * 12)

cat("=== Target 3: Annual mortality from DCC ===\n")
cat("  Model predicted:  ", round(dcc_mort_annual * 100, 1), "% /year\n")
cat("  Literature range: ", target_dcc_mortality$annual_low * 100, "-",
    target_dcc_mortality$annual_high * 100, "% /year\n")
cat("  Source:", target_dcc_mortality$source, "\n\n")

# ---- Validation summary table ----
df_val_summary <- data.frame(
  Target = c("Annual HCC incidence (F4)", "Annual DCC rate (F4)", "Annual mortality from DCC"),
  Literature_Low = c(
    paste0(target_hcc_incidence$annual_low  * 100, "%"),
    paste0(target_dcc_rate$annual_low       * 100, "%"),
    paste0(target_dcc_mortality$annual_low  * 100, "%")
  ),
  Literature_High = c(
    paste0(target_hcc_incidence$annual_high * 100, "%"),
    paste0(target_dcc_rate$annual_high      * 100, "%"),
    paste0(target_dcc_mortality$annual_high * 100, "%")
  ),
  Model_Predicted = c(
    paste0(round(hcc_incidence_annual * 100, 2), "%"),
    paste0(round(dcc_rate_annual      * 100, 2), "%"),
    paste0(round(dcc_mort_annual      * 100, 1), "%")
  ),
  Pass = c(
    ifelse(hcc_incidence_annual >= target_hcc_incidence$annual_low &
           hcc_incidence_annual <= target_hcc_incidence$annual_high, "PASS", "FAIL"),
    ifelse(dcc_rate_annual      >= target_dcc_rate$annual_low &
           dcc_rate_annual      <= target_dcc_rate$annual_high, "PASS", "FAIL"),
    ifelse(dcc_mort_annual      >= target_dcc_mortality$annual_low &
           dcc_mort_annual      <= target_dcc_mortality$annual_high, "PASS", "FAIL")
  ),
  Source = c(target_hcc_incidence$source, target_dcc_rate$source, target_dcc_mortality$source),
  stringsAsFactors = FALSE
)

cat("\n============================================================\n")
cat(" EXTERNAL VALIDATION SUMMARY\n")
cat("============================================================\n")
print(df_val_summary, row.names = FALSE)
cat("============================================================\n\n")

# ---- Plot A: State occupancy over time ----
df_trace_lsm <- as.data.frame(trace_lsm)
df_trace_lsm$Age   <- ages_trace
df_trace_lsm$Cycle <- 0:n_cycles

fibrosis_states <- c("F0","F1","F2","F3","F4_CC","DCC","HCC","Dead")

df_trace_long <- df_trace_lsm %>%
  select(Age, all_of(fibrosis_states)) %>%
  pivot_longer(-Age, names_to = "State", values_to = "Proportion") %>%
  mutate(State = factor(State, levels = fibrosis_states))

state_colors <- c(
  F0    = "#56B4E9", F1    = "#0072B2", F2  = "#009E73",
  F3    = "#E69F00", F4_CC = "#D55E00", DCC = "#CC79A7",
  HCC   = "#F0E442", Dead  = "gray40"
)

ggplot(df_trace_long, aes(x = Age, y = Proportion, fill = State)) +
  geom_area(alpha = 0.85) +
  scale_fill_manual(values = state_colors) +
  scale_x_continuous("Age (years)", breaks = seq(10, 90, by = 10)) +
  scale_y_continuous("Proportion of cohort", labels = percent_format()) +
  labs(title    = "Model-Predicted State Occupancy Over Time (LSM / Natural History)",
       subtitle = "Starting cohort: F2/F3 adolescents at age 12; 80-year horizon") +
  theme_bw(base_size = txtsize) +
  theme(legend.position = "right")

# ---- Plot B: Cumulative mortality under LSM ----
df_dead <- df_trace_lsm %>%
  mutate(pct_dead = Dead) %>%
  select(Age, pct_dead)

ggplot(df_dead, aes(x = Age, y = pct_dead)) +
  geom_line(color = "gray30", linewidth = 1.1) +
  scale_x_continuous("Age (years)", breaks = seq(10, 90, by = 10)) +
  scale_y_continuous("Cumulative proportion dead", labels = percent_format()) +
  labs(title    = "Model-Predicted Cumulative Mortality (LSM cohort)",
       subtitle = "Includes background + disease-specific mortality") +
  theme_bw(base_size = txtsize)

cat("04_base_case.R complete.\n")
