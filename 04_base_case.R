# 04_base_case.R
# Requires: source("00_parameters.R"), source("00b_le_transitions.R"), source("02_calibration.R"), source("03_validation.R")
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

trace_lsm  <- traces_imm[["LSM"]]
ages_trace <- age_start + (0:n_cycles) * cycle_length
txtsize    <- 13

# ---- Plot A: State occupancy over time ----
df_trace_lsm <- as.data.frame(trace_lsm)
df_trace_lsm$Age   <- ages_trace
df_trace_lsm$Cycle <- 0:n_cycles

fibrosis_states <- c("F0","F1","F2","F3","F4_CC","DCC","HCC","LT_Y1","Post_LT","Dead")

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

row_sums <- rowSums(trace_lsm)
cat("Min row sum:", min(row_sums), "\n")
cat("Max row sum:", max(row_sums), "\n")
cat("Any below 0.999?", any(row_sums < 0.999), "\n")
