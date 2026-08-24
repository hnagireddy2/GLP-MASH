# 05b_elasticity.R
# Requires: source("00_parameters.R"), source("00b_le_transitions.R"),
#           source("01_model_functions.R"), source("02_calibration.R")
#
# Elasticity check: nudge each parameter +/-25% one at a time (all else
# held at base case), report the resulting ICER and an elasticity score.
#
# Unlike OWSA (05_owsa.R), which sweeps each parameter over its own
# literature-derived CI, every parameter here gets the SAME relative
# nudge. That isolates how sensitive the model's ICER is to each input
# structurally, independent of how uncertain that input actually is --
# useful for deciding which parameters are worth a closer OWSA/PSA look.
#
# Elasticity = (% change in ICER from -25% to +25%) / (% change in the
# parameter, i.e. 50%). Larger |elasticity| = ICER moves more for a given
# relative change in that parameter.

pct_shift   <- 0.25
base_icer_v <- run_owsa_icer()

elasticity_row <- function(name, base_val, run_fn) {
  icer_lo <- run_fn(base_val * (1 - pct_shift))
  icer_hi <- run_fn(base_val * (1 + pct_shift))
  data.frame(
    Parameter    = name,
    Base_Value   = signif(base_val, 4),
    ICER_minus25 = round(icer_lo),
    ICER_plus25  = round(icer_hi),
    Elasticity   = round(((icer_hi - icer_lo) / base_icer_v) / (2 * pct_shift), 3)
  )
}

rows <- list()

###### Treatment effects (RR) ######

rows[["RR_regress"]] <- elasticity_row("Semaglutide RR: regression", RR_regress, function(x) {
  rv <- rr_regress; rv["Semaglutide"] <- x
  run_owsa_icer(rr_regress_vec = rv)
})
rows[["RR_progress"]] <- elasticity_row("Semaglutide RR: progression", RR_progress, function(x) {
  rv <- rr_progress; rv["Semaglutide"] <- x
  run_owsa_icer(rr_progress_vec = rv)
})

###### Fibrosis transitions + advanced-disease hazards (monthly probs) ######

trans_names <- c("F0_F1", "F1_F0", "F1_F2", "F2_F1", "F2_F3", "F3_F2", "F3_F4", "F4_F3",
                 "F3_HCC", "F4_HCC", "F4_DCC", "DCC_HCC", "DCC_LT", "DCC_Death",
                 "HCC_LT", "HCC_Death", "LT_Death", "PostLT_Death")

for (nm in trans_names) {
  base_p <- p_prog_month[[nm]]
  rows[[paste0("trans_", nm)]] <- elasticity_row(paste0("Transition: ", nm), base_p, function(x) {
    pm <- p_prog_month; pm[[nm]] <- x
    run_owsa_icer(p_prog_month_local = pm)
  })
}

###### Costs ######
## F0/F1/F2 share one literature-sourced cost estimate, so they're nudged
## together as a single parameter rather than three identical rows.

cost_groups <- list(F0_F2 = c("F0", "F1", "F2"), F3 = "F3", F4_CC = "F4_CC",
                    DCC = "DCC", HCC = "HCC", LT = "LT", Post_LT = "Post_LT")

for (grp in names(cost_groups)) {
  states <- cost_groups[[grp]]
  base_c <- costs_base[states[1]]
  rows[[paste0("cost_", grp)]] <- elasticity_row(paste0("Cost: ", grp), base_c, function(x) {
    cv <- costs_base; cv[states] <- x
    run_owsa_icer(cost_vector = cv)
  })
}

rows[["drug_sema"]] <- elasticity_row("Semaglutide annual drug cost", cost_sema_base, function(x) {
  dv <- drug_cost; dv["Semaglutide"] <- x
  run_owsa_icer(drug_cost_vec = dv)
})

###### Health state utility decrements ######

for (grp in names(cost_groups)) {          # same state groupings as costs
  states <- cost_groups[[grp]]
  base_q <- qaly_dec_base[states[1]]
  rows[[paste0("qdec_", grp)]] <- elasticity_row(paste0("Util decrement: ", grp), base_q, function(x) {
    qd <- qaly_dec_base; qd[states] <- x
    um <- build_util_matrix(v_util_age_base, qd)
    run_owsa_icer(util_matrix = um)
  })
}

###### Assemble + report ######

df_elasticity <- do.call(rbind, rows)
rownames(df_elasticity) <- NULL
df_elasticity <- df_elasticity[order(-abs(df_elasticity$Elasticity)), ]

cat("\n=================================================================\n")
cat(" ELASTICITY CHECK: +/-25% one-at-a-time, all parameters\n")
cat(" Base-case ICER: $", format(round(base_icer_v), big.mark = ","), "/QALY\n", sep = "")
cat("=================================================================\n\n")
print(df_elasticity, row.names = FALSE)

cat("\n05b_elasticity.R complete.\n")
