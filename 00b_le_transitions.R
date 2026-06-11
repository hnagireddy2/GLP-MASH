# 00b_le_transitions.R
# Requires: source("00_parameters.R")
# ---------------------------------------------------------------------------
# Reproducible derivation of fibrosis transition probabilities from the
# REPORTED pooled incidence rates in Le et al. 2023, Clin Gastroenterol
# Hepatol 2023;21:1154-1168.
#
# Source tables (NAFLD rows = all-locations pooled estimate):
#   Table 2 -> RCTs              (cases/100 PY, 95% CI, by baseline stage)
#   Table 3 -> Observational     (cases/100 PY, 95% CI, by baseline stage)
# Both pooled with a random-effects Poisson model. We use cases/100 PY (one
# event = a patient advancing/regressing >= 1 stage) to
# match the single-step structure of a monthly Markov cycle.
#
# ---------------------------------------------------------------------------

## ---- Conversion: cases per 100 person-years -> annual probability ----------
# Rate r (per 100 PY) implies annual event probability 1 - exp(-r/100).
# Point estimate AND each 95% CI bound are converted separately (nonlinear).
rate100_to_p <- function(r_per_100py) 1 - exp(-r_per_100py / 100)

## ---- Le et al. Table 2: RCT pooled cases/100 PY (NAFLD) ---------------------
# Progression rows: F0,F1,F2,F3 ; Regression rows: F1,F2,F3,F4.
le_trial_rates <- data.frame(
  transition = c("F0_F1", "F1_F2", "F2_F3", "F3_F4",
                 "F1_F0", "F2_F1", "F3_F2", "F4_F3"),
  rate = c(19.1, 25.4, 18.5,  8.0,
           12.7, 21.2, 17.2,  9.0),
  lo   = c( 9.4, 13.4, 12.7,  4.5,
            7.5, 14.9, 12.9,  4.6),
  hi   = c(39.1, 48.1, 26.9, 14.0,
           21.5, 30.1, 22.9, 17.9),
  source = "Table 2 (RCT, NAFLD)",
  stringsAsFactors = FALSE
)

## ---- Le et al. Table 3: Observational pooled cases/100 PY (NAFLD) -----------
le_obs_rates <- data.frame(
  transition = c("F0_F1", "F1_F2", "F2_F3", "F3_F4",
                 "F1_F0", "F2_F1", "F3_F2", "F4_F3"),
  rate = c(6.5, 6.9, 5.8, 4.5,
           2.5, 4.6, 6.0, 4.7),
  lo   = c(4.6, 4.8, 4.6, 3.3,
           1.4, 2.4, 3.7, 2.1),
  hi   = c(9.0, 10.0, 7.3, 6.0,
           4.6, 8.8, 9.6, 10.8),
  source = "Table 3 (Observational, NAFLD)",
  stringsAsFactors = FALSE
)

## ---- Build named annual-probability vectors --------------------------------
.named <- function(df, col) setNames(rate100_to_p(df[[col]]), df$transition)

le_trial_low  <- .named(le_trial_rates, "lo")
le_trial_mean <- .named(le_trial_rates, "rate")
le_trial_high <- .named(le_trial_rates, "hi")

le_obs_low  <- .named(le_obs_rates, "lo")
le_obs_mean <- .named(le_obs_rates, "rate")
le_obs_high <- .named(le_obs_rates, "hi")

## ---- Report ----------------------------------------------------------------
fmt <- function(df, mean_v) data.frame(
  transition    = df$transition,
  rate_per100PY = df$rate,
  CI_per100PY   = sprintf("%.1f-%.1f", df$lo, df$hi),
  annual_prob   = round(mean_v[df$transition], 4),
  row.names = NULL
)

cat("\n=== RCT (trial) set: Le Table 2 cases/100 PY -> annual probability ===\n")
print(fmt(le_trial_rates, le_trial_mean), row.names = FALSE)
cat("\n=== Observational set: Le Table 3 cases/100 PY -> annual probability ===\n")
print(fmt(le_obs_rates, le_obs_mean), row.names = FALSE)
cat("\n00b_le_transitions.R complete:",
    "le_trial_{low,mean,high} and le_obs_{low,mean,high} defined.\n")