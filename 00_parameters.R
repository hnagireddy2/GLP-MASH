# 00_parameters.R  
# All model inputs: libraries, helper functions, costs, utilities,
# transition probabilities, treatment effects, scenarios.
# No model execution happens here.

setwd("~/GitHub/GLP-MASH") 

rm(list = ls())

## Libraries
library(dplyr)
library(tidyr)
library(ggplot2) 
library(ggrepel)
library(scales)      
library(dampack)
library(truncnorm)
library(matrixStats) 
library(reshape2)
library(mgcv)
library(gtools)

## Colorblind-safe palette
okabe_ito_strat <- c(
  LSM         = "#0072B2",
  Semaglutide = "#009E73"
)

############################################################
######################   Functions  ########################
############################################################

prob_to_rate <- function(p, t = 1) {
  -log(1 - p) / t
}
rate_to_prob <- function(r, t = 1) {
  1 - exp(-r * t)
}
row_normalize <- function(m) sweep(m, 1, rowSums(m), FUN = "/")

clip01        <- function(x) pmin(pmax(x, 0), 0.999)

apply_rr_at_trial_timescale <- function(p_cycle_baseline, RR,
                                        trial_weeks = 72,
                                        weeks_per_year = 52.1429,
                                        cycle_length_local = cycle_length) {
  # Step 1: Cycle probability -> cycle rate
  r_cycle <- -log(1 - p_cycle_baseline)
  
  # Step 2: Scale cycle rate to 72-week rate
  cycles_in_trial <- trial_weeks / (weeks_per_year * cycle_length_local)
  r_72wk <- r_cycle * cycles_in_trial
  
  # Step 3: Convert to 72-week probability
  p_72wk <- 1 - exp(-r_72wk)
  
  # Step 4: Apply RR at the 72-week timescale
  p_72wk_treated <- p_72wk * RR
  
  # Safety check
  if (p_72wk_treated >= 1) {
    warning("Treated 72-week probability >= 1. Clamping to 0.999.")
    p_72wk_treated <- 0.999
  }
  
  # Step 5: Convert treated 72-week probability back to rate
  r_72wk_treated <- -log(1 - p_72wk_treated)
  
  # Step 6: Scale treated rate back to cycle length
  r_cycle_treated <- r_72wk_treated / cycles_in_trial
  
  # Step 7: Convert to cycle probability
  p_cycle_treated <- 1 - exp(-r_cycle_treated)
  
  return(p_cycle_treated)
}

############################################################
########### Load Mortality + Background Costs ##############
############################################################

mort_cost_df <- read.csv("age_mort_background_costs.csv")

#############################################################
######################## Model Specs ########################
#############################################################

cycle_length <- 4 / 52.1429  # 4-week cycles 
age_start    <- 12
age_max      <- max(mort_cost_df$age)
time_horizon <- age_max - age_start
n_cycles     <- ceiling(time_horizon / cycle_length)

#############################################################
################## Health States + Initial ##################
#############################################################

v_states <- c("F0","F1","F2","F3","F4_CC","DCC","HCC","LT_Y1","Post_LT","Dead")

n_states <- length(v_states)

treatments <- c("LSM", "Semaglutide")

# Initial F2/F3 distribution based on normalizing NHANES-based F2/F3 distribution 
v_init <- c(F0 = 0,    F1 = 0,    F2 = 0.686, F3 = 0.314,
            F4_CC = 0, DCC = 0,   HCC = 0, LT_Y1 = 0, Post_LT = 0, Dead = 0)

#############################################################
############################ Costs ##########################
#############################################################

# Annual state costs 
costs_base <- c(
  F0=8698, F1=8698, F2=8698, F3=10372, F4_CC=42207,
  DCC=195156, HCC=141615, LT_Y1=262900,
  Post_LT=2344, Dead=0
)

costs_cycle <- costs_base * cycle_length 

### Drug costs (annual)
cost_lsm        <- 0
cost_sema_base  <- 14072     
cost_sema_low   <- 4188    
cost_sema_high  <- 16188    

drug_cost <- c(
  LSM         = cost_lsm,
  Semaglutide = cost_sema_base   
)

#############################################################
############### Discounting + Half-cycle ####################
#############################################################

disc_cost <- 0.03
disc_qaly <- 0.03

v_dwc <- 1/(1+disc_cost)^((0:n_cycles)*cycle_length)
v_dwu <- 1/(1+disc_qaly)^((0:n_cycles)*cycle_length)
v_wcc <- c(0.5, rep(1, n_cycles-1), 0.5)

#############################################################
######## Age-specific Mortality + Background Costs ##########
#############################################################

ages_cycles <- age_start + (0:(n_cycles - 1)) * cycle_length   # length n_cycles
ages_states <- age_start + (0:n_cycles)       * cycle_length   # length n_cycles + 1

## Overall background mortality (annual) -> monthly
## prob_to_rate -> rate_to_prob 
v_p_bg_annual <- approx(mort_cost_df$age,
                        mort_cost_df$background_mortality_prob,
                        ages_cycles, rule=2)$y

v_r_bg_annual <- prob_to_rate(v_p_bg_annual)               # hazard per year
v_p_bg_month  <- rate_to_prob(v_r_bg_annual, cycle_length) # monthly prob

## Background costs
v_background_cost_annual <-
  approx(mort_cost_df$age,
         mort_cost_df$mean_background_cost_2025,
         ages_cycles, rule = 2)$y   

v_background_cost_cycle <- v_background_cost_annual * cycle_length

#############################################################
########################## Utilities ########################
#############################################################

age_vec <- c(12,25,35,45,55,65,75)

util_age_base <- c(0.919,0.911,0.841,0.816,0.815,0.824,0.811)
qaly_dec_base <- c(
  Healthy=0, F0=0.016, F1=0.016, F2=0.016, F3=0.145, F4_CC=0.145,
  HCC=0.165, DCC=0.155, LT_Y1=0.286,
  Post_LT=0.036, Dead=1
)

## Linear interpolation across ages 
v_util_age_base <- approx(age_vec, util_age_base, ages_states, rule=2)$y

state_order_for_util <- c(
  "F0","F1","F2","F3","F4_CC","HCC",
  "DCC","LT_Y1","Post_LT","Dead"
)

#Builds age x state utility matrix (961 x 11) 
build_util_matrix <- function(v_util_age, qdec){
  m <- sapply(state_order_for_util, function(st) {
    if (st=="Dead") rep(0,length(v_util_age))
    else pmax(v_util_age - qdec[st], 0)
  })
  colnames(m) <- state_order_for_util
  m
}

#Call function and store results 
m_util_base <- build_util_matrix(v_util_age_base, qaly_dec_base)

#############################################################
################### TRANSITION PROBABILITIES ################
#############################################################

# Function for annual probability -> monthly probability
annual_to_month <- function(p_annual) {
  rate_to_prob(prob_to_rate(p_annual), cycle_length)
}

# ---- Fixed annual transition probabilities ----
nonfib_annual <- list(
  F3_HCC        = 0.0034,   # Kim 2025 Suppl. Table 1 (0.34%) -> Sanyal 2021
  F4_HCC        = 0.0378,   # Kim 2025 Suppl. Table 1 (3.78%) -> Orci 2022
  F4_DCC        = 0.0659,   # Kim 2025 Suppl. Table 1 (6.59%) -> Younossi 2020
  DCC_HCC       = 0.0378,   # Kim 2025 Suppl. Table 1 (3.78%) -> Orci 2022
  DCC_LT        = 0.023,    # Rustgi 2022 Table 1 (DCC -> LT)
  DCC_Death     = 0.1549,   # Rustgi 2022 Table 1 (DCC -> Death)
  HCC_LT        = 0.03,     # Rustgi 2022 Table 1 (HCC -> LT)
  HCC_Death     = 0.1305,   # Kim 2025 Suppl. Table 1 (13.05%) 
  F4_LT         = 0.0,      # structural
  DCC_RegressF4 = 0.0,      # structural
  HCC_RegressF4 = 0.0,      # structural
  LT_Y1_Post_LT = 1.0,       # structural (deterministic)
  LT_Y1_Death    = 0.0909,     # Rustgi 2022 Table 1, liver-related mortality only (LRM)
  Post_LT_Death      = 0.0909   # Rustgi 2022 Table 1, liver-related mortality only (LRM)
)

# Convert to monthly probabilities (except the deterministic post-LT ones)
p_prog_month <- list(
  
  # ---- Fibrosis transitions: pre-calibration placeholder zeros ----
  F0_F1 = 0, F1_F0 = 0, F1_F2 = 0, F2_F1 = 0,
  F2_F3 = 0, F3_F2 = 0, F3_F4 = 0, F4_F3 = 0,

  # ---- Fixed transitions ----
  F3_HCC        = annual_to_month(nonfib_annual$F3_HCC),
  F4_HCC        = annual_to_month(nonfib_annual$F4_HCC),
  F4_DCC        = annual_to_month(nonfib_annual$F4_DCC),
  DCC_HCC       = annual_to_month(nonfib_annual$DCC_HCC),
  DCC_LT        = annual_to_month(nonfib_annual$DCC_LT),
  DCC_Death     = annual_to_month(nonfib_annual$DCC_Death),
  HCC_LT        = annual_to_month(nonfib_annual$HCC_LT),
  HCC_Death     = annual_to_month(nonfib_annual$HCC_Death),
  F4_LT         = annual_to_month(nonfib_annual$F4_LT),
  DCC_RegressF4 = annual_to_month(nonfib_annual$DCC_RegressF4),
  HCC_RegressF4 = annual_to_month(nonfib_annual$HCC_RegressF4),
  LT_Y1_Post_LT   = nonfib_annual$LT_Y1_Post_LT,
  LT_Y1_Death    = annual_to_month(nonfib_annual$LT_Y1_Death),
  Post_LT_Death  = annual_to_month(nonfib_annual$Post_LT_Death)
)

#############################################################
###################    TREATMENT EFFECTS   ##################
#############################################################

all_strat_labels <- c(
  "LSM",
  "Sema 72w (Age 12)",
  "Sema 72w (Age 18)"
)

#############################################################
## RR Derivation from ESSENCE Trial (Sanyal et al. 2025)   ##
#############################################################

# --- Regression RR (fibrosis improvement, no worsening of steatohepatitis) ---
# Source: Supplement page 27, composite table
a_reg <- 197   # sema responders (regressed)
b_reg <- 337   # sema non-responders
c_reg <- 59    # placebo responders (regressed)
d_reg <- 207   # placebo non-responders

N_tx_reg      <- a_reg + b_reg   # 534
N_placebo_reg <- c_reg + d_reg   # 266

p_tx_reg      <- a_reg / N_tx_reg        # 0.3690
p_placebo_reg <- c_reg / N_placebo_reg   # 0.2218

RR_regress <- p_tx_reg / p_placebo_reg   # 1.6637

# SE on the log scale
SE_ln_RR_regress <- sqrt(1/a_reg - 1/N_tx_reg + 1/c_reg - 1/N_placebo_reg)
ln_RR_regress    <- log(RR_regress)

# 95% CI
RR_reg_lo <- exp(ln_RR_regress - 1.96 * SE_ln_RR_regress)  # 1.295
RR_reg_hi <- exp(ln_RR_regress + 1.96 * SE_ln_RR_regress)  # 2.139

cat("=== Regression RR ===\n")
cat("RR:", round(RR_regress, 4), "\n")
cat("SE(ln RR):", round(SE_ln_RR_regress, 4), "\n")
cat("95% CI:", round(RR_reg_lo, 3), "-", round(RR_reg_hi, 3), "\n\n")

# --- Progression RR (F2 worsening) ---
# Source: Table S6 / Figure S13 (F2 subgroup only)
a_prog <- 35   # sema F2 patients who progressed
b_prog <- 134  # sema F2 patients who did not progress
c_prog <- 29   # placebo F2 patients who progressed
d_prog <- 52   # placebo F2 patients who did not progress

N_tx_prog      <- a_prog + b_prog   # 169
N_placebo_prog <- c_prog + d_prog   # 81

p_tx_prog      <- a_prog / N_tx_prog        # 0.2071
p_placebo_prog <- c_prog / N_placebo_prog   # 0.3580

RR_progress <- p_tx_prog / p_placebo_prog   # 0.5785

# SE on the log scale
SE_ln_RR_progress <- sqrt(1/a_prog - 1/N_tx_prog + 1/c_prog - 1/N_placebo_prog)
ln_RR_progress    <- log(RR_progress)

# 95% CI
RR_prog_lo <- exp(ln_RR_progress - 1.96 * SE_ln_RR_progress)  # 0.382
RR_prog_hi <- exp(ln_RR_progress + 1.96 * SE_ln_RR_progress)  # 0.876

cat("=== Progression RR ===\n")
cat("RR:", round(RR_progress, 4), "\n")
cat("SE(ln RR):", round(SE_ln_RR_progress, 4), "\n")
cat("95% CI:", round(RR_prog_lo, 3), "-", round(RR_prog_hi, 3), "\n\n")

#############################################################
## Treatment Effect RR Vectors (used by build_a_P)         ##
#############################################################

rr_regress <- c(
  LSM         = 1.0,
  Semaglutide = RR_regress
)

rr_progress <- c(
  LSM         = 1.0,
  Semaglutide = RR_progress
)

##############################################################
########################### SCENARIOS ########################
##############################################################

# Only one duration: 72 weeks (ESSENCE trial duration)
treat_dur_72w_years <- 72/52   # ~1.385 years

treat_dur_72w_cycles <- c(LSM = 0L,
                          Semaglutide = round(treat_dur_72w_years / cycle_length))
# ~17 monthly cycles

base_treat_dur_cycles <- treat_dur_72w_cycles

# Delay from age 12 to age 18 (6 years = 72 cycles)
delay_to_18_cycles <- round((18 - 12) / cycle_length)

treat_start_immediate <- c(LSM = 1L, Semaglutide = 1L)
treat_start_age18     <- c(LSM = 1L, Semaglutide = delay_to_18_cycles + 1L)

cat("00_parameters.R loaded.\n")
