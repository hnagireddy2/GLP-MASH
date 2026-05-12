# 00_parameters.R
# All model inputs: libraries, helper functions, costs, utilities,
# transition probabilities, treatment effects, scenarios.
# No model execution happens here.

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

prob_to_rate  <- function(p, t = 1) -(1/t) * log(1 - p)
rate_to_prob  <- function(r, t = 1) 1 - exp(-r * t)
row_normalize <- function(m) sweep(m, 1, rowSums(m), FUN = "/")
clip01        <- function(x) pmin(pmax(x, 0), 0.999)

############################################################
########### Load Mortality + Background Costs ##############
############################################################

mort_cost_df <- read.csv("~/GitHub/GLP-MASH/age_mort_background_costs.csv")

#############################################################
######################## Model Specs ########################
#############################################################

cycle_length <- 1/12         # Monthly
age_start    <- 12
time_horizon <- 80           # 80 years of cycles
age_end      <- age_start + time_horizon
n_cycles     <- time_horizon / cycle_length   # 960 cycles

#############################################################
################## Health States + Initial ##################
#############################################################

v_states <- c(
  "F0","F1","F2","F3","F4_CC","DCC",
  "HCC","LT_Y1","LT_Y1_P","Post_LT","Dead"
)
n_states <- length(v_states)

# Initial F2/F3 distribution based on normalizing NHANES-based F2/F3 distribution 
v_init <- c(F0 = 0,    F1 = 0,    F2 = 0.686, F3 = 0.314,
            F4_CC = 0, DCC = 0,   HCC = 0,
            LT_Y1_P = 0, LT_Y1 = 0, Post_LT = 0, Dead = 0)

#############################################################
############################ Costs ##########################
#############################################################

# Annual state costs 
costs_base <- c(
  F0=8698, F1=8698, F2=8698, F3=10372, F4_CC=42207,
  DCC=195156, HCC=141615,
  LT_Y1_P=185743, LT_Y1=262900,
  Post_LT=49851, Dead=0
)

### LT Complications 

# Vector of probabilities for each LT complication 
lt_comp_probs <- c(
  acr=0.2, biliary_comp=0.105, HAT=0.0545, skin_infection=0.21,
  pneumonia=0.155, bloodstream_inf=0.29, peritonitis=0.0765,
  uti=0.17, cdiff=0.0385, other_infection=0.55,
  VTE=0.031, reoperation=0.125, primary_nonfxn=0.226,
  HVS=0.035, renal_failure=0.1
)

# Lower Bounds Vector 
lt_comp_probs_low <- c(
  acr=0.15, biliary_comp=0.02, HAT=0.019,
  skin_infection=0.16, pneumonia=0.08, bloodstream_inf=0.19,
  peritonitis=0.063, uti=0.16, cdiff=0.027,
  other_infection=0.41, VTE=0.02, reoperation=0.08,
  primary_nonfxn=0.052, HVS=0.01, renal_failure=0.05
)

lt_comp_probs_high <- c(
  acr=0.25, biliary_comp=0.19, HAT=0.09, skin_infection=0.26,
  pneumonia=0.23, bloodstream_inf=0.40, peritonitis=0.09,
  uti=0.18, cdiff=0.05, other_infection=0.69, VTE=0.04,
  reoperation=0.22, primary_nonfxn=0.40, HVS=0.06, renal_failure=0.20
)

# Vector of costs for each LT complication 
lt_comp_cost <- c(
  acr=28950, biliary_comp=54943, HAT=112834, skin_infection=3915,
  pneumonia=80291, bloodstream_inf=102690, peritonitis=119762,
  uti=68730, cdiff=46091, other_infection=68063,
  VTE=53165, reoperation=111674, primary_nonfxn=107031,
  HVS=73838, renal_failure=82524
)

# Lower Bounds Vector
lt_comp_cost_low <- c(
  acr=14476, biliary_comp=27472, HAT=56418, skin_infection=1958,
  pneumonia=40145, bloodstream_inf=51345, peritonitis=59882,
  uti=34366, cdiff=23046, other_infection=34032,
  VTE=26583, reoperation=55838, primary_nonfxn=53516,
  HVS=36919, renal_failure=41262
)

# Upper Bounds Vector
lt_comp_cost_high <- c(
  acr=43425, biliary_comp=82415, HAT=169252, skin_infection=5874,
  pneumonia=120436, bloodstream_inf=154036, peritonitis=179645,
  uti=103095, cdiff=69137, other_infection=102095,
  VTE=79816, reoperation=167512, primary_nonfxn=160547,
  HVS=110756, renal_failure=123785
)

### Expected added cost from LT complications
lt_add_cost_base  <- sum(lt_comp_cost       * lt_comp_probs)
lt_add_cost_low   <- sum(lt_comp_cost_low   * lt_comp_probs_low)
lt_add_cost_high  <- sum(lt_comp_cost_high  * lt_comp_probs_high)

costs_base["LT_Y1"]   <- costs_base["LT_Y1"]   + lt_add_cost_base
costs_base["LT_Y1_P"] <- costs_base["LT_Y1_P"] + lt_add_cost_base

### Drug costs (annual)
cost_lsm        <- 0
cost_sema_base  <- 6829     
cost_sema_low   <- 2940     
cost_sema_high  <- 13658    

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

ages_states <- seq(age_start, age_end, by=cycle_length)
ages_cycles <- ages_states[-length(ages_states)]

## Overall background mortality (annual) -> monthly
## prob_to_rate -> rate_to_prob 
v_p_bg_annual <- approx(mort_cost_df$Age,
                        mort_cost_df$overall_mortality_avg,
                        ages_cycles, rule=2)$y

v_r_bg_annual <- prob_to_rate(v_p_bg_annual)               # hazard per year
v_p_bg_month  <- rate_to_prob(v_r_bg_annual, cycle_length) # monthly prob

## Background costs
v_background_cost_annual <-
  approx(mort_cost_df$Age,
         mort_cost_df$background_cost_2025,
         ages_cycles, rule = 2)$y   

v_background_cost_month <- v_background_cost_annual / 12

#############################################################
########################## Utilities ########################
#############################################################

age_vec <- c(12,25,35,45,55,65,75)

util_age_base <- c(0.919,0.911,0.841,0.816,0.815,0.824,0.811)
qaly_dec_base <- c(
  Healthy=0, F0=0.016, F1=0.016, F2=0.016, F3=0.145, F4_CC=0.145,
  HCC=0.165, DCC=0.155, LT_Y1=0.286, LT_Y1_P=0.286,
  Post_LT=0.036, Dead=1
)

## Linear interpolation across ages 
v_util_age_base <- approx(age_vec, util_age_base, ages_states, rule=2)$y

state_order_for_util <- c(
  "F0","F1","F2","F3","F4_CC","HCC",
  "DCC","LT_Y1","LT_Y1_P","Post_LT","Dead"
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
  F3_HCC        = 0.0034,   
  F4_HCC        = 0.0378,   
  F4_DCC        = 0.0659,   
  DCC_HCC       = 0.0378,   
  DCC_LT        = 0.023,   
  DCC_Death     = 0.1549,     
  HCC_LT        = 0.03,   
  HCC_Death     = 0.3488,     
  F4_LT         = 0.0,      
  DCC_RegressF4 = 0.0,      
  HCC_RegressF4 = 0.0,      
  LT1_to_PostLT = 1.0,      
  LT1P_to_PostLT= 1.0,
  LT_Death      = 0.0909     
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
  LT1_to_PostLT = nonfib_annual$LT1_to_PostLT,
  LT1P_to_PostLT= nonfib_annual$LT1P_to_PostLT,
  LT_Death      = annual_to_month(nonfib_annual$LT_Death),
  PostLT_Death  = rate_to_prob(0.036, cycle_length)
)

#############################################################
###################    TREATMENT EFFECTS   ##################
#############################################################

treatments <- c("LSM","Semaglutide")
all_strat_labels <- c(
  "LSM",
  "Sema 72w (Age 12)",
  "Sema 72w (Age 18)"
)

## Semaglutide progression effect (vs LSM)
rr_sema_progress <- 0.6097

## Semaglutide regression effect (vs LSM) 
rr_sema_regress  <- 1.5625

rr_regress <- c(
  LSM         = 1.0,
  Semaglutide = rr_sema_regress
)

rr_progress <- c(
  LSM         = 1.0,
  Semaglutide = rr_sema_progress
)

##############################################################
########################### SCENARIOS ########################
##############################################################
# Moved here from later in the Rmd so that base_treat_dur_cycles
# is available when build_a_P uses it as a default argument.

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
