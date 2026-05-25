# 01_model_functions.R
# Requires: source("00_parameters.R")

#############################################################
######################## CYCLE TRANSITION ###################
#############################################################

apply_rr <- function(p, rr) {
  if (rr != 1.0) apply_rr_at_trial_timescale(p, rr) else p
}

build_cycle_P <- function(rr_regress_now, rr_progress_now, age_t,
                          p_prog_month_local = p_prog_month) {
  
  P <- matrix(0, n_states, n_states, dimnames=list(v_states,v_states))
  
### Fibrosis progression — progression RR applies only to F2->F3
  P["F0","F1"]    <- clip01(p_prog_month_local$F0_F1)
  P["F1","F2"]    <- clip01(p_prog_month_local$F1_F2)
  P["F2","F3"]    <- clip01(apply_rr(p_prog_month_local$F2_F3, rr_progress_now))
  P["F3","F4_CC"] <- clip01(p_prog_month_local$F3_F4)                   

### Fibrosis regression — regression RR applies to F2->F1 and F3->F2
  P["F1","F0"]    <- clip01(p_prog_month_local$F1_F0)
  P["F2","F1"]    <- clip01(apply_rr(p_prog_month_local$F2_F1, rr_regress_now))
  P["F3","F2"]    <- clip01(apply_rr(p_prog_month_local$F3_F2, rr_regress_now))
  P["F4_CC","F3"] <- clip01(p_prog_month_local$F4_F3)                               
  
  ### Advanced liver routes 
  P["F3","HCC"]    <- clip01(p_prog_month_local$F3_HCC)
  P["F4_CC","HCC"] <- clip01(p_prog_month_local$F4_HCC)
  P["F4_CC","DCC"] <- clip01(p_prog_month_local$F4_DCC)
  P["DCC","HCC"]   <- clip01(p_prog_month_local$DCC_HCC)
  
  ### Death from DC / HCC / LT / Post-LT
  P["DCC","Dead"]     <- clip01(p_prog_month_local$DCC_Death)
  P["HCC","Dead"]     <- clip01(p_prog_month_local$HCC_Death)
  P["LT_Y1","Dead"]   <- clip01(p_prog_month_local$LT_Death)
  P["LT_Y1_P","Dead"] <- clip01(p_prog_month_local$LT_Death)
  P["Post_LT","Dead"] <- clip01(p_prog_month_local$PostLT_Death)  
  
  ### Liver transplant
  P["F4_CC","LT_Y1"] <- clip01(p_prog_month_local$F4_LT)
  P["DCC","LT_Y1"]   <- clip01(p_prog_month_local$DCC_LT)
  P["HCC","LT_Y1"]   <- clip01(p_prog_month_local$HCC_LT)
  
  ### LT -> Post-LT
  P["LT_Y1","Post_LT"]   <- clip01(p_prog_month_local$LT1_to_PostLT)
  P["LT_Y1_P","Post_LT"] <- clip01(p_prog_month_local$LT1P_to_PostLT)
  P
}

##############################################################
################### BUILD 4D TRANSITION ARRAY ################
##############################################################

build_a_P <- function(rr_reg, rr_prog,
                      p_prog_month_local = p_prog_month,
                      treat_dur_cycles   = NULL,
                      treat_start_cycles = NULL) {

  if(is.null(treat_dur_cycles)){
    treat_dur_cycles <- base_treat_dur_cycles
  }
  if(is.null(treat_start_cycles)){
    treat_start_cycles <- setNames(rep(1L, length(treatments)), treatments)
  }

  a_P <- array(0,
               dim=c(n_states,n_states,n_cycles,length(treatments)),
               dimnames=list(v_states,v_states,1:n_cycles,treatments))

  for(stg in treatments){
    
    dur_stg   <- treat_dur_cycles[stg]
    start_stg <- treat_start_cycles[stg]

    for(t in 1:n_cycles){

      age_t <- age_start + (t-1)*cycle_length

      if(dur_stg > 0 && t >= start_stg && t < (start_stg + dur_stg)){
        rrreg <- rr_reg[stg]
        rrprog<- rr_prog[stg]
      } else {
        rrreg <- 1.0
        rrprog<- 1.0
      }

      P <- build_cycle_P(rrreg, rrprog, age_t, p_prog_month_local)

      ### Add age-specific background mortality 
      p_bg_t <- v_p_bg_month[t]

      for (s in v_states[v_states != "Dead"]) {
        base_p_death   <- P[s, "Dead"]
        p_total_death  <- 1 - (1 - base_p_death) * (1 - p_bg_t)
        P[s, "Dead"]   <- clip01(p_total_death)
      }

      P["Dead",]       <- 0
      P["Dead","Dead"] <- 1

      ### Stay probabilities
      for(s in v_states){
        if(s != "Dead"){
          row_sum <- sum(P[s,])
          if(row_sum >= 1){
            P[s,] <- P[s,] / row_sum
          } else {
            P[s,s] <- 1 - row_sum
          }
        }
      }

      P <- row_normalize(P)
      a_P[,,t,stg] <- P
    }
  }

  a_P
}

#############################################################
########################  Markov Model  #####################
#############################################################

run_markov <- function(P4d_single, v_init){

  m_M <- matrix(NA, nrow=n_cycles+1, ncol=n_states,
                dimnames=list(0:n_cycles, v_states))

  m_M[1,] <- v_init[v_states]

  for(t in 1:n_cycles){
    m_M[t+1,] <- m_M[t,] %*% P4d_single[,,t]
  }

  m_M
}

#############################################################
######################## OUTCOMES ###########################
#############################################################

summarize_outcomes <- function(trace, drug_cost_per_year,
                               util_matrix, bg_cost_cycle,
                               cost_vector,
                               treat_dur_cycles) {

  util_mat <- util_matrix[, v_states]

  # State costs
  v_cost_state_annual <- trace %*% matrix(cost_vector[v_states], ncol=1)
  v_cost_state_cycle  <- v_cost_state_annual * cycle_length

  # Drug cost while on treatment
    on_vec <- rep(0, n_cycles + 1)
    if (treat_dur_cycles > 0) {
    idx_end <- min(treat_dur_cycles, n_cycles + 1)
    on_vec[1:idx_end] <- 1
     }

    treat_year_marker <- ceiling(cumsum(on_vec) / round(1/cycle_length))  
    treat_year_marker[on_vec == 0] <- 0

    years_treated <- length(unique(treat_year_marker[treat_year_marker > 0]))
    total_drug_cost <- drug_cost_per_year * years_treated

    v_cost_drug_cycle <- rep(0, n_cycles + 1)
    if (sum(on_vec) > 0) {
     v_cost_drug_cycle[on_vec == 1] <- total_drug_cost / sum(on_vec)
      }

  # Background cost 
  alive <- 1 - trace[, "Dead"]
  bg_cost_extended <- c(bg_cost_cycle[1], bg_cost_cycle)
  v_bg_cost_cycle <- bg_cost_extended * alive

  # Total cost
  v_cost_total <- v_cost_state_cycle + v_cost_drug_cycle + v_bg_cost_cycle

  # QALYs
  v_qaly_state <- rowSums(trace * util_mat)

  # Life-years (1 - Dead)
  LY <- rowSums(trace[, v_states != "Dead"])

  # Discount + half-cycle
  tot_LY   <- sum(LY * v_wcc * v_dwu) * cycle_length
  tot_QALY <- sum(v_qaly_state * v_wcc * v_dwu) * cycle_length
  tot_cost <- sum(v_cost_total * v_wcc * v_dwc)

  c(LY=tot_LY, QALY=tot_QALY, Cost=tot_cost)
}

#############################################################
################# SUMMARIZE STRATEGIES ######################
#############################################################

summarize_strategies <- function(traces_list, scenario_name,
                                 util_matrix, bg_cost_cycle,
                                 cost_vector, drug_cost_vec,
                                 treat_dur_cycles_vec) {

  res_mat <- sapply(treatments, function(stg){

    summarize_outcomes(
      trace              = traces_list[[stg]],
      drug_cost_per_year = drug_cost_vec[stg],
      util_matrix        = util_matrix,
      bg_cost_cycle      = bg_cost_cycle,
      cost_vector        = cost_vector,
      treat_dur_cycles   = treat_dur_cycles_vec[stg]
    )
  })

  df <- as.data.frame(t(res_mat))
  df$Strategy <- rownames(df)
  df$Scenario <- scenario_name
  df[, c("Scenario","Strategy","Cost","QALY","LY")]
}

#############################################################
################ ICER vs LSM function (OWSA) ################
#############################################################

icer_vs_lsm <- function(res_df){
  lsm  <- res_df[res_df$Strategy == "LSM", ]
  sema <- res_df[res_df$Strategy == "Semaglutide", ]
  if (nrow(lsm) == 0 || nrow(sema) == 0) return(NA_real_)
  dCost <- sema$Cost - lsm$Cost
  dQALY <- sema$QALY - lsm$QALY
  if (dQALY == 0) return(NA_real_)
  dCost / dQALY
}

#############################################################
################ GENERIC ICER FUNCTION ######################
#############################################################

run_owsa_icer <- function(rr_regress_vec      = rr_regress,
                          rr_progress_vec     = rr_progress,
                          util_matrix         = m_util_base,
                          cost_vector         = costs_base,
                          drug_cost_vec       = drug_cost,
                          bg_cost_cycle_vec   = v_background_cost_cycle,
                          p_prog_month_local  = p_prog_month,
                          treat_dur_cycles_vec = base_treat_dur_cycles,
                          treat_start_cycles   = treat_start_immediate,
                          v_init_vec          = v_init) {

  a_P_loc <- build_a_P(rr_regress_vec, rr_progress_vec,
                       p_prog_month_local,
                       treat_dur_cycles   = treat_dur_cycles_vec,
                       treat_start_cycles = treat_start_cycles)

  traces_loc <- lapply(treatments, function(stg) {
    run_markov(a_P_loc[,,,stg], v_init_vec)
  })
  names(traces_loc) <- treatments

  res_loc <- summarize_strategies(
    traces_loc, "SA",
    util_matrix,
    bg_cost_cycle_vec,
    cost_vector,
    drug_cost_vec,
    treat_dur_cycles_vec
  )

  icer_vs_lsm(res_loc)
}

#############################################################
################ run_strategy helper ########################
#############################################################

run_strategy <- function(treat_start, label) {
  aP_sc <- build_a_P(rr_regress, rr_progress,
                     p_prog_month_local = p_prog_month,
                     treat_dur_cycles   = treat_dur_72w_cycles,
                     treat_start_cycles = treat_start)

  traces_sc <- lapply(treatments, function(stg)
    run_markov(aP_sc[,,,stg], v_init))
  names(traces_sc) <- treatments

  list(
    summary = summarize_strategies(traces_sc, label, m_util_base,
                                   v_background_cost_cycle, costs_base,
                                   drug_cost, treat_dur_cycles_vec = treat_dur_72w_cycles),
    traces  = traces_sc
  )
}

cat("01_model_functions.R loaded.\n")

