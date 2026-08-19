# 08_voi.R
# Requires: source("06_psa.R"), source("07_psa_postprocess.R")
# Value of Information: EVPI, EVPPI, EVSI, ENBS
  
# ============================================================
# Value of Information (VOI) Analysis — Strategies
# ============================================================

# ---- 0. Setup ----------------------------------------------
source("GA_functions.R")
wtp        <- 150000
txtsize    <- 14

# Use all-strategy PSA outputs
df_c_voi   <- df_c_all
df_e_voi   <- df_e_all
strategies <- all_strat_labels
n_sim      <- nrow(df_c_voi)
n_strat    <- length(strategies)

# ---- 1. Net Monetary Benefit (NMB) -------------------------
df_nmb <- as.data.frame(
  sapply(strategies, function(s) wtp * df_e_voi[[s]] - df_c_voi[[s]])
)
colnames(df_nmb) <- strategies

cat("===== NMB Summary (WTP = $150,000/QALY) =====\n")
print(colMeans(df_nmb))

## -- Plot NMB distributions --
df_nmb_long <- reshape2::melt(df_nmb,
                              variable.name = "Strategy",
                              value.name    = "NMB")
ggplot(df_nmb_long, aes(x = NMB / 1000)) +
  geom_histogram(aes(y = after_stat(density)), col = "black", fill = "gray80", bins = 60) +
  geom_density(color = "#E69F00", linewidth = 1) +
  facet_wrap(~ Strategy, scales = "free_y", ncol = 3) +
  labs(x = "Net Monetary Benefit (thousands $)",
       title = "PSA: NMB Distributions — All Strategies") +
  scale_x_continuous(n.breaks = 5) +
  scale_y_continuous("", breaks = NULL) +
  theme_bw(base_size = txtsize)

## -- Incremental NMB vs LSM --
lsm_col <- df_nmb[["LSM"]]
df_inmb <- as.data.frame(
  sapply(strategies[strategies != "LSM"], function(s) df_nmb[[s]] - lsm_col)
)
df_inmb_long <- reshape2::melt(df_inmb,
                               variable.name = "Comparison",
                               value.name    = "INMB")
ggplot(df_inmb_long, aes(x = INMB / 1000)) +
  geom_histogram(aes(y = after_stat(density)), col = "black", fill = "gray80", bins = 60) +
  geom_density(color = "#E69F00", linewidth = 1) +
  geom_vline(xintercept = 0, col = "#0072B2", linewidth = 1.2, linetype = "dashed") +
  facet_wrap(~ Comparison, scales = "free_y", ncol = 4) +
  labs(x = "Incremental NMB vs LSM (thousands $)",
       title = "PSA: Incremental NMB — All Strategies vs LSM") +
  scale_y_continuous("", breaks = NULL) +
  theme_bw(base_size = txtsize)

# ---- 2. Loss Matrix ----------------------------------------
d_star <- which.max(colMeans(df_nmb))
cat("\nOptimal strategy (highest mean NMB):", strategies[d_star], "\n")

# Loss = NMB of optimal - NMB of each strategy (per simulation)
m_loss <- as.matrix(df_nmb - df_nmb[, d_star])
colnames(m_loss) <- strategies

# ---- 3. EVPI -----------------------------------------------
v_max_loss_i <- rowMaxs(m_loss)
evpi         <- mean(v_max_loss_i)
cat("\n===== EVPI =====\n")
cat("EVPI (per patient): $", format(round(evpi), big.mark = ","), "\n")

## EVPI vs WTP
v_wtp  <- seq(0, 600000, by = 10000)
v_evpi <- numeric(length(v_wtp))
for (k in seq_along(v_wtp)) {
  nmb_k   <- as.data.frame(
    sapply(strategies, function(s) v_wtp[k] * df_e_voi[[s]] - df_c_voi[[s]])
  )
  d_k       <- which.max(colMeans(nmb_k))
  loss_k    <- as.matrix(nmb_k - nmb_k[, d_k])
  v_evpi[k] <- mean(rowMaxs(loss_k))
}

df_evpi_wtp <- data.frame(WTP = v_wtp, EVPI = v_evpi)
ggplot(df_evpi_wtp, aes(x = WTP / 1000, y = EVPI / 1000)) +
  geom_line(linewidth = 1, color = "#0072B2") +
  geom_vline(xintercept = wtp / 1000, linetype = "dashed", color = "gray40") +
  annotate("text", x = wtp / 1000 + 8, y = max(v_evpi) * 0.9 / 1000,
           label = paste0("WTP = $", wtp / 1000, "K"), hjust = 0, size = 4) +
  labs(x = "Willingness-to-Pay (thousands $/QALY)",
       y = "EVPI (thousands $)",
       title = "EVPI Across WTP Thresholds — Strategies") +
  scale_x_continuous(labels = comma) +
  scale_y_continuous(labels = comma) +
  theme_bw(base_size = txtsize)

# ---- 4. EVPPI (individual parameters) ----------------------
v_names_params <- c(
  "rr_sema_regress", "rr_sema_progress",
  "h_F0_F1", "h_F1_F2", "h_F2_F3", "h_F3_F4", "h_F4_DCC",
  "h_F1_F0", "h_F2_F1", "h_F3_F2", "h_F4_F3",
  "h_DCC_Death", "h_HCC_Death", "h_LT_Death", "h_PostLT_Death",
  "cost_F0_F2", "cost_F3", "cost_F4_CC", "cost_DCC", "cost_HCC",
  "qdec_F0_F2", "qdec_F3", "qdec_F4_CC", "qdec_DCC", "qdec_HCC", "qdec_PostLT"
)

param_labels_voi <- c(
  rr_sema_regress  = "Sema RR: Regression",
  rr_sema_progress = "Sema RR: Progression",
  h_F0_F1 = "Hazard F0\u2192F1", h_F1_F2 = "Hazard F1\u2192F2",
  h_F2_F3 = "Hazard F2\u2192F3", h_F3_F4 = "Hazard F3\u2192F4",
  h_F4_DCC = "Hazard F4\u2192DCC",
  h_F1_F0 = "Hazard F1\u2192F0 (reg)", h_F2_F1 = "Hazard F2\u2192F1 (reg)",
  h_F3_F2 = "Hazard F3\u2192F2 (reg)", h_F4_F3 = "Hazard F4\u2192F3 (reg)",
  h_DCC_Death = "DCC\u2192Death", h_HCC_Death = "HCC\u2192Death",
  h_LT_Death = "LT\u2192Death (yr 1)", h_PostLT_Death = "Post-LT\u2192Death",
  cost_F0_F2 = "Cost: F0\u2013F2", cost_F3 = "Cost: F3",
  cost_F4_CC = "Cost: F4/CC", cost_DCC = "Cost: DCC", cost_HCC = "Cost: HCC",
  qdec_F0_F2 = "Util dec: F0\u2013F2", qdec_F3 = "Util dec: F3",
  qdec_F4_CC = "Util dec: F4/CC", qdec_DCC = "Util dec: DCC",
  qdec_HCC = "Util dec: HCC", qdec_PostLT = "Util dec: Post-LT"
)

# Use df_psa_input_all since this VOI is based on the all-strategy PSA
df_params <- df_psa_input_all[, v_names_params]
n_params  <- ncol(df_params)

cat("\n===== Computing EVPPI for individual parameters =====\n")

v_evppi  <- numeric(n_params)
lmm_list <- vector("list", n_params)   # one list of n_strat GAMs per parameter

for (p in seq_len(n_params)) {
  cat("  EVPPI:", v_names_params[p], "\n")
  # Fit one GAM per strategy
  gams_p <- lapply(seq_len(n_strat), function(j)
    gam(m_loss[, j] ~ s(df_params[, p]))
  )
  lmm_list[[p]] <- gams_p
  m_Lhat        <- do.call(cbind, lapply(gams_p, `[[`, "fitted.values"))
  v_evppi[p]    <- mean(rowMaxs(m_Lhat))
}

df_evppi <- data.frame(
  Parameter = factor(
    unname(param_labels_voi[v_names_params]),
    levels = unname(param_labels_voi[v_names_params])[order(v_evppi, decreasing = TRUE)]
  ),
  EVPPI = v_evppi
)

cat("\n===== EVPPI Results =====\n")
print(df_evppi[order(df_evppi$EVPPI, decreasing = TRUE), ])

ggplot(df_evppi, aes(x = Parameter, y = EVPPI)) +
  geom_bar(stat = "identity", fill = "#0072B2") +
  geom_hline(yintercept = evpi, linetype = "dashed", color = "red", linewidth = 0.9) +
  annotate("text", x = n_params - 1, y = evpi * 1.05,
           label = "EVPI", color = "red", size = 4) +
  coord_flip() +
  labs(x = "", y = "EVPPI ($ per patient)",
       title = "EVPPI — Strategies",
       subtitle = paste0("WTP = $", format(wtp, big.mark = ","),
                         "/QALY | EVPI = $", format(round(evpi), big.mark = ","))) +
  scale_y_continuous(labels = comma, n.breaks = 6) +
  theme_bw(base_size = txtsize)

# ---- 5. EVPPI for Grouped Parameters -----------------------
param_groups <- list(
  "Sema Treatment Effects" = c("rr_sema_regress", "rr_sema_progress"),
  "Fibrosis Progression"   = c("h_F0_F1","h_F1_F2","h_F2_F3","h_F3_F4","h_F4_DCC"),
  "Fibrosis Regression"    = c("h_F1_F0","h_F2_F1","h_F3_F2","h_F4_F3"),
  "Mortality Hazards"      = c("h_DCC_Death","h_HCC_Death","h_LT_Death","h_PostLT_Death"),
  "State Costs"            = c("cost_F0_F2","cost_F3","cost_F4_CC","cost_DCC","cost_HCC"),
  "Health State Utilities" = c("qdec_F0_F2","qdec_F3","qdec_F4_CC",
                                "qdec_DCC","qdec_HCC","qdec_PostLT")
)

cat("\n===== Computing EVPPI for parameter groups =====\n")
v_evppi_groups <- numeric(length(param_groups))
names(v_evppi_groups) <- names(param_groups)

for (g in seq_along(param_groups)) {
  grp_name   <- names(param_groups)[g]
  grp_params <- param_groups[[g]]
  grp_idx    <- which(v_names_params %in% grp_params)
  cat("  Group:", grp_name, "(", length(grp_params), "params)\n")

  smooth_terms  <- paste0("s(df_params[,", grp_idx, "])", collapse = " + ")
  if (length(grp_idx) > 1) {
    interact_term <- paste0("ti(df_params[,", grp_idx[1], "], df_params[,",
                            grp_idx[2], "])")
    f_str <- paste("~", smooth_terms, "+", interact_term)
  } else {
    f_str <- paste("~", smooth_terms)
  }

  gams_g <- lapply(seq_len(n_strat), function(j)
    gam(as.formula(paste("m_loss[,", j, "]", f_str)))
  )
  m_Lhat_grp      <- do.call(cbind, lapply(gams_g, `[[`, "fitted.values"))
  v_evppi_groups[g] <- mean(rowMaxs(m_Lhat_grp))
}

df_evppi_groups <- data.frame(
  Group = factor(names(v_evppi_groups),
                 levels = names(v_evppi_groups)[order(v_evppi_groups, decreasing = TRUE)]),
  EVPPI = v_evppi_groups
)

cat("\n===== EVPPI by Parameter Group =====\n")
print(df_evppi_groups[order(df_evppi_groups$EVPPI, decreasing = TRUE), ])

ggplot(df_evppi_groups, aes(x = Group, y = EVPPI)) +
  geom_bar(stat = "identity", fill = "#009E73") +
  geom_hline(yintercept = evpi, linetype = "dashed", color = "red", linewidth = 0.9) +
  annotate("text", x = 1, y = evpi * 1.05, label = "EVPI", color = "red", size = 4) +
  coord_flip() +
  labs(x = "", y = "EVPPI ($ per patient)",
       title = "EVPPI by Parameter Group — Strategies",
       subtitle = paste0("WTP = $", format(wtp, big.mark = ","), "/QALY")) +
  scale_y_continuous(labels = comma, n.breaks = 6) +
  theme_bw(base_size = txtsize)

# ---- 6. EVSI -----------------------------------------------
n0_defaults <- c(
  rr_sema_regress = 400, rr_sema_progress = 400,
  h_F0_F1 = 2000, h_F1_F2 = 2000, h_F2_F3 = 2000,
  h_F3_F4 = 2000, h_F4_DCC = 1500,
  h_F1_F0 = 2000, h_F2_F1 = 2000, h_F3_F2 = 2000, h_F4_F3 = 1500,
  h_DCC_Death = 800, h_HCC_Death = 800,
  h_LT_Death = 500, h_PostLT_Death = 500,
  cost_F0_F2 = 300, cost_F3 = 300, cost_F4_CC = 300,
  cost_DCC = 300, cost_HCC = 300,
  qdec_F0_F2 = 500, qdec_F3 = 500, qdec_F4_CC = 500,
  qdec_DCC = 300, qdec_HCC = 300, qdec_PostLT = 200
)

v_n <- c(0, 10, 25, 50, 100, 200, 300, 500, 750, 1000, 1500, 2000, 3000, 5000)
n_samples <- length(v_n)

cat("\n===== Computing EVSI for individual parameters =====\n")

df_evsi <- data.frame(N = v_n,
                      matrix(0, nrow = n_samples, ncol = n_params))
colnames(df_evsi)[-1] <- v_names_params

for (p in seq_len(n_params)) {
  cat("  EVSI:", v_names_params[p], "\n")
  gams_p <- lmm_list[[p]]   # reuse GAMs fitted in EVPPI step
  for (s in seq_len(n_samples)) {
    Ltilde_list <- lapply(gams_p, function(gm)
      predict.ga(gm, n = v_n[s], n0 = n0_defaults[v_names_params[p]], verbose = FALSE)
    )
    m_Ltilde          <- do.call(cbind, Ltilde_list)
    df_evsi[s, p + 1] <- mean(rowMaxs(m_Ltilde))
  }
}

# Reshape + label
df_evsi_long <- reshape2::melt(df_evsi, id.vars = "N",
                               variable.name = "Parameter", value.name = "EVSI")
df_evsi_long$Label <- param_labels_voi[as.character(df_evsi_long$Parameter)]

df_evppi_ref <- data.frame(
  Parameter = v_names_params,
  Label     = unname(param_labels_voi[v_names_params]),
  EVPPI     = v_evppi
)

## Top 8 by EVPPI
top8         <- df_evppi_ref$Parameter[order(df_evppi_ref$EVPPI, decreasing = TRUE)][1:8]
df_evsi_top  <- df_evsi_long[df_evsi_long$Parameter %in% top8, ]
df_evppi_top <- df_evppi_ref[df_evppi_ref$Parameter %in% top8, ]

lvls <- df_evppi_ref$Label[order(df_evppi_ref$EVPPI, decreasing = TRUE)]
df_evsi_top$Label  <- factor(df_evsi_top$Label,  levels = lvls)
df_evppi_top$Label <- factor(df_evppi_top$Label, levels = lvls)

ggplot(df_evsi_top, aes(x = N, y = EVSI)) +
  geom_line(aes(linetype = "EVSI"), linewidth = 0.9) +
  geom_point(size = 1.5) +
  geom_hline(aes(yintercept = EVPPI, linetype = "EVPPI"),
             data = df_evppi_top, color = "red", linewidth = 0.8) +
  facet_wrap(~ Label, scales = "free_y", ncol = 4) +
  scale_linetype_manual(name = "",
                        values = c("EVSI" = "solid", "EVPPI" = "dashed")) +
  labs(x = "Future study sample size (n)", y = "EVSI ($ per patient)",
       title = "EVSI — Top 8 Parameters (Strategies)",
       subtitle = "Dashed = EVPPI ceiling") +
  scale_x_continuous(n.breaks = 5) +
  scale_y_continuous(labels = dollar, n.breaks = 5) +
  theme_bw(base_size = txtsize) +
  theme(legend.position = "bottom")

# ---- 7. EVSI Study Designs & ENBS --------------------------
disc_pop   <- 0.03
LT_years   <- 10
v_time_pop <- seq(0, LT_years)

# Define functions BEFORE they are called
TotPop <- function(time, prev, incid, disc = 0) {
  pop_time <- c(prev, rep(incid, length(time) - 1))
  sum(pop_time / (1 + disc)^time)
}
CostRes <- function(fixed.cost = 0, samp.size, cost.per.patient,
                    INMB = 0, clin.trial = TRUE, n.arms = 2) {
  if (clin.trial) fixed.cost + n.arms * samp.size * cost.per.patient + samp.size * INMB
  else fixed.cost + samp.size * cost.per.patient
}

# Population parameters for ENBS (millions)
mash_prev  <- 0.300
mash_incid <- 0.040

# Run once per scenario and print all three
for (pop_scenario in list(
  c(prev=0.120, inc=0.015, label="Pediatric only"),
  c(prev=0.660, inc=0.080, label="All MASH"),
  c(prev=0.300, inc=0.040, label="Conservative mid")
)) {
  tot_pop_s <- TotPop(v_time_pop, as.numeric(pop_scenario["prev"]),
                      as.numeric(pop_scenario["inc"]), disc_pop)
  cat(pop_scenario["label"], "— Pop EVSI: $", 
      round(evpi * tot_pop_s, 2), "M\n")
}

tot_pop <- TotPop(v_time_pop, mash_prev, mash_incid, disc_pop)

cat("\n===== Population EVPI =====\n")
cat("Total affected population (discounted, millions):", round(tot_pop, 4), "\n")
cat("Population EVPI ($M):", round(evpi * tot_pop, 2), "\n")

## -- Study A: RCT of Semaglutide --
v_sel_rct  <- which(v_names_params %in% c("rr_sema_regress","rr_sema_progress"))
n0_rct     <- n0_defaults[v_names_params[v_sel_rct]]
v_n_rct <- c(0, 25, 50, 100, 200, 400, 600, 800, 1000, 1500, 2000, 3000, 5000)

cat("\n===== Fitting GAMs for RCT EVSI =====\n")
gams_rct <- lapply(seq_len(n_strat), function(j)
  gam(m_loss[, j] ~
        s(df_params[, v_sel_rct[1]]) +
        s(df_params[, v_sel_rct[2]]) +
        ti(df_params[, v_sel_rct[1]], df_params[, v_sel_rct[2]]))
)
evppi_rct <- mean(rowMaxs(do.call(cbind, lapply(gams_rct, `[[`, "fitted.values"))))

df_evsi_rct <- data.frame(Study = "RCT (Sema efficacy)", N = v_n_rct, EVSI = 0)
for (s in seq_along(v_n_rct)) {
  Ltilde_list <- lapply(gams_rct, function(gm)
    predict.ga(gm, n = v_n_rct[s], n0 = n0_rct, verbose = FALSE))
  df_evsi_rct$EVSI[s] <- mean(rowMaxs(do.call(cbind, Ltilde_list)))
}

## -- Study B: Natural History Cohort --
v_sel_nh <- which(v_names_params %in%
                    c("h_F0_F1","h_F1_F2","h_F2_F3","h_F3_F4","h_F4_DCC",
                      "h_F1_F0","h_F2_F1","h_F3_F2","h_F4_F3"))
n0_nh    <- n0_defaults[v_names_params[v_sel_nh]]
v_n_nh  <- c(0, 50, 100, 250, 500, 750, 1000, 1500, 2000, 3000, 5000, 7500, 10000)

cat("===== Fitting GAMs for Natural History Cohort EVSI =====\n")
smooth_nh   <- paste0("s(df_params[,", v_sel_nh, "])", collapse = " + ")
interact_nh <- paste0("ti(df_params[,", v_sel_nh[1], "], df_params[,", v_sel_nh[2], "])")
gams_nh <- lapply(seq_len(n_strat), function(j)
  gam(as.formula(paste("m_loss[,", j, "] ~", smooth_nh, "+", interact_nh)))
)
evppi_nh <- mean(rowMaxs(do.call(cbind, lapply(gams_nh, `[[`, "fitted.values"))))

df_evsi_nh <- data.frame(Study = "Cohort (fibrosis transitions)", N = v_n_nh, EVSI = 0)
for (s in seq_along(v_n_nh)) {
  Ltilde_list <- lapply(gams_nh, function(gm)
    predict.ga(gm, n = v_n_nh[s], n0 = n0_nh, verbose = FALSE))
  df_evsi_nh$EVSI[s] <- mean(rowMaxs(do.call(cbind, Ltilde_list)))
}

## -- EVSI plot --
df_evppi_designs <- data.frame(
  Study = c("RCT (Sema efficacy)", "Cohort (fibrosis transitions)"),
  EVPPI = c(evppi_rct, evppi_nh)
)
ggplot(bind_rows(df_evsi_rct, df_evsi_nh), aes(x = N, y = EVSI)) +
  geom_line(aes(linetype = "EVSI"), linewidth = 1) +
  geom_point() +
  geom_hline(aes(yintercept = EVPPI, linetype = "EVPPI"),
             data = df_evppi_designs, color = "red", linewidth = 0.9) +
  facet_wrap(~ Study, scales = "free_x") +
  scale_linetype_manual(name = "",
                        values = c("EVSI" = "solid", "EVPPI" = "dashed")) +
  labs(x = "Sample size (n)", y = "EVSI ($ per patient)",
       title = "EVSI by Study Design — Strategies vs. LSM",
       subtitle = "Dashed = EVPPI ceiling for that parameter group") +
  scale_x_continuous(labels = comma, n.breaks = 5) +
  scale_y_continuous(labels = dollar, n.breaks = 6) +
  theme_bw(base_size = txtsize) +
  theme(legend.position = c(0.85, 0.2))

## -- ENBS --
df_evsi_rct$popEVSI <- df_evsi_rct$EVSI * tot_pop
df_evsi_nh$popEVSI  <- df_evsi_nh$EVSI  * tot_pop

v_cost_rct <- CostRes(fixed.cost = 10e6 * 1e-6, samp.size = v_n_rct,
                      cost.per.patient = 30e3 * 1e-6, clin.trial = TRUE, n.arms = 2)
v_cost_nh  <- CostRes(fixed.cost = 500e3 * 1e-6, samp.size = v_n_nh,
                      cost.per.patient = 5e3 * 1e-6, clin.trial = FALSE)

# Cost breakdown table 
cat("\n===== Research Cost Breakdown =====\n")
cat("RCT — Fixed: $10M | Variable: $30K/patient/arm × 2 arms\n")
cat(sprintf("  At n=100:  $%.1fM\n", 10 + 2*100*0.03))
cat(sprintf("  At n=500:  $%.1fM\n", 10 + 2*500*0.03))
cat(sprintf("  At n=1000: $%.1fM\n", 10 + 2*1000*0.03))
cat(sprintf("  At n=5000: $%.1fM\n", 10 + 2*5000*0.03))
cat("Cohort — Fixed: $0.5M | Variable: $5K/patient\n")
cat(sprintf("  At n=1000:  $%.1fM\n",  0.5 + 1000*0.005))
cat(sprintf("  At n=5000:  $%.1fM\n",  0.5 + 5000*0.005))
cat(sprintf("  At n=10000: $%.1fM\n",  0.5 + 10000*0.005))

df_enbs_rct <- merge(df_evsi_rct, data.frame(N = v_n_rct, CS = v_cost_rct), by = "N")
df_enbs_nh  <- merge(df_evsi_nh,  data.frame(N = v_n_nh,  CS = v_cost_nh),  by = "N")
df_enbs_rct$ENBS <- df_enbs_rct$popEVSI - df_enbs_rct$CS
df_enbs_nh$ENBS  <- df_enbs_nh$popEVSI  - df_enbs_nh$CS

oss_rct      <- df_enbs_rct$N[which.max(df_enbs_rct$ENBS)]
oss_nh       <- df_enbs_nh$N[which.max(df_enbs_nh$ENBS)]
max_enbs_rct <- max(df_enbs_rct$ENBS)
max_enbs_nh  <- max(df_enbs_nh$ENBS)

cat("\n===== ENBS =====\n")
cat("RCT  — n* =", oss_rct, "| Max ENBS: $",
    format(round(max_enbs_rct * 1e6), big.mark = ","), "\n")
cat("Cohort — n* =", oss_nh,  "| Max ENBS: $",
    format(round(max_enbs_nh  * 1e6), big.mark = ","), "\n")

df_enbs_rct$nstar <- oss_rct
df_enbs_nh$nstar  <- oss_nh

enbs_all_long <- bind_rows(
  reshape2::melt(df_enbs_rct[, c("Study","N","nstar","popEVSI","CS","ENBS")],
                 id.vars = c("Study","N","nstar"), value.name = "Million"),
  reshape2::melt(df_enbs_nh[,  c("Study","N","nstar","popEVSI","CS","ENBS")],
                 id.vars = c("Study","N","nstar"), value.name = "Million")
) %>% mutate(Study_label = paste0(Study, "\n(n* = ", comma(nstar), ")"))

df_enbs_labels <- enbs_all_long %>%
  group_by(Study_label, variable) %>%
  slice_max(N, n = 1) %>%
  ungroup() %>%
  mutate(line_label = case_when(
    variable == "popEVSI" ~ "Pop. EVSI",
    variable == "CS"      ~ "Research cost",
    variable == "ENBS"    ~ "ENBS"
  ))

ggplot(enbs_all_long, aes(x = N, y = Million, colour = variable, group = variable)) +
  facet_wrap(~ Study_label, scales = "free_x") +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "gray50") +
  geom_vline(aes(xintercept = nstar), linetype = "dashed", colour = "gray50") +
  geom_line(linewidth = 1) +
  geom_point() +
  geom_text(
    data    = df_enbs_labels,
    aes(label = line_label),
    hjust   = -0.08,
    size    = 3.2,
    fontface = "bold",
    show.legend = FALSE
  ) +
  scale_x_continuous(
    "Sample size (N)",
    labels = comma,
    n.breaks = 5,
    expand = expansion(mult = c(0.02, 0.18))  
  ) +
  scale_y_continuous(
    "Value (millions USD)\nNote: popEVSI and Research cost both in $M",
    labels = dollar,
    n.breaks = 6
  ) +
  scale_colour_manual(
    name   = "",
    values = c("popEVSI" = "#0072B2", "CS" = "#D55E00", "ENBS" = "#009E73"),
    labels = c("popEVSI" = "Population EVSI ($M)",
               "CS"      = "Cost of research ($M)",
               "ENBS"    = "ENBS ($M)")
  ) +
  labs(
    title    = "Expected Net Benefit of Sampling — Strategies vs LSM",
    subtitle = "Population: US adolescents with F2/F3 MASH | Disc. rate 3% | All values in millions USD",
    caption  = "Dashed vertical line = optimal sample size (n*). Research cost = fixed + per-patient variable cost."
  ) +
  theme_bw(base_size = txtsize) +
  theme(legend.position = "bottom", panel.spacing = unit(2, "lines"))

# ---- 8. Summary Table --------------------------------------
cat("\n=================================================\n")
cat("  VOI SUMMARY — Strategies vs. LSM\n")
cat("  WTP = $", format(wtp, big.mark = ","), "/QALY\n")
cat("=================================================\n")
cat("Optimal strategy:               ", strategies[d_star], "\n")
cat("EVPI (per patient):             $", format(round(evpi), big.mark = ","), "\n")
cat("Population EVPI (discounted):   $",
    format(round(evpi * tot_pop * 1e6), big.mark = ","), "\n")
cat("\nTop 5 parameters by EVPPI:\n")
top5 <- head(df_evppi[order(df_evppi$EVPPI, decreasing = TRUE), ], 5)
for (i in seq_len(nrow(top5))) {
  cat("  ", i, ". ", as.character(top5$Parameter[i]),
      "  EVPPI = $", format(round(top5$EVPPI[i]), big.mark = ","), "\n", sep = "")
}
cat("\nENBS Optimal Sample Sizes:\n")
cat("  RCT (sema efficacy):       n* =", oss_rct, "\n")
cat("  Natural history cohort:    n* =", oss_nh,  "\n")
cat("=================================================\n")

