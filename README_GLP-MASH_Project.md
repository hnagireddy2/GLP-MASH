# HRP 263: GLP-1 MASH Progression Model

Markov cost-effectiveness analysis comparing semaglutide vs lifestyle modification (LSM) for adolescent MASH (metabolic-associated steatohepatitis).

---

## File Structure

### Primary Model File

| File | Description |
|------|-------------|
| `HRP_263_GLP-1_MASH_Progression_Model_UPDATED_3.11.26.Rmd` | Complete R Markdown model. Contains all code from calibration through VOI analysis. Knits to PDF via XeLaTeX. |

### Required Input Files

All three should be placed in `~/Downloads/` before running.

| File | Description |
|------|-------------|
| `age_mort_background_costs.csv` | Age-specific overall mortality rates and background healthcare costs |
| `baseline_transition_probs.monthly.csv` | Monthly baseline transition probabilities between MASH health states |
| `GA_functions.R` | GAM-based prediction functions for EVSI computation (`predict.ga`, `Predict.smooth.ga`, `Predict.matrix.tensor.smooth.ga`) |

### Cache Files (Auto-Generated)

On first run, these `.rds` files are created automatically in the working directory. On subsequent runs, they are loaded from disk to skip expensive computations. Cache file must be deleted to force that section to re-run.

| File | First-Run Time | Contents |
|------|---------------|----------|
| `calibration_results.rds` | ~30–60 min | IMIS posterior samples, MAP/mean estimates, calibrated `v_init` and `ped_mult_calib` |
| `psa_results_base.rds` | ~45–90 min | Base-case PSA (3,000 sims): cost and QALY matrices for LSM vs Semaglutide 1.5y |
| `psa_results_all.rds` | ~60–120 min | All-strategy PSA (1,000 sims): cost and QALY matrices for all 9 strategies |

---

## How to Run

### Dependencies

```r
install.packages(c(
  "dplyr", "tidyr", "ggplot2", "ggrepel", "scales", "reshape2",    # core
  "dampack",                                                       # CEA
  "truncnorm", "matrixStats", "lhs", "gtools",                     # PSA sampling
  "IMIS",                                                          # calibration
  "mgcv",                                                          # VOI (GAMs)
  "GGally"                                                         # diagnostics
))
```

XeLaTeX is required for PDF output (specified in the YAML header).

### Steps

1. Place all three input files in `~/Downloads/`.
2. Open the `.Rmd` in RStudio.
3. Click **Knit to PDF**.
4. **First run:** expect 2–4 hours total as calibration and both PSA loops run and cache.
5. **Subsequent runs:** with cache files present, completes in ~15–30 minutes (OWSA, threshold SA, and VOI still run each time).

### Expected Runtime Breakdown

| Section | First Run | With Cache |
|---------|-----------|------------|
| IMIS Calibration | 30–60 min | Instant |
| Base Case + OWSA | 5–10 min | 5–10 min |
| Dampack NMB OWSA (100 samples) | 10–15 min | 10–15 min |
| Threshold SA (bisection) | 5–10 min | 5–10 min |
| Base PSA (3,000 sims) | 45–90 min | Instant |
| All-Strategy PSA (1,000 sims) | 60–120 min | Instant |
| VOI (EVPPI + EVSI) | 10–20 min | 10–20 min |
| Validation + Plots | 5–10 min | 5–10 min |

---

## Key Outputs

### Figures

| # | Figure |
|---|--------|
| 1 | Calibration: Model-Predicted vs Literature Targets |
| 2 | Posterior Predictive Check (violin + boxplot) |
| 3 | Prior vs Posterior distributions for calibrated parameters |
| 4 | Posterior correlation matrix (pairs plot) |
| 5 | CE Plane: LSM vs Semaglutide (1.5y base case) |
| 6 | NMB OWSA: All 9 strategies by parameter |
| 7 | Optimal strategy heatmap by parameter value |
| 8 | Threshold Sensitivity Analysis: Non-dominated strategies |
| 9 | State occupancy traces: LSM and Semaglutide by duration |
| 10 | Combined CEA Frontier: All 9 strategies (deterministic) |
| 11 | PSA parameter distributions: All sampled parameters |
| 12 | PSA scatter: All 9 strategies |
| 13 | Post-PSA Combined CEA Frontier: All strategies |
| 14 | CEAC: All 9 strategies |
| 15 | Expected Loss Curve: All 9 strategies |
| 16 | NMB distributions: All 9 strategies |
| 17 | Incremental NMB vs LSM: All strategies |
| 18 | EVPI across WTP thresholds |
| 19 | EVPPI: Individual parameters |
| 20 | EVPPI: Parameter groups |
| 21 | EVSI curves: Top 8 parameters by sample size |
| 22 | EVSI by study design (RCT vs Cohort) |
| 23 | ENBS: Expected Net Benefit of Sampling |
| 24 | Validation: State occupancy over time (LSM) |
| 25 | Validation: % F3+ trajectory with posterior credible interval |
| 26 | Validation: Cumulative mortality |

### Tables

| Table | Description |
|-------|-------------|
| IMIS Posterior Summary | Mean, 95% CrI, and MAP for all 5 calibrated parameters |
| Combined CEA Frontier | All 9 strategies with incremental costs, QALYs, ICERs, and dominance status |
| Post-PSA CEA Frontier | Probabilistic frontier with dominance classification |
| VOI Summary | EVPI, top 5 EVPPI parameters, optimal sample sizes |
| External Validation | 8 epidemiological targets with model predictions and pass/fail |
