# main.R — Sources all scripts in dependency order
# Run this file to execute the full analysis

source("00_parameters.R")
source("00b_le_transitions.R")
source("01_model_functions.R")
source("02_calibration.R")
source("03_validation.R")
source("04_base_case.R")
source("05_owsa.R")
source("05b_elasticity.R")
source("06_psa.R")
source("07_psa_postprocess.R")
source("08_voi.R")
