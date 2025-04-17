library(CINsim)
library(dplyr)
library(tidyr)

# sets wd to location of the script - this is dirty, and I do not 
# recommend it for anyone else following in my footsteps.
setwd(Sys.getenv("TMPDIR"))

if (!dir.exists("plots/figure_3E")){
  dir.create("plots/figure_3E", recursive = TRUE)
}

if (!dir.exists("results/T_ALL_params_3E")){
  dir.create("results/T_ALL_params_3E", recursive = TRUE)
}

# make plots comparing the CnFS, # of generations, aneuploidy, and heterogeneity between
# simulated with 100% p_surv, simulated karyotype dependant p_surv, and real metrics
optimal_sims_surv <- readRDS()
optimal_sims_div <- readRDS()
optimal_sims_full <- readRDS()

# Starts with the CnFS barplot
# makes a data frame with all the required entries
