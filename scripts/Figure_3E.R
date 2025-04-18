library(CINsim)
library(dplyr)
library(tidyr)

# sets wd to location of the script - this is dirty, and I do not 
# recommend it for anyone else following in my footsteps.
setwd(Sys.getenv("TMPDIR"))

if (!dir.exists("plots/figure_3E")){
  dir.create("plots/figure_3E", recursive = TRUE)
}

# source some utils functions for CnFS and viability
source("scripts/utils/CnFS.R")
source("scripts/utils/check_viability.R")
source("scripts/utils/cinsim_theme.R")

# MAGIC VARIABLES
ITERATIONS <- 100

# these paths are hardcoded! Sorry yet again.
# the files contain the list of best sims (red star) for plots 2B, 3A & 3C
optimal_sims_surv <- readRDS("/scratch/p319788/CINsim/best_sims/survival.Rds")
optimal_sims_div <- readRDS("/scratch/p319788/CINsim/best_sims/division.Rds")
optimal_sims_full <- readRDS("/scratch/p319788/CINsim/best_sims/full.Rds")

# makes a data frame with all the required entries
plot_df <- data.frame(
  CnFS = NA,
  CnFS_error = NA,
  generations = NA,
  generations_error = NA,
  aneuploidy = NA,
  aneuploidy_error = NA,
  heterogeneity = NA,
  heterogeneity_error = NA,
  sim_type = c("Survival", "Division", "Full", "Mps1"))

best_sims <- c(optimal_sims_surv, optimal_sims_div, optimal_sims_full)
names(best_sims) <- c("Survival", "Division", "Full")

# loops over the conditions fills in their plotting info
for (condition in names(best_sims)){
  row <- plot_df[plot_df$sim_type == condition , ]
  sims <- best_sims[condition]

  #TODO
  MAX_CELLS <- 
  GENERATIONS <- 


  surviving_sims <- unlist(map(1:ITERATIONS, ~ check_viability(sims, .x, 
                                              threshold_value = MAX_CELLS, max_g = GENERATIONS)))

  # calc the CnFS for each individual sim within the list of sims
  CnFS <- unlist(lapply(sims, get_CnFS))
  row$
}