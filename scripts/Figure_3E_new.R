library(CINsim)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(ggsignif)

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
# for plotting colours
condition_colors <- c(
  "Survival" = "#1b9e77",
  "Division" = "#d95f02",
  "Full"     = "#7570b3",
  "Mps1"     = "#e7298a"
)

# these paths are hardcoded! Sorry yet again.
# the files contain the list of best sims (red star) for plots 2B, 3A & 3C
optimal_sims_surv <- readRDS("/scratch/p319788/CINsim/best_sims/survival.Rds")
optimal_sims_div <- readRDS("/scratch/p319788/CINsim/best_sims/division.Rds")
optimal_sims_full <- readRDS("/scratch/p319788/CINsim/best_sims/full.Rds")

# concatenates the sims for ease of access
best_sims <- list(optimal_sims_surv, optimal_sims_div, optimal_sims_full)
names(best_sims) <- c("Survival", "Division", "Full")

# cleans up the memory
rm(optimal_sims_surv, optimal_sims_div, optimal_sims_full)

# makes a data frame with all the required entries
plot_df <- data.frame(
  CnFS = NA,
  generations = NA,
  aneuploidy = NA,
  heterogeneity = NA,
  sim_type = rep(c("Survival", "Division", "Full", "Mps1"), each = 100)
)

i <- 1

# need to extract all the information from the simulation data
for (condition in names(best_sims)){
    # indicates which slice of the dataframe we're filling atm
    slice <- (1 + (i-1) * 100):(i * 100)

    sims <- best_sims[[condition]]
    CnFS <- lapply(sims, get_CnFS)
    plot_df[slice, "CnFS"] <- CnFS


    # iterate counter for proper slicing
    i <- i + 1
}