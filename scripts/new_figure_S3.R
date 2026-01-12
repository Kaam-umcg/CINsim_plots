library(CINsim)
library(tidyr)
library(reshape2)
library(ggplot2)
library(purrr)
library(dplyr)

# sets wd to location of the script - this is dirty, and I do not 
# recommend it for anyone else following in my footsteps.
# Since this script runs a lot of simulations, it's made for HPC using slurm.
setwd(Sys.getenv("TMPDIR"))

# MAGIC NUMBERS
# sets variables for the simulations
ITERATIONS <- 100
GENERATIONS <- 250
MAX_CELLS <- 10e10

set.seed(42)

# source some utils functions for CnFS and viability
source("scripts/utils/CnFS.R")
source("scripts/utils/check_viability.R")
source("scripts/utils/cinsim_theme.R")

# plotting theme and scale colours for copy numbers
# copy number colors
copy_num_cols <- c("gray90", "darkorchid3", "springgreen2", "red3", "gold2", "navy",
                   "lemonchiffon", "dodgerblue", "chartreuse4")
names(copy_num_cols) <- c("0", "1", "2", "3", "4", "5", "6", "7", "8")

# set the seed for the file, as we do simulations.
# this makes future reruns reproducible!
set.seed(42)

full_results_sweep_results <- readRDS("/scratch/p319788/CINsim/single_sweep_T_ALL/combined_results.Rds")

optimal_params_sim_params <- full_results_sweep_results %>%
    filter(new_CnFS == max(new_CnFS))

optimal_params_sim_fn <- file.path("/scratch/p319788/CINsim/single_sweep_T_ALL",
                                   optimal_params_sim_params$surv_FC,
                                   "tmp", 
                                   "results",
                                   "single_sweep_T_ALL_",
                                   paste0("misseg_", round(optimal_params_sim_params$pMisseg, 7), "_div_FC_", optimal_params_sim_params$division_FCs, ".Rds"))
print(optimal_params_sim_fn)
optimal_params_sim <- readRDS(optimal_params_sim_fn)
# Figure 2d can be recreated by running the code below:
# we keep randomly sampling simulations untill we find a viable one
# most simulations will be viable, but need to be sure by checking
# against the max_g as that is the usual exit we get for good params
# iters is just so that this code doesn't go infinite if something
# went horribly wrong earlier in the script, sloppy but it works
viable_sim <- FALSE
iters <- 0
while (!viable_sim & (iters < length(optimal_params_sim))) {
  selected_sim <- optimal_params_sim[[sample.int(length(optimal_params_sim), size = 1)]]
  # checks whether the true cell count exceeded the max cell count at any point
  # if it was, the simulation was viable
  if (max(selected_sim$gen_measures$true_cell_count) > as.numeric(selected_sim$sim_info[["max_num_cells"]])) {
    viable_sim <- TRUE
  }
  iters <- iters + 1
}

p <- cnvHeatmap(selected_sim, subset_size = 1000)
p <- p + theme(axis.text.x = element_text(size = 15),
          axis.title = element_text(size = 15), aspect.ratio = 1,
          plot.title = element_text(hjust = 0.5, size = 18),
          plot.subtitle = element_text(hjust = 0.5, size = 14)) +
  labs(title = "Karyotype landscape", subtitle = "(optimal p_misseg)") +
  scale_x_discrete(guide = guide_axis(check.overlap = TRUE, n.dodge = 2))
ggsave("plots/new_Fig_S3/fig_S3_D.pdf", plot = p, create.dir = TRUE)

rm(list = ls())
# we also do the same for the other "cube" simulations of T302_p3
full_results_sweep_results <- readRDS("/scratch/p319788/CINsim/single_sweep_T302_p3/combined_results.Rds")

optimal_params_sim_params <- full_results_sweep_results %>%
    filter(new_CnFS == max(new_CnFS))

optimal_params_sim_fn <- file.path("/scratch/p319788/CINsim/single_sweep_T302_p3",
                                   optimal_params_sim_params$surv_FC, 
                                   "tmp", 
                                   "results",
                                   "single_sweep_T302_p3_",
                                   paste0("misseg_", round(optimal_params_sim_params$pMisseg, 7), "_div_FC_", optimal_params_sim_params$division_FCs, ".Rds"))
print(optimal_params_sim_fn)
optimal_params_sim <- readRDS(optimal_params_sim_fn)
# Figure 2d can be recreated by running the code below:
# we keep randomly sampling simulations untill we find a viable one
# most simulations will be viable, but need to be sure by checking
# against the max_g as that is the usual exit we get for good params
# iters is just so that this code doesn't go infinite if something
# went horribly wrong earlier in the script, sloppy but it works
viable_sim <- FALSE
iters <- 0
while (!viable_sim & (iters < length(optimal_params_sim))) {
  selected_sim <- optimal_params_sim[[sample.int(length(optimal_params_sim), size = 1)]]
  # checks whether the true cell count exceeded the max cell count at any point
  # if it was, the simulation was viable
  if (max(selected_sim$gen_measures$true_cell_count) > as.numeric(selected_sim$sim_info[["max_num_cells"]])) {
    viable_sim <- TRUE
  }
  iters <- iters + 1
}

p <- cnvHeatmap(selected_sim, subset_size = 1000)
p <- p + theme(axis.text.x = element_text(size = 15),
          axis.title = element_text(size = 15), aspect.ratio = 1,
          plot.title = element_text(hjust = 0.5, size = 18),
          plot.subtitle = element_text(hjust = 0.5, size = 14)) +
  labs(title = "Karyotype landscape", subtitle = "(optimal p_misseg)") +
  scale_x_discrete(guide = guide_axis(check.overlap = TRUE, n.dodge = 2))
ggsave("plots/new_Fig_S5/fig_S5_B.pdf", plot = p, create.dir = TRUE)