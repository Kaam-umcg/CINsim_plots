library(CINsim)
library(tidyr)
library(reshape2)
library(ggplot2)
library(purrr)
library(dplyr)

# sets wd to location of the script - this is dirty, and I do not 
# recommend it for anyone else following in my footsteps.
setwd(Sys.getenv("TMPDIR"))

if (!dir.exists("plots/figure_3C")){
  dir.create("plots/figure_3C", recursive = TRUE)
}

if (!dir.exists("results/T_ALL_params_3C")){
  dir.create("results/T_ALL_params_3C", recursive = TRUE)
}

cores_avail <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))

# MAGIC NUMBERS
# sets variables for the simulations
ITERATIONS <- 100
GENERATIONS <- 250 # because we vary the div_FC we need more generations
MAX_CELLS <- 10e10

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

sim_df <- readRDS("/scratch/p319788/CINsim/single_sweep_T_ALL/1.11/tmp/results/single_sweep_T_ALL_/full_sweep_T_ALL_results.Rds")

# convert some vars to factor for the heatmap
sim_df$division_FCs <- round(sim_df$division_FCs, 2)
sim_df$division_FCs <- as.factor(sim_df$division_FCs)

# presets p_misseg with subscript
p_lab <- bquote(italic(p)[misseg]) 

# gets the location of highest CnFS
red_star <- sim_df %>%
  filter(new_CnFS == max(new_CnFS)) %>%
  select(division_FCs, pMisseg)

# plotting the heatmap of p_misseg, div_FC, viability and CnFS
ggplot(sim_df, aes(x = division_FCs, y = pMisseg)) +
  geom_tile(aes(fill = new_CnFS, alpha = viability), color = "white", lwd = 0.2, linetype = 1) +
  scale_fill_gradient(low = "blue", high = "yellow", name = "CnFS") +
  scale_alpha_continuous(range = c(0.1, 1), name = "Viability") +
  scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x)),
                breaks = 10^(-6:-1)) +
  scale_x_discrete(name = "Division FC", breaks = c(10, 5, 3.33, 2.5, 2, 1.67, 1.43, 1.25, 1.11),
                   limits = unique(sim_df$division_FCs)) +
  geom_point(data = red_star, aes(x = division_FCs, y = pMisseg), color = "red",
             size = 3, shape = 4) +  
  labs(x = "Division FC", y = p_lab, 
       title = "Karyotype similarity", subtitle = expression("optimal p_misseg")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 15),
        aspect.ratio = 1, panel.grid = element_blank())
ggsave("plots/figure_3C/fig_3c.pdf")
