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
BASE_DATA_PATH <- "/scratch/p319788/CINsim/Figure_2_backup/tmp/results/T_ALL_params"


if (!dir.exists("plots/n_generations")){
  dir.create("plots/n_generations", recursive = TRUE)
}
if (!dir.exists("results/n_generations")){
  dir.create("results/n_generations", recursive = TRUE)
}

# MAGIC NUMBERS
# sets variables for the simulations
ITERATIONS <- 100
GENERATIONS <- 100
MAX_CELLS <- 2e9

# source some utils functions for CnFS and viability
source("scripts/utils/CnFS.R")
source("scripts/utils/check_viability.R")
source("scripts/utils/cinsim_theme.R")

# plotting theme and scale colours for copy numbers
# copy number colors
copy_num_cols <- c("gray90", "darkorchid3", "springgreen2", "red3", "gold2", "navy",
                   "lemonchiffon", "dodgerblue", "chartreuse4")
names(copy_num_cols) <- c("0", "1", "2", "3", "4", "5", "6", "7", "8")


# first we recreate all the survival_FCs
# because these aren't directly given to CINsim but used as a fit, we need
# to determine the slope and intercept for a specific survival FC
survival_FCs <- make_cinsim_coeffcients(selection_metric = Mps1,
                                        euploid_copy = 2,
                                        min_survival_euploid = 0.1, #survival_FC = 10
                                        max_survival_euploid =  0.9, #survival_FC = 1.111.. 
                                        max_survival = 1,
                                        interval = 0.8 / 24,         # takes 25 samples
                                        probability_types = c("pSurvival"))
survival_FC_vals = (1 / seq(0.1, 0.9, 0.8 / 24))
# sets the names for easy indexing later
names(survival_FCs) <- survival_FC_vals

# creates the p_missegs we will iterate over
pMissegs <- 10^(seq(-6, -1, length.out = 30))

# sets the base CnFS - we keep track of the highest for plotting and loading
# back some simulation data later
highest_CnFS <- 0
best_sim <- "no best sim yet"

# we create an object to store the results in and init some NA values
sim_df <- data.frame(
  pMisseg = rep(pMissegs, each = length(survival_FC_vals)),
  survival_FCs = rep(survival_FC_vals, times = length(pMissegs)),
  viability = NA, 
  viable_sims = NA,
  avg_generations = NA)

# I am unsure how foreach works with code that uses snow,
# so I opt for doing the simulations by row.
# this also allows me to save the intermediary result sometimes
# which is recommended for such a chunky piece of code.
# not very R standard, but such is life (sorry!)
for (i in 1:nrow(sim_df)){
  row <- sim_df[i,]

  # because some of the survival_FCs are ints we need
  # to cast them to char so that we don't get the wrong value
  idx_surv_FC <- as.character(row$survival_FCs)

  cat("Log: loading for p_misseg:", row$pMisseg,
  "\nand survival_FC:", row$survival_FCs)

  # we don't simulate, we load the already done simulation based on the metrics in sim_df
  sim_name <- paste0("misseg_", round(row$pMisseg, 7),
                     "_surv_FC_", round(row$survival_FCs, 2), ".Rds")

  sim_list <- readRDS(file.path(BASE_DATA_PATH, sim_name))

  # viability of the combination of p_misseg and survival_FC is defined
  # as the amount of simulations that get above the MAX_CELLS threshold or reach GENERATIONS generations
  surviving_sims <- unlist(map(1:ITERATIONS, ~ check_viability(sim_list, .x, 
                                                threshold_value = MAX_CELLS, max_g = GENERATIONS)))
  row$viability <- sum(surviving_sims)/ ITERATIONS
  row$viable_sims <- sum(surviving_sims)

  row$avg_generations <- sim_list %>%
    map_dbl(~ max(.x$gen_measures$g)) %>%
    mean()

  # the matching score of the combination of p_misseg and survival_FC is defined
  # as the mean CnFS (Copy number Frequency Score) over all the viable sims +
  # 0 for all non-viable sims

  # finally, writes the relevant data to out plotting object
  sim_df[i, ] <- row
  print(row)
}

# convert some vars to factor for the heatmap
sim_df$survival_FCs <- round(sim_df$survival_FCs, 2)
sim_df$survival_FCs <- as.factor(sim_df$survival_FCs)

# presets p_misseg with subscript
p_lab <- bquote(italic(p)[misseg]) 

# plotting the heatmap of p_misseg, surv_FC, viability and CnFS
ggplot(sim_df, aes(x = survival_FCs, y = pMisseg)) +
  geom_tile(aes(fill = avg_generations, alpha = viability), color = "white", lwd = 0.2, linetype = 1) +
  scale_fill_gradient(low = "blue", high = "yellow", name = "Generations") +
  scale_alpha_continuous(range = c(0.1, 1), name = "Viability") +
  scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x)),
                breaks = 10^(-6:-1)) +
  scale_x_discrete(name = "Survival FC", breaks = c(10, 5, 3.33, 2.5, 2, 1.67, 1.43, 1.25, 1.11),
                   limits = unique(sim_df$survival_FCs)) +
  labs(x = "Survival FC", y = p_lab, 
       title = "Karyotype similarity", subtitle = expression("optimal p_misseg")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 15),
        aspect.ratio = 1, panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))
ggsave("plots/n_generations/n_generations_viability_heatmap.pdf")


ggplot(sim_df, aes(x = survival_FCs, y = pMisseg)) +
  geom_tile(aes(fill = avg_generations), color = "white", lwd = 0.2, linetype = 1) +
  scale_fill_gradient(low = "blue", high = "yellow", name = "Generations") +
  scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x)),
                breaks = 10^(-6:-1)) +
  scale_x_discrete(name = "Survival FC", breaks = c(10, 5, 3.33, 2.5, 2, 1.67, 1.43, 1.25, 1.11),
                   limits = unique(sim_df$survival_FCs)) +
  labs(x = "Survival FC", y = p_lab, 
       title = "Karyotype similarity", subtitle = expression("optimal p_misseg")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 15),
        aspect.ratio = 1, panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))
ggsave("plots/n_generations/n_generations_heatmap.pdf")

saveRDS(sim_df, "results/sim_df_generations.Rds")