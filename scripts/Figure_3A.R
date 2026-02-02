library(CINsim)
library(tidyr)
library(reshape2)
library(ggplot2)
library(purrr)
library(dplyr)

# sets wd to location of the script - this is dirty, and I do not 
# recommend it for anyone else following in my footsteps.
setwd(Sys.getenv("TMPDIR"))

if (!dir.exists("plots/figure_3A")){
  dir.create("plots/figure_3A", recursive = TRUE)
}

if (!dir.exists("results/T_ALL_params_div")){
  dir.create("results/T_ALL_params_div", recursive = TRUE)
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

# Figure 3A can be reproduced using the code below: 

# first we recreate all the division_FCs
# because these aren't directly given to CINsim but used as a fit, we need
# to determine the slope and intercept for a specific survival FC
division_FCs <- make_cinsim_coeffcients(selection_metric = Mps1,
                                        euploid_copy = 2,
                                        min_survival_euploid = 0.1, #  division_FC = 10
                                        max_survival_euploid =  0.9, # division_FC = 1.111.. 
                                        max_survival = 1,
                                        interval = 0.8 / 24,         # takes 25 samples
                                        probability_types = c("pDivision"))
division_FC_vals = (1 / seq(0.1, 0.9, 0.8 / 24))
# sets the names for easy indexing later
names(division_FCs) <- division_FC_vals

# creates the p_missegs we will iterate over
pMissegs <- 10^(seq(-6, -1, length.out = 30))

# sets the base CnFS - we keep track of the highest for plotting and loading
# back some simulation data later
highest_CnFS <- 0
best_sim <- "no best sim yet"

# sets the a and b values for the pSurvival - since we want this to be karyotype
# independant, we set a = 0 and b > 1000 so that every cell always survives, effectively
# making this simulation more about clonal expansion than selection
coeff_struct <- list(NULL, NULL)
names(coeff_struct) <- c("pDivision", "pSurvival")
coeff_struct$pSurvival <- c(0, 1000)
names(coeff_struct$pSurvival) <- c("a", "b")

# we create an object to store the results in and init some NA values
sim_df <- data.frame(
  pMisseg = rep(pMissegs, each = length(division_FC_vals)),
  division_FCs = rep(division_FC_vals, times = length(pMissegs)),
  viability = NA, 
  CnFS = NA,
  viable_sims = NA,
  non_zero_CnFS = NA)

# I am unsure how foreach works with code that uses snow,
# so I opt for doing the simulations by row.
# this also allows me to save the intermediary result sometimes
# which is recommended for such a chunky piece of code.
# not very R standard, but such is life (sorry!)
for (i in 1:nrow(sim_df)){
  row <- sim_df[i,]

  # because some of the division_FCs are ints we need
  # to cast them to char so that we don't get the wrong value
  idx_div_FC <- as.character(row$division_FCs)

  cat("Log: simulating for p_misseg:", row$pMisseg,
  "\nand division_FC:", row$division_FCs)

  # collects the pDivision that changed with the consistent pSurvival
  coeff_struct$pDivision <- division_FCs[[idx_div_FC]]$pDivision
  coeff_struct$pMisseg <- division_FCs[[idx_div_FC]]$pMisseg

  # runs the simulations for 1 entry in the heatmap. We only vary
  # the pDivision, as pMisseg is iterated in the loop.
  sim_list <- parallelCinsim(iterations = ITERATIONS,                   
                   cores = cores_avail,
                   karyotypes = NULL,
                   euploid_ref = 2,
                   g = GENERATIONS,
                   max_num_cells = MAX_CELLS,
                   pMisseg = row$pMisseg,
                   selection_mode = "cn_based",
                   selection_metric = Mps1,
                   coef = coeff_struct,
                   fit_division = TRUE,
                   collect_fitness_score = TRUE,
                   pDivision = 0) # funny quirk in the CINsim code, but pDivision needs to be 0
                                  # otherwise there will be no division at all.
  
  # viability of the combination of p_misseg and division_FC is defined
  # as the amount of simulations that get above the max_cells threshold
  # of 10^10 before 250 generations
  surviving_sims <- unlist(map(1:ITERATIONS, ~ check_viability(sim_list, .x, 
                                                threshold_value = MAX_CELLS, max_g = GENERATIONS)))
  row$viability <- sum(surviving_sims)/ ITERATIONS
  row$viable_sims <- sum(surviving_sims)

  # the matching score of the combination of p_misseg and division_FC is defined
  # as the mean CnFS (Copy number Frequency Score) over all the viable sims +
  # 0 for all non-viable sims
  CnFS <- unlist(lapply(sim_list[surviving_sims], get_CnFS))
  row$non_zero_CnFS <- mean(CnFS)

  row$CnFS <- mean(append(CnFS, rep(0, sum(!surviving_sims))))
  
  # saves the sim results to disk
  sim_name <- paste0("misseg_", round(row$pMisseg, 7),
                     "_div_FC_", round(row$division_FCs, 2), ".Rds")
  saveRDS(sim_list, file.path("results/T_ALL_params_div", sim_name))
  
  # updates the best scoring simulation
  if (row$CnFS > highest_CnFS){
    highest_CnFS <- row$CnFS
    best_sim <- sim_name
  }
  # finally, writes the relevant data to out plotting object
  sim_df[i, ] <- row
}

# saving the metrics of the simulations for later use (if needed)
saveRDS(sim_df, file.path("results/sim_div_df.Rds"))

# convert some vars to factor for the heatmap
sim_df$division_FCs <- round(sim_df$division_FCs, 2)
sim_df$division_FCs <- as.factor(sim_df$division_FCs)

# presets p_misseg with subscript
p_lab <- bquote(italic(p)[misseg]) 

# gets the location of highest CnFS
red_star <- sim_df %>%
  filter(CnFS == max(CnFS)) %>%
  select(division_FCs, pMisseg)

# plotting the heatmap of p_misseg, div_FC, viability and CnFS
ggplot(sim_df, aes(x = division_FCs, y = pMisseg)) +
  geom_tile(aes(fill = CnFS, alpha = viability), color = "white", lwd = 0.2, linetype = 1) +
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
        axis.text.y = element_text(vjust = 1, size = 15),
        axis.title = element_text(size = 15),
        aspect.ratio = 1, panel.grid = element_blank())
ggsave("plots/figure_3A/fig_3a.pdf")

# we then create all the other relevant plots for this section of Figure 3, specifically 3B & C
# as they're based on the results of this simulation
# For figure 3B:
optimal_params_sim <- readRDS(file.path("results/T_ALL_params_div", best_sim))
saveRDS(optimal_params_sim, "/scratch/p319788/CINsim/best_sims/division.Rds")

# uses the build-in plot_cn function from CINsim to get our initial plot
p <- plot_cn(optimal_params_sim) 

# then we do some manual adding of extra elements, mainly styling text
p + theme(axis.text.x = element_text(hjust = 1, size = 15),
        axis.text.y = element_text(vjust = 1, size = 15),
        axis.title = element_text(size = 15), aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5, size = 18),
        plot.subtitle = element_text(hjust = 0.5, size = 14)) +
    labs(title = "Copy number frequency", subtitle = "(optimal Division FC and p_misseg)") +
    scale_x_discrete(guide = guide_axis(check.overlap = TRUE, n.dodge = 2))
ggsave("plots/figure_3A/fig_3b1.pdf", plot = p)

# Figure 3C can be recreated by running the code below:
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
          axis.title = element_text(hjust = 0.5, size = 15), aspect.ratio = 1,
          plot.title = element_text(hjust = 0.5, size = 18),
          plot.subtitle = element_text(hjust = 0.5, size = 14)) +
  labs(title = "Karyotype landscape", subtitle = "(optimal Division FC and p_misseg)") +
  scale_x_discrete(guide = guide_axis(check.overlap = TRUE, n.dodge = 2))
ggsave("plots/figure_3A/fig_3b2.pdf", plot = p)