library(CINsim)
library(tidyr)
library(reshape2)
library(ggplot2)

# sets wd to location of the script - this is dirty, and I do not 
# recommend it for anyone else following in my footsteps.
# Only works in Rstudio
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

if (!dir.exists("../plots/figure_2")){
  dir.create("../plots/figure_2")
}

# plotting theme and scale colours for copy numbers
# copy number colors
copy_num_cols <- c("gray90", "darkorchid3", "springgreen2", "red3", "gold2", "navy",
                   "lemonchiffon", "dodgerblue", "chartreuse4")
names(copy_num_cols) <- c("0", "1", "2", "3", "4", "5", "6", "7", "8")

# plotting theme
cinsim_theme <- function() {
  theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text = element_text(colour = "black"),
          aspect.ratio = 1)
}
# set the seed for the file, as we do simulations.
# this makes future reruns reproducable!
set.seed(42)

# Figure 2a can be recreated by running the code below:
# loads the observed (target) karyotype
# More specifically, loads the observed frequencies of chromosomes
#TODO need to reconstruct this directly from the scWGS data
load("../CINsim/data/Mps1.RData")

# reshapes for ggplotting
melt_Mps1 <- melt(Mps1)
colnames(melt_Mps1) <- c("Copy_number", "Chromosome", "Fraction")
melt_Mps1$Copy_number <- factor(melt_Mps1$Copy_number)

# plots the observed CN states per chromosome in the CINsim style
ggplot(melt_Mps1, aes(x = Chromosome, y = Fraction, 
                      fill = factor(Copy_number, levels = rev(levels(Copy_number))))) + # have to reverse the levels
  geom_bar(stat = "identity", position = "stack") +  # Stacked bars
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  scale_fill_manual(values = copy_num_cols) +  # Unique colors for each row
  labs(x = "Chromosome", y = "Frequency", fill = "Copies") +
  coord_cartesian(ylim = c(0, 1)) +
  cinsim_theme()

# Figure 2b can be recreated by running the code below:
# first we recreate all the survival_FCs
# because these aren't directly given to CINsim but used as a fit, we need
# to determine the slope and intercept for a specific survival FC
survival_FCs <- make_cinsim_coeffcients(selection_metric = Mps1_X,
                                        euploid_copy = 2,
                                        min_survival_euploid = 0.1, #survival_FC = 10
                                        max_survival_euploid =  0.9, #survival_FC = 1.111. 
                                        max_survival = 1,
                                        interval = 0.8 / 24,         # takes 25 samples
                                        probability_types = c("pSurvival"))

# sets variables for the simulations
n_reps <- 6
pMissegs <- 10^(seq(-6, -1, length.out = 30))
generations <- 100

test_sim <- parallelCinsim(iterations = 6,
                           cores = 6,
                           karyotypes = NULL,
                           euploid_ref = 2,
                           g = generations,
                           pMisseg = pMissegs[20],
                           selection_mode = "cn_based",
                           selection_metric = Mps1,
                           probability_types = c("pSurvival"),
                           coef = survival_FCs[[20]],
                           collect_fitness_score = TRUE)

# a simulation is viable if the 
# max($gen_measures$true_cell_count) > 2*10^9 (max_num_cells) BOOL
# then take the ratio of surviving sims for that pMisseg and survival_FC
# the pMisseg can be found in $sim_info["pMisseg"] float
# survival_FC = index of (1 / seq(0.1, 0.9, 0.8 / 24))
# 



# goes over all the pMissegs and survival_FC coefficients and
# simulates 6 replicates for each combination.
# This is a total of 30 * 25 * 6 simulations, with a running time
# of ~330(!) hours for a desktop computer on 6 cores.
# core efficiency diminishes with more cores than iterations (sims),
# see the benchmarking in Supplementary for more information.

for (i in 1:n_reps) {
  sim_data <- parallelCinsim(iterations = 6,
                 cores = 6,
                 karyotypes = NULL,
                 euploid_ref = 2,
                 g = generations,
                 pMisseg = pMissegs,
                 selection_mode = "cn_based",
                 selection_metric = Mps1)
  save(sim_data, file = file.path("outputs"), paste0("replicate_", i, ".Rds"))
  
  # collect simulation info
  sim_info <- tibble(sim_label = paste0("sim_", 1:length(sim_data)),
                     p_misseg = sim_data %>% map_dbl(function(x) {
                       p_misseg <- x$sim_info[["pMisseg"]]
                       return(as.numeric(p_misseg))
                     }))
  
}


# sample for plotting the heatmap:
df <- data.frame(
  p_test = rep(c("A", "B", "C", "D", "E"), each = 4),
  surf = rep(c("V", "W", "X", "Y", "Z"), 4),
  similarity = c(0.8, 0.7, 0.9, 0.4, 0.3, 0.5, 0.8, 0.1, 0.6, 0.3, 0.1, 0.9, 0.8, 0.2, 0.7, 0.5, 0.8, 0.9, 0.2, 0.3),
  vibes = c(0.6, 0.3, 0.1, 0.9, 0.8, 0.2, 0.7, 0.5, 0.6, 0.3, 0.1, 0.9, 0.8, 0.2, 0.7, 0.5, 0.8, 0.9, 0.2, 0.3)
)

# Plotting the heatmap
ggplot(df, aes(x = surf, y = p_test)) +
  geom_tile(aes(fill = similarity, alpha = vibes), color = "white", width = 0.9, height = 0.9) + 
  scale_fill_gradient(low = "white", high = "blue", name = "Similarity") +  # Color for similarity
  scale_alpha_continuous(range = c(0.2, 1), guide = "none") +  # Transparency for vibes
  labs(x = "Surface", y = "Test", title = "Heatmap of Similarity and Vibes") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))