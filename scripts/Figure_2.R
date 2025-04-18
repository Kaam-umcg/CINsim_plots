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

if (!dir.exists("plots/figure_2")){
  dir.create("plots/figure_2", recursive = TRUE)
}

if (!dir.exists("results/T_ALL_params")){
  dir.create("results/T_ALL_params", recursive = TRUE)
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

# set the seed for the file, as we do simulations.
# this makes future reruns reproducible!
set.seed(42)

# Figure 2a can be recreated by running the code below:
#TODO need to reconstruct this directly from the scWGS data
# that will require some bash scripting, and maybe coding the wget
# for the available data as well. Discuss with Floris whether to do that.
Mps1 <- CINsim::Mps1 # import the observed Mps1 copy number frequencies

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
  labs(x = "Chromosome", y = "Frequency", fill = "Copies", title = "382 T-ALL cells",
    subtitle = "(Mps1DK/DK; p53f/f/; Lck-Cre+)") +
  coord_cartesian(ylim = c(0, 1)) +
  cinsim_theme() +
  theme(axis.text.x = element_text(hjust = 1, size = 15),
        axis.text.y = element_text(vjust = 1, size = 15),
        axis.title = element_text(size = 15), aspect.ratio = 1)
ggsave(file = "plots/figure_2/figure_2a.pdf")

# Figure 2b can be recreated by running the code below:

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

  # because some of the survival_FCs are ints we need
  # to cast them to char so that we don't get the wrong value
  idx_surv_FC <- as.character(row$survival_FCs)

  cat("Log: simulating for p_misseg:", row$pMisseg,
  "\nand survival_FC:", row$survival_FCs)
  # runs the simulations for 1 entry in the heatmap. We only vary
  # the pSurvival, as pMisseg is iterated in the loop.
  sim_list <- parallelCinsim(iterations = ITERATIONS,                   
                   cores = 10,
                   karyotypes = NULL,
                   euploid_ref = 2,
                   g = GENERATIONS,
                   max_num_cells = MAX_CELLS,
                   pMisseg = row$pMisseg,
                   selection_mode = "cn_based",
                   selection_metric = Mps1,
                   probability_types = c("pSurvival"),
                   coef = survival_FCs[[idx_surv_FC]],
                   collect_fitness_score = TRUE)
  
  # viability of the combination of p_misseg and survival_FC is defined
  # as the amount of simulations that get above the MAX_CELLS threshold or reach GENERATIONS generations
  surviving_sims <- unlist(map(1:ITERATIONS, ~ check_viability(sim_list, .x, 
                                                threshold_value = MAX_CELLS, max_g = GENERATIONS)))
  row$viability <- sum(surviving_sims)/ ITERATIONS
  row$viable_sims <- sum(surviving_sims)

  # the matching score of the combination of p_misseg and survival_FC is defined
  # as the mean CnFS (Copy number Frequency Score) over all the viable sims +
  # 0 for all non-viable sims
  CnFS <- unlist(lapply(sim_list[surviving_sims], get_CnFS))
  row$non_zero_CnFS <- mean(CnFS)

  row$CnFS <- mean(append(CnFS, rep(0, sum(!surviving_sims))))
  
  # saves the sim results to disk
  sim_name <- paste0("misseg_", round(row$pMisseg, 7),
                     "_surv_FC_", round(row$survival_FCs, 2), ".Rds")
  saveRDS(sim_list, file.path("results/T_ALL_params", sim_name))
  
  # updates the best scoring simulation
  if (row$CnFS > highest_CnFS){
    highest_CnFS <- row$CnFS
    best_sim <- sim_name
  }
  # finally, writes the relevant data to out plotting object
  sim_df[i, ] <- row
}

# saving the metrics of the simulations for later use (if needed)
saveRDS(sim_df, file.path("results/sim_df.Rds"))

# convert some vars to factor for the heatmap
sim_df$survival_FCs <- round(sim_df$survival_FCs, 2)
sim_df$survival_FCs <- as.factor(sim_df$survival_FCs)

# presets p_misseg with subscript
p_lab <- bquote(italic(p)[misseg]) 

# gets the location of highest CnFS
red_star <- sim_df %>%
  filter(CnFS == max(CnFS)) %>%
  select(survival_FCs, pMisseg)

# plotting the heatmap of p_misseg, surv_FC, viability and CnFS
ggplot(sim_df, aes(x = survival_FCs, y = pMisseg)) +
  geom_tile(aes(fill = CnFS, alpha = viability), color = "white", lwd = 0.2, linetype = 1) +
  scale_fill_gradient(low = "blue", high = "yellow", name = "CnFS") +
  scale_alpha_continuous(range = c(0.1, 1), name = "Viability") +
  scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x)),
                breaks = 10^(-6:-1)) +
  scale_x_discrete(name = "Survival FC", breaks = c(10, 5, 3.33, 2.5, 2, 1.67, 1.43, 1.25, 1.11),
                   limits = unique(sim_df$survival_FCs)) +
  geom_point(data = red_star, aes(x = survival_FCs, y = pMisseg), color = "red",
             size = 3, shape = 4) +  
  labs(x = "Survival FC", y = p_lab, 
       title = "Karyotype similarity", subtitle = expression("optimal p_misseg")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15),
        axis.text.y = element_text(vjust = 1, size = 15),
        axis.title = element_text(size = 15),
        aspect.ratio = 1, panel.grid = element_blank())
ggsave("plots/figure_2/figure_2b.pdf")

# Figure 2c can be recreated by running the code below:
# We want to plot the final karyotypes of the last generation for one of the simulations
# at optimal p_misseg and survival_FC.
# This is effectively the location of the red star in plot 2B, so we can take those
# survival_FC and p_misseg and use those for loading the relevant simulations
optimal_params_sim <- readRDS(file.path("results/T_ALL_params", best_sim))

# TODO: hardcoded filepath for saving the best sim for Fig 3E
saveRDS(optimal_params_sim, "/scratch/p319788/CINsim/best_sims/survival.Rds")

# uses the build-in plot_cn function from CINsim to get our initial plot
p <- plot_cn(optimal_params_sim) 

# then we do some manual adding of extra elements, mainly styling text
p + theme(axis.text.x = element_text(hjust = 1, size = 15),
        axis.text.y = element_text(vjust = 1, size = 15),
        axis.title = element_text(size = 15), aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5, size = 18),
        plot.subtitle = element_text(hjust = 0.5, size = 14)) +
    labs(title = "Copy number frequency", subtitle = "(optimal p_misseg)") +
    scale_x_discrete(guide = guide_axis(check.overlap = TRUE, n.dodge = 2))
ggsave("plots/figure_2/figure_2c.pdf", plot = p)

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
p + theme(axis.text.x = element_text(hjust = 1, size = 15),
          axis.title = element_text(size = 15), aspect.ratio = 1,
          plot.title = element_text(hjust = 0.5, size = 18),
          plot.subtitle = element_text(hjust = 0.5, size = 14)) +
  labs(title = "Karyotype landscape", subtitle = "(optimal p_misseg)") +
  scale_x_discrete(guide = guide_axis(check.overlap = TRUE, n.dodge = 2))
ggsave("plots/figure_2/figure_2d.pdf", plot = p)


# Figure 2e can be recreated by running the code below:
# function for extracting the aneu/het scores from a karyosim object
get_karyoscores <- function(karyo_sim){
  last_g <- max(karyo_sim$gen_measures$g)
  last_g_measures <- karyo_sim$pop_measures[karyo_sim$pop_measures$g == last_g, ]$measures[[1]]
  
  # we need to strip the X chromosome from the calculation
  last_g_measures <- last_g_measures[last_g_measures$chromosome != "X", ]
  
  # returns the relevant values
  aneu_score <- mean(last_g_measures$aneuploidy)
  het_score <- mean(last_g_measures$heterogeneity)
  
  return(c(aneu_score, het_score))
}

# gets the het/aneu scores from all optimal param simulations
pop_scores <- map(optimal_params_sim, get_karyoscores)

# and uses whacky R indexing to get the specific scores per type
pop_scores <- unlist(pop_scores)
sim_aneu_scores <- pop_scores[c(TRUE, FALSE)]
sim_het_scores <- pop_scores[c(FALSE, TRUE)]

# Extract the heterogeneity and aneuploidy scores of the observed, which are stored in the CINsim package
observed_karyoscores <- CINsim::karyotype_measures

# creates a dataframe to do the ggplotting
bar_plot_frame <- data.frame(Group = rep(c("Observed", "Simulated"), each = 2),
                             Score = rep(c("Aneuploidy", "Heterogeneity"), times = 2),
                             Value = NA, 
                             Error = NA)

# fills in the dataframe using a ton of conditionals - this is more Python like than R, 
# but it works I suppose
for (i in 1:nrow(bar_plot_frame)){
  row <- bar_plot_frame[i, ]
  if (row$Group == "Simulated"){
    if (row$Score == "Aneuploidy"){
      row$Value <- mean(sim_aneu_scores)
      row$Error <- sd(sim_aneu_scores) / sqrt(length(sim_aneu_scores))
    }
    else {
      row$Value <- mean(sim_het_scores)
      row$Error <- sd(sim_het_scores) / sqrt(length(sim_het_scores))
    }
  }
  else {
    if (row$Score == "Aneuploidy"){
      row$Value <- mean(observed_karyoscores[observed_karyoscores$parameter == "aneuploidy", ]$score)
      # we get the error by getting the error per chromosome - this doesn't feel like the right metric!
      row$Error <- sd(observed_karyoscores[observed_karyoscores$parameter == "aneuploidy", ]$score) / sqrt(length(observed_karyoscores[observed_karyoscores$parameter == "aneuploidy", ]$score))
    }
    else {
      row$Value <- mean(observed_karyoscores[observed_karyoscores$parameter == "heterogeneity", ]$score)
      row$Error <- sd(observed_karyoscores[observed_karyoscores$parameter == "heterogeneity", ]$score) /sqrt(length(observed_karyoscores[observed_karyoscores$parameter == "heterogeneity", ]$score))     
    }
  }
  bar_plot_frame[i, ] <- row
}

# makes the actual plot
p <- ggplot(bar_plot_frame, aes(x = Score, y = Value, fill = Group, group = Group)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  geom_errorbar(
    aes(ymin = Value - Error, ymax = Value + Error), 
    position = position_dodge(width = 0.7), 
    width = 0.25 
  ) +
  labs(title = "Karyotype measures",
       y = "Score") +
  theme_minimal()
ggsave("plots/figure_2/figure_2e.pdf", plot = p)