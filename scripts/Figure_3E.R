library(CINsim)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)

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

best_sims <- list(optimal_sims_surv, optimal_sims_div, optimal_sims_full)
names(best_sims) <- c("Survival", "Division", "Full")

# cleans up some memory space - can this be done now that best_sims exists?
rm(optimal_sims_surv, optimal_sims_div, optimal_sims_full)

# loops over the conditions and fills in their plotting info
# this could be the gnarliest map() statement if properly functioned out
for (condition in names(best_sims)){
  row <- plot_df[plot_df$sim_type == condition , ]
  sims <- best_sims[[condition]]
  # gets generation and cell cap metrics
  MAX_CELLS <- as.numeric(sims[[1]]$sim_info["max_num_cells"])
  MAX_GENERATIONS <- as.numeric(sims[[1]]$sim_info["max_g"])

  surviving_sims <- unlist(map(1:ITERATIONS, ~ check_viability(sims, .x, 
                                              threshold_value = MAX_CELLS, max_g = MAX_GENERATIONS)))

  # calc the CnFS for each individual sim within the list of sims
  #TODO get the error (SEM) for this as well, need to know whether to ignore non-viable
  CnFS <- unlist(lapply(sims, get_CnFS))
  row$CnFS <- mean(CnFS)
  row$CnFS_error <- sd(CnFS) / sqrt(length(CnFS))

  # gets the average number of generations within simulations
  generations <- unlist(map(sims, function(x) as.integer(x$sim_info[["g"]])))
  row$generations <- mean(generations)
  row$generations_error <- sd(generations) / sqrt(length(generations))
  
  # then we move on to the karyotype measures of heterogeneity & aneuploidy
  aneu_scores_per_chrom <- map(sims, function(x) CINsim:::qAneuploidy(karyoMat = x$karyotypes))
  aneu_scores_per_sim <- unlist(map(aneu_scores_per_chrom, function(x) mean(x)))

  row$aneuploidy <- mean(aneu_scores_per_sim)
  row$aneuploidy_error <- sd(aneu_scores_per_sim) / sqrt(length(aneu_scores_per_sim))

  # then we move on to the karyotype measures of heterogeneity & aneuploidy
  het_scores_per_chrom <- map(sims, function(x) CINsim:::qHeterogeneity(karyoMat = x$karyotypes))
  het_scores_per_sim <- unlist(map(het_scores_per_chrom, function(x) mean(x)))

  row$heterogeneity <- mean(het_scores_per_sim)
  row$heterogeneity_error <- sd(het_scores_per_sim) / sqrt(length(het_scores_per_sim))

  plot_df[plot_df$sim_type == condition , ] <- row
}

# now the plot_df is complete, with the exception of the Mps1 statistics
# we extract those from the karyotype_measures object included in CINsim
Mps1_measures <- CINsim::karyotype_measures
row <- plot_df[plot_df$sim_type == "Mps1", ]

# need to skip X chromosome in this calculation - currently it's taken with the result
row$aneuploidy <- mean(unlist(Mps1_measures[Mps1_measures$parameter == "aneuploidy", "score"]))
row$heterogeneity <- mean(unlist(Mps1_measures[Mps1_measures$parameter == "heterogeneity", "score"]))

plot_df[plot_df$sim_type == "Mps1", ] <- row
df_3cond <- filter(plot_df, sim_type %in% c("Survival", "Division", "Full"))

plot_df$sim_type <- factor(plot_df$sim_type)
# now we can actually start making the plots

p1 <- ggplot(df_3cond, aes(x = sim_type, y = CnFS, fill = sim_type)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7, color = "black", show.legend = FALSE) +
  geom_errorbar(aes(ymin = CnFS - CnFS_error, ymax = CnFS + CnFS_error), width = 0.2, color = "black") +
  scale_fill_manual(values = condition_colors) +
  theme_minimal() +
  labs(title = "Karyotype similarity", y = "CnFS")
ggsave("plots/figure_3E/fig_3e.pdf", plot = p1)

p2 <- ggplot(df_3cond, aes(x = sim_type, y = generations, fill = sim_type)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7, color = "black", show.legend = FALSE) +
  geom_errorbar(aes(ymin = generations - generations_error, ymax = generations + generations_error),
   width = 0.2, color = "black") +
  scale_fill_manual(values = condition_colors) +
  theme_minimal() +
  labs(title = "Time", y = "Generations")
ggsave("plots/figure_3E/fig_3f.pdf", plot = p2)

p3 <- ggplot(plot_df, aes(x = sim_type, y = aneuploidy, fill = sim_type)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7, color = "black", show.legend = FALSE) +
  geom_errorbar(aes(ymin = aneuploidy - aneuploidy_error, ymax = aneuploidy + aneuploidy_error),
   width = 0.2, color = "black") +
  scale_fill_manual(values = condition_colors) +
  theme_minimal() +
  labs(title = "Aneuploidy", y = "Score")
ggsave("plots/figure_3E/fig_3g.pdf", plot = p3)

p4 <- ggplot(plot_df, aes(x = sim_type, y = heterogeneity, fill = sim_type)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7, color = "black", show.legend = FALSE) +
  geom_errorbar(aes(ymin = heterogeneity - heterogeneity_error, ymax = heterogeneity + heterogeneity_error),
   width = 0.2, color = "black") +
  scale_fill_manual(values = condition_colors) +
  theme_minimal() +
  labs(title = "Heterogeneity", y = "Score")
ggsave("plots/figure_3E/fig_3h.pdf", plot = p4)