library(CINsim)
library(tidyr)
library(reshape2)
library(ggplot2)
library(purrr)
library(dplyr)

setwd(Sys.getenv("TMPDIR"))

# this pattern matches the string for the organoid names
# pattern <- "^[0-9]+[A-Za-z]+\\.Rds$"

# creates dir to store our figures in later
if (!dir.exists("plots/figure_5")){
  dir.create("plots/figure_5", recursive = TRUE)
}

if (!dir.exists("results")){
  dir.create("results", recursive = TRUE)
}

# MAGIC NUMBERS
# sets variables for the simulations
GENERATIONS <- 250
ITERATIONS <- 100
MAX_CELLS <- 5e10

# and some names of organoids and locations of sim results
# the paths might be different for you if you're running this later!
ORGANOID_NAMES <- c("14T", "16T", "9T", "24TB")
DIPLOID_RESULTS <- "/scratch/p319788/CINsim/Figure_S4a_diploid/tmp/results"
WGD_RESULTS <- "/scratch/p319788/CINsim/Figure_S4a_WGD/tmp/results"

# gets the observed CN freqs from scWGS
organoid_cn_freqs <- readRDS("data/organoid_cn_freqs.Rds")

# iteration variable so we can store and make the patchwork from a single object later
sim_counter <- 1

# makes objects to store the sims and plots into
sim_object <- vector("list", length = length(ORGANOID_NAMES) * 2)
plot_object <- vector("list", length = length(ORGANOID_NAMES) * 2)

# defines the names for the plots/sims
plot_names <- unlist(lapply(ORGANOID_NAMES, function(x) paste0(x, c("_diploid", "_WGD"))))
names(sim_object) <- plot_names
names(plot_object) <- plot_names

plot_df <- data.frame(
    CnFS = NA, 
    sim_type = rep(c("diploid", "WGD"), times = 4),
    organoid = rep(ORGANOID_NAMES, each = 2)
)

for (sample_name in ORGANOID_NAMES){
    # reconstructs the filenames for the results.Rds file
    organoid_cn <- organoid_cn_freqs[[sample_name]]
    
    # get the correct hyperparams for the sims from the S4 results, starting with diploid
    sim_results_diploid <- readRDS(file.path(DIPLOID_RESULTS, paste0(sample_name ".Rds")))
    best_params_diploid <- sim_results_diploid[which.max(sim_results_diploid$CnFS), ]

    # if we load the old simulation, we can immediately get the coeffs
    diploid_sims <- readRDS("diploid", "_misseg_", round(best_params_diploid$p_misseg, 7),
                   "_div_FC_", best_params_diploid$division_FC,
                   "_surv_FC_", best_params_diploid$survival_FC, ".Rds")
    sim_object[[sim_counter]] <- diploid_sims

    # saves the CnFS for plot B
    plot_df[sim_counter, "CnFS"] <- best_params_diploid$CnFS

    # select a the first sim for the cnvHeatmap
    p <- cnvHeatmap(diploid_sims[[1]])
    ggsave(filename = paste0("results/", sample_name, "_diploid.pdf"), plot = p)
    plot_object[[sim_counter]] <- p

    # iters the sim counter
    sim_counter <- sim_counter + 1

    sim_results_WGD <- readRDS(file.path(WGD_RESULTS, paste0(sample_name ".Rds")))
    best_params_WGD <- sim_results_WGD[which.max(sim_results_WGD$CnFS), ]

    # if we load the old simulation, we can immediately get the coeffs
    WGD_sims <- readRDS("WGD", "_misseg_", round(best_params_diploid$p_misseg, 7),
                   "_div_FC_", best_params_diploid$division_FC,
                   "_surv_FC_", best_params_diploid$survival_FC, ".Rds")
    sim_object[[sim_counter]] <- WGD_sims
    
    plot_df[sim_counter, "CnFS"] <- best_params_WGD$CnFS

    # select a the first sim for the cnvHeatmap
    p <- cnvHeatmap(WGD_sims[[1]])
    ggsave(filename = paste0("results/", sample_name, "_WGD.pdf"), plot = p)
    plot_object[[sim_counter]] <- p

    # iters the sim counter
    sim_counter <- sim_counter + 1
}

# dumps the sim and plots objects for patchwork later
saveRDS(sim_object, "results/sims_figure_5.Rds")
saveRDS(plot_object, "results/plots_figure_5.Rds")

# makes plot B
ggplot(plot_df, aes(x = organoid, y = CnFS, fill = sim_type)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Organoid", y = "CnFS", fill = "sim_type") +
  theme_minimal()
ggsave("plot/figure_5/fig_5b.pdf")