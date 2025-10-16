library(CINsim)
library(tidyr)
library(reshape2)
library(ggplot2)
library(purrr)
library(dplyr)

setwd(Sys.getenv("TMPDIR"))

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
ORGANOID_NAMES <- c("14T", "16T", "9T", "24TB")
# the paths might be different for you if you're running this later!
DIPLOID_RESULTS <- "/scratch/p319788/CINsim/Figure_S4a_diploid/tmp/results"
WGD_RESULTS <- "/scratch/p319788/CINsim/Figure_S4a_WGD/tmp/results"

# we just want to generate plots if the simulation was viable, otherwise
# we output a "non-viable" signal and put that in the paper. That is the best representation
# of the data I can make where we're not just presenting artefacts as the real deal.
# quiver before messed up path management - sorry if you're debugging this
plot_df <- data.frame(
  CnFS = NA, 
  sim_type = rep(c("diploid", "WGD"), times = length(ORGANOID_NAMES)),
  organoid = rep(ORGANOID_NAMES, each = 2),
  load_path = file.path(rep(c(DIPLOID_RESULTS, WGD_RESULTS), times = length(ORGANOID_NAMES)), 
                        paste0(rep(ORGANOID_NAMES, each = 2), ".Rds"))
)

# now we sequentially load the relevant files and fill in the CnFS for each row
for (i in seq_along(1:nrow(plot_df))){
  row <- plot_df[i, ]
  sim_results <- readRDS(row$load_path)
  row$CnFS <- max(sim_results$CnFS, na.rm = TRUE)
  plot_df[i, ] <- row

  # gets the parameters for the best performing sim (CnFS)
  best_sim_row <- sim_results[which.max(sim_results$CnFS), ]
  sim_fn <- paste0(row$sim_type, "_", row$organoid, "_misseg_", 
                    round(best_sim_row$p_misseg, 7), "_div_FC_",
                    best_sim_row$division_FC, "_surv_FC_",
                    best_sim_row$survival_FC, ".Rds")
  
  if (row$sim_type == "WGD"){
    best_sims <- readRDS(file.path(WGD_RESULTS, sim_fn))
  }else{
    best_sims <- readRDS(file.path(DIPLOID_RESULTS, sim_fn))
  }

  # we loop over all the sims that were done and plot them all
  # we select one of them as a representative sample for the article
  for (i in seq(1, 6)){
    best_sim_selected <- best_sims[[i]]

    plot_fn <- paste0(row$organoid, "_sim_", i, "_", row$sim_type, ".pdf")

    p1 <- cnvHeatmap(karyoSim = best_sim_selected, subset_size = 1000)

    p1 <- p1 + labs(
      title = paste0("CN of best performing parameters in organoid ", row$organoid),
      subtitle = element_blank())
    ggsave(file.path("plots/figure_5", plot_fn), plot = p1)
  }
}
