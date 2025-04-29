library(ggplot2)
library(CINsim)

# sets wd to node specific folder
setwd(Sys.getenv("TMPDIR"))

# extracts the cores assigned in .sh script - doesn't work for multinode jobs!
cores_avail <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))

if (!dir.exists("plots/figure_S2")){
  dir.create("plots/figure_S2", recursive = TRUE)
}

# figure S2a

# sets parameters for the sims
p_misseg <- 0.0025
GENERATIONS <- 250

