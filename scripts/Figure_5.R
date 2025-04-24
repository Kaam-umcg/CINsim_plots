library(CINsim)
library(tidyr)
library(reshape2)
library(ggplot2)
library(purrr)
library(dplyr)

setwd(Sys.getenv("TMPDIR"))

cores_avail <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))

# this pattern matches the string for the organoid names
pattern <- "^[0-9]+[A-Za-z]+\\.Rds$"

