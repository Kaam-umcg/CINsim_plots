library(CINsim)
library(tidyr)
library(reshape2)
library(ggplot2)
library(purrr)
library(dplyr)
# sets wd to location of the script - this is dirty, and I do not 
# recommend it for anyone else following in my footsteps.
setwd(Sys.getenv("TMPDIR"))

if (!dir.exists("plots/figure_3")){
  dir.create("plots/figure_3", recursive = TRUE)
}

if (!dir.exists("results/T_ALL_params")){
  dir.create("results/T_ALL_params", recursive = TRUE)
}

# MAGIC NUMBERS
# sets variables for the simulations
ITERATIONS <- 100
GENERATIONS <- 200 # because we vary the div_FC we need more generations

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
# this makes future reruns reproducible!
set.seed(42)

# define function for checking viability later
check_viability <- function(sim_list, sim_number = 1:100, threshold_value = 2e9) {
  # Access the last value of the true_cell_count
  last_value <- sim_list[[paste0("sim_", sim_number)]]$gen_measures$true_cell_count[length(sim_list[[paste0("sim_", sim_number)]]$gen_measures$true_cell_count)]
  
  # Check if it's greater than the threshold
  return(last_value > threshold_value)
}

# define function to calculate the CnFS, as it's unclear to me 
# if the package does this at any point
get_CnFS <- function(sim_results) {
    if (!class(sim_results) == "karyoSim"){
    message("Incorrect type for given object, needs to be karyoSim")
    return(FALSE)
    }
  # inits the inverse score to 0, smaller = better
  inverse_CnFS <- 0
  
  # gets the observed (target) karyotype frequencies
  karyo_obs <- sim_results$selection_metric

  # and the copy number frequencies of the final generation
  final_g <- max(sim_results$pop_measures$g)
  karyo_sim_final <- sim_results$pop_measures[sim_results$pop_measures$g == final_g, ]
  
  # calculates the CnFS (Copy number Frequency Score) for the final generation
  for (chrom in colnames(karyo_obs)){
    # early exit for X, since we skip that chrom for the CnFS
    if (chrom == "X"){
      break
    }

    # cast to int as we need it for comparisons
    chrom <- as.integer(chrom)
    chrom_df <- karyo_sim_final[karyo_sim_final$chromosome == chrom, ]
    
    for (cn in rownames(karyo_obs)){
      freq_obs <- karyo_obs[as.character(cn), as.character(chrom)]
      
      # test whether there is a value for cn in chrom_df,
      # if not we assign it as 0.
      if (!(as.integer(cn) %in% chrom_df$copy)){
        freq_sim <- 0
      }else{
        freq_sim <- chrom_df[chrom_df$copy == as.integer(cn), ]$cn_freq[[1]][["freq"]]
      }
      
      chrom_cn_CnFS <- (freq_sim - freq_obs)^2
      inverse_CnFS <- inverse_CnFS + chrom_cn_CnFS
    }
  }
  CnFS <- 1 / inverse_CnFS
  return(CnFS)
}

# Figure 3A can be reproduced using the code below: 
