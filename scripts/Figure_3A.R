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

# MAGIC NUMBERS
# sets variables for the simulations
ITERATIONS <- 10
GENERATIONS <- 250 # because we vary the div_FC we need more generations
MAX_CELLS <- 10e10

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
check_viability <- function(sim_list, sim_number = 1:100, threshold_value = 10e10) {
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

# first we recreate all the division_FCs
# because these aren't directly given to CINsim but used as a fit, we need
# to determine the slope and intercept for a specific survival FC
division_FCs <- make_cinsim_coeffcients(selection_metric = Mps1_X,
                                        euploid_copy = 2,
                                        min_survival_euploid = 0.1, #division_FC = 10
                                        max_survival_euploid =  0.9, #division_FC = 1.111.. 
                                        max_survival = 1,
                                        interval = 0.8 / 24,         # takes 25 samples
                                        probability_types = c("pDivision"))
division_FC_vals = (1 / seq(0.1, 0.9, 0.8 / 24))
# sets the names for easy indexing later
names(division_FCs) <- division_FC_vals

# creates the p_missegs we will iterate over
pMissegs <- 10^(seq(-6, -1, length.out = 30))

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
  # runs the simulations for 1 entry in the heatmap. We only vary
  # the pDivision, as pMisseg is iterated in the loop.
  sim_list <- parallelCinsim(iterations = ITERATIONS,                   
                   cores = 10,
                   karyotypes = NULL,
                   euploid_ref = 2,
                   g = GENERATIONS,
                   max_num_cells = MAX_CELLS,
                   pMisseg = row$pMisseg,
                   selection_mode = "cn_based",
                   selection_metric = Mps1,
                   probability_types = c("pDivision"),
                   coef = division_FCs[[idx_div_FC]],
                   fit_division = TRUE,
                   collect_fitness_score = TRUE)
  
  # viability of the combination of p_misseg and division_FC is defined
  # as the amount of simulations that get above the max_cells threshold
  # of 10^10 before 250 generations
  surviving_sims <- unlist(map(1:ITERATIONS, ~ check_viability(sim_list, .x)))
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
ggsave("plots/figure_3/figure_3a.pdf")