library(CINsim)
library(tidyr)
library(reshape2)
library(ggplot2)
library(purrr)

# sets wd to location of the script - this is dirty, and I do not 
# recommend it for anyone else following in my footsteps.
# Only works in Rstudio
#setwd(dirname(rstudioapi::getSourceEditorContext()$path))
setwd(Sys.getenv("TMPDIR"))
cat(getwd())

if (!dir.exists("plots/figure_2")){
  dir.create("plots/figure_2", recursive = TRUE)
}
if (!dir.exists("results/T_ALL_params")){
  dir.create("results/T_ALL_params", recursive = TRUE)
}

# MAGIC NUMBERS
# sets variables for the simulations
ITERATIONS <- 10
GENERATIONS <- 100

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
  labs(x = "Chromosome", y = "Frequency", fill = "Copies") +
  coord_cartesian(ylim = c(0, 1)) +
  cinsim_theme()
ggsave(file = "plots/figure_2/figure_2a.pdf")


# Figure 2b can be recreated by running the code below:
########### This code is effectively hyperparameter tuning
# WARNING # Read: it takes very long (80 hours for 10 cores) to run.
########### I would recommend running this particular script on HPC

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
  # inits the score to 0, smaller = better
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

# first we recreate all the survival_FCs
# because these aren't directly given to CINsim but used as a fit, we need
# to determine the slope and intercept for a specific survival FC
survival_FCs <- make_cinsim_coeffcients(selection_metric = Mps1_X,
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

# we create an object to store the results in and init some NA values
sim_df <- data.frame(
  pMisseg = rep(pMissegs, each = length(survival_FC_vals)),
  survival_FCs = rep(survival_FC_vals, times = length(pMissegs)),
  viability = NA, 
  CnFS = NA)

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
                   pMisseg = row$pMisseg,
                   selection_mode = "cn_based",
                   selection_metric = Mps1,
                   probability_types = c("pSurvival"),
                   coef = survival_FCs[[idx_surv_FC]],
                   collect_fitness_score = TRUE)
  
  # viability of the combination of p_misseg and survival_FC is defined
  # as the amount of simulations that get above the max_cells threshold
  # of 2*10^9 before 100 generations
  surviving_sims <- unlist(map(1:ITERATIONS, ~ check_viability(sim_list, .x)))
  row$viability <- sum(surviving_sims)/ ITERATIONS
  
  # the matching score of the combination of p_misseg and survival_FC is defined
  # as the mean CnFS (Copy number Frequency Score) over all the viable sims +
  # 0 for all non-viable sims
  CnFS <- unlist(lapply(sim_list[surviving_sims], get_CnFS))
  row$CnFS <- mean(append(CnFS, rep(0, sum(!surviving_sims))))
  
  # saves the sim results to disk
  sim_name <- paste0("misseg_", round(row$pMisseg, 7),
                     "_surv_FC_", round(row$survival_FCs, 2), ".Rds")
  saveRDS(sim_list, file.path("results/T_ALL_params", sim_name))
  
  # finally, writes the relevant data to out plotting object
  sim_df[i, ] <- row
}

# saving the metrics of the simulations for later use (if needed)
saveRDS(sim_df, file.path("results/sim_df.Rds"))

# Plotting the heatmap
ggplot(sim_df, aes(x = survival_FC, y = )) +
  geom_title(aes(fill = CnFS, alpha = viability), color = "white", width = 0.9, height = 0.9) + 
  scale_fill_gradient(low = "white", high = "blue", name = "Similarity") +
  scale_alpha_continuous(range = c(0.2, 1), guide = "none") +
  labs(x = "Survival FC", y = "p_misseg", title = "Karyotype similarity") +
  cinsim_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("plots/figure_2/figure_2b.pdf")