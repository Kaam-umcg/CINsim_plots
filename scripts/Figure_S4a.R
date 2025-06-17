library(CINsim)
library(ggplot2)
library(tidyr)
library(purrr)

# MAGIC NUMBERS
GENERATIONS <- 250
MAX_CELLS <- 5e10
ITERATIONS <- 6
START_CELLS <- as.integer(Sys.getenv("START_CELLS"))
CORES_AVAIL <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
results_dir <- file.path("results", as.character(START_CELLS))

# changes wd to node specific folder
setwd(Sys.getenv("TMPDIR"))

cat("Using number of starting cells:\n")
print(START_CELLS)

if (!dir.exists("plots/figure_S4a")){
  dir.create("plots/figure_S4a", recursive = TRUE)
}

if (!dir.exists(results_dir)){
  dir.create(results_dir, recursive = TRUE)
}

# sources all required util scripts - part of CINsim package in full release
source("scripts/utils/check_viability.R")

# sets seed for reproducibility
set.seed(42)

# gets the organoid CNs, as we need those to determine our parameters
organoid_cns <- readRDS("data/organoid_cn_freqs.Rds")

# determines WGD or diploid starting karyotype
sim_type <- Sys.getenv("SIM_TYPE")
if (sim_type == "WGD"){
  aneuploid_base <- 2 # used to be 4, but has huge effect on coefficients
  start_karyotype <- makeKaryotypes(numCell = START_CELLS, species = "human", copies = 4)
  cat("Starting for a standard CN of 4")
}else if (sim_type == "diploid"){
  aneuploid_base <- 2
  start_karyotype <- makeKaryotypes(numCell = START_CELLS, species = "human", copies = 2)
  cat("Starting for a standard CN of 2")
}else{
  cat("Incorrect sim_type - needs to be 'WGD' or 'diploid'")
}

#TODO - due to a bug somewhere in CINsim (issue#9) I need to set the rownames
# of all the cells to the same to prevent crashing at the end of a sim
rownames(start_karyotype) <- rep("cell_1", START_CELLS)

# loop over all the unique organoid cn_frequencies
for (organoid_name in names(organoid_cns)){
    organoid_cn <- organoid_cns[[organoid_name]]

    # gets all the required simulation parameters
    p_missegs <- 10^seq(-6, -1, length.out = 10)
    division_FCs <- make_cinsim_coeffcients(selection_metric = organoid_cn,
                                            euploid_copy = aneuploid_base,
                                            min_survival_euploid = 0.1,
                                            max_survival_euploid = 0.9,
                                            probability_types = c("pDivision"),
                                            interval = 0.1)
    division_FC_vals <- round((1 / seq(0.1, 0.9, 0.1)), 2)
    names(division_FCs) <- division_FC_vals
    survival_FCs <- make_cinsim_coeffcients(selection_metric = organoid_cn,
                                            euploid_copy = aneuploid_base,
                                            min_survival_euploid = 0.6,
                                            max_survival_euploid = 0.9,
                                            probability_types = c("pSurvival"),
                                            interval = 0.1)
    survival_FC_vals <- round((1 / seq(0.6, 0.9, 0.1)), 2)
    names(survival_FCs) <- survival_FC_vals

    # creates a dataframe with all the simulations that need to be done
    # for this specific organoid
    sim_data_frame <- data.frame(
    CnFS = NA,
    viability = NA,
    p_misseg = rep(p_missegs, each = length(division_FC_vals) * length(survival_FC_vals)),
    division_FC = rep(division_FC_vals, times = length(survival_FC_vals)), 
    survival_FC = rep(survival_FC_vals, each = length(division_FC_vals)),
    sim_type = sim_type)

    cat("\nLog: need to do a total of", dim(sim_data_frame)[1], "simulation for organoid", 
          organoid_name, "for condition", sim_type, "\n")

    # does every simulation 1 by 1, doing 6 repeats of each sim
    for (i in 1:nrow(sim_data_frame)){
        row <- sim_data_frame[i, ]

        # compiles the coefficients with the specified div/surv_FCs
        # we need the as.character, as otherwise ints will try to index out of bounds
        coeff_struct <- survival_FCs[[as.character(row$survival_FC)]]
        coeff_struct$pDivision <- division_FCs[[as.character(row$division_FC)]]$pDivision

        # does the simulation for the given parameters in row
        sim_list <- parallelCinsim(selection_metric = organoid_cn,
                                selection_mode = "cn_based",
                                iterations = ITERATIONS, cores = CORES_AVAIL, 
                                karyotypes = start_karyotype, 
                                euploid_ref = aneuploid_base,
                                g = GENERATIONS, max_num_cells = MAX_CELLS,
                                coef = coeff_struct, pMisseg = row$p_misseg,
                                fit_misseg = FALSE, fit_division = TRUE,
                                pDivision = 0, CnFS = TRUE, 
                                collect_fitness_score = TRUE)

        # saves the sim results to disk
        sim_name <- paste0(sim_type, "_misseg_", round(row$p_misseg, 7),
                        "_div_FC_", row$division_FC,
                        "_surv_FC_", row$survival_FC, ".Rds")
        saveRDS(sim_list, file.path("results", sim_name))

        # uses a util function to determine which simulations were viable (max_cells or max_generations)
        surviving_sims <- unlist(map(1:ITERATIONS, ~ check_viability(sim_list, .x, 
                                                                    threshold_value = MAX_CELLS, 
                                                                    max_g = GENERATIONS)))
        row$viability <- sum(surviving_sims)/ ITERATIONS

        # gets both the CnFS of all viable sims and the score is non_surviving sims are assigned 0
        # CnFS <- unlist(lapply(sim_list[surviving_sims], get_CnFS))
        # row$CnFS <- mean(CnFS)
        # row$zeroed_CnFS <- mean(append(CnFS, rep(0, sum(!surviving_sims))))
        
        final_CnFS <- sapply(1:ITERATIONS, function(i) {
          sim_list[[i]]$gen_measures$CnFS[length(sim_list[[i]]$gen_measures$CnFS)]
          })
        row$CnFS <- mean(final_CnFS)

        # prints some summary statistics for checking in the logs
        cat("Results for this row are: \n")
        print(row)
        cat("With final_CnFS per simulation values of: \n")
        print(final_CnFS)
        cat("And the following coefficients: \n")
        print(coeff_struct)

        # resets the row back into the main dataframe
        sim_data_frame[i, ] <- row
    }

    # at the end of the loop, we save the sim_data_frame with viability and CnFS for later plotting
    saveRDS(sim_data_frame, paste0(results_dir, "/", organoid_name, ".Rds"))
}
