library(CINsim)
library(ggplot2)
library(tidyr)
library(purrr)

# MAGIC NUMBERS
GENERATIONS <- 250
MAX_CELLS <- 5e10
START_CELLS <- 10
CORES_AVAIL <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
ITERATIONS <- CORES_AVAIL * 2

sim_type <- Sys.getenv("SIM_TYPE")

# changes wd to node specific folder
TMP_DIR <- Sys.getenv("TMPDIR")
setwd(TMP_DIR)

plot_dir <- file.path(
  TMP_DIR,
  sim_type,
  "plots"
)

dir.create(
  plot_dir,
  recursive = TRUE)

output_dir <- file.path(
  TMP_DIR,
  sim_type,
  "results"
)

sim_dir <- file.path(
  output_dir,
  "sims"
)

dir.create(
  sim_dir, 
  showWarnings = FALSE,
  recursive = TRUE)

# sets seed for reproducibility
set.seed(42)

# gets the organoid CNs, as we need those to determine our parameters
organoid_cns <- readRDS("~/CINsim_plots/data/organoid_cn_freqs.Rds")

# determines WGD or diploid starting karyotype
if (sim_type == "WGD"){
  aneuploid_base <- 2
  start_karyotype <- makeKaryotypes(
    numCell = START_CELLS, 
    species = "human", 
    copies = 4)
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
  p_missegs <- 10^seq(-4, -1, length.out = 10)
  division_FCs <- make_cinsim_coeffcients(
    selection_metric = organoid_cn,
    euploid_copy = aneuploid_base,
    min_survival_euploid = 0.1,
    max_survival_euploid = 0.9,
    probability_types = c("pDivision"),
    interval = 0.1)

  division_FC_vals <- round((1 / seq(0.1, 0.9, 0.1)), 2)
  names(division_FCs) <- division_FC_vals
  survival_FCs <- make_cinsim_coeffcients(
    selection_metric = organoid_cn,
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

  # does every simulation 1 by 1, doing ITERATIONS repeats of each set of params
  for (i in 1:nrow(sim_data_frame)){
    row <- sim_data_frame[i, ]

    # compiles the coefficients with the specified div/surv_FCs
    # we need the as.character, as otherwise ints will try to index out of bounds
    coeff_struct <- survival_FCs[[as.character(row$survival_FC)]]
    coeff_struct$pDivision <- division_FCs[[as.character(row$division_FC)]]$pDivision

    # does the simulation for the given parameters in row
    sim_list <- parallelCinsim(selection_metric = organoid_cn,
                            selection_mode = "cn_based",
                            iterations = ITERATIONS, 
                            cores = CORES_AVAIL, 
                            karyotypes = start_karyotype, 
                            euploid_ref = aneuploid_base,
                            g = GENERATIONS, 
                            max_num_cells = MAX_CELLS,
                            coef = coeff_struct, 
                            pMisseg = row$p_misseg,
                            fit_misseg = FALSE, 
                            fit_division = TRUE,
                            pDivision = 0, 
                            CnFS = TRUE,
                            monosomy_penalty = TRUE, 
                            penalty_fraction = 0.1, 
                            collect_fitness_score = TRUE)

    # saves the sim results to disk
    sim_name <- paste0(sim_type, "_", organoid_name, "_misseg_", round(row$p_misseg, 7),
                    "_div_FC_", row$division_FC,
                    "_surv_FC_", row$survival_FC, ".Rds")

    saveRDS(
      sim_list, 
      file.path(
        sim_dir, 
        sim_name))

    # uses a util function to determine which simulations were viable (max_cells or max_generations)
    surviving_sims <- unlist(map(1:ITERATIONS, ~ check_viability(sim_list, .x, 
                                                                threshold_value = MAX_CELLS, 
                                                                max_g = GENERATIONS)))
    row$viability <- sum(surviving_sims)/ ITERATIONS
    
    CnFS <- unlist(lapply(sim_list[surviving_sims], get_final_CnFS))
    row$CnFS <- mean(append(CnFS, rep(0, sum(!surviving_sims))))

    # prints some summary statistics for checking in the logs
    cat("final_CnFS per simulation values of: \n")
    print(CnFS)

    # resets the row back into the main dataframe
    sim_data_frame[i, ] <- row
  }

  # at the end of the loop, we save the sim_data_frame with viability and CnFS for later plotting
  saveRDS(
    object = sim_data_frame, 
    file =  file.path(
      output_dir,
      paste0(
        organoid_name,
        "_",
        sim_type,
        "_sim_df.Rds"
      )
    )
    )
}
