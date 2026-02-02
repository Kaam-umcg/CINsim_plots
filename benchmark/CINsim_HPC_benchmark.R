# gets required variables from SLURM
cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))

# and sets magic variables for the vars we test
TEST_RUNS <- 3
max_g_vars <- c(24, 36, 48, 60, 72, 84, 96, 108)
max_cells_vars <- c(2.5e8, 5e8, 1e9, 2e9, 4e9, 8e9, 1.6e9)
iterations_vars <- c(1, 2, 4, 6, 8, 10, 12, 14, 16)

# loads libraries
library(CINsim)

# wrapper for parallelCINsim
# uses sensible defaults as we got from the Mps1 karyotype data
benchmark_CINsim <- function(cores, generations, iterations, max_cells) {
  # we set sensible defaults for surv_FC and p_misseg based on our
  # parameter search in the main section of the paper - for more info, 
  # check figure 2B. These are effectively magic numbers
  surv_FC <- get_coefficients(selection_metric = Mps1_X, euploid_copy = 2, euploid_survival = 0.7352941, max_survival = 1)
  p_misseg <- 6.21e-3

  time0 <- proc.time()
  parallelCinsim(karyotypes = NULL, 
         g = generations,
         cores = cores,
         iterations = iterations,
         max_num_cells = max_cells,
         selection_metric = Mps1_X,
         probability_types = c("pSurvival"),
         coef = surv_FC,
         pMisseg = p_misseg)

  final_time <- proc.time() - time0
  return(final_time)
}


# benchmark performance for max_g
cat("Starting benchmarking for max generations")
timing_df <- data.frame(
  run = rep(1:TEST_RUNS, each = length(max_g_vars)),
  cores = cores,
  generations = rep(max_g_vars, times = TEST_RUNS),
  time_taken = NA)

# don't want to vectorise this as that might influence performance towards sims
# that get more cores assigned through SNOW
for (i in 1:nrow(timing_df)){
    max_g <- timing_df[i, "generations"]    

    # we set iterations at 8 for this to see the speed up for more cores till 8
    # and then the plateau if we add even more cores (maybe doubles still speeds up?)
    # max_cells is the CINsim default
    time_taken <- benchmark_CINsim(cores = cores, 
                                  generations = max_g, 
                                  iterations = 8,
                                  max_cells = 2e9)

    timing_df[i, "time_taken"] <- time_taken[3]
}

# this path is hardcoded because I want to get this done 
# if you're running this yourself later, I'm so sorry
# if I had more time, I would've written a better script - Mark GPTwain
save_fn <- file.path("~/CINsim_plots/benchmark/results/max_g", paste0(cores, ".Rds"))
saveRDS(timing_df, save_fn)

# benchmark performance for max_cells
cat("Starting benchmarking for max cells")
timing_df <- data.frame(
  run = rep(1:TEST_RUNS, each = length(max_cells_vars)),
  cores = cores,
  real_max_cells = rep(max_cells_vars, times = TEST_RUNS),
  time_taken = NA)

# don't want to vectorise this as that might influence performance towards sims
# that get more cores assigned through SNOW
for (i in 1:nrow(timing_df)){
    max_cells <- timing_df[i, "real_max_cells"]    

    # we set iterations at 8 for this to see the speed up for more cores till 8
    # and then the plateau if we add even more cores (maybe doubles still speeds up?)
    # generations are at 96 due to that being (normally) enough time
    # to reach the max population for good pMisseg and surv_FC
    time_taken <- benchmark_CINsim(cores = cores, 
                                  generations = 96, 
                                  iterations = 8,
                                  max_cells = max_cells)

    timing_df[i, "time_taken"] <- time_taken[3]
}

# this path is hardcoded because I want to get this done 
# if you're running this yourself later, I'm so sorry
# if I had more time, I would've written a better script - Mark GPTwain
save_fn <- file.path("~/CINsim_plots/benchmark/results/max_cells", paste0(cores, ".Rds"))
saveRDS(timing_df, save_fn)


# benchmark performance for iterations
cat("Starting benchmarking for max iterations")
timing_df <- data.frame(
  run = rep(1:TEST_RUNS, each = length(iterations_vars)),
  cores = cores,
  iterations = rep(iterations_vars, times = TEST_RUNS),
  time_taken = NA)

# don't want to vectorise this as that might influence performance towards sims
# that get more cores assigned through SNOW
for (i in 1:nrow(timing_df)){
    iterations <- timing_df[i, "iterations"]    

    # generations are at 96 due to that being (normally) enough time
    # to reach the max population for good pMisseg and surv_FC
    # max_cells is CINsim default
    time_taken <- benchmark_CINsim(cores = cores, 
                                  generations = 96, 
                                  iterations = iterations,
                                  max_cells = 2e9)

    timing_df[i, "time_taken"] <- time_taken[3]
}

# this path is hardcoded because I want to get this done 
# if you're running this yourself later, I'm so sorry
# if I had more time, I would've written a better script - Mark GPTwain
save_fn <- file.path("~/CINsim_plots/benchmark/results/iterations", paste0(cores, ".Rds"))
saveRDS(timing_df, save_fn)