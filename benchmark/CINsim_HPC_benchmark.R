# add the CINsim benchmarking code here, make it plyable 
# to run with different amount of cores by doing it Sys.getenv
# for the amount of cores. Then use parallelCinsim to get it all going
# need to install CINsim into this version of R for that to work as well.

# gets required variables from SLURM
cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))

# and the magic variable of replicates
TEST_RUNS <- 3

# loads libraries
library(CINsim)

# starts with running parallelCINsim with some basic parameters
benchmark_CINsim <- function(cores, generations, iterations = 24) {
  parallelCinsim(karyotypes = NULL, 
         g = generations,
         cores = cores,
         iterations = iterations)
}

timing <- vector(length = TEST_RUNS)

for (run in 1:TEST_RUNS){
    time0 <- proc.time()
    benchmark_CINsim(cores, generations = 100, iterations = 24)

    final_time <- proc.time() - time0
    timing[run] <- final_time[3]
}


saveRDS(timing, paste0("timing_", cores, "_cores", ".Rds"))