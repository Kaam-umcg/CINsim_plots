library(CINsim)
library(ggplot2)
library(tidyr)
library(dplyr)

# Processor	Intel(R) Core(TM) Ultra 7 155U   1.70 GHz
# Installed RAM	32,0 GB (31,5 GB usable)
# System type	64-bit operating system, x64-based processor

# Does some benchmarking computational demands of CINsim
# I'll start with local, we'll do single core first

benchmark_CINsim <- function(cores, generations, iterations = 24) {
  parallelCinsim(karyotypes = NULL, 
         g = generations,
         cores = cores,
         iterations = iterations)
}

cores_to_test <- c(1, 2, 4, 6, 8, 10, 12)
generations_to_test <- c(12, 24, 48, 72, 96, 144, 192, 288, 384)
test_runs <- 3

df_list <- vector("list", length = test_runs)
for (i in 1:test_runs) {
  df_list[[i]] <- data.frame(matrix(0, nrow = length(cores_to_test), ncol = length(generations_to_test)))
  colnames(df_list[[i]]) <- paste0(generations_to_test,"_gens")
  rownames(df_list[[i]]) <- paste0(cores_to_test, "_cores")
}

mean_df <- data.frame(matrix(0, nrow = length(cores_to_test), ncol = length(generations_to_test)))
sd_df <- data.frame(matrix(0, nrow = length(cores_to_test), ncol = length(generations_to_test)))
colnames(mean_df) <- paste0(generations_to_test,"_gens")
rownames(mean_df) <- paste0(cores_to_test, "_cores")
colnames(sd_df) <- paste0(generations_to_test,"_gens")
rownames(sd_df) <- paste0(cores_to_test, "_cores")

# Benchmarks CINsim: check whether the suite is correctly set up atm
# Mostly getting good performance on 6 cores because
# 6 iterations are being run in parallel so getting perfect 1:1
# Want to see the performance when there are more cores than sims
for (run in 1:test_runs){
  cat("Starting benchmarking suite", run, "\n")
  for (gens in generations_to_test) {
    cat("Testing for", gens, "generations of CINsim\n")
    for (cores in cores_to_test) {
      cat("Testing for", cores, "cores\n")
      time0 <- proc.time()
      benchmark_CINsim(cores, gens)
      final_time <- proc.time() - time0
      
      df_list[[run]][paste0(cores, "_cores"), paste0(gens, "_gens")] <- final_time[3]
    }
  }
}

# Summarizes the data by calcing the mean and sd for plotting later
for (gens in colnames(df_list[[1]])){
  for (cores in rownames(df_list[[1]])){
    time_avg <- vector("numeric", test_runs)
    time_sd <- 0
    for (run in 1:test_runs) {
      time_avg[[run]] <- df_list[[run]][cores, gens]
    }
    time_sd <- sd(time_avg)
    time_avg <- sum(time_avg) / test_runs
    
    mean_df[cores, gens] <- time_avg
    sd_df[cores, gens] <- time_sd
  }
}

# Makes plots for visualising the benchmarking performances
# starts by reshaping the data to long: I need to get used to 
# doing this from the start
mean_df$cores <- rownames(mean_df)

# ChatGPT code
# Reshape the mean data from wide to long format using reshape
long_df_mean <- reshape(mean_df, 
                   varying = names(mean_df)[1:8], # Selecting the generation columns (12_gens, 24_gens, ...)
                   v.names = "mean_time",         # The values in the columns will be stored in 'mean_time'
                   times = gsub("_gens", "", names(mean_df)[1:8]),  # Extract the number of generations
                   timevar = "generations",       # Name of the new column for generations
                   direction = "long")            # Specify long format

# Clean up the columns
long_df_mean$cores <- gsub("_cores", "", long_df_mean$cores)
long_df_mean$`384_gens` <- NULL
long_df_mean$id <- NULL
rownames(long_df_mean <- 1:length(rownames(long_df_mean)))

# repeat for the sd dataframe
sd_df$cores <- rownames(sd_df)

# ChatGPT code
# Reshape the mean data from wide to long format using reshape
long_df_sd <- reshape(sd_df, 
                   varying = names(sd_df)[1:8], # Selecting the generation columns (12_gens, 24_gens, ...)
                   v.names = "mean_time",         # The values in the columns will be stored in 'mean_time'
                   times = gsub("_gens", "", names(sd_df)[1:8]),  # Extract the number of generations
                   timevar = "generations",       # Name of the new column for generations
                   direction = "long")            # Specify long format

# Clean up the columns
long_df_sd$cores <- gsub("_cores", "", long_df_sd$cores)
long_df_sd$`384_gens` <- NULL
long_df_sd$id <- NULL
rownames(long_df_sd) <- 1:length(rownames(long_df_sd))

# merge the dataframes
data <- merge(long_df_mean, long_df_sd, by = c("cores", "generations"))
colnames(data) <- c("cores", "generations", "mean_time", "sd_time")

# now we need to do some sorting because I was sloppy with the data types
data$cores <- as.integer(data$cores)
data$generations <- as.integer(data$generations)

# Need to fix the x axis labelling in this plot
ggplot(data, aes(x = factor(cores), y = mean_time, color = factor(generations), group = factor(generations))) +
  geom_line() +  # Plot lines for each generation
  geom_point() +  # Add points at each (core, mean_time) intersection
  geom_errorbar(aes(ymin = mean_time - sd_time, ymax = mean_time + sd_time), width = 0.2) +  # Add error bars
  scale_x_discrete(labels = factor(cores_to_test)) +  # Set x-axis to have the cores values
  labs(x = "Number of Cores", y = "Time (s)", color = "Generations") +  # Axis labels
  theme_minimal() +  # Clean theme
  theme(legend.position = "top")  # Place legend at the top

