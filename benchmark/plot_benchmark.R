library(CINsim)
library(ggplot2)
library(dplyr)

# hardcoded paths because I love making this someone else their problem I suppose
base_data_path <- "~/CINsim_plots/benchmark/results"
base_plot_path <- "~/CINsim_plots/benchmark/plots"

if (!(dir.exists(base_plot_path))){
    dir.create(base_plot_path)
}
benchmark_vars <- basename(list.dirs(base_data_path, recursive = FALSE))

# collect all the data together
for(var in benchmark_vars){
  data_folder <- file.path(base_data_path, var)
  # init null
  merge_frame <- NULL
  for(file in list.files(data_folder)){
    loaded_file <- readRDS(file.path(data_folder, file))
    
    # sets the first entry to merging frame
    if (is.null(merge_frame)){
      merge_frame <- loaded_file
      next
    }
    # merges every next frame into the merging frame
    merge_frame <- rbind(merge_frame, loaded_file)
  }
  # gets the param name from the dataframe - not consistent with var
  param <- setdiff(colnames(merge_frame), c("run", "cores", "time_taken"))[1]
  
  # summarises for plotting with ggplot
  df_summary <- merge_frame %>%
    group_by(cores, !!sym(param)) %>%
    summarise(
      mean_time = mean(time_taken),
      se_time = sd(time_taken) / sqrt(n()),
      .groups = "drop")
  
  # create line plot with benchmark_var on the x axis, time taken on y and error bars
  p <- ggplot(df_summary, aes(x = !!sym(param), y = mean_time, color = factor(cores), group = cores)) +
    geom_line(size = 1) +
    geom_point() +
    geom_errorbar(aes(ymin = mean_time - se_time, ymax = mean_time + se_time), width = 0.2) +
    labs(
      x = var,
      y = "Mean Time Taken",
      color = "Cores",
      title = "Simulation Time by Iterations and Number of Cores"
    ) +
    theme_minimal()
  
  # sets the saving path for the plot
  plot_fn <- file.path(base_plot_path, paste0(param ,".pdf"))
  # saves the plot and then we move on to the next var in the next iteration of the loop
  ggsave(plot_fn, plot = p)
}

