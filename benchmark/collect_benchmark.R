library(CINsim)
library(ggplot2)

# we want to collect the information per variable we perturbed
# and make line plots from all that information
base_data_path <- "~/CINsim_plots/benchmark/results"
benchmark_vars <- list.dirs("results", recursive = FALSE)

# collect all the data together
for(var in benchmark_vars){
    data_folder <- file.path(base_data_path, benchmark_vars)
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

    # makes a plot 
    p <- ggplot() +
            geom_errorbar() +
            theme() +
            

    # sets the saving path for the plot
    plot_fn <- file.path()

    # saves the plot and then we move on to the next var
}